import os
from math import sqrt
import numpy as np
import torch
from functions.ckpt_util import download
from functions.denoising import efficient_generalized_steps
from guided_diffusion.script_util import create_model
import mat73
from scipy.io import savemat


def get_beta_schedule(beta_schedule, *, beta_start, beta_end, num_diffusion_timesteps):
    def sigmoid(x):
        return 1 / (np.exp(-x) + 1)

    if beta_schedule == "quad":
        betas = (
                np.linspace(
                    beta_start ** 0.5,
                    beta_end ** 0.5,
                    num_diffusion_timesteps,
                    dtype=np.float64,
                )
                ** 2
        )
    elif beta_schedule == "linear":
        betas = np.linspace(
            beta_start, beta_end, num_diffusion_timesteps, dtype=np.float64
        )
    elif beta_schedule == "const":
        betas = beta_end * np.ones(num_diffusion_timesteps, dtype=np.float64)
    elif beta_schedule == "jsd":  # 1/T, 1/(T-1), 1/(T-2), ..., 1
        betas = 1.0 / np.linspace(
            num_diffusion_timesteps, 1, num_diffusion_timesteps, dtype=np.float64
        )
    elif beta_schedule == "sigmoid":
        betas = np.linspace(-6, 6, num_diffusion_timesteps)
        betas = sigmoid(betas) * (beta_end - beta_start) + beta_start
    else:
        raise NotImplementedError(beta_schedule)
    assert betas.shape == (num_diffusion_timesteps,)
    return betas


class Diffusion(object):
    def __init__(self, args, config, device=None):
        self.args = args
        self.config = config
        if device is None:
            device = (
                torch.device("cuda")
                if torch.cuda.is_available()
                else torch.device("cpu")
            )
        self.device = device

        self.model_var_type = config.model.var_type
        betas = get_beta_schedule(
            beta_schedule=config.diffusion.beta_schedule,
            beta_start=config.diffusion.beta_start,
            beta_end=config.diffusion.beta_end,
            num_diffusion_timesteps=config.diffusion.num_diffusion_timesteps,
        )
        betas = self.betas = torch.from_numpy(betas).float().to(self.device)
        self.num_timesteps = betas.shape[0]

        alphas = 1.0 - betas
        alphas_cumprod = alphas.cumprod(dim=0)
        alphas_cumprod_prev = torch.cat(
            [torch.ones(1).to(device), alphas_cumprod[:-1]], dim=0
        )
        self.alphas_cumprod_prev = alphas_cumprod_prev
        posterior_variance = (
                betas * (1.0 - alphas_cumprod_prev) / (1.0 - alphas_cumprod)
        )
        if self.model_var_type == "fixedlarge":
            self.logvar = betas.log()
            # torch.cat(
            # [posterior_variance[1:2], betas[1:]], dim=0).log()
        elif self.model_var_type == "fixedsmall":
            self.logvar = posterior_variance.clamp(min=1e-20).log()

    def sample(self):
        cls_fn = None

        config_dict = vars(self.config.model)
        model = create_model(**config_dict)
        if self.config.model.use_fp16:
            model.convert_to_fp16()
        if self.config.model.in_channels == 3:
            ckpt = os.path.join(self.args.log_path, self.config.model.name)
            print('Path of the current ckpt: ' + ckpt)
            if not os.path.exists(ckpt):
                print('The model does not exist, downloading...')
                download(
                    'https://openaipublic.blob.core.windows.net/diffusion/jul-2021/%dx%d_diffusion_uncond.pt' % (
                        self.config.data.image_size, self.config.data.image_size), ckpt)
        elif self.config.model.in_channels == 1:  # todo: for the model with out_channel == 1 (won't happen for now)
            ckpt = os.path.join(self.args.log_path, self.config.model.name)
            print('Path of the current ckpt: ' + ckpt)

        else:
            raise ValueError
        model.load_state_dict(torch.load(ckpt, map_location=self.device))
        model.to(self.device)
        model.eval()
        model = torch.nn.DataParallel(model)

        self.sample_sequence(model, cls_fn)

    def sample_sequence(self, model, cls_fn=None):
        args, config = self.args, self.config
        print("data channels : " + str(self.config.data.channels))
        print("model in_channels : " + str(self.config.model.in_channels))
        print('The corresponding MATLAB path: ' + self.args.matlab_path)

        # ** define the dasLst and the gammaLst **
        if self.config.model.problem_model in ["HtH", "CHtH"]:
            # =============================== simuData (use model000000.pt (Imagenet) in DGM4MICCAI) ===============
            phantomIdies = ['1', '2', '3']
            phantomIdx = phantomIdies[2]  # select from [0 - kidney | | 1 - fetus | | 2 - simuComplex]
            gammaLst = [0.3, 0.7, 1, 1.5, 2, 2.5, 3, 3.5]  # the sequence of the additive noise level

            if self.config.model.problem_model == "HtH":
                dasLst = ['simulation' + str(gammaLst[i]) + '_Hty_' + phantomIdx for i in range(len(gammaLst))]
                gammaLst = [gammaLst[i] * 13.6 for i in range(len(gammaLst))]  # Htn has a larger std than n
            elif self.config.model.problem_model == "CHtH":
                dasLst = ['simulation' + str(gammaLst[i]) + '_CHty_' + phantomIdx for i in range(len(gammaLst))]
            else:
                raise ValueError

        elif self.config.model.problem_model in ["DRUS", "WDRUS"]:
            # =============================== picmus (use model006000.pt in DGM4MICCAI) ============================
            dasLst = ['simu_reso', 'simu_cont', 'expe_reso', 'expe_cont']
            if self.config.model.problem_model == "DRUS":
                gammaLst = [25, 50, 75, 20]  # gamma is a hyperparameter, the values here are set by grid-searching
            elif self.config.model.problem_model == "WDRUS":
                gammaLst = [4.6, 7.1, 2.2, 5.5]  # gamma is a hyperparameter, the values here are set by grid-searching
            else:
                raise ValueError
        else:
            raise ValueError

        # ** get SVD results of the model matrix **
        print(f'Loading the SVD of the degradation matrix (' + self.config.model.problem_model + ')')
        if self.config.model.problem_model in ["DRUS", "WDRUS", "HtH"]:
            from functions.svd_replacement import ultrasound0

            if self.config.model.problem_model == "DRUS":
                svdPath = self.args.matlab_path + 'SVD/02_picmus/DRUS/svd/'
                Up = torch.from_numpy(mat73.loadmat(svdPath + 'Ud.mat')['Ud'])
                lbdp = torch.from_numpy(mat73.loadmat(svdPath + 'Sigma.mat')['Sigma'])
                Vp = torch.from_numpy(mat73.loadmat(svdPath + 'Vd.mat')['Vd'])

            elif self.config.model.problem_model == "WDRUS":
                svdPath = self.args.matlab_path + 'SVD/02_picmus/WDRUS/svd/'
                Up = torch.from_numpy(mat73.loadmat(svdPath + 'Ud.mat')['Ud'])
                lbdp = torch.from_numpy(mat73.loadmat(svdPath + 'Sigma.mat')['Sigma'])
                Vp = torch.from_numpy(mat73.loadmat(svdPath + 'Vd.mat')['Vd'])

            elif self.config.model.problem_model == "HtH":
                svdPath = self.args.matlab_path + 'SVD/01_simulation/svd/'
                Up = torch.from_numpy(mat73.loadmat(svdPath + 'V.mat')['V'])
                lbdp = torch.from_numpy(mat73.loadmat(svdPath + 'LBD.mat')['LBD'])
                Vp = torch.from_numpy(mat73.loadmat(svdPath + 'Vp.mat')['Vp'])
            else:
                raise ValueError
            H_funcs = ultrasound0(config.data.channels, Up, lbdp, Vp, self.device)

        elif self.config.model.problem_model == "CHtH":
            from functions.svd_replacement import ultrasound1
            svdPath = self.args.matlab_path + 'SVD/01_simulation/svd/'
            lbdp = torch.from_numpy(mat73.loadmat(svdPath + 'lbd.mat')['lbd'])
            Vp = torch.from_numpy(mat73.loadmat(svdPath + 'Vp.mat')['Vp'])
            H_funcs = ultrasound1(config.data.channels, lbdp, Vp, self.device)

        else:
            print("ERROR: problem_model type not supported")
            quit()

        idx_so_far = 0
        print(f'Start restoration')
        for _ in range(len(dasLst)):
            dasSaveName = dasLst[idx_so_far] + '.mat'
            # ** load the observation y_0 **
            if self.config.model.problem_model in ("DRUS", "WDRUS"):
                if self.config.model.problem_model == "DRUS":
                    y_0 = torch.from_numpy(mat73.loadmat(self.args.matlab_path + 'Results/02_picmus/BH/yd/' + dasSaveName)['By'])
                elif self.config.model.problem_model == "WDRUS":
                    y_0 = torch.from_numpy(mat73.loadmat(self.args.matlab_path + 'Results/02_picmus/CBH/yd/' + dasSaveName)['CBy'])
                else:
                    raise ValueError
                y_0 = (y_0.view(1, -1)).repeat(1, config.data.channels).to(self.device)

                gamma = gammaLst[idx_so_far] * 1
                if gamma <= 0:
                    gamma = 0.1
                if self.config.model.in_channels == 3:
                    gamma = gamma * sqrt(3)

            elif self.config.model.problem_model in ("HtH", "CHtH"):
                dataPath = self.args.matlab_path + 'Results/01_simulation/SimulationResults/' + phantomIdx + '/yd/'
                if self.config.model.problem_model == "HtH":
                    y_0 = torch.from_numpy(mat73.loadmat(dataPath + dasSaveName)['o_Hty'])
                elif self.config.model.problem_model == "CHtH":
                    y_0 = torch.from_numpy(mat73.loadmat(dataPath + dasSaveName)['o_CHty'])
                else:
                    raise ValueError
                y_0 = (y_0.view(1, -1)).to(self.device)
                gamma = gammaLst[idx_so_far]

            else:
                raise ValueError

            # ===========================================================================================================

            # ** Begin DDRM **
            x = torch.randn(
                y_0.shape[0],
                config.data.channels,
                config.data.image_size,
                config.data.image_size,
                device=self.device,
            )

            # NOTE: This means that we are producing each predicted x0, not x_{t-1} at timestep t.
            with torch.no_grad():
                x, _ = self.sample_image(x, model, H_funcs, y_0, gamma, last=False, cls_fn=cls_fn)
            for i in [-1]:  # range(len(x)):
                for j in range(x[i].size(0)):
                    # ** save the DDRM restored image as .mat **
                    savemat(os.path.join(self.args.image_folder, f"{idx_so_far + j + 1}_{i}.mat"),
                            {'x': x[i][j].detach().cpu().numpy()})

            idx_so_far += y_0.shape[0]  # iterate multiple images
            print(f'Finish {idx_so_far}')

    def sample_image(self, x, model, H_funcs, y_0, gamma, last=False, cls_fn=None, classes=None):
        skip = self.num_timesteps // self.args.timesteps
        seq = range(0, self.num_timesteps, skip)
        x = efficient_generalized_steps(x, seq, model, self.betas, H_funcs, y_0, gamma, etaB=self.args.etaB,
                                        etaA=self.args.eta, etaC=self.args.eta, cls_fn=cls_fn,
                                        classes=classes)
        if last:
            x = x[0][-1]
        return x
