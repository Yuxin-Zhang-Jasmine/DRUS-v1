# Ultrasound Image Reconstruction with Denoising Diffusion Restoration Models 

The supplement materials for this repo is on [GoogleDrive](https://drive.google.com/drive/folders/14snEBvvgLAAqhFtH6gRbbpqty_rvY1Ra?usp=sharing), 
including 2 folders MATLAB and exp. It is neccessary to download them for reproducing the results.

After downloading the supplement materials, the folder's structure is supposed to be:
```bash
<MATLAB> 
├── PICMUS  # (download from drive) contains some tools required for running other scripts 
│   ├── ... 
├── Results # work in this directory for pre(post)-processing or showing the results in paper
│   ├── 01_simulation  
│   ├── 02_spicmus
│   ├── ...  # some functions required for running other scripts 
├── SVD  # (download from drive) contains the SVD results of model matrices 
│   ├── 01_simulation  
│   ├── 02_spicmus
<exp>  # (download from drive)
├── logs 
│   ├── imagenet  # a folder named by the argument `--doc` given to main.py
│   │   ├──model000000.pt  # checkpoint (2.2GB)  before fine-tuning
│   │   ├──model006000.pt  # checkpoint (2.2GB)  after fine-tuning
├── image_samples 
│   ├── us # contains the restored data
<runners> 
├── ... 
<gided_diffusion> 
├── ... 
<functions> 
├── ... 
<configs>
├── ... 
main.py
environment.yml
```

### For simply displaying the figures and the table values in the paper
Only `MATLAB > PICMUS` is required to be downloaded from drive, and then you can simply run the scripts
- `MATLAB > Results > 01_simulation > ResultShower_simulation_images.m` for Fig.2 (a-b)
- `MATLAB > Results > 01_simulation > ResultShower_simulation_metrics_6figs.m` for Fig.2 (c)
- `MATLAB > Results > 02_picmus > ResultShower_picmus_images.m` for Fig.3
- `MATLAB > Results > 02_picmus > ResultShower_picmus_metrics.m` for Tab.1


### Inputs of DDRM
The inputs used in the paper are saved in the folders `MATLAB/Results/.../yd`. These inputs were calculated by using the scripts: 
- `MATLAB > Results > 01_simulation > compute_Hty_CHty.m` (for synthetic datasets with both DRUS and WDRUS models)
- `MATLAB > Results > 02_picmus > BH > compute_By.m` (for picmus datasets with DRUS model)
- `MATLAB > Results > 02_picmus > CBH > compute_CBy.m` (for picmus datasets with WDRUS model)

### Sampling from the model
The general command to do the restoration is as follows:
```
python main.py --ni --config {CONFIG.yml} --doc {MODELFOLDER} --timesteps {STEPS} --matlab_path {MATLABPATH} 
```

where the following are options
- `CONFIG` is the name of the config file (see `configs/`), including hyperparameters such as batch size and network architectures.
- `MODELFOLDER` is the name of the folder saving the diffusion model checkpoints 
- `STEPS` controls how many timesteps (in [1,1000]) used in the process. (e.g. 50)
- `MATLABPATH` is the path of the folder `<MATLAB>`.  

For example
```
python main.py --ni --config imagenet_256.yml --doc imagenet --timesteps 50 --matlab_path /home/user/Documents/MATLAB/ 
```

## References and Acknowledgements
```
@inproceedings{kawar2022denoising,
    title={Denoising Diffusion Restoration Models},
    author={Bahjat Kawar and Michael Elad and Stefano Ermon and Jiaming Song},
    booktitle={Advances in Neural Information Processing Systems},
    year={2022}
}
```

This implementation is based on / inspired by:
- [https://ddrm-ml.github.io/](https://ddrm-ml.github.io/) (the DDRM repo)
