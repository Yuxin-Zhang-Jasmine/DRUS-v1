# Ultrasound Image Reconstruction with Denoising Diffusion Restoration Models 

The code has been tested on PyTorch 1.10. 
The supplement materials for this repo is on [GoogleDrive](https://drive.google.com/drive/folders/14snEBvvgLAAqhFtH6gRbbpqty_rvY1Ra?usp=share_link), 
including the SVD results of DRUS (BH, HtH) and WDRUS (CBH, CHtH), and the picmus database. It is neccessary to download them for reproducing the results.
### Inputs of DDRM
The ultrasound RF channel data `y` is pre-processed in the scripts: 
- `compute_By.m` (for picmus datasets with DRUS model)
- `compute_CBy.m` (for picmus datasets with WDRUS model)
- `compute_Hty_CHty.m` (for simulated datasets with both DRUS and WDRUS models)

And the generated files which are the inputs of DDRM would be saved in folders `yd`. 

### Pretrained models
- For the 3-channel (RGB) input case, we use pretrained models from [https://github.com/openai/guided-diffusion](https://github.com/openai/guided-diffusion)
- For the 1-channel (Gray) input case, the pretrained models can be downloaded from [GoogleDrive](https://drive.google.com/drive/folders/1ANNCX53r7LHC76tIo-y_g4QmJEfOOUV5?usp=share_link)

The models and datasets are placed in the `exp/` folder as follows:
```bash
<exp> # a folder named by the argument `--exp` given to main.py
├── datasets # all dataset files 
├── logs # contains checkpoints produced during training
│   ├── imagenet # ImageNet checkpoint files (2.2GB for each)
│   │   ├── Gray
│   │   │   ├──256x256_diffusion_uncond.pt
│   │   ├── RGB
│   │   │   ├──256x256_diffusion_uncond.pt
├── image_samples # contains generated samples
```

### Sampling from the model
It is necessary to modify the paths for SVD results (line120-130) and for the inputs `yd` (line159-161) in the script `runners > diffusion.py` for executing the DDRM process.
And the general command to sample from the model is as follows:
```
python main.py --ni --config {CONFIG}.yml --doc {DATASET} --timesteps {STEPS} --eta {ETA} --etaB {ETA_B} --deg {DEGRADATION} --sigma_0 {SIGMA_0} -i {IMAGE_FOLDER}
```
where the following are options
- `ETA` is the eta hyperparameter in the paper. (default: `0.85`)
- `ETA_B` is the eta_b hyperparameter in the paper. (default: `1`)
- `STEPS` controls how many timesteps used in the process.
- `DEGREDATION` is the type of degredation allowed. (default: `us0`)
- `SIGMA_0` is the noise observed in y.
- `CONFIG` is the name of the config file (see `configs/` for a list), including hyperparameters such as batch size and network architectures.
- `DATASET` is the name of the dataset used, to determine where the checkpoint file is found.
- `IMAGE_FOLDER` is the name of the folder the resulting images will be placed in (default: `images`)

The restored data is placed in the `<exp>/image_samples/{IMAGE_FOLDER}` folder


### Example

```
python main.py --ni --config imagenet_256.yml --doc imagenet --timesteps 100 --eta 0.85 --etaB 1 --deg us0 --sigma_0 0.5 -i us
```
and the restored data would be saved in a folder `us` under the directory `image_samples`. you can then move these results to 
`MATLAB > workInthisDirectory > picmus > BH/CBH > Results > [foldername]`

You can also run the scripts 
- `ResultShower_simulation_images.m`
- `ResultShower_simulation_metrics.m`
- `ResultShower_picmus_images.m`
- `script_evaluation_scores.m`

to have a look of the saved results, and all of the codes in MATLAB is preferred to be run by `Run Section` under the folder `workInthisDirectory` to avoid the `directory/file not found` error.
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
