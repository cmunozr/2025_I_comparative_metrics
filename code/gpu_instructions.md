# GPU instructions

## 0. Pre-Execution Setup
Establish what server you will use and how many GPUs are available to work with, you will need a string and a number to use in step 1, ex: "kingma", 2 (gpus)

## 0.1 Environment Preparation & Driver Fixes
*(Note: If the lab updates the host Ubuntu version, delete and completely recreate this virtual environment).*

#### A. Initialize and upgrade virtual environment
```bash
python -m venv hmsc-venv
source hmsc-venv/bin/activate
python -m pip install --upgrade pip
```

#### B. Install Hmsc-HPC framework and requirements
```bash
python -m pip install git+https://github.com/hmsc-r/hmsc-hpc.git@main
python -m pip install tensorflow-probability==0.24.0
```

#### C. FIX: Install native CUDA capabilities for TensorFlow
```bash
python -m pip install tensorflow[and-cuda]
export TF_CPP_MIN_LOG_LEVEL=0
python -c "import tensorflow as tf; print('Available GPUs:', tf.config.list_physical_devices('GPU'))"
```

It should show something like:

```bash
Available GPUs: [PhysicalDevice(name='/physical_device:GPU:0', device_type='GPU')]
```

If not proceed to D and E otherwise go to F

#### D. FIX: If TensorFlow fails to find the GPU (empty list) or throws "Cannot dlopen" errors:
```bash
nano ~/Documents/hmsc-venv/bin/activate
```

Scroll to the bottom and paste at the end the following dynamic path injector


```nano
# -------------------------------------------------------------------
VIRTUAL_ENV_BASE_PATH="$HOME/Documents/hmsc-venv"
PY_LIB_PATH=$(dirname $(find "$VIRTUAL_ENV_BASE_PATH/lib" -maxdepth 2 -name "site-packages" | head -n 1))

# Inject nested nvidia library folders into the system linker path
if [ -d "$PY_LIB_PATH/site-packages/nvidia" ]; then
    export LD_LIBRARY_PATH="$PY_LIB_PATH/site-packages/nvidia/cuda_nvcc/lib:$PY_LIB_PATH/site-packages/nvidia/cuda_runtime/lib:$PY_LIB_PATH/site-packages/nvidia/cublas/lib:$PY_LIB_PATH/site-packages/nvidia/cudnn/lib:$PY_LIB_PATH/site-packages/nvidia/cufft/lib:$PY_LIB_PATH/site-packages/nvidia/curand/lib:$PY_LIB_PATH/site-packages/nvidia/cusolver/lib:$PY_LIB_PATH/site-packages/nvidia/cusparse/lib:$PY_LIB_PATH/site-packages/nvidia/nccl/lib:$PY_LIB_PATH/site-packages/nvidia/nvtx/lib:$LD_LIBRARY_PATH"
fi
# -------------------------------------------------------------------
```

(Press Ctrl+O to save, Enter to confirm, Ctrl+X to close nano)

#### E. Reload the environment to load the updated linker script and Verify GPU visibility
```bash
deactivate
source hmsc-venv/bin/activate
python -c "import tensorflow as tf; print('Available GPUs:', tf.config.list_physical_devices('GPU'))"
```

#### F. Verify Hmsc-HPC framework
```bash
export TF_CPP_MIN_LOG_LEVEL=0
python -c "import hmsc; print('Hmsc installed successfully.')"
pip show hmsc
```

## 1. Local Initialization

Initiate the RDS and create the GPU call shell scripts locally by running the R preparation script:
`S02_fit_models_gpu.R`

## 2. Server Data Transfer

Using FileZilla, copy the folder containing the generated model scripts and `.rds` targets to the server. Ensure they land in the target directory configured in R script.

## 3. Make them executable in linux

Move to your target directory on the server and make the execution scripts executable:
```bash
chmod +x ~/Documents/*_run_gpu_*.sh
# Example: chmod +x ~/Documents/fbs_M003_run_gpu_*.sh
```

## 4. Run 

Run the python command using something like:

```bash
nohup /home/avesta/munozcs/Documents/*_run_gpu_0.sh & nohup /home/avesta/munozcs/Documents/*_run_gpu_1.sh &
```

## 5. Wait

Wait until finishing, revise nohup_chains log. transfer to local pc. Import the chains to the Hmsc model used, using:
`S02B_import_computed_posterior_gpu.R`
