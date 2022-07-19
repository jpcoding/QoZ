QoZ: Dynamic Quality Metric Oriented Error Bounded Lossy Compression for Scientific Datasets

## Introduction

This is the source code of QoZ: Dynamic Quality Metric Oriented Error Bounded Lossy Compression for Scientific Datasets

## Dependencies

Please Installing the following dependencies before running the artiact evaluation experiments:

* Python >= 3.6
* numpy >= 1.21.2
* pandas >= 1.4.1
* qcat (from https://github.com/Meso272/qcat, check its readme for installation guides. Make sure the following executables are successfully installed: calculateSSIM and computeErrAutoCorrelation)

## 3rd party libraries/tools

* Zstd >= 1.3.5 (https://facebook.github.io/zstd/). Not mandatory to be mannually installed as Zstandard v1.4.5 is included and will be used if libzstd can not be found by
  pkg-config.

## Installation

* mkdir build && cd build
* cmake -DCMAKE_INSTALL_PREFIX:PATH=[INSTALL_DIR] ..
* make
* make install

Then, you'll find all the executables in [INSTALL_DIR]/bin and header files in [INSTALL_DIR]/include. A Cmake version >= 3.13.0 is needed and we recommend to use gcc version 9.x to compile the code. 

## Single compression/decompression testing Examples

You can use the executable 'qoz' command to do the compression/decompression. Just run "qoz" command to check the instructions for its arguments.
Currently you need to add a configuration file to the argument line (-c) for activating the new features of QoZ. 
The corresponding cofiguration files for each test dataset can be generated by generate_config.py (details shown on following).

## Evaluation guides

Step 1: Download the dataset from the following links,then unzip them:

* CESM-ATM: https://g-8d6b0.fd635.8443.data.globus.org/ds131.2/Data-Reduction-Repo/raw-data/CESM-ATM/SDRBENCH-CESM-ATM-cleared-1800x3600.tar.gz 
* Miranda: https://g-8d6b0.fd635.8443.data.globus.org/ds131.2/Data-Reduction-Repo/raw-data/Miranda/SDRBENCH-Miranda-256x384x384.tar.gz
* Hurricane-ISABEL: https://g-8d6b0.fd635.8443.data.globus.org/ds131.2/Data-Reduction-Repo/raw-data/Hurricane-ISABEL/SDRBENCH-Hurricane-ISABEL-100x500x500_log.tar.gz
* NYX: https://g-8d6b0.fd635.8443.data.globus.org/ds131.2/Data-Reduction-Repo/raw-data/EXASKY/NYX/SDRBENCH-EXASKY-NYX-512x512x512_log.tar.gz
* SCALE-LETKF: https://g-8d6b0.fd635.8443.data.globus.org/ds131.2/Data-Reduction-Repo/raw-data/SCALE_LETKF/SDRBENCH-SCALE_98x1200x1200_log.tar.gz
* Data download command: wget {DataLink} --no-check-certificate (DataLink: The link of data)
* Data unzip command: tar -zxvf {DataFile} (DataFile: the .tar.gz file of downloaded data)

Step 2: Preprocess the downloaded Miranda data with the preprocess_data.py:

* Usage: python preprocess_data -m {MirandaPath} (MirandaPath is the folder of the Miranda Dataset)

Step 3: Run the generate_config.py to generate the configuration files for QoZ:

* python generate_config.py

After that, the configuration files for QoZ test will generated in configs/ folder in the name format of {DatasetName}\_{Target}.config. 

* DatasetName: The name of dataset (cesm, miranda, hurricane, nyx, scale)
* Target: The optimization target, includes cr (compression ratio), psnr (PSNR), ssim (SSIM), ac (Autocorrelation). For each mode of the QoZ, use the configuration file of the corresponding target.

Step 4: Run test_qoz.py to generate the test results.

* Command: test_qoz.py -i {Datapath} -o {OutputPath} -d {DatasetName} -t {Target}
* Datapath: the folder path of the dataset
* OutputPath: the output data file prefix. The output files will be in format of Outputpath_{Metric}.tsv
* DatasetName: See step 3
* Target: See step 3

The output will contain:
* overall_cr: the overall compression ratio under different vr-rel error bounds.
* overall_psnr: the overall compression PSNR under different vr-rel error bounds.
* cspeed: the compression speed under different vr-rel error bounds (the speed results in the paper are generated with psnr mode).
* dspeed: the decompression speed under different vr-rel error bounds (the speed results in the paper are generated with psnr mode).
* overall_ssim: the overall compression SSIM under different vr-rel error bounds (ssim target only).
* overall_ac: the overall compression Autocorrelation under different vr-rel error bounds (ac target only).

Output examples:

* For fast execution only part of the data points (error bounds) are covered in the output.
* The example outputs are in the results folder. Notice that the results are different from the results in the paper as different machine and nodes are used.

Tips for plotting and comparing the results genereted with the results in the paper:

* Bit rate = 32 / compression ratio.
* The compression ratios generated using mode cr correspond to the table 3 in the paper.
* The cr and PSNR generated using mode psnr correspond to the figure 8 in the paper.
* The cr and SSIM generated using mode ssim correspond to the figure 9 in the paper.
* The cr and AC generated using mode ac correspond to the figure 10 in the paper.
* The compression/decompression speeds generated using mode psnr correspond to the table 4 in the paper (different results may achieved if not run on the same nodes listed in the paper).


