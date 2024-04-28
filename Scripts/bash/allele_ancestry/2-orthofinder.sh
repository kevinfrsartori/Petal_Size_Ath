#!/bin/bash

#SBATCH -A naiss2023-22-1378
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH -J 447_OF
#SBATCH --mail-user kevin.sartori@slu.se
#SBATCH --mail-type=ALL

ml bioinfo-tools OrthoFinder/2.5.2

orthofinder.py -f /crex/proj/snic2020-16-182/proteomes/A.tha447.hal.lyr.013024/ -y