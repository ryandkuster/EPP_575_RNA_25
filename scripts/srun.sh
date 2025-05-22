#!/usr/bin/env bash
srun --account acf-utk0011 \
     --partition=short \
     --qos=short \
     --nodes=1 \
     --cpus-per-task=10 \
     --mem=10G \
     --time=0-3:00:00 \
     --pty /bin/bash
