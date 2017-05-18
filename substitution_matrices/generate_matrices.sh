#!/bin/bash


# VTML240 60 (MUSCLE)
perl vt_scores.pl VTML 240 vtml240 60
python3 serialize_matrices.py vtml240
# VTML 200 3
perl vt_scores.pl VTML 200 vtml200
python3 serialize_matrices.py vtml200
