conda create --name capSACIN python=3.9 -y
conda activate capSACIN
conda install -c conda-forge mdanalysis -y
pip install -r requirements.txt
pip install -r requirements-dev.txt
