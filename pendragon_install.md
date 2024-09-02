```bash
wget https://sourceforge.net/projects/proper-library/files/proper_v3.0b_python_6mar18.tar.gz/download
mv download proper_v3.0b_python_6mar18.tar.gz
tar -zxf proper_v3.0b_python_6mar18.tar.gz 
cd proper_v3.0b_python_6mar18
python setup.py install --user
pip install --user pyfftw
cd 
python -c 'import proper; proper.prop_use_fftw()'
python -c 'import proper; proper.prop_fftw_wisdom(2048)'
python -c 'import proper; proper.prop_fftw_wisdom(1024)'
python -c 'import proper; proper.prop_fftw_wisdom(512)'
```

```bash
mkdir toliman
git clone git@github.com:beldaz/toliman.git toliman
```