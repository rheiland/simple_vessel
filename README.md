# Model to test handshake between PhysiCell and SimVascular

Clone this repo (or use a release). From the root directory:

```
cd src
make -j2
cp ../data/multiple_vessels.xml .
simple_vessel multiple_vessels.xml

# kill (ctl-c) the sim after the first 1 or 2 time steps

# analyze, visualize the output files
cd output/
cp ../../examples/info_cells3D.py .
cp ../../examples/fury_cells.py .
cp ../../examples/pyMCDS_cells.py .

python info_cells3D.py

python fury_cells.py   # assuming you have done "pip install fury"
```