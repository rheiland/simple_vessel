# Model to test handshake between PhysiCell and SimVascular

Clone this repo (or use a release). From the root directory:

Build the model's executable:
```
cd src
make -j2
```

To test in /src:
```
cp ../data/simple_vessel.xml .
simple_vessel simple_vessel.xml

# kill (ctl-c) the sim after the first 1 or 2 time steps

# analyze, visualize the output files
cd output/
cp ../../examples/pyMCDS_cells.py .
cp ../../examples/info_cells3D.py .
cp ../../examples/fury_cells.py .

python info_cells3D.py

python fury_cells.py   # assuming you have done "pip install fury"
```

To use the Model Builder (+ Studio) GUI:

From the root directory (if you are in /src, just `cd ..`)
```
mv src/simple_vessel .    # don't just copy it, move it!
python bin/pmb.py --studio
```
Click the `Run` tab then click `Run Simulation` (might take a few secs to start).
Click the `Plot` tab and `Play` output results.
