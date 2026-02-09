# F1F209
Fit model from Peter Bosted.


## Installation
```bash
make clean && make -j<no-of-cpu-cores>
```
It will create an executable `./Fit`

## Usage
```bash
./Fit --input <params.txt> --output <desired-outputfilename.txt> [--p beam_enegy=<value> --p outgoing_lepton_angle=<value>]
```
The structure of the params.txt is as follows:

```bash
number_of_events      = 200
number_of_prtons      = 6
number_of_neutrons    = 6
beam_enegy            = 0.68  // In GeV
outgoing_lepton_angle = 36    // In degrees
```

One can run in batch mode for a specific nucleus using the helper script `BatchMode.sh`
### Usage
```bash
./BatchMode.sh <config.txt> [--verbose]
```

The `config.txt` flie contains different kinematics for which the fit should be run. The contents of the file looks as follows:

```
# E (GeV) theta (deg)
0.1597  90
0.1597  145
...
```




