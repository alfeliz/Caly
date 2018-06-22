# Caly

ALEX probes calibration program.


## Usage

This program needs to be in the same folder than the shots information, directly copied. It assumes that you only take RAW information of 4 channels:
1. The Rogowsky probe channel,
2. 2Resistive divider channel,
3. 3Resistive divider channel, and 
4. RC integrator from the Rogowsky probe.

The channel order is the listed here, but it can be changed without problem.

It requests the folder "peakutlis" here stored as a ZIP file. 
To use it, just execute the program with the command "python caly.py" and you will obtain a text file with the calibration of both, the Rogowsky and RC integrator signals. It assumes that you know already the calibration of hte voltages probes 2resistivie divider and 3Resistive divider.
