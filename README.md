# Caly

ALEX probes calibration program.


## Usage

This program needs to be in the same folder than the shots information, directly copied. It assumes that you only take RAW information of 4 channels:
1. The Rogowsky probe channel,
2. 4Resistive divider channel,and
3. RC integrator from the Rogowsky probe.

The channel order is the listed here, but it can be changed without problem.

Folder "peakutlis", here stored as a ZIP file, was requested by the first version of the progrma, but now is not ncessary at all.
To use it, just execute the program with the command "python3 caly.py" and you will obtain a text file with the calibration of the Rogowsky and RC integrator signals with some information on the ALEX circuit. Namely, its period and inductance. It assumes that you know already the calibration of the voltages probe 4resistivie divider used to calculate the voltage.
