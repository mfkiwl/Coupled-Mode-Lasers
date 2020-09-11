# Overlap Factor Model
Matlab routines for laser dynamics.  

## INSTALLATION NOTES:
For these routines to work properly, both the 'Coupled Mode Twin' directory
and the 'Numeric' directory should be in a directory called 
'Overlap-Factor-Model' on the default MATLAB path. This is usually something
like:

```
    C:\Users\<user_name>\Documents\MATLAB
```

Type 'userpath' at the MATLAB command line to check.

## USER GUIDE:
The 'Coupled Mode Twin' should have a sub-directory called 'html'. From 
within MATLAB, double-click on the 'CoupledMode.html' file. This opens a 
MATLAB browser with a help file for the routines.

You may also type: 

```
    help <matlab_filename>
``` 

at the command-line if you are in the 'Coupled Mode Twin' directory.

To run the help file examples, type:

```
    CoupledMode
```

at the command-line.

## IN DEVELOPMENT
### Additional help files
Help files have been written for the following routines:

```
    coupled1D.m
    coupled1DS.m
    loadParams.m
    runCoupled1D.m
    singleSlab.m
```

The script 

```
    CoupledMode.m
```

is essentially a help file and may be run directly.

The files above may be used for time simulations. Other files in the 
directory may be used for finding steady state solutions. These have 
some commented documentation but are not yet written up as proper help 
files.

### Stability mapping
A stability mapping GUI app is currently being developed. The installation 
file may be found in the ```Stability Map GUI``` folder

```
    StabilityMap.mlappinstall
```

Double clicking on this installs the app.


