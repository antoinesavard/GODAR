# GODAR
## Granular Discrete Arctic Rheology: DEM for sea ice modeling
---
DEM model. To compile, run scons-3 with the appropriate parameters.

scons-3 -j [n]        ! compile the code on [n] cores
scons-3 -c            ! clear the build
scons-3 debug=1       ! debug the code
./godar < input_file  ! run the executable

The code executes automatically a python script to output a video file of the collisions. This is saved in collisionXY.mp4 file. To suppress, just comment the line 

```fortran
call execute_command_line("/aos/home/asavard/anaconda3/bin/python /storage/asavard/DEM/plots/video.py")
```

To use it, change the path to your python interpreter and the location of the python code.
