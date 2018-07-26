# Code Base Purpose
Rigid alignment of two non-overlapping scans of a common 3d object, using learning techniques on stitched surfaces.

# What each codebase component represents
- <path.h> :: file  - specifies common router logic to directorise, configurations
- <BUILD> :: dir - where build results end up.
- <BENCHMARKING> :: dir - where performance metrics are outputted. 
	- all large pipeline function timing statistics
- <FIGURES> :: dir - where all end result meshes are ouputted
- <MESHES> :: dir - where input meshes are stored. Receieved from dataset of interest.
- <NOTES> :: dir - collection of useful notes
- <PROCESSING> :: dir - where temporary meshes are outputted
- <PROGRAMS> :: dir - where project programs are stored.
- <TEST> :: dir - unsure what this accomplishes.

#TODO :: Describe how to test out stuff, as a user
Note that code coverage is limited, and tests will not encompass the entire codebase.
Note that testing was introduced only recently into the codebase, and for specific files.

#TODO :: Describe future expansions that might prove possible


#TODO :: Set up cron job to update libraries. Useful for external libs [ igl, cgal ]
