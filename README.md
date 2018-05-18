# THESIS PROJECT

#TODO :: Describe what this code base is for
Rigid alignment of two non-overlapping scans of a common 3d object, using learning techniques on stitching surfaces.

#TODO :: Describe what each component represents
- <path.h> :: file serves as common router logic to different directorise and config areas of the project.
- <FIGURES> :: dir - where all end result meshes are ouputted
- <BENCHMARK> :: dir - where performance metrics are outputted. Includes
	- timing for all large pipeline functions
- <PROCESSING> :: dir - where temporary meshes are outputted
	- example-1 :: meshes across pipeline stages

#TODO :: Describe how to test out stuff, as a user
Note that code coverage is limited, and tests will not encompass the entire codebase.
Note that testing was introduced only recently into the codebase, and for specific files.

#TODO :: Describe future expansions that might prove possible
#TODO :: rename folder from <shared> to <meshes>
#TODO :: Set up an automated script to update libraries, based on current time. This is useful for external libs [ e.g. igl, cgal ]


