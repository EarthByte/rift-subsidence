#### The easiest way to run the Jupyter notebooks in this folder is to use a docker container.

```docker run -p 8889:8888 -it --rm -v`pwd`:/workspace earthbyte/rift_subsidence_docker```

Open http://localhost:8889 in a web browser

#### build "rift_subsidence_docker" docker container

* Step 1: download folder https://github.sydney.edu.au/EarthByte/TeachingToolkit/tree/master/GEOS3104_3804/docker
* Step 2: go inside the "docker" folder and run `docker build -f Dockerfile -t my_rift_subsidence_docker .`

The RiftSubsidence.py is converted from Fortran program gforisos.f. Read the comments inside RiftSubsidence.py about the authors and change history.

#### The 4 columns in the RiftSubsidence output file are:
* Time (my) after initial rifting
* tectonic water-loaded basement subsidence (km)
* heat flow (mW/m2) 
* strain rate (per billion year)

