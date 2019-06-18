## Using iEBE-MUSIC via Docker

Docker is a software tool that allows one to deploy an application in a portable environment. A docker "image" can be created for the application, allowing any user to run a docker "container" from this image.


### 1. Build a new Docker image
We can build a docker image for the iEBE-MUSIC package using the following command,

	`docker build -t iebe-music .`
	
### 2. Run iEBE-MUSIC
The docker container has ready compiled all the software packages for iEBE-MUSIC.

`docker run -it -v ~/Desktop/coding/iEBE-MUSIC:/home/iEBE-MUSIC/ --name myiEBE iebe-music`

### 3. To delete all the Docker images in your laptop

`docker system prune -a`
