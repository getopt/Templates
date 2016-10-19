# check whether docker daemon is running
sudo status docker

# see all imanges (currently running and stopped)
sudo docker ps -a

# search docker hub
sudo docker search bowtie

# start in foreground
sudo docker run -i -t ubuntu /bin/bash

# start in dettached mode as deamon
sudo docker run --name daemon_bob -d ubuntu /bin/sh -c "while true; do echo hello world; sleep 1; done"

# attach to the container
sudo docker attach bob_the_container

# see what's inside the container 
sudo docker logs -ft daemon_bob
sudo docker top daemon_bob
sudo docker insepct daemon_bob

# run process inside the container
sudo docker exec -d daemon_bob touch /etc/new_config_file

# stopping daemonized container
sudo docker stop daemon_bob

# chech if you are logged in
sudo docker login

#  build the image from Dockerfile (in . )
sudo docker build -t="siarheimanakov/bowtie:1.1.2" .

#  build the image from Dockerfile in Git repo
sudo docker build -t="siarheimanakov/bowtie:1.1.2" git@github.com:getopt/Templates/Docker/bowtie/1.1.2/

# launching container from your own image
sudo docker run -d --name bowtie siarheimanakov/bowtie:1.1.2

# push image to Docker Hub
sudo docker push siarheimanakov/bowtie:1.1.2
