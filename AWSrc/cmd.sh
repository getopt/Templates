
sudo apt-get install git
mkdir Git
cd Git/
git clone https://github.com/getopt/Templates.git
cd Templates/HomeSettings/
cp vimrc           ~/.vimrc
cp profile         ~/.profile         
cp inputrc         ~/.inputrc
cp bashrc          ~/.bashrc
cp bash_profile    ~/.bash_profile
cp bash_path       ~/.bash_path
cp bash_generic    ~/.bash_generic
cp bash_apparix    ~/.bash_apparix


cd
sudo sh -c "echo deb https://apt.dockerproject.org/repo ubuntu-trusty main > /etc/apt/sources.list.d/docker.list"
sudo apt-key adv --keyserver hkp://p80.pool.sks-keyservers.net:80 --recv-keys 58118E89F3A912897C070ADBF76221572C52609D
sudo apt-get update
sudo apt-get install docker-engine
