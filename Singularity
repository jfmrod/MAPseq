Bootstrap: docker
From: ubuntu:20.04
IncludeCmd: yes

%environment
#R_VERSION=4.1
#export R_VERSION
#R_CONFIG_DIR=/etc/R/
#export R_CONFIG_DIR
export LC_ALL=C
export PATH=$PATH:/opt/software/miniconda3/bin

%post
  apt-get update
  apt-get install -y apt-transport-https apt-utils software-properties-common
  apt-get install -y add-apt-key
  export DEBIAN_FRONTEND=noninteractive
  ln -fs /usr/share/zoneinfo/America/New_York /etc/localtime
  apt-get install -y tzdata
  dpkg-reconfigure --frontend noninteractive tzdata

  apt-get install -y wget

  mkdir -p /opt/software

  # MAPseq installation
  cd /opt/software
  wget -q https://github.com/jfmrod/MAPseq/releases/download/2.0.1alpha/mapseq-2.0.1alpha-linux.tar.gz
  tar xvzf mapseq-2.0.1alpha-linux.tar.gz
  rm mapseq-2.0.1alpha-linux.tar.gz
  mv mapseq-2.0.1alpha-linux mapseq

  ln -s /opt/software/mapseq/mapseq /usr/bin/
  ln -s /opt/softwaremapseq/share /usr/bin/
