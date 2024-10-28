# Registration and Installations
This is the start of this project...

## Registrations
The nature of this project involves handling sensitive patient data. In addition to working on a Unix platform, all data must be analysed in a secure high-computing cluster (HPC). COSMOS-SENS, a cluster managed by LUNARC was used. To ensure access to needed resources and platforms for data analysis, the following initial steps were necesaary:
- Register on LUNARC to gain access to the project.
- Follow instructions from LUNARC and create login password for pocket pass and login username for COSMOS-SENS.
- create a Static VPN account to ensure access to the cluster

## Installations
Local computer OS: macOS on 14.6.1
HPC OS: Ubuntu Rocky Linux 9.3 (Blue Onyx), rhel Centos Fedora

To manage pacakges on the MacBook Unix terminal, install Homebrew:
```bash
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```
File transfer to/from CSENS is enabled via SFTP only against the ```cs-diode.lunarc.lu.se``` server.   
The  macFUSE system often doesn't allow mac users to use sftp. So, the solution was to replace macFUSE with FUSE-T, so sftp could work.

- First uninstall the current package going to System Preferences -> FUSE -> Remove FUSE
- Then install FUSE-T and FUSE-T-sshfs from a terminal.

```bash
$ brew tap macos-fuse-t/homebrew-cask
$ brew install fuse-t
$ brew install fuse-t-sshfs
# retart the computer and then the terminal
```

## Connecting and Moving Files to the Cluster
```bash 
sftp elugoh@cs-diode.lunarc.lu.se # requires password and Pocket Pass OTP
get -r Elijah/ /home/elugoh/Documents/DataFolder # move all data files from local directory to cluster
exit # quit sftp
rm -r Elijah/ # remove all data files from local directory
```

## Loggin in to the Computer Cluster 
- Turn on VPN
- Open ```https://sens12.lunarc.lu.se/``` on a web browser, and
- Login using username, password, and Pocket Pass OTP

## Software Programmes
[Plink 1.9](https://www.cog-genomics.org/plink/)
[GCTA 1.94.1](https://yanglab.westlake.edu.cn/software/gcta/#Download)
[BCFTolls 1.18](http://www.htslib.org/download/)
[R 4.3.2](https://cran.r-project.org/bin/windows/base/old/4.3.2/)
[Python 3.9.18](https://www.python.org/downloads/release/python-3918/)
