# Registrations and Installations

## Registrations
This project involves handling sensitive patient data. In addition to working on a Unix platform, all data were analysed in a secure high-computing cluster (HPC). COSMOS-SENS (CSENS), a cluster managed by LUNARC was used. To ensure access to needed resources and platforms for data analysis, the following initial steps were necesaary:
- Register a project account on LUNARC and request needed compute and storage resources (Principal Investigator).
- Follow instructions from LUNARC and create a login account for COSMOS-SENS and a 2FA password for Pocket Pass.
- Create a Static VPN account to ensure access to the cluster

## Installations
The work was carried out on:
Local computer OS: macOS on 14.6.1 and Linux Ubuntu 24.04
HPC OS: Ubuntu Rocky Linux 9.3 (Blue Onyx), rhel Centos Fedora

To manage pacakges on the MacBook Unix terminal, install Homebrew:
```bash
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```
File transfer to/from CSENS is enabled via SFTP only against the ```cs-diode.lunarc.lu.se``` server.   
The  macFUSE system often doesn't allow mac users to use sftp. So, the solution was to replace macFUSE with FUSE-T, so sftp could work.

- First uninstall the current macFUSE package by going to System Preferences -> FUSE -> Remove FUSE
- Then install FUSE-T and FUSE-T-sshfs from a terminal.

```bash
$ brew tap macos-fuse-t/homebrew-cask
$ brew install fuse-t
$ brew install fuse-t-sshfs
# retart the computer and then the terminal
```

## Connecting and Moving Files to CSENS
```bash 
sftp elugoh@cs-diode.lunarc.lu.se # requires LUNARC password and Pocket Pass OTP
sftp> put -r Elijah/ /home/elugoh/Documents/DataFolder # move all data files from local directory to cluster
sftp> put -r supplementary_data/ /home/elugoh/Documents/DataFolder 
exit # quit sftp
rm -r Elijah/ # remove all data files from local directory
```

## Loggin in to the Compute Node on CSENS 
- Turn on VPN
- Open ```https://sens12.lunarc.lu.se/``` on a web browser, and
- Login using username, password, and Pocket Pass OTP

## Software Programmes Used
- [Plink v1.9 and v2.0](https://www.cog-genomics.org/plink/)
- [PRSice v2.3.5](https://choishingwan.github.io/PRSice/)
- [GCTA v1.94.1](https://yanglab.westlake.edu.cn/software/gcta/#Download)
- [BCFTools v1.18](http://www.htslib.org/download/)
- [R/4.3.2](https://cran.r-project.org/bin/windows/base/old/4.3.2/)
- [Python 3.9.18](https://www.python.org/downloads/release/python-3918/)
- [Perl v5.32.1](https://www.perl.org/get.html)
