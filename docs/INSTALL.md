Installation
============

Quickie Installation
--------------------

1. Run `make install` while you are in the Lyve-SET directory. This step probably takes 10-20 minutes.
2. Update the path to include the scripts subdirectory. You can do this yourself if you are comfortable or run `make env`.
3. Update your current session's path: `source ~/.bashrc`

Requirements
------------
* **Perl, multithreaded**
* **BLAST+**
* **GIT**, **SVN** (for installation and updating)

Installation
------------

    make

It doesn't get easier than this.

###Testing your installation

    set_test.pl lambda lambda

Upgrading
---------
### By stable releases
Unfortunately the best way to get the next stable release is to download the full version like usual, followed by `make`.  If successful, then delete the directory containing the older version.

    cd ~/tmp
    wget http://latest/release.tar.gz
    tar zxvf release.tar.gz
    cd Lyve-SET
    make # takes 10-20 minutes to download packages on broadband; install
    cd ~/bin
    rm -r Lyve-SET && mv ~/tmp/Lyve-SET .

### By `git`
    git pull -u origin master
    make
