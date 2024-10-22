## RootFreqs

[BEAST 2](http://beast2.org) package for specifying equilibrium (root) frequencies 
that differ from frequencies provided by the site model.
A sequence can be specified to define a different distribution at each site. 
This can be useful for example when a root sequence is known.

The package should work with the following features:

* gamma rate heterogeneity
* proportion invariable sites
* BEAGLE
* ascertainment correction
* ambiguous sites at tips and/or at root
* uncertain sites at tips and/or at root
* threaded likelihood calculation
* with/without an error model at the tips


## Installation

## Install through BEAUti

# Installing the package

Rootfreqs is a [BEAST2](http://beast2.org) package that requires BEAST 2 v2.7.
If you have not already done so, you can get BEAST 2 from [here](http://beast2.org).

To install rootfreqs, it is easiest to start BEAUti (a program that is part of BEAST), and select the menu `File/Manage packages`. A package manager dialog pops up, that looks something like this:

![Package Manager](https://github.com/rbouckaert/rootfreqs/blob/main/doc/package_repos.png?raw=true)


If the rootfreqs package is listed, just click on it to select it, and hit the `Install/Upgrade` button.

If the rootfreqs package is not listed, you may need to add a package repository by clicking the `Package repositories` button. A window pops up where you can click `Add URL` and add `https://raw.githubusercontent.com/CompEvol/CBAN/master/packages-extra-2.7.xml` in the entry. After clicking OK, the dialog should look something like this:

![Package Repositories](https://github.com/CompEvol/CCD/raw/master/doc/package_repos0.png)

Click OK and now rootfreqs should be listed in the package manager (as in the first dialog above). Select and click Install/Upgrade to install.



## Install by hand

* Download the package from [here](https://github.com/rbouckaert/rootfreqs/releases/download/v0.0.2/rootfreqs.package.v0.0.2.zip)
* Create rootfreqs directory inside BEAST package directory
  * for Windows in Users\<YourName>\BEAST\2.X\rootfreqs
  * for Mac in /Users/<YourName>\/Library/Application Support/BEAST/2.X/rootfreqs
  * for Linux /home/<YourName>/.beast/2.X/rootfreqs
  Here <YourName> is the username you use, and in “2.X” the X refers to the major version of BEAST, so 2.X=2.7 for version 2.7.6.
* Unzip the file `rootfreqs.package.v0.0.2.zip` inside the rootfreqs directory

## Build from code

* Get code for beast2, BeastFX and CCD repositories:
  * git clone https://github.com/CompEvol/beast2.git
  * git clone https://github.com/CompEvol/BeastFX.git
  * git clone https://github.com/rbouckaert/rootfreqs.git
* Run `ant install` from the rootfreqs directory
  
## Usage

There is no BEAUti support yet, so you have to edit the XML by hand: see [ example](https://github.com/rbouckaert/rootfreqs/blob/main/examples/testRootFreqSequence.xml) file for details.
