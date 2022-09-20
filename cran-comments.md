## Resubmission
This is a resubmission. In this version I have:

* Modified the DESCRIPTION to not start with title or "This package".
* Added a reference to the paper for the package in DESCRIPTION.
* Added \value tags to geny.Rd, MixtureModel.Rd, qmultinom.Rd.
* Removed print()/cat() from mixture.R and sim.R, except for in summary().
* Updated all examples to use no more than 2 cores.

## Test Environments
- ubuntu-gcc-release
- macos-highsierra-release
- windows-x86_64-devel
- windows-x86_64-release
- fedora-clang-devel

## R CMD check results

There were no ERRORs or WARNINGs.

There are two NOTES on Windows Server 2022, R-devel, 64 bit:

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Phuc H. Nguyen <phuc.nguyen.rcran@gmail.com>'

New submission

Uses the superseded packages: 'doSNOW', 'snow'
```

I'm using 'doSNOW' and 'snow' because they support the display of progress bar inside 'foreach' loop, while 'doParallel' do not.

```
* checking for detritus in the temp directory ... NOTE
Found the following files/directories:
  'lastMiKTeXException'
```
