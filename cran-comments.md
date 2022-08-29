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
