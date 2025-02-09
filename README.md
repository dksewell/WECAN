
# WECAN

## Overview

## Installation

Package aLSEC needs to be installed first.

``` r
install.packages("devtools")
```

    ## Error in install.packages : Updating loaded packages

``` r
# To install the aLSEC
devtools::install_github('hanhtdpham/aLSEC')
```

    ## Using GitHub PAT from the git credential store.

    ## Downloading GitHub repo hanhtdpham/aLSEC@HEAD

    ## These packages have more recent versions available.
    ## It is recommended to update all of them.
    ## Which would you like to update?
    ## 
    ##  1: All                                         
    ##  2: CRAN packages only                          
    ##  3: None                                        
    ##  4: clue         (0.3-65    -> 0.3-66   ) [CRAN]
    ##  5: slam         (0.1-50    -> 0.1-55   ) [CRAN]
    ##  6: RcppEigen    (0.3.4.0.0 -> 0.3.4.0.2) [CRAN]
    ##  7: Rcpp         (1.0.13    -> 1.0.14   ) [CRAN]
    ##  8: rlang        (1.1.4     -> 1.1.5    ) [CRAN]
    ##  9: glue         (1.7.0     -> 1.8.0    ) [CRAN]
    ## 10: skmeans      (0.2-16    -> 0.2-18   ) [CRAN]
    ## 11: RSpectra     (0.16-1    -> 0.16-2   ) [CRAN]
    ## 12: cpp11        (0.5.0     -> 0.5.1    ) [CRAN]
    ## 13: RcppArmad... (14.0.0-1  -> 14.2.3-1 ) [CRAN]
    ## 14: mvtnorm      (1.2-4     -> 1.3-3    ) [CRAN]
    ## 15: movMF        (0.2-8     -> 0.2-9    ) [CRAN]
    ## 16: igraph       (2.0.3     -> 2.1.4    ) [CRAN]
    ## 
    ## clue         (0.3-65    -> 0.3-66   ) [CRAN]
    ## slam         (0.1-50    -> 0.1-55   ) [CRAN]
    ## RcppEigen    (0.3.4.0.0 -> 0.3.4.0.2) [CRAN]
    ## Rcpp         (1.0.13    -> 1.0.14   ) [CRAN]
    ## rlang        (1.1.4     -> 1.1.5    ) [CRAN]
    ## glue         (1.7.0     -> 1.8.0    ) [CRAN]
    ## skmeans      (0.2-16    -> 0.2-18   ) [CRAN]
    ## RSpectra     (0.16-1    -> 0.16-2   ) [CRAN]
    ## cpp11        (0.5.0     -> 0.5.1    ) [CRAN]
    ## RcppArmad... (14.0.0-1  -> 14.2.3-1 ) [CRAN]
    ## mvtnorm      (1.2-4     -> 1.3-3    ) [CRAN]
    ## movMF        (0.2-8     -> 0.2-9    ) [CRAN]
    ## igraph       (2.0.3     -> 2.1.4    ) [CRAN]

    ## Installing 13 packages: clue, slam, RcppEigen, Rcpp, rlang, glue, skmeans, RSpectra, cpp11, RcppArmadillo, mvtnorm, movMF, igraph

    ## 将程序包安装入'C:/Users/tamia/AppData/Local/R/win-library/4.4'
    ## (因为'lib'没有被指定)

    ## 程序包'clue'打开成功，MD5和检查也通过
    ## 程序包'slam'打开成功，MD5和检查也通过
    ## 程序包'RcppEigen'打开成功，MD5和检查也通过
    ## 程序包'Rcpp'打开成功，MD5和检查也通过

    ## Warning: 无法删除软件包 'Rcpp' 的先前安装

    ## Warning in file.copy(savedcopy, lib, recursive = TRUE): problem copying
    ## C:\Users\tamia\AppData\Local\R\win-library\4.4\00LOCK\Rcpp\libs\x64\Rcpp.dll to
    ## C:\Users\tamia\AppData\Local\R\win-library\4.4\Rcpp\libs\x64\Rcpp.dll: Permission denied

    ## Warning: 回复了'Rcpp'

    ## 程序包'rlang'打开成功，MD5和检查也通过

    ## Warning: 无法删除软件包 'rlang' 的先前安装

    ## Warning in file.copy(savedcopy, lib, recursive = TRUE): problem copying
    ## C:\Users\tamia\AppData\Local\R\win-library\4.4\00LOCK\rlang\libs\x64\rlang.dll to
    ## C:\Users\tamia\AppData\Local\R\win-library\4.4\rlang\libs\x64\rlang.dll: Permission denied

    ## Warning: 回复了'rlang'

    ## 程序包'glue'打开成功，MD5和检查也通过
    ## 程序包'skmeans'打开成功，MD5和检查也通过
    ## 程序包'RSpectra'打开成功，MD5和检查也通过
    ## 程序包'cpp11'打开成功，MD5和检查也通过
    ## 程序包'RcppArmadillo'打开成功，MD5和检查也通过
    ## 程序包'mvtnorm'打开成功，MD5和检查也通过
    ## 程序包'movMF'打开成功，MD5和检查也通过
    ## 程序包'igraph'打开成功，MD5和检查也通过
    ## 
    ## 下载的二进制程序包在
    ##  C:\Users\tamia\AppData\Local\Temp\RtmpEPaAn8\downloaded_packages里
    ## ── R CMD build ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
    ##          checking for file 'C:\Users\tamia\AppData\Local\Temp\RtmpEPaAn8\remotes5a8c2e017cd2\hanhtdpham-aLSEC-d13d954/DESCRIPTION' ... OK  ✔  checking for file 'C:\Users\tamia\AppData\Local\Temp\RtmpEPaAn8\remotes5a8c2e017cd2\hanhtdpham-aLSEC-d13d954/DESCRIPTION' (389ms)
    ##       ─  preparing 'aLSEC': (1.2s)
    ##    checking DESCRIPTION meta-information ...  ✔  checking DESCRIPTION meta-information
    ## ─  cleaning src
    ##       ─  checking for LF line-endings in source and make files and shell scripts (505ms)
    ##       ─  checking for empty or unneeded directories
    ##   ─  building 'aLSEC_0.1.0.tar.gz'
    ##      
    ## 

    ## 将程序包安装入'C:/Users/tamia/AppData/Local/R/win-library/4.4'
    ## (因为'lib'没有被指定)

``` r
# To install the WECAN
devtools::install_github("HaominLi7/WECAN")
```

    ## Using GitHub PAT from the git credential store.

    ## Downloading GitHub repo HaominLi7/WECAN@HEAD

    ## These packages have more recent versions available.
    ## It is recommended to update all of them.
    ## Which would you like to update?
    ## 
    ## 1: All                             
    ## 2: CRAN packages only              
    ## 3: None                            
    ## 4: rlang  (1.1.4  -> 1.1.5 ) [CRAN]
    ## 5: withr  (3.0.1  -> 3.0.2 ) [CRAN]
    ## 6: pillar (1.9.0  -> 1.10.1) [CRAN]
    ## 7: Rcpp   (1.0.13 -> 1.0.14) [CRAN]
    ## 
    ## ── R CMD build ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
    ##          checking for file 'C:\Users\tamia\AppData\Local\Temp\RtmpEPaAn8\remotes5a8c3f4e4cc8\HaominLi7-WECAN-70e76cf/DESCRIPTION' ...  ✔  checking for file 'C:\Users\tamia\AppData\Local\Temp\RtmpEPaAn8\remotes5a8c3f4e4cc8\HaominLi7-WECAN-70e76cf/DESCRIPTION' (389ms)
    ##       ─  preparing 'WECAN': (708ms)
    ##      checking DESCRIPTION meta-information ...     checking DESCRIPTION meta-information ...   ✔  checking DESCRIPTION meta-information
    ## ─  cleaning src
    ##       ─  checking for LF line-endings in source and make files and shell scripts (365ms)
    ##       ─  checking for empty or unneeded directories
    ##       ─  building 'WECAN_0.0.0.9000.tar.gz'
    ##      
    ## 

    ## 将程序包安装入'C:/Users/tamia/AppData/Local/R/win-library/4.4'
    ## (因为'lib'没有被指定)

    ## Warning in i.p(...):
    ## 安装程序包'C:/Users/tamia/AppData/Local/Temp/RtmpEPaAn8/file5a8c530f7489/WECAN_0.0.0.9000.tar.gz'时退出狀態的值不是0

## Example of running the algorithm
