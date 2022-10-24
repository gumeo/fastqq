## Test environments
* rhub with devtools::check_rhub
* local OS X install, R 4.1.2
* win builder with devtools::check_win_devel

## R CMD check results
There were no ERRORs or WARNINGs. But I am making an update because of an
email I received from Kurt Hornik about a warning on debian and fedora with
clang 

(see here: https://cran.r-project.org/web/checks/check_results_fastqq.html).

I have replaced the bitwise or (which was ok here, so everything should
work the same).

I have also replaces old http links with https in the readme, (comment from
last submission).
