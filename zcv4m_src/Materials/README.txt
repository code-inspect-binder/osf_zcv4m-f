README.txt 26 Sep 2017
Jeffrey R. Stevens (jeffrey.r.stevens@gmail.com)

This file describes how to reproduce the results from the following paper: Winke, T. & Stevens, J.R. (2017). Is cooperative memory special? The role of costly errors, context, and social network size when remembering cooperative actions. Frontiers in Robotics and AI. doi:10.3389/frobt.2017.00052.

All materials presented here are released under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International license (CC BY-NC-SA 4.0). You may share and adapt this content with attribution, for non-commercial purposes if you ShareAlike (distribute any contributions under the same license).

To reproduce these results, first unzip winke_stevens_2017_rr.zip into a folder.  Then, ensure that a subfolder named "figures" is in the folder with all other files.  Next, open winke_stevens_2017_rcode.R and ensure that all packages mentioned at the top of the script are installed.  Once all packages are install, run the script in R using "source("winke_stevens_2017_rcode.R")".  

Once the script runs without errors, you can compile the R Markdown document winke_stevens_2017.Rmd.  Open this file in RStudio and ensure that you have packages 'knitr' and 'rmarkdown' installed.  Once installed, use knitr to compile the document (control-shift-k).  Use the same process to compile winke_stevens_2017_SM.Rmd.
