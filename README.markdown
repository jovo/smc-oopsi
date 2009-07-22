This is a repository containing the code to implement the algorithms employed in Vogelstein et al, 2009, and subsequent related advances.  Publication version of manuscript is now available from [Cell Press] [4] with subscription, or [here] [5] for free


This work may be cited as:

Vogelstein, Joshua T.; Watson, Brendon O.; Packer, Adam M.; Yuste, Rafael; Jedynak, Bruno; Paninski, Liam. Spike Inference from Calcium Imaging Using Sequential Monte Carlo Methods. Biophysical Journal, Volume 97, Issue 2, pp.636 - 655, 2009.


Repository Organization
=======================

* text: contains tex and compiled other files of submitted manuscript (and some minor revisions) 
* functions: contains functions required to use our smc algorithm to infer spikes (and a subfolder called for_figs containing older versions used in creating some of the figures)
* scripts: each script calls various functions, saves the data (simulated or analyzed), and generates figures
* data: raw and simulated data used to make figures
* figs: pdf files for each final figure
* general_background: a document containing some general background information that might be useful
* reviews&proofs: reviews and proofs from BJ to get this paper published

Quick Tips
==========

* final version of proof may be found [here] [1] 
* to learn how to use the algorithm, i suggest playing with scripts in the [script folder] [3]
* any questions or issues, please use the [issues tab] [2]
* any problems with that, or questions inappropriate for that forum, please contact me at: joshuav@jhu.edu
* please respect our policy of openness, and operate accordingly.

[1]: http://github.com/jovo/smc-oopsi/raw/master/reviews&proofs/proof_v6.pdf "here"
[2]: http://github.com/jovo/smc-oopsi/issues "issues tab"
[3]: http://github.com/jovo/smc-oopsi/tree/master/scripts "scripts folder"
[4]: http://www.cell.com/biophysj/abstract/S0006-3495(09)00311-7 "Cell Press"
[5]: http://github.com/jovo/smc-oopsi/raw/master/Vogelstein2009.pdf "here"

License
=======

License stuff (apparently, this is useful for some people):
Copyright (c) 2008 Joshua Vogelstein

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.