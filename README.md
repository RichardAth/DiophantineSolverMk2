# DiophantineSolverMk2
Solves diophantine equations of the form Ax² + Bxy + Cy² +Dx + Ey +F = 0

This is based entirely on Dario Alpert's well-known solver.

His version was written as a web page using Java around 2003, and is still available. However modern browsers don't like Java. Firefox no 
longer support Java at all. I just could not get it to work with Internet Explorer, Edge or Chrome so I downloaded the source and converted 
it to a console program, which does work.

To make it work as a console program it was necessary to remove all the HTML tags, which means that the output is not as 'pretty' and the 
layout is a bit different.

Many other technical changes were necessary:

In an effort to tidy it up I created a C++ equivalent program. In this program the functions were restructured radically and many 
comments added to make it a bit easier to follow. Also, I discoverd that there were two completely different types of home-made 
bigintegers used, whch needed completely different functions to manipulate them. I replaced both types with GMP/MPIR extended precision 
functions, whch allowed me to simply remove a significant amount of code, and we can be confident that the extended-precision functions 
are well-documented and reliable. The division in particular was troublesome. It turned out that in some cases it is essential to use 
floor' division and in other cases 'truncation' division must be used.

I also grabbed Dario Alpert's web page describing his methods and converted it to word-processor document.
