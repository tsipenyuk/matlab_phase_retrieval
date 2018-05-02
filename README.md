##Matlab Phase Retrieval Sandbox
----------------------------------------------------------------------------

This project is no longer supported.

Matlab Phase Retrieval Sandbox is a collection of scripts and utilities for
developement and testing of phase retrieval algorithms. 

There are no stable releases available yet.

###Installation
1. Download the repository.
2. Unpack the zip-file.
3. From the folder, launch MATLAB and type `phase_retrieval` for a quick 
tutorial on how to use this package.

###Licensing
All the code in the folder 'lib' and folders therein is developed by 
other parties.
The corresponding author credits, licensing information, and links to the original source code are stored in files 'lib\*\README.md'

###Structure  
phase_retrieval  
|-- cls  - contains definitions of classes that describe molecules and algorithms.  
|-- fun  - contains functions (class-independent), such as projections and algorithm update steps.  
|-- inst - contains scripts that instantiate various classes; simplifies matlab-shell use of the code.  
|-- lib  - contains third-party software.  
  
The code is split into procedural part ('fun') and object-oriented part ('cls', 'inst'). The procedural part
does not depend on any objects --- it contains algorithm updates (such as ER update) and auxilary functions 
(such as projections onto constraint sets). The object-oriented part is introduced to ease instantiation
and plotting of case studies.

----------------------------------------------------------------------------
(c) 2016 Arseniy Tsipenyuk (TUM M7)
