## GA Jigsaw puzzle solver
---------------
Project for rearranging fragments of scrambled images in order to obtain the original ones.

It's assumed that each scrambled image was obtained by dividing the original into square pieces of a certain size and shuffling them randomly.
## Requierements
---------------
The app is **multi-thread, cross-platform**, to run it only **Qt** is needed (i used 5.12 version).
## Setting up
--------------
To process the image do following: 

1. Configure paths variables and algorithm's constants (their meaning is described in code).
2. Run parser to convert scrambled images ('cause pure **C++** doesn't work with images).
3. Run algorithm (it'll write the result in specified place).
4. Run to-image to convert algorithm's output to image.

## Algorithm
------------
My solution is based on **Genetic algorithm**. It's made up of a few repeating steps:

1. Initial *population* of image-candidates is generated randomly.
2. Each of them is evaluated via *fitness function*.
3. Best candidates are *selected* and retained for the further phases.
4. *Crossing* of solutions occurs till a certain population's size.
5. If the termination condition hasn't been reached, go to step 2.,
    otherwise print the best found image-candidate.

## Aspects
-----------------------
The algorithm implements and augments the ideas of article *An Automatic Solver for Very Large Jigsaw Puzzles
Using Genetic Algorithms* prepared by *Dror Sholomon, Eli (Omid) David, Nathan S. Netanyahu*.

Information about genetic operators:

* *Mutation operator* cycle-shifts the image-candidate (to prevent dissarrangement of already matched pieces).
* *Crossover operator* consists of three phases prioritized in such a way: *2-agreement*, *best-buddy*, *greedy*.
* *Elitist selection* with a certain part of *randoms* (to enhance the diversity of genes) is provided.  

Algorithm allows you to customize some of its components to change the balance between run-time and quality of the output:

* Set *number of threads* used
* Change *dissimilarity metric* of image fragments (e.g. **L2,3 - norms, Cosine similarity**, etc.)
* Configure:
    * *number of runs*, *number of generations* in each run.
    * *population's size*, *retain ratio*.
    * *mutation rate*, *percentage of randoms* in selection.
    * when to start using *2-agreement* phase (to avert early convergence).

## Experimental results & observations
------------------------
The algorithm was tested on about **3600** images split into **32x32**, **16x16**, **8x8** fragments.

The observation is that the image's resolution practically doesn't impact the processing time as well as quality of output image, it more depends on the number of fragments.

After fitting constants the *neighbor comparison* measurement of packs was:

- **8x8** fragments - **99.4%** accuracy
- **16x16** fragments - **97.6%** accuracy
- **32x32** fragments - **76.2%** accuracy

It's seen that genetic algorithm gives near-perfect result on less than about **500** fragments, but then accuracy gradually falls.
## Author
---------
The whole project was contributed by me, [**Kiryl Shyrma**](https://github.com/ruct), **HSE Moscow**.

This **Jigsaw-solver** was developed for the following contest: [**HONORCUP**](https://codeforces.com/honorcup), where more than two hundreds of individuals and teams participated, i took the [**31st place**](https://codeforces.com/contest/1235/standings).
## Future plans
---------------
Project's still in development, the main goals are:

* Improve accuracy of algorithm.
* Add handling of non-square images split into non-square fragments.
* Optimise code to increase the perfomance of algorithm.
* Add **scrambler** so that user can test algo having original images.
* Simplify launch to one console command and without **Qt**.

## Contact
--------------
For any queries/issues regarding the project, you can contact me at **pufbsd@gmail.com** or at **@flyce32** (Telegram).
