# brownlab-workflow-mcmctree
Tutorial for using MCMCTree in the PAML package as implemented by the Brown Lab at SIUC. For divergence time estimation on UCE datasets, most geared toward data from dendrobatid frogs.

Here I will give a brief tutorial into using MCMCTree for divergence time estimation in dendrobatid frogs. If you haven't yet obtained a maximum likelihood tree for your dataset, you can learn how to do so using Wilson's tutorial on the other parts of the Brown Lab pipeline [here](https://github.com/wxguillo/brownlab-workflow). We'll use example data to walk through prepping files for analysis and running the actual analysis. PAML is set up on the lab Mac, so we'll walk through as such, though everything should apply to any system. This also isn't a comprehensive deep dive into all of MCMCTree's features, so if you want to adjust more parameters or need help with something not covered in this tutorial, here are some great resources I used to learn the program:
[PAML website](http://abacus.gene.ucl.ac.uk/software/paml.html)
[PAML user manual](http://abacus.gene.ucl.ac.uk/software/pamlDOC.pdf)
[Jun Inoue tutorial](http://www.fish-evol.org/mcmctreeExampleVert6/text1Eng.html)
[Mario dos Reis tutorial](http://abacus.gene.ucl.ac.uk/software/MCMCtree.Tutorials.pdf)
[PAML Google group](https://groups.google.com/forum/#!forum/pamlsoftware)

For this tutorial, we will be doing everything in PAML 4.8, which I will assume you already have installed (if not, there are easy instructions on the PAML website). This is not the most current version, which is 4.9j, but I was unable to get 4.9j to work because I found it would abort the MCMCTree run just before the end and fail to produce a treefile with confidence intervals. If a new version comes out I would highly encourage you to give that a try to see if this issue has been debugged, since it's always better to use the most up-to-date packages, or even give 4.9j a try on your own (all aspects of this tutorial will remain the same). You will also need R handy for parts of this tutorial.

## Before the analysis
Before you begin, you will need three files.
1) sequence file in phylip format
2) treefile #output from ML analysis
3) control file 

Your sequence file should be a concatenated matrix of all the UCEs you captured. This is the same sequence file you inputted into IQ-TREE for your maximum likelihood analysis.

The treefile is the output treefile of your maximum likelihood analysis. One thing that MCMCTree does to save analysis time is constrain the analysis to a predetermined topology, which you will provide here, in the form of an ML tree you have already run.

The control file is what MCMCTree uses to specify model parameters for the analysis.

Each of these requires some prep prior to starting the analysis.

### Prepping the sequence file
Here I have a sequence file for 9 samples of Ranitomeya and one Andinobates outgroup sample, comprised of a little less than 1200 loci and ~10,000 PIS.

Your sequence file should be in phylip format (I don't believe PAML will accept Nexus files, though I could be hallucinating on that one). Your sequence file should already be in phylip format anyway if you're taking the sequence you used for IQ-TREE. 

The second thing you have to do with your sequence file is add spaces between the name of the sample and the beginning of the sequence. PAML requires sequence files to have at least 2 spaces in between the sample name and beginning of the sequence, which is not the default formatting for a phylip file. Your sequence file is probably too large to edit manually, and you've got better things to do with your time than that, so you can add these extra spaces with a simple line of code from the terminal.

```
sed -r -i 's/(.*) (.*)/\1  \2/g' example-seq.phylip
```

This tells the `sed` command to add two spaces in the file (which are above in between `\1` and `\2`) between the sample name and the beginning of the sequence.

Now, your sequence file should be prepared. As always, do be 100% sure the names in your sequence file exactly match those in the treefile, but that should be done already if your treefile was produced using these sequence data.

## Prepping the treefile
Here I have a rooted tree. This is the output of a maximum likelihood analysis that has been rooted in FigTree and exported from Figtree to a Newick file maintaining the rooted form of the tree. We will constrain the MCMCTree analysis to this topology, which saves a lot of computation time.

To prepare your treefile you need to do the following in R:

1) remove branchlengths
2) add calibration point(s) 

#### Step 1: removing branch lengths
We need to remove the branchlengths from our maximum likelihood tree. MCMCTree requires input topologies not to have any branch lengths in the file, because we will calculate branch lengths prior to the MCMCTree analysis.

We will be removing branch lengths and calibrating the tree in R. Make sure to change the working directory to wherever your data is.

```
setwd("~/Desktop/PAML_tutorial")
library(MCMCtreeR)
```
First, we'll load in the tree.
```
frogtree <- read.tree("uncalibrated-example.treefile")
```
Now we'll remove the branch lengths.
```
frogtree$edge.length<-NULL
frogtree
```
Now, when you type `frogtree`, R should tell you you have a rooted tree with no branch lengths.

#### Step 2: Calibrating treefile
Okay, here's what's arguably the most important step of divergence time estimation in any program. We have to give MCMCTree a calibration node to use as a reference for determining node ages. Ideally this would include multiple calibration points based on fossil evidence or simulations tied to paleogeographic evidence, however we are not so lucky. There are no poison frog fossils, making the calibration process quite difficult. One study, Santos et al. (2009), has successfully determined divergence time estimates using fossil evidence from the rest of Amphibia to determine estimates for Dendrobatidae. They also used three different paleogeographic scenarios to calibrate their analysis, and thus have calculated sets of divergence times (sets including mean, standard dev, etc.) for each node in their tree. All their estimates are in the supplement of their paper, which can be found [here](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1000056).

We will only calibrate one node, at the common ancestor of *Ranitomeya* and its sister genus, *Andinobates*. Typically in past analyses for the lab we have chosen to spend our one calibration at a deeper node, though obviously multiple calibrations at both shallow and deep timescales would have been ideal. I will use the average of the three means and the average of the three standard deviations derived for the common ancestor of *Ranitomeya* and *Andinobates* by Santos et al. (2009) to set the calibration point. However, if you have the time, I would encourage you to do three different MCMCTree runs, calibrate them using these three different numerical values, and compare the difference in results. For now, we will just use the average. 

If you average the three mean divergence times from Santos et al. (2009) for the node representing the common ancestor between *Ranitomeya* and *Andinobates*, you get an average mean of 12.651 million years ago. If you average the three standard deviations, you get an average standard deviation of 2.576 million years. Now, we can use those two values to calculate the boundaries of a 95% confidence interval, which come out to a minimum of 7.601 million years, and a maximum of 17.701 million years. We will use these values to set our calibration.

The package `MCMCTreeR` will visualize our calibrations and write them to a new treefile. You could write in the calibration manually since we only have one calibration point, but it's easy to get confused about where to put it, so we're going to do it in R. There are lots of options for distribution shapes with this package (and distribution shapes are very important in calibration for determining the likelihood of an average straying from a prior), but I will just show you how to do the uniform distribution, which I chose to use for a couple reasons you can feel free to inquire with me about. That being said, [the guide] (https://github.com/PuttickMacroevolution/MCMCtreeR) to MCMCTreeR by Mark Puttick is very nice and has clear examples that show you how to look at other distribution shapes.

First, let's look at our input tree.
```
plot(frogtree, cex=0.5) #plot your tree with a smaller font to see labels
nodelabels() #view node label numbers
```
Execute both at once to produce a tree with node labels. Take note of the number that appears at our calibration node of interest. In this case, it is 11.

[put the plot in here]

Next we'll make a few objects related to our calibration.
```
MinTime <- .07601 #minimum time for Andinobates-Ranitomeya node calculated before
names(MinTime) <- "node_11" #label with name of node we saw in above tree
MaxTime <- .17701 # maximum time for Andinobates-Ranitomeya calculated before
names(MaxTime) <- "node_11" #do the same with labeling
MonoGroup <- tipDes(frogtree, 11) #denotes monophyletic group at Andinobates-Ranitomeya LCA
```
Like I mentioned in the annotation above, change the name assignment to whatever node number on the tree represents your calibration node. Make sure to change the value of the number in the monophyletic group section based on how many samples are in your tree.

Alright, now let's use that to calculate the results of a uniform distribution.
```
uniform_results <- estimateBound(minAge = MinTime, maxAge = MaxTime, 
                                 monoGroups = MonoGroup, phy = frogtree,
                                 plot = FALSE)

uniform_results$parameters 
```
[put in results of uniform parameters]

In the uniform distribution, tL is the lower bound and tU is the upper bound. pL and pU are the probabilities of the lower and upper bounds being violated, with default values of 0.025.

Now, let's plot those results to make sure they match what we have in mind for the calibration.
```
plotMCMCtree(uniform_results$parameters[1, ],
             method = "bound", title = "Uniform Andino-Ranito Node", 
             upperTime = MaxTime + 0.5)
```
Remember that with brackets in R, we are denoting `[rows,columns]`. This plots the calibration node for the first row of parameters. We only have one row because we only have one calibration node; if we had multiple calibration nodes, we would change this.

Now we can see the shape of the uniform distribution. There is no difference in weight between the bounds (both are hard boundaries, not soft boundaries), and the upper and lower bounds are in place with equal probability of any value in between (this is visually a bit coarse, but they are the appropriate boundaries of the plot).

[put in plot result]

If you're happy with the shape of your distribution, run the uniform distribution command again, and also tell it to write the tree to a text file.
```
uniform_results <- estimateBound(minAge = MinTime, maxAge = MaxTime, monoGroups = MonoGroup, phy = frogtree, plot = FALSE, 
                                 writeMCMCtree=TRUE, MCMCtreeName='calibrated_example.treefile')
```
Now, our tree has been calibrated and has branchlengths removed. You can now find it in your working directory.

### Prepping the control file
The control file has a lot of parts. We will only be adjusting some of them, so I will only describe those below. I would encourage you to understand the others using the PAML manual in case you think changing other parameters may aid in your analysis, though.
```
          seed = -12345 
       seqfile = ../muscle-nexus-clean-70p_mk2_top75.phylip
      treefile = ../calibrated_finalML_1176loci.treefile
       outfile = out

         ndata = 0
       seqtype = 0  
       usedata = 2 in.BV
         clock = 2    
       RootAge = <.20  

         model = 7  
         alpha = 1    
         ncatG = 4    

     cleandata = 0   

       BDparas = 1 1 0.1    
   kappa_gamma = 6 2      
   alpha_gamma = 1 1 

   rgene_gamma = 2 2   
  sigma2_gamma = 1 10    

      finetune = 0:  0.04 0.2 0.3 0.1 0.3

         print = 1
        burnin = 2000000
      sampfreq = 1000
       nsample = 20000
 ```
Of course, change the directories for `seqfile` and `treefile` so that they exactly match the names and locations of the files, and name your `outfile` whatever you like. 

`ndata` tells MCMCTree what kind of data the sequence data is. `0` is for nucleotide data, so we will always use it. 

`seqtype` will always be `0` for us, as other codes deal with protein sequences.

`usedata` is very important and we will change it with different steps of the analysis. Setting it to `3` will have the program run branch length approximation using BASEML, and setting it to `2` will have the program run the MCMC chain. When setting it to `2`, you always have to remember to put in.BV after the name so it knows to use that file. Setting it to `0` tell it not to use sequence data. More on that later.

`clock` is set to `2` for an independent rates clock model. You can also select a strict clock or a couple other clock options detailed in the tutorials and manual.

`RootAge` sets the boundary for the age of the root. I chosen 20 million years as a liberal boundary for the root in comparison to estimates from Santos et al. (2009).

`model` specifies the nucleotide substitution models. As I understand it, there are currently 7 models available in MCMCTree: JC69 (0), K80 (1), F81 (2), F84 (3), HKY (4), T92 (5), TN93 (6), and GTR (7). I would choose whichever most closely matches the model you selected for on sequence data in IQ-TREE (using -m MFP when running IQ-TREE). In this case, we will use GTR, which is the closest matching model to what ModelFinder selected when I ran a maximum likelihood analysis on this dataset. Not all models in IQ-TREE are available in MCMCTree, so we picked the next best model selected by IQ-TREE that was available in MCMCTree.

The remaining parameters change the shape of the gamma prior on other parameters, but I mostly left them untouched. I would definitely encourage you to tinker with them as you are able to, if that is in your interest. 

I read in a guide that the `finetune` function is deprecated in recent versions of PAML, so I did not worry about it.

The bottom 3 parameters are incredibly important. You need to set up your analysis so that your `burnin` number matches 10% of the total number of generations of your analysis, just like any other Bayesian program. 

To set these last three parameters, set `nsample` somewhere between 10,000 and 20,000 (I choose 20k usually). Determine how many iterations you want to run. Say I want to run 20,000,000 total iterations. Set the sampling frequency `sampfreq` equal to 20,000,000/`nsample`, which is in this case 1,000. Then, you set your `burnin` to 2,000,000, which is 10% of your total iterations. Now your analysis is set up to have the appropriate burn-in percentage and collect the right number of samples to match your total desired iterations. If your analysis does not converge, simply increase the total number of iterations, and adjust your sampfreq and burn-in values accordingly.

A common reason analyses do not converge is not running them for long enough. The solution to this is to extend the run so that the `burnin` is longer. You DON'T want to do this by increasing `nsample`, because you will just burn up your computer and not even progress much further into treespace. Don't do it. INSTEAD, increase the `sampfreq` and adjust your `burnin` value accordingly to 10% of the new total iterations that comes from changing your `sampfreq`  value.

## Running the Analysis
Once you have your sequence file, treefile, and control file set up, running the analysis is quite simple and comes in three main steps, plus two quality control steps. I am a little bit type A, so I typically run the following commands to set up my file structure prior to running the analysis.
```
mkdir run0
mkdir run1
mkdir run2
cp run.ctl run0/run0.ctl
cp run.ctl run1/run1.ctl
cp run.ctl run2/run2.ctl
```
The second thing to do before you start is to make sure the bin with the PAML executables is in your path. This isn't necessary for MCMCTree on its own, but in order to run branch length approximation, BASEML (another program in the PAML package) has to be callable in your path when you run the command. Here's how to do it in your open terminal window if you haven't done it before.

First, check to see what your current path looks like.
```
$PATH
```
If the path to the folder to the PAML folder does not appear, add it like this (changing the path based on where YOUR PAML executables are located, of course).
```
PATH=$PATH:~/Desktop/paml4.8/bin
```
Now, check your path again to see if it was successfully added.
```
$PATH 
```
You should see the path to the PAML bin added to the end of the list of directories in your path.

### Step 1: Running analysis without sequence data
In this step of the analysis we're going to run MCMCTree without sequence data, only using model parameters and our calibration point. This will tell us whether our model priors are reasonable, or should be adjusted, and will tell us once we run the analysis whether or sequence data made much of a difference from our priors.

We're going to make one small adjustment to the control file for this part. Modify the `usedata=2 in.BV` argument to `usedata=0` to tell the program not to use the sequence data. Don't keep `in.BV` after the 0, and don't modify the sequence data argument.
```
cd run0

open run0.ctl
```
Make sure to save the file!

Then, run the mcmctree command:
```
mcmctree run0.ctl
```
This should not take very long.

### Step 2: Branch-length approximation
Another major step MCMCTree takes to cut down on analysis time is approximating branch length values prior to running the MCMC chains. Here, the program will call on BASEML to calculate branch length approximations, which will then be inputted into your MCMCTree analyses alongside your sequence data and calibrated constraint tree. This step is the reason you put BASEML in your path.

Once again, we will make a minor adjustment to the control file. This time, we will change `usedata=2 in.BV` to  `usedata=3`.
```
cd ../run1

open run1.ctl
```
Then, run the command!
```
mcmctree run1.ctl
```
This step will take awhile for a real dataset. For example, my Masters dataset with 1,176 loci and around 72k PIS for 67 taxa took about 4 days to run at 60 million iterations.

### Step 3: MCMCTree
Once step 2 is complete, you should have a file in the run1 folder called `out.BV`. Rename the file to `in.BV`.
```
mv out.BV in.BV

open run1.ctl
```
Then, change your ctl file to modify `usedata=3` back to `usedata=2 in.BV`. You MUST put the file name `in.BV` after the `2`, or it will not run.

Then, run it again!
```
mcmctree run1.ctl
```
This will also take a few days for a larger dataset.

Once the analysis is complete, you should have a number of files, most notably of which are your outfile (whatever you named it), `mcmc.txt`, and a file called `Figtree.tre`. Your file `mcmc.txt` shows you the numerical values the MCMC chain sampled, and you will use data stored in this file to assess convergence. You do not need to wait to do run #2 to assess ESS values for convergence, in fact, you should check ESS values before potentially wasting time on run #2 to see if your run was long enough for run #1. See below. As for your `Figtree.tre` file, as suggested by the name, if inputted into Figtree you will see mean divergence time values at each node and 95% confidence interval error bars on each node. Raw values for mean node values and values for confidence intervals will also be at the bottom of the outfile. Your tree will look something like this.

[figtree photo]

OK, now do it again! Run number 2. You shouldn't need to make any adjustments in this control file; the same branch length estimates you calculated for run #1 also apply to run #2 because you are using the same sequence data.
```
cp in.BV ../run2

cd ../run2

mcmctree run2.ctl
```
You should see the same set of output files as with run #1, though the folder will not contain as many files because you didn't run the branch length approximation in this one.

### Checking for convergence
This is quality control step #1. Here, we will take two measures to assess whether our values converged between the two MCMC chains. This is where you will find out how large a grain of salt you need to take your results with. Some results, like a single ESS value of 199 on one of 50 nodes, are pretty easy to take. Others, such as 20 out of 50 nodes with ESS values at 49 and two runs with high dissimilar mean node values, are so large that they resemble an icy asteroid of poor nodal support hurtling toward Bayesian Earth, thereby you cannot take them. Regardless of my understanding of how the take-with-a-grain-of-salt analogy works, the take-home message is that convergence is very important!! 

First, we will look at both runs in Tracer. Open up Tracer and select `File -> Import trace file`. Then, select the file `mcmc.txt` in the run1 folder. You can visualize the values at which the mcmc chain sampled each node. Ideally, you want a "fuzzy caterpillar"  distribution, which shows your analysis was converging around a value. I prefer to use the fourth tab to visualize.

[input picture of tracer]

On the left panel you will see a different line for each node in the analysis and an ESS value (estimated sample size). This tells you how many times you MCMC chain sampled a certain value. Generally, ESS values are regarded as good if they are higher than 200, meaning that your analysis sampled that divergence time estimate for that node many times, giving you confidence the analysis converged on a value for that node. Scroll through each node to view the distribution, and check to see if they are above 200. Tracer will color code the values if your ESS at that node did not reach 200.

Repeat this step for the mcmc.txt file in run2.

If all your nodes had ESS values above 200 for both runs, that's fantastic! Then you can move to the second step to assess convergence described below. If they did not, extend the length of your analysis, dramatically if you must. Remember, change the sampling frequency, NOT the number of samples.

If you did have high enough ESS values, let's go to step 2 of assessing convergence. I have also included an Excel document I already set up in the tutorial folder called convergence-check.xlsx. Open that up and there should be three tabs. Use control-a (or command-a) to select all of the values in `mcmc.txt` for run1 and paste them into the first tab. Repeat the same thing for the second tab with  mcmc.txt from run 2.

Now go to the third tab. Here, you will see that the average node value for each node in (each column in the first two tabs) for run 1 and run 2. You can visually assess how similar those look in the scatter plot. If the plot shows a straight line when plotting the two columns against each other, you have achieved convergence because the mean values for run 1 and run 2 at each node (x and y values respectively in the plot) are approximately the same. I've also included a function to tell you the slope of the line to give you a different numerical feel for the similarity. The further the slope value is from 1, the more different your values were from each other, and the probable culprits for those differences can be identified by seeing which nodes did not have high enough ESS values on Tracer. It is okay if the values are not EXACTLY the same; decimal points will be off slightly. The important thing is how overall close they are.

You will have to adjust tab 3 with the averages for the number of nodes in your tree, but this sheet should give you an easy-to-follow skeleton for how to set that up and save you some time. The calculations are very simple.

I made the Excel step the second step because I didn't want hopes getting up with the spreadsheet just to see cruddy values in Tracer. You should do both things to assess convergence, though their results really go hand in hand. If you are struggling to achieve full convergence, again, EXTEND YOUR ANALYSIS, particularly your burn-in. If you feel like you've tried a number of very long analyses with no luck, I would suggest revisiting your data filtering schemes and testing the waters for different subsets of loci. I have had instances where slightly smaller datasets achieve convergence, but adding 100 more loci with varying numbers of PIS completely throw it off and make it much more difficult to achieve convergence. It's a matter of trial and error.

### One last step!

Okay, so you've got your output treefiles, your analyses converged, and you are victorious! Here is one last thing to do.

MCMCTree's default units entail that .12 means 12 million years, and outputs values as such. To change the values in your output treefiles to shifted over decimal points that better fit the timescale of the data we are working with, use the following in the search and replace function when you open text files of the treefiles:

[put here]

I hope this was helpful! If you have specific questions you think I might be able to answer, feel free to contact me at mrm0161@auburn.edu. Otherwise post in the PAML google group I linked at the top of the page, where a number of brilliant folks in the community, including Dr. Yang himself (the writer of PAML), may help you out faster and more concisely.

## Citing PAML and Santos et al. (2009)

Putting this here just to be extra sure this program gets its credit. Here are recommended citations by the author:

[put them in]

And here is the citation for Santos et al. (2009) for our source of calibrations:

[put it here]

### Tips for Troubleshooting: A list of errors that commonly ruined my day

--If you are attempting to run Step 1 of the analysis with approximate branch length estimation, be 100% sure you set `usedata=3`, NOT `usedata=1`. No less than 3 times, I put 1 instead of 3 because my brain got confused that step 1 of the analysis meant inputting a 3. This results in a memory related error message (OOm), at least for large datasets. Don't be like me, don't do that.

--ABSOLUTELY DO NOT run the last step without branch length approximation. Your analysis will never finish

--It's okay to use the full path name for the program if you're new to all this and struggling to put the program in your path, e.g., `~/Desktop/paml4.8/bin/mcmctree run1.ctl`, however it is still necessary to put BASEML in the path no matter what to run step 1
