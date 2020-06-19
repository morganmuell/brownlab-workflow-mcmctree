# brownlab-workflow-mcmctree
Tutorial for using MCMCTree as implemented by Brown Lab at SIUC. For divergence time estimation on UCE datasets from dendrobatid frogs.

Here I will give a brief tutorial into using MCMCTree for divergence time estimation in dendrobatid frogs. We'll use example data to walk through prepping files for analysis and running the actual analysis. PAML is set up on the lab Mac, so we'll walk through as such, though everything should apply to any system. This also isn't a comprehensive deep dive into all of MCMCTree's features, so if you want to adjust more parameters and such, here are some great resources I used to learn the program: #get links
PAML website (here is where you can install the program if you don't have it yet) 
PAML manual
Jun Inuoue tutorial
Mario dos Reis tutorial
google group 

For this tutorial, we will be doing everything in PAML 4.8. This is not the most current version, which is 4.9j, but I was unable to get 4.9j to work because I found it would abort the MCMCTree run just before the end and fail to produce a treefile, which is the whole reason we want to run MCMCTree in the first place. If a new version comes out I would highly encourage you to give that a shot, since it's always better to use the most up-to-date packages. You will also need R handy and available for prepping calibration points.

##Before the analysis
Before you begin, you will need three files.
1) sequence file in phylip format
2) treefile #output from ML analysis
3) control file 

Your sequence file should be a concatenated matrix of all the UCEs you captured. This is the same sequence file you inputted into IQ-TREE for your maximum likelihood analysis.

The treefile is the output treefile of your maximum likelihood analysis. One thing that MCMCTree does to save analysis time is constrain the analysis to a predetermined topology, which you will provide here, in the form of an ML tree you have already run.

The control file is what MCMCTree uses to specify model parameters for the analysis.

Each of these requires some prep prior to starting the analysis.

###Prepping the sequence file
Here I have a sequence file for 9 samples of Ranitomeya and one Andinobates outgroup sample, comprised of a little less than 1200 loci and ~10,000 PIS.

Your sequence file should be in phylip format (I don't believe PAML will accept Nexus files, though I could be hallucinating on that one). Your sequence file should already be in phylip format anyway if you're taking the sequence you used for IQ-TREE. 

The second thing you have to do with your sequence file is add spaces between the name of the sample and the beginning of the sequence. PAML requires sequence files to have at least 2 spaces in between, which is not the default setting for a phylip file. Your sequence file is probably too large to edit manually (and you've got better things to do with your time than that, king/queen/non-binary monarch) so you can add these extra spaces with a simple line of code from the terminal.

```
sed -r -i 's/(.*) (.*)/\1  \2/g' example-seq.phylip
```

This tells the sed command to add two spaces in the file (which are above in between `\1` and `\2`) between the sample name and the beginning of the sequence. Just be sure to replace "sequence.phylip" with your own sequence name, obviously.

Now, your sequence file should be prepared. As always, do be 100% sure the names in your sequence file exactly match those in the treefile, but that should be done already if your treefile was produced using these sequence data.

#Prepping the treefile (heading)
Here I have a rooted tree. This is the output of a maximum likelihood analysis that has been rooted in FigTree and exported to a Newick file in the rooted form.

This step is important and multifaceted. To prepare your treefile you need to do the following in R:

1) remove branchlengths
2) add calibration point 

####Step 1: removing branch lengths
We need to remove the branchlengths from our maximum likelihood tree. MCMCTree requires input topologies not to have any branch lengths in the file, because we willcalculate branch lengths prior to the MCMCTree analysis.

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
(show result)

We still need to write the tree to a file to use, but before we do that we need to do step 2, calibration, in R as well.

####Step 2: calibrating treefile
Okay, here's what's arguably the most important step of divergence time estimation in any program. We have to give MCMCTree a reference for node ages so it knows where tostart in calculating the rest of ages. Ideally, for a group of organisms, this would include multiple calibration points based on fossil evidence of simulations tied to paleogeographic scenarios, such as the uplift of the Andes or other events that we know influenced when lineages could come about. We will use estimates calculated by Santos et al. (2009) to calibrate our analysis.

Below I describe how I got to the values I used to calibrate the node; skip over this section if you already have your node values and just need to know how to implement them.

As you probably already know, there are no poison frog fossils, making this process quite difficult. One study, Santos et al. (2009), has successfully determined divergence time estimates, by using fossil evidence from the rest of Amphibia to calibrate a tree for all of Amphibia, and in turn determining estimates for dendrobatid frogs. They also used three different paleogeographic scenarios to calibrate their analysis, and thus have calculated sets of divergence times (mean, standard dev, etc.) for each node in their analysis. All their estimates are in the supplement of their paper. This paper is an excellent resource, and we will use their estimates to calibrate our analysis.

We will only calibrate one node. For this tree, I decided to calibrate the node denoting the common ancestor of Ranitomeya and its sister genus, Andinobates. Typically in past analyses for the lab we have chosen to calibrate at a deeper node like this. In ideal scenarios you would calibrate multiple nodes at deep and shallow time scales, but since we only have this one resource we will only calibrate one deep node. There is a lot of literature out there about the effects of which nodes you select to calibrate and how that node's age affects calculation of divergence times, if you're more interested in this decision.

As I mentioned, Santos et al. (2009) calculated three sets of estimates based on three different paleogeographic scenarios. Here, I will average the means and standard deviations of those three estimates and use that to set the calibration point. However, if you have the time, I would encourage you to do three different MCMCTree runs, calibrate them using these three different numerical values, and compare the difference in results. For now, we will just use the average.

OK. So, if you average the 3 mean divergence times from Santos et al. (2009) for the node representing the common ancestor between Ranitomeya and Andinobates, you get an average mean of 12.651 million years ago. If you average the 3 standard deviations, you get an average standard deviation of 2.576 million years. We will use these two values to set the calibration.

Alright, so we know our values, and we know where we want to put them on the tree. However, we still need to decide on the shape of the distribution. This is important, because you're telling MCMCTree how likely it is based on your prior information that the value it calculates based on the sequence today gets closer to or further from this average. The MCMCTreeR package in R allows you to plot different shaped distributions and visualize what they look like.

I have used uniform distributions on my calibration nodes in MCMCTree for two reasons. One, by default, most other distribution shapes for calibration nodes in MCMCTree are automatically treated as fossils. This means that the tail for the older end of the distribution will be much longer than the tail at the younger end of the distribution (a "soft" boundary), which shoots down very abruptly (a "hard" boundary). The reason for this shape is that we have much more certainty about the minimum age of the fossil because that is the age that fossil we found was dated to and we know THAT fossil exists, but we have
much less certainty about the maximum age of the fossil. Because our calibration is already based on results of previous analyses, I do not want to attribute that sort of certainty to my priors. This is a good segue to my second reason for using the uniform distribution, which is that I wanted a diffuse prior on account of not really knowing much concrete information about when Ranitomeya's ancestors may have diverged, outside of estimates calculated by Santos et al. Uniform distributions have lead to problems with convergence in other programs like BEAST because your priors are too broad, but I have yet
to experience that problem with MCMCTree. This is NOT to say that that problem is not a possibility here; you should be aware of it. Nevertheless, I have still found the uniform distribution easier to implement and sufficient for my analyses thus far.

OK. So, if you average the 3 mean divergence times from Santos et al. (2009) for the node representing the common ancestor between Ranitomeya and Andinobates, you get an average mean of 12.651 million years ago. If you average the 3 standard deviations, you get an average standard deviation of 2.576 million years. Now, we can use those two values to calculate the boundaries of a 95% confidence interval, which rounds out to a minimum of
7.601 million years, and a maximum of 17.701 million years. We will use these values as the boundaries of our uniform distribution.

You can add calibrations manually by entering them into the treefile, or MCMCTreeR can do it for you. You can do it manually since we only have one, but it's easy to get confused about where to put it, so we're going to do it in R. There are lots of options for distribution shapes in R, but I will just show you how to do the uniform distribution. The R documentation for this package by Mark Puttick is very nice and has clear examples that show how to try out other distribution shapes at nodes.

First, let's look at our input tree.

```
plot(frogtree, cex=0.5) #plot your tree with a smaller font to see labels
nodelabels() #view node label numbers
```

Execute both at once to produce a tree with nodelabels. Take note of the
number that appears at our calibration node of interest. In this case, it
is 11.

[put the plot in here]

Next we'll make a few objects related to our calibration.
```
MinTime <- .07601 #minimum time for Andino-Ranito node calculated before
names(MinTime) <- "node_11" #label with name of node we saw in above tree
MaxTime <- .17701 # maximum time for Anino-Ranito calculated before
names(MaxTime) <- "node_11" #do the same with labelling
MonoGroup <- tipDes(frogtree, 11) #denotes monophyletic group at andino_ranito LCA
```
Like I mentioned in the annotation above, change the name assignment to 
whatever node number on the tree represents your calibration node.

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
Remember that with brackets in R, we are denoting `[rows,columns]`. This plots
the calibration node for the first row of parameters. We only have one row because we only have one calibration node, but if you ever use multiple calibration nodes, you may want to change this.

Now we can see the shape of the uniform distributions. There is no difference in weight between the bounds, and the upper and lower bounds are in place with equal probability of any value in between (this is visually a bit course, but
they are the boundaries of the plot).

[put in plot result]

If you're happy with the shape of your distribution again, run the uniform distribution command again, and also tell it to write the tree to a text file.
```
uniform_results <- estimateBound(minAge = MinTime, maxAge = MaxTime,
                                 monoGroups = MonoGroup, phy =
                                   frogtree, plot = FALSE, 
                                   writeMCMCtree=TRUE,                                   MCMCtreeName='calibrated_example.treefile')
```
Now, our tree has been calibrated and has branchlengths removed. You can now find it in your working directory.

###Prepping the control file
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
Of course, change the directories for seqfile and treefile so that they exactly match the names and locations of the files, and name your outfile whatever you like. 

`ndata` tells MCMCTree what kind of data the sequence data is. `0` is for nucleotide data,so we will always use it. 

`seqtype` will always be `0` for us.

`usedata` is very important and we will change it with different steps of the analysis. Setting it to `3` will have it run branch length approximation using BASEML, and setting it to `2` will have it run the MCMC chain. When setting it to `2`, you always have to remember to put in.BV after the name so it knows to use that file. Setting it to `0` tell it not to use sequence data. More on that later.

`clock` is set to `2` for an independent rates clock model. You can also select a strict clock or a couple other clock options detailed in the tutorials and manual.

`RootAge` sets the boundary for the age of the root. I chosen 20 million years as a liberal boundary for the root in comparison to estimates from Santos et al. (2009).

`model` specifies the nucleotide substitution models. As I understand it, there are currently 7 models available in MCMCTree: JC69 (0), K80 (1), F81 (2), F84 (3), HKY (4), T92 (5), TN93 (6), and GTR (7). I would choose whichever most closely matches the model you selected for on sequence data in IQ-TREE (using -m MFP when running IQ-TREE). In this case, we will use GTR, which is the closest matching model to what ModelFinder selected when I ran a maximum likelihood analysis on this dataset. Not all models in IQ-TREE are available in MCMCTree, so we picked the next best model selected by IQ-TREE that was available in MCMCTree.

The remaining parameters change the shape of the gamma prior on other parameters, but I mostly left them untouched. I would definitely encourage you to tinker with them as you are able to, if that is in your interest. As I understand it, the finetune function is deprecated.

The bottom 3 parameters are incredibly important. You need to set up your analysis so that your burn-in number matches 10% of the total number of generations of your analysis, just like any other Bayesian program. 

To set these last three parameters, set `nsample` somewhere between 10,000 and 20,000 (I choose 20k usually). Determine how many iterations you want to run. Say I want to run 20,000,000 total iterations. Set the sampling frequency `sampfreq` equal to 20,000,000/`nsample`, which is in this case 1,000. Then, you set your `burnin` to 2,000,000, which is 10% of your total iterations. Now your analysis is set up to have the appropriate burn-in percentage and collect the right number of samples to match your total desired iterations. If your analysis does not converge, simply increase the total number of iterations, and adjust your sampfreq and burn-in values accordingly.

A common reason analyses do not converge is not running them for long enough. The solution to this is to extend the run so that the `burnin` is longer. You DON'T want to do this by increasing `nsample`, because you will just burn up your computer and not even progress much further into treespace. Don't do it. INSTEAD, increase the `sampfreq` and adjust your `burnin` value accordingly to 10% of the new total iterations that comes from changing your `sampfreq`  value.

##Running the Analysis
Once you have your sequence file, treefile, and control file set up, running the analysis is quite simple and comes in three steps, plus two quality control steps. I am a little bit type A, so I typically run the following commands to set up my file structure prior to running the analysis.
```
mkdir run0
mkdir run1
mkdir run2
cp run.ctl run0/run0.ctl
cp run.ctl run1/run1.ctl
cp run.ctl run2/run2.ctl
```
The second thing to do before you start is to make sure the bin with the paml executables is in your path. This isn't necessary for MCMCTree, but in order to run baseml for branch length approximation, it has to be callable in your path when you run the command. Here's how to do it in your open terminal window if you haven't done it before.
```
$PATH # check to see what your current path looks like

PATH=$PATH:~/Desktop/paml4.8/bin # input path to where your executables are

$PATH # check to see if it was successfully added
```
###Step 1: Running analysis without sequence data
This step of the analysis we're going to run MCMCTree without sequence data, and only using model parameters and our calibration point. This will tell us whether our model priors are reasonable, or should be adjusted, and will tell us once we run the analysis whether or sequence data made much of a difference from our priors.

Now we're going to make one small adjustment to the ctl file for this part. Modify the usedata argument to 0 to tell the program not to use the sequence data. Don't modify the sequence data argument.
```
cd run0

open run0.ctl
```
Make sure to save the file!!

Then, run the mcmctree command:
```
mcmctree run0.ctl
```
This should not take that long.

###Step 2: Branch-length approximation
Another major step MCMCTree takes to cut down on analysis time is approximating branch length values prior to running the MCMC chains. Here, the program will call on BASEML to calculate branch length approximations, which will then be inputted into your MCMCTree analyses alongside your sequence data and calibrated constraint tree. This call on BASEML is the reason you had to put PAML in your path. 

Once again, we will make a minor adjustment to the control file. This time, we will change `usedata` to  `3`.
```
cd ../run1

open run1.ctl
```
Then, run the command!
```
mcmctree run1.ctl
```
This step will take awhile for a real dataset. My dataset with 1,176 and around 72k PIS for 67 taxa took about 4 days to run at 60 million iterations.

###Step 3: MCMCTree
Once step 2 is complete, you should have a file in the run1 file called out.BV. Rename the file to in.BV.
```
mv out.BV in.BV

open run1.ctl
```
Then, change your ctl file back to `usedata=2 in.BV`. You MUST put the file name `in.BV` after the `2`, or it will not run.

Then, run it again!
```
mcmctree run1.ctl
```
This will also take a few days for a larger dataset.

Once the analysis is complete, you should have a number of files, most notably of which are your outfile (whatever you named it), mcmc.txt, and a file called Figtree.tre. As suggested by the name, if you input this file into Figtree it will have divergence times and 95% confidence interval error bars. Raw values for mean node values and values for confidence intervals will be at the bottom of the outfile. It'll look something like this.

[figtree photo]

OK, now do it again! Run number 2. You shouldn't need to make any adjustments in the control file; the same branch length estimates you calculated for run #1 also apply to run #2 because you are using the same sequence data.
```
cp in.BV ../run2

cd ../run2

mcmctree run2.ctl
```
You should see the same sort of output as with run #1 with your outfiles, though the folder will not contain as many files because you didn't run the branch length approximation in this one.

###Checking for convergence
This is quality control step #1. Here, we will take two measures to assess whether our values converged between the two MCMC chains. This is where you will find out how large a grain of salt you need to take your results with. Some results, like a single ESS value of 199 on one of 50 nodes, are pretty easy to take. Others, such as 20 out of 50 nodes with ESS values at 49 and two runs with high dissimilar mean node values, are so large that they resemble an icy asteroid of poor nodal support hurtling toward Bayesian Earth, thereby you cannot take them. Regardless of my understanding of how the take-with-a-grain-of-salt analogy works, the take-home message is that convergence is very important!! 

First, we will look at both runs in Tracer. Open up Tracer and select `File -> Impot trace file`. Then, select the file mcmc.txt in the run1 folder. You can visualize the values at which the mcmc chain sampled each node. Ideally, you want a "fuzzy caterpillar"  distribution, which shows your analysis was centering around a value. I prefer to use the fourth tab to visualize.

[input picture of tracer]

On the left panel you will see a different line for each node in the analysis and an
ESS value (estimated sample size). This tells you how many times you MCMC chain sampled
a certain value. Generally, ESS values are regarded as good if they are higher than 200,
meaning that your analysis sampled that divergence time estimate for that node many times,
giving you confidence the analysis converged on a value for that node. Scroll through each
node to view the distribution, and check to see if they are above 200. Tracer will color
code the values if your ESS at that node did not reach 200.

Repeat this step for the mcmc.txt file in run2.

If all your nodes had ESS values above 200 for both runs, that's fantastic! Then you can
move to the second step to assess convergence described below. If they did not, extend
the length of your analysis, dramatically if you must. Extending the burn-in will give
your analysis more space to converge on a value. Remember, change the sampling frequency,
NOT the number of samples.

If you did have high enough ESS values, let's go to step 2 of assessing convergence. I 
have also included an Excel document I already set up in the tutorial folder called
convergence-check.xlsx (or whatever the file extension is). Open that up and there should
be three tabs. Use control-a (or command-a) to select all of the values in mcmc.txt for
run1 and paste them into the first tab. Repeat the same thing for the second tab with 
mcmc.txt from run 2.

Now go to the third tab. Here, you will see that the average node value for each node in
(each column in the first two tabs) for run 1 and run 2. You can visually assess how
similar those look in the line graph. If the line graph is a straight line when plotting
the two columns against each other, you have achieved convergence because the mean values
are approximately the same. I've included a function to tell you the slope of the line
(look at the m value in y=mx+b, CLASSIC) to give you a numerical feel for the similarity.
The further the slope value is from 1, the more different your values were from each 
other, and the culprits for those differences can be identified by seeing which nodes
did not have high enough ESS values on Tracer.

You will have to adjust  tab 3 with the averages for the number of nodes in your tree,
but this sheet should give you an easy-to-follow skeleton for how to set that up and save
you some time. The calculations are very simple.

I made the Excel step the second step because I didn't want hopes getting up with the
spreadsheet just to see cruddy values in Tracer. You should do both things to assess
convergence, though their results really go hand in hand. If you are struggling to achieve
full convergence, again, EXTEND YOUR ANALYSIS, particularly your burn-in. If you feel 
like you've tried a number of very long analyses with no luck, I would suggest revisiting
your data filtering schemes and testing the waters for different subsets of loci. I have
had instances where slightly smaller datasets achieving convergence, but adding 100 more
loci with varying numbers of PIS completely throw it off and make it impossible to achieve
convergence. It's a matter of trial and error.

#One last step!

Okay, so you've got your treefiles, your analyses converged, and you are victorious! Here
is one last thing to do.

MCMCTree's default units entail that .12 means 12 million years, and outputs values as
such. To change the values in your treefiles to shifted over decimal points that better
fit the timescale of the data we are working with, use the following in the search and
replace function when you open text files of the treefiles:

I hope this was helpful! If you have specific questions you think I might be able to
answer, feel free to contact me at mrm0161@auburn.edu. Otherwise post in the google group
where a number of brilliant folks including Dr. Yang himself, the writer of PAML, may
help you out faster and more concisely.

Tips for Troubleshooting: A list of errors that commonly upset me
-- if you are attempting to run Step 1 of the analysis with approximate branch length
estimation, be 100% sure you set usedata=3, NOT usedata=1. No less than 3 times, I put
1 instead of 3 because my brain got confused that step 1 of the analysis meant inputting
a 3. This results in a memory related error message. Don't be like me, don't do that.

--ABSOLUTELY DO NOT run the last step without branch length approximation. I never did 
this but your analysis will never finish

--if you're bugged out about the path thing, just use the full file name to get to the
program and that's fine, e.g., ~/Desktop/paml4.8/bin/mcmctree run1.ctl
