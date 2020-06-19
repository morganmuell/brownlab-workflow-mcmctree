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

