##################################################
####         RNAseq data SIMULATION           ####
##################################################


# By Sonia Tarazona
# Created: 4-oct-2012



### Function to simulate RNA-seq data

simula.counts = function (counts0 = NULL, ngenes = 20000, nrepl = 5,  
                          depth = 30, propdeg = 0.05, noise = 0.1, beta = 6) {
  
  # counts0: Counts per gene used to derive counts in condition 1. 
  #          If NULL, random numbers are generated from a power-law distribution.
  # ngenes: Number of genes in the dataset when counts0 is NULL.
  # nrepl: Number of replicates per condition.
  # depth: Minimum estimation for the total number of counts (in millions) for each sample
  # propdeg: Proportion of genes expected to be differentially expressed.
  # noise: Variability for mu_0 between conditions.
  # beta: Beta distribution parameter (shape2)
 
  
  depth = depth*10^6
  ncond = 2

  
  if (length(nrepl) == 1) {
    nrepl = rep(nrepl, ncond)
  }
  
  
  if (!is.null(counts0)) {   # counts0 are given
    ngenes = length(counts0)
    
    if (is.null(names(counts0))) {
      mygenes = paste("g", formatC(1:ngenes, width = nchar(ngenes), 
                                   format = "d", flag = "0"), sep ="")
    } else {
      mygenes = names(counts0)
    }
        
  } else {   # counts0 are randomly generated
    potencia = 0.5
    x = 1:round(depth/1000,0)
    myprob = x^(-potencia)    
    counts0 = sample(x, size = ngenes, replace = TRUE, prob = myprob)
    
    mygenes = paste("g", formatC(1:ngenes, width = nchar(ngenes), 
                                 format = "d", flag = "0"), sep ="")
  }
  
  names(counts0) = mygenes
  
  # DEG
  ndeg = round(propdeg*ngenes, 0)  # number of DEG
  
  deg = sample(mygenes, ndeg)  # names of DEG
  
  p.up = runif(1, 0.25, 0.75)  # percentage of DEG up-regulated
  
  deg1 = sample(deg, round(p.up*ndeg,0))  # names of DEG up-regulated
  
  deg2 = setdiff(deg, deg1)  # names of DEG down-regulated
  
  k1 = k2 = 5  # to avoid multiply foldchange by 0 or very low counts


  forma1 = 1.5
  forma2 = beta
  
  
  # DEG-up        
  cond1.change = rbeta(length(deg1), shape1 = forma1, shape2 = forma2)*100 + 1.5
  
  counts1 = counts0
  counts1[deg1] = (counts1[deg1] + k1) * cond1.change
  
  
  # DEG-down
  cond2.change = rbeta(length(deg2), shape1 = forma1, shape2 = forma2)*100 + 1.5
  
  counts2 = counts0
  counts2[deg2] = (counts2[deg2] + k2) * cond2.change
  
    
  
  # Adjusting counts for the desired depth  
  counts1 = round(depth * counts1/sum(counts1),0)
  names(counts1) = mygenes
  
  counts2 = round(depth * counts2/sum(counts2),0)
  names(counts2) = mygenes
    

  
  # Condition 1    
  cond1 = gener.cond(counts = counts1, noise = noise, nrepl = nrepl[1])
  colnames(cond1) = paste("cond1", 1:nrepl[1], sep = ".")
  rownames(cond1) = mygenes
      
  
  # Condition 2    
  cond2 = gener.cond(counts = counts2, noise = noise, nrepl = nrepl[2])
  colnames(cond2) = paste("cond", 2, ".", 1:nrepl[2], sep = "")
  rownames(cond2) = mygenes

          
    
  # Results
  
  if (propdeg > 0) {
    
    changes = data.frame("gene" = deg1, "regu" = 1, "foldchange" = cond1.change)
    changes = rbind(changes,
                    data.frame("gene" = deg2, "regu" = 2, "foldchange" = cond2.change))
    changes$gene = as.character(changes$gene)
        
  } else { changes = NULL }  
  
  
  list("simucounts" = cbind(cond1,cond2), "change" = changes, 
       "changedepth" = sum(counts1)-sum(counts2))
  
  
}





##-------------------------------------------------------------------##
##-------------------------------------------------------------------##
##-------------------------------------------------------------------##
##-------------------------------------------------------------------##



### Function phi(mu) obtained from real datasets

# (mu,phi) taken from real data sets + interpolation function
lasmus = c(0.388888888888889,0.708333333333333,0.881038647342995,1.03381642512077,1.19444444444444,
           1.3780193236715,1.58913043478261,1.85,2.15217391304348,2.4993961352657,2.90277777777778,
           3.36111111111111,3.90485185185186,4.55581803542673,5.31763285024155,6.16111111111111,
           7.12608695652174,8.26358695652174,9.64583333333333,11.2518115942029,13.0646376811594,
           15.1733333333333,17.6209239130435,20.4375,23.6354166666667,27.3089814814815,
           31.3814452495974,35.8328804347826,40.8017451690821,46.3384299516908,52.6304347826087,
           59.5,66.6265217391304,74.3904106280193,82.9713768115942,92.4259118357488,
           102.620053542673,113.546009661836,125.34107568438,137.798917874396,151.36309178744,
           165.45229468599,180.591425120773,196.93825,213.954101449275,232.221954508857,
           251.616769726248,272.179304347826,293.808748792271,317.25170531401,342.512927536232,
           369.294992753623,397.492036231884,428.298376811594,462.24537037037,498.544742351047,
           538.499621980676,580.50025,625.7,675.786956521739,732.746376811594,798.775333333334,
           874.794357487924,959.35825925926,1053.70944444444,1163.35739774557,1294.77174879227,
           1456.64353623188,1660.2097294686,1919.08176811594,2264.31477938808,2813.65251851852,
           3755.676,5910.27933333333,369699.430057971)

lasphi = c(1.4921875,2.1625,2.65,2.4703125,2.4375,2.1375,1.94375,1.7828125,1.5890625,1.45078125,
           1.2921875,1.203515625,1.1140625,1.0265625,0.923046875000001,0.7978515625,0.7607421875,
           0.718164062500001,0.625,0.5740234375,0.5375,0.5091796875,0.486328125,0.458203125,
           0.457421875,0.38359375,0.384375,0.3875,0.36318359375,0.33466796875,0.333203125,
           0.323828125,0.3205078125,0.2845703125,0.27119140625,0.2830078125,0.2716796875,
           0.265234375,0.24609375,0.2521484375,0.2572265625,0.244921875,0.22734375,0.22275390625,
           0.240625,0.210546875,0.20947265625,0.2203125,0.228125,0.21796875,0.1970703125,0.1859375,
           0.18828125,0.1763671875,0.178173828125,0.175390625,0.1888671875,0.178125,0.1591796875,
           0.169140625,0.1595703125,0.158203125,0.149609375,0.1580078125,0.15556640625,0.14453125,
           0.1390625,0.14453125,0.13046875,0.139453125,0.1353515625,0.140625,0.15078125,
           0.1775390625,0.2140625)


phimean = approxfun(lasmus, lasphi, rule = 2)




##-------------------------------------------------------------------##




### Function to generate the replicates for a given condition from Binomial distribution + noise

gener.cond = function (counts, noise, nrepl) {
 
  # Allowing for some noise in the NB mean
  counts.noise = t(sapply(counts, function(x) { x + c(-1,1)*noise*x }))
  counts.noise[which(counts.noise < 0)] = 0
  
  mu.noise = apply(counts.noise, 1, function (x) { runif(1, x[1], x[2]) })
  #mu.noise = counts + runif(length(counts), -noise, noise)
  #mu.noise = sapply(mu.noise, function(x) { max(x,0) })
  
  # Computing phi
  phi.mean = phimean(mu.noise)
  phi.sd = 1 / (1 + mu.noise^0.25)
  
  phi = apply(cbind(phi.mean, phi.sd), 1, function (x) { max(0.000001,rnorm(1, x[1], x[2])) } )
    
  # Results
  t(apply(cbind(mu.noise, phi), 1, function (x) { rnbinom(n = nrepl, size = 1/x[2], 
                                                   mu = max(x[1], 0.1)) }))    
} 





##-------------------------------------------------------------------##



