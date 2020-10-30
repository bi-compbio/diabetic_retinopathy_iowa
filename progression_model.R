library(caret)
library(parallel)
library(doParallel)
library(foreach)
library(data.table)
library(MEIGOR)
source('./funcs.R')


### load data
# corrected data
iowa.full = fread('./results/data_corrected.RData') %>%
  data.frame %>%
  tibble::column_to_rownames('ensemblID')
# metadata
iowa.metadata = fread('./data/TableS1.csv') %>%
  dplyr::filter(sampleID %in% colnames(iowa.full)) %>%
  data.frame %>%
  magrittr::set_rownames(.$sampleID)

### global settings
# samples for run
samples = rownames(iowa.metadata[iowa.metadata$sample_site == 'Macula', ])
# features to include
features = rownames(iowa.full)
# threads
n.threads = 20
# k-fold split
n.fold = 5
# k-fold repeats
n.repeats = 20       # number of k-fold repeats
# seed
seed = 420            # global seed
# max number of feval for optimisation
n.feval = 100
# max number of features
max.n = 1500
# output file
file.out = './results/df.cost.macula.txt'

### subset data, numeric disease variable
# subset samples
meta = iowa.metadata %>%
  dplyr::filter(disease_group_detailed != 'PDR + DME') %>%
  dplyr::filter(sampleID %in% samples) %>%
  mutate(progAB = factor(disease_group_detailed, 
    levels=c('Control', 'Diabetic', 'NPDR', 'NPDR + DME')) %>% 
    as.numeric) %>%
  mutate(progBA = factor(disease_group_detailed, 
    levels=c('Control', 'NPDR', 'Diabetic', 'NPDR + DME')) %>% 
    as.numeric) %>%
  dplyr::select(sampleID, disease_group_detailed, progAB, progBA) %>%
  data.frame %>%
  tibble::column_to_rownames('sampleID')
# align data and metadata
data = t(iowa.full[features, rownames(meta)])

### prepare cross validation and optimizatin
# grid
splsGrid = expand.grid(
  eta = seq(from = 0.1, to = 0.9, by = 0.05), 
  K = 1:15, 
  kappa = 0.5)
set.seed(seed)
seeds = lapply(1:(n.fold*n.repeats), function(x) sample.int(n=1000, nrow(splsGrid)))
seeds[[(n.fold*n.repeats) + 1]] = sample.int(1000, 1)
fitControl = trainControl(
  method = 'repeatedcv',
  search = 'grid',
  number = n.fold,
  repeats = n.repeats,
  allowParallel = T,
  seeds = seeds)
# register parallel backend
cl = makePSOCKcluster(n.threads)
registerDoParallel(cl)
# set up optimization problem
x_L =  c(1.1, 1.1, 0.1, 1, 50)
x_U =  c(3.9, 3.9, 0.9, 15, max.n)
problem.settings = list(
  f = getRMSE_filter, 
  x_L = x_L, 
  x_U = x_U,
  int_var = 2)
# optimization options
opts = list(maxeval=n.feval, local_solver=0)

### ranking of features
df.rank = lapply(c('progAB', 'progBA'), function(x)
  cor(data, meta[, x], method='spearman') %>%
    data.frame %>%
    magrittr::set_names('correlation') %>%
    tibble::rownames_to_column('ensemblID') %>%
    arrange(-abs(correlation)) %>%
    mutate(rank = 1:nrow(.))) %>%
  magrittr::set_names(c('progAB', 'progBA')) %>%
  bind_rows(.id = 'prog') %>%
  data.table


### MAIN
# estimate parameters
result = MEIGO(
  problem = problem.settings, 
  opts = opts, 
  algorithm = 'eSS', 
  data, meta, df.rank, fitControl, file.out)

# get features
d1.hat = result$xbest[1]
d2.hat = result$xbest[2]
eta.hat = result$xbest[3]
K.hat = result$xbest[4]
n.hat = result$xbest[5]
y.hat = meta %>%
  mutate(tmp = disease_group_detailed) %>%
  mutate(tmp = gsub('Control', 1, tmp)) %>%
  mutate(tmp = gsub('Diabetic', d1.hat, tmp)) %>%
  mutate(tmp = gsub('NPDR [+] DME', 4, tmp)) %>%
  mutate(tmp = gsub('NPDR', d2.hat, tmp)) %>%
  .$tmp %>%
  as.numeric
if (d1.hat < d2.hat){
  top.feats = df.rank[prog == 'progAB'][order(rank)][1:n.hat, ensemblID]
} else {
  top.feats = df.rank[prog == 'progBA'][order(rank)][1:n.hat, ensemblID]
}
x.top = data[, top.feats]
set.seed(946)
spls.fit = caret::train(
  x.top,
  y.hat,
  method = "spls",
  preProcess = c("center", "scale"),
  trControl = fitControl,
  tuneGrid = data.frame(
    eta = eta.hat, K = K.hat, kappa = 0.5))
df.beta = data.table(
  ensemblID = rownames(spls.fit$finalModel$betahat),
  beta = as.numeric(spls.fit$finalModel$betahat))

# stop parallel
stopCluster(cl)

