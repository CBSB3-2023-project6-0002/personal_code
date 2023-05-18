# CBSB3 Mini Project 6: Modelling Gene Regulatory Circuitry

#### Import ####
library(tidyverse)
library(deSolve)
library(data.table)

#### Activation & Repression Functions ####
activation_regulator <- function(u, Ku, v=0, Kv=1, H=2){
  return((u/Ku)^H / (1 + (u/Ku)^H + (v/Kv)^H))
}

repression_regulator <- function(u, Ku, v=0, Kv=1, H=2){
  return(1 / (1 + (u/Ku)^H + (v/Kv)^H))
}

#### And & Or Gates ####
and_gate <- function(X, Kx, X_sign, Y, Ky, Y_sign, H=2){
  # get X effect
  if(X_sign == 1){
    fx <- activation_regulator(X, Kx, H=H)
  }else{
    fx <- repression_regulator(X, Kx, H=H)
  }
  # get Y effect
  if(Y_sign == 1){
    fy <- activation_regulator(Y, Ky, H=H)
  }else{
    fy <- repression_regulator(Y, Ky, H=H)
  }
  return(fx * fy)
}

or_gate <- function(X, Kx, X_sign, Y, Ky, Y_sign, H=2){
  # get X effect
  if(X_sign == 1){
    fx <- activation_regulator(X, Kx, Y, Ky, H=H)
  }else{
    fx <- repression_regulator(X, Kx, Y, Ky, H=H)
  }
  # get Y effect
  if(Y_sign == 1){
    fy <- activation_regulator(Y, Ky, X, Kx, H=H)
  }else{
    fy <- repression_regulator(Y, Ky, X, Kx, H=H)
  }
  return(fx + fy)
}

#### Generalized Feed-forward Motif Model ####
ffl_model <- function(t, vars, params){
  # extract variables
  X <- vars[1]
  Y <- vars[2]
  Z <- vars[3]
  # extract params
  By=params["By"]
  Bz=params["Bz"]
  alphay=params["alphay"]
  alphaz=params["alphaz"]
  betay=params["betay"]
  betaz=params["betaz"]
  Kxy=params["Kxy"]
  Kxz=params["Kxz"]
  Kyz=params["Kyz"]
  H=params["H"]
  sign.xy=params["sign.xy"]
  sign.xz=params["sign.xz"]
  sign.yz=params["sign.yz"]
  and.or=params["and.or"]
  
  # dYdt
  if(sign.xy == 1){
    y.regulation <- activation_regulator(X, Kxy, H=H)
  }else{
    y.regulation <- repression_regulator(X, Kxy, H=H)
  }
  dYdt <- By + betay * y.regulation - alphay * Y
  
  # dZdt
  if(and.or == 1){
    z.regulation <- and_gate(X, Kxz, sign.xz, Y, Kyz, sign.yz, H)
  }else{
    z.regulation <- or_gate(X, Kxz, sign.xz, Y, Kyz, sign.yz, H)
  }
  dZdt <- Bz + betaz * z.regulation - alphaz * Z
  
  return(list(c(0, dYdt, dZdt)))
}

####Simulation & Plot ####
simulation_model <- function(model, vars, times, params, X.times, X.vals, 
                             noise.mean=0, noise.sd=NA,
                             mode='df', pheno=NA, noise.thres=0){
  # X input events
  events <- data.frame(var="X", time=X.times, value=X.vals,
                       method="add")
  # phenotypic noise
  if(!is.na(noise.sd)){
    X.noise <- data.frame(var="X", time=times, 
                          value=rnorm(length(times), noise.mean, noise.sd),
                          method="add")
    Y.noise <- data.frame(var="Y", time=times, 
                          value=rnorm(length(times), noise.mean, noise.sd),
                          method="add")
    Z.noise <- data.frame(var="Z", time=times, 
                          value=rnorm(length(times), noise.mean, noise.sd),
                          method="add")
    events <- rbind(events, X.noise, Y.noise, Z.noise) %>% arrange(time)
  }
  
  # simulation
  sim.res <- ode(
    func=model,
    y=vars,
    times=times,
    parms=params,
    events=list(data = events)
  ) %>% as.data.frame()
  # remove negative points
  sim.res[sim.res < 0] <- 0
  
  # if we filter by noise sd:
  if(noise.thres > 0){
    true.noise.sd <- 0
    while(true.noise.sd < noise.thres){
      # X input events
      events <- data.frame(var="X", time=X.times, value=X.vals,
                           method="add")
      X.noise <- data.frame(var="X", time=times, 
                            value=rnorm(length(times), noise.mean, noise.sd),
                            method="add")
      Y.noise <- data.frame(var="Y", time=times, 
                            value=rnorm(length(times), noise.mean, noise.sd),
                            method="add")
      Z.noise <- data.frame(var="Z", time=times, 
                            value=rnorm(length(times), noise.mean, noise.sd),
                            method="add")
      events <- rbind(events, X.noise, Y.noise, Z.noise) %>% arrange(time)
      # simulation
      sim.res <- ode(
        func=model,
        y=vars,
        times=times,
        parms=params,
        events=list(data = events)
      ) %>% as.data.frame()
      # remove negative points
      sim.res[sim.res < 0] <- 0
      
      true.noise.sd <- sd(sim.res$Z[sim.res$time < X.times[1]])
      #print(true.noise.sd)
    }
  }
  
  if(mode=='df'){
    return(sim.res)
  }else if(mode=='pheno'){
    # get baseline
    baseline.rows <- sim.res %>%
      filter(time < X.times[1])
    Z.mean <- mean(baseline.rows$Z)
    Z.sd <- sd(baseline.rows$Z)
    # get after onset
    onset.rows <- sim.res %>%
      filter(time >= X.times[1])
    
    if(pheno=='response.time.on'){
      sig.rows <- onset.rows %>%
        filter(Z > Z.mean + 3*Z.sd)
      res <- sig.rows$time[1] - X.times[1]
      
    }else if(pheno=='pulse'){
      sig.rows <- sim.res %>% filter(Z > Z.mean + 3*Z.sd)
      res <- max(sig.rows$time) - min(sig.rows$time)
      
    }else if(pheno=='response.time.off'){
      sig.rows <- onset.rows %>%
        filter((Z < Z.mean + 3*Z.sd) & time >= X.times[2])
      res <- sig.rows$time[1] - X.times[2]
      
    }else if(pheno=='FCD'){
      first.fcd <- sim.res %>% filter(time > 5 & time < 10)
      fcd1 <- max(first.fcd$Z)
      second.fcd <- sim.res %>% filter(time > 10 & time < 15)
      fcd2 <- max(second.fcd$Z)
      res <- fcd1 / fcd2
      
    }else if(pheno=='persistence'){
      onset.rows <- sim.res %>% filter(time > 1 & time < 2)
      res <- mean(baseline.rows$Z) - mean(onset.rows$Z)
    }
    return(res)
  }
}

plot_time_course <- function(sim.res, x.time=1){
  # get baseline
  baseline.rows <- sim.res %>%
    filter(time < x.time)
  Z.mean <- mean(baseline.rows$Z)
  Z.sd <- sd(baseline.rows$Z)
  # plot
  ggplot(sim.res %>% gather(variable, value, -time)) +
    geom_line(aes(x=time, y=value)) +
    geom_hline(data=sim.res %>% gather(variable, value, -time) %>% filter(variable=="Z"),
               linetype="dashed", aes(yintercept = Z.mean + 3*Z.sd), color="red") +
    theme_classic() +
    ylab("Concentration (arbitrary unit)") +
    xlab("Time (s)") +
    facet_wrap(~variable, ncol=1)
}

#### Test Cases ####

## CFFL1

# AND gate
default.vars <- c(
  X=0,
  Y=0,
  Z=0
)
global.params.cffl1.and <- c(
  # reaction related
  By=0.0,
  Bz=0.0,
  alphay=1,
  alphaz=1,
  betay=1,
  betaz=1,
  # regulation related
  Kxy=0.5,
  Kxz=0.5,
  Kyz=0.5,
  H=2
)
global.params.cffl1.and.ctrl <- c(
  # reaction related
  By=0.0,
  Bz=0.0,
  alphay=1,
  alphaz=1,
  betay=1,
  betaz=1,
  # regulation related
  Kxy=0.5,
  Kxz=0.5,
  Kyz=0.00001,
  H=2
)
global.params.cffl1.and.mut <- c(
  # reaction related
  By=0.0,
  Bz=0.0,
  alphay=1,
  alphaz=1,
  betay=1,
  betaz=1,
  # regulation related
  Kxy=0.5,
  Kxz=0.5,
  Kyz=0.25,
  H=2
)
cffl1.and.params <- c(
  sign.xy=1,
  sign.xz=1,
  sign.yz=1,
  and.or=1
)

# facet wrap plot for mutant vs WT
combine_mut_plot <- function(df, dfmut, x.time=1){
  # get baselines
  baseline.rows <- df %>%
    filter(time < x.time)
  Z.mean <- mean(baseline.rows$Z)
  Z.sd <- sd(baseline.rows$Z)
  
  baseline.rows.mut <- dfmut %>%
    filter(time < x.time)
  Z.mean.mut <- mean(baseline.rows.mut$Z)
  Z.sd.mut <- sd(baseline.rows.mut$Z)
  
  # name rows
  df$name <- "ctrl"
  dfmut$name <- "mut"
  sim.res <- rbind(df, dfmut)
  
  # plot
  ggplot(sim.res %>% gather(variable, value, -c(time, name))) +
    geom_line(aes(x=time, y=value, color=name)) +
    geom_hline(data=sim.res %>% gather(variable, value, -time) %>% filter(variable=="Z"),
               linetype="dashed", aes(yintercept = Z.mean + 3*Z.sd), color="#050505") +
    geom_hline(data=sim.res %>% gather(variable, value, -time) %>% filter(variable=="Z"),
               linetype="dashed", aes(yintercept = Z.mean.mut + 3*Z.sd.mut), color="#fc6100") +
    theme_classic() +
    ylab("Concentration (arbitrary unit)") +
    scale_color_manual(values = c("#050505", "#fc6100")) +
    xlab("Time (s)") +
    facet_wrap(~variable, ncol=1)
}


# persistent X
cffl1.and.pers <- simulation_model(ffl_model, default.vars, 
                                   seq(0, 5, 0.005), c(global.params.cffl1.and, cffl1.and.params),
                                   c(1, 3), c(2, -2), noise.sd = 0.002, noise.thres = 0.005)
plot_time_course(cffl1.and.pers)

# mut
cffl1.and.pers.mut <- simulation_model(ffl_model, default.vars, 
                                       seq(0, 5, 0.005), c(global.params.cffl1.and.mut, cffl1.and.params),
                                       c(1, 3), c(2, -2), noise.sd = 0.002, noise.thres = 0.005)
plot_time_course(cffl1.and.pers.mut)

# put together
combine_mut_plot(cffl1.and.pers, cffl1.and.pers.mut)



# persistent X ctrl simple regulation
cffl1.and.pers.ctrl <- simulation_model(ffl_model, default.vars, 
                                        seq(0, 5, 0.005), c(global.params.cffl1.and.ctrl, 
                                                            cffl1.and.params),
                                        c(1, 3), c(2, -2), noise.sd = 0.002, noise.thres = 0.005)
plot_time_course(cffl1.and.pers.ctrl)

# transient X
cffl1.and.trans <- simulation_model(ffl_model, default.vars, 
                                    seq(0, 5, 0.005), c(global.params.cffl1.and, cffl1.and.params),
                                    c(1, 1.41), c(2, -2), noise.sd = 0.002, noise.thres = 0.005)
plot_time_course(cffl1.and.trans)

cffl1.and.trans.mut <- simulation_model(ffl_model, default.vars, 
                                        seq(0, 5, 0.005), c(global.params.cffl1.and.mut, cffl1.and.params),
                                        c(1, 1.12), c(2, -2), noise.sd = 0.002, noise.thres = 0.005)
plot_time_course(cffl1.and.trans.mut)

# plot together
combine_mut_plot(cffl1.and.trans, cffl1.and.trans.mut)



## IFFL1

# pulse
global.params.pulse <- c(
  # reaction related
  By=0,
  Bz=0,
  alphay=1,
  alphaz=2,
  betay=1,
  betaz=2,
  # regulation related
  Kxy=0.5,
  Kxz=0.5,
  Kyz=0.5,
  H=2
)
global.params.pulse.noy <- c(
  # reaction related
  By=0,
  Bz=0,
  alphay=1,
  alphaz=2,
  betay=0,
  betaz=2,
  # regulation related
  Kxy=0.5,
  Kxz=0.5,
  Kyz=0.5,
  H=2
)
global.params.pulse.mut <- c(
  # reaction related
  By=0,
  Bz=0,
  alphay=1,
  alphaz=2,
  betay=1,
  betaz=2,
  # regulation related
  Kxy=0.5,
  Kxz=0.5,
  Kyz=1,
  H=2
)
iffl1.and.params <- c(
  sign.xy=1,
  sign.xz=1,
  sign.yz=-1,
  and.or=1
)

# WT
iffl1.pulse <- simulation_model(ffl_model, default.vars, 
                                seq(0, 5, 0.005), c(global.params.pulse, iffl1.and.params),
                                c(1, 3), c(2, -2), noise.sd = 0.002, noise.thres = 0.005)
plot_time_course(iffl1.pulse)

# mutant
iffl1.pulse.mut <- simulation_model(ffl_model, default.vars, 
                                    seq(0, 5, 0.005), c(global.params.pulse.mut, iffl1.and.params),
                                    c(1, 3), c(2, -2), noise.sd = 0.002, noise.thres = 0.005)
plot_time_course(iffl1.pulse.mut)

# plot together
combine_mut_plot(iffl1.pulse, iffl1.pulse.mut)

# simple regulation
iffl1.pulse.ctrl <- simulation_model(ffl_model, default.vars, 
                                seq(0, 5, 0.005), 
                                c(global.params.pulse.noy, iffl1.and.params),
                                c(1, 3), c(2, -2), 
                                noise.sd = 0.002, noise.thres = 0.005)
plot_time_course(iffl1.pulse.ctrl)



# FCD
ffl_fcd_model <- function(t, vars, params){
  # extract variables
  X <- vars[1]
  Y <- vars[2]
  Z <- vars[3]
  # extract params
  By=params["By"]
  Bz=params["Bz"]
  alphay=params["alphay"]
  alphaz=params["alphaz"]
  betay=params["betay"]
  betaz=params["betaz"]
  Kxy=params["Kxy"]
  Kxz=params["Kxz"]
  Kyz=params["Kyz"]
  H=params["H"]
  
  # dYdt
  dYdt <- By + betay * X - alphay * Y
  
  # dZdt
  dZdt <- Bz + betaz * X / Y - alphaz * Z
  
  return(list(c(0, dYdt, dZdt)))
}

default.vars.fcd <- c(
  X=2,
  Y=2,
  Z=5
)
global.params.fcd <- c(
  # reaction related
  By=0,
  Bz=0,
  alphay=1,
  alphaz=1,
  betay=1,
  betaz=5,
  # # regulation related
  # Kxy=0.5,
  # Kxz=0.5,
  # Kyz=0.5,
  H=2
)
iffl1.fcd <- simulation_model(ffl_fcd_model, default.vars.fcd, 
                              seq(0, 20, 0.01), c(global.params.fcd, iffl1.and.params),
                              c(5, 10), c(4, 8), noise.sd = 0.005)
plot_time_course(iffl1.fcd)

mean.fcd <- mean_time_course(ffl_fcd_model, default.vars.fcd, 
                             seq(0, 20, 0.01), c(global.params.fcd, iffl1.and.params),
                             c(5, 10), c(4, 8), noise.sd = 0.005,
                             n=10)
plot_time_course(mean.fcd)


#### Robustness Related####

## Add Noise to time course sim (done built-in)

## Define Measurable Phenotype from SE of baseline level

## generalized phenotype simulation function
batch_simulation_pheno <- function(model, default.vars, times, global.params, model.params,
                                   X.times, X.vals,
                                   pheno, n=1000, noise.sd=0.002, noise.thres=0.01){
  phenos <- c()
  # start simulation
  for(i in 1:n){
    res <- simulation_model(model, default.vars, 
                            times, c(global.params, model.params),
                            X.times, X.vals, noise.sd = noise.sd,
                            mode = 'pheno', pheno = pheno, noise.thres=noise.thres)
    phenos <- c(phenos, res)
    print(i)
  }
  return(phenos)
}

## delay
and.response.time.phenos <- batch_simulation_pheno(ffl_model, default.vars,
                                                   seq(0, 5, 0.005),
                                                   global.params.cffl1.and,
                                                   cffl1.and.params,
                                                   c(1, 2), c(2, -2),
                                                   'response.time.on',
                                                   n=100, noise.sd=0.002,
                                                   noise.thres=0.005)
hist(and.response.time.phenos)
mean(and.response.time.phenos)

# ctrl
and.response.time.phenos.ctrl <- batch_simulation_pheno(ffl_model, default.vars,
                                                        seq(0, 5, 0.005),
                                                        global.params.cffl1.and.ctrl,
                                                        cffl1.and.params,
                                                        c(1, 2), c(2, -2),
                                                        'response.time.on',
                                                        n=100, noise.sd=0.002,
                                                        noise.thres=0.005)
hist(and.response.time.phenos.ctrl)


# before adding mutations, first quantify phenotype
mean(and.response.time.phenos)
t.test(and.response.time.phenos, and.response.time.phenos.ctrl)

# mutate
and.response.time.phenos.mut <- batch_simulation_pheno(ffl_model, default.vars,
                                                       seq(0, 5, 0.005),
                                                       global.params.cffl1.and.mut,
                                                       cffl1.and.params,
                                                       c(1, 2), c(2, -2),
                                                       'response.time.on',
                                                       n=100, noise.sd=0.002,
                                                       noise.thres=0.005)

t.test(and.response.time.phenos, and.response.time.phenos.mut)

## persistence using binary search
find_persistence <- function(mode="sequential", interval=seq(0.1,1,0.1), n=10, precision=0.01,
                             start=0, end=1, params=global.params.cffl1.and){
  if(mode=="sequential"){
    for(int in interval){
      diff <- tryCatch(
        expr = {
          batch_simulation_pheno(ffl_model, default.vars,
                                 seq(0, 5, 0.005),
                                 global.params.cffl1.and,
                                 cffl1.and.params,
                                 c(1, 1 + int), c(2, -2),
                                 'persistence',
                                 n=n, noise.sd=0.002,
                                 noise.thres=0.005)
        },
        error = function(e){
          batch_simulation_pheno(ffl_model, default.vars,
                                 seq(0, 5, 0.005),
                                 global.params.cffl1.and,
                                 cffl1.and.params,
                                 c(1, 1 + int + 0.001), c(2, -2),
                                 'persistence',
                                 n=n, noise.sd=0.002,
                                 noise.thres=0.005)
        }
      )
      
      ttest.res <- t.test(diff, mu=0)
      if(ttest.res$p.value <= 0.05){
        return(int)
      }
      print(1+int)
    }
  }else if(mode=="binary"){
    while(end - start > precision){
      curr <- (end + start)/2
      
      print(paste("testing:", curr))
      
      curr.diff <- tryCatch(
        expr = {
          batch_simulation_pheno(ffl_model, default.vars,
                                 seq(0, 5, 0.005),
                                 params,
                                 cffl1.and.params,
                                 c(1, 1 + curr), c(2, -2),
                                 'persistence',
                                 n=n, noise.sd=0.002,
                                 noise.thres=0.005)
        },
        error = function(e){
          batch_simulation_pheno(ffl_model, default.vars,
                                 seq(0, 5, 0.005),
                                 params,
                                 cffl1.and.params,
                                 c(1, 1 + curr + 0.001), c(2, -2),
                                 'persistence',
                                 n=n, noise.sd=0.002,
                                 noise.thres=0.005)
        }
      )
      # if sig, move end to curr, else, move start to curr
      if(t.test(curr.diff, mu=0, alternative="less")$p.value <= 0.05){
        end <- curr
        print(paste("moving end to:", curr))
      }else{
        start <- curr
        print(paste("moving start to:", curr))
      }
      
    }
    return(curr)
  }
  
  return(NA)
}

pers.base.onset.diff <- batch_simulation_pheno(ffl_model, default.vars,
                                               seq(0, 5, 0.005),
                                               global.params.cffl1.and,
                                               cffl1.and.params,
                                               c(1, 1.33), c(2, -2),
                                               'persistence',
                                               n=5, noise.sd=0.002,
                                               noise.thres=0.005)


hist(pers.base.onset.diff)
t.test(pers.base.onset.diff, mu=0)


# find persistence
find_persistence(mode="sequential", seq(0.05,0.2,0.01), n=5)

find_persistence(mode="binary", n=10)

find_persistence(mode="binary", n=10, params = global.params.cffl1.and.mut)

find_persistence(mode="binary", n=10, params = global.params.cffl1.and.ctrl, precision = 0.001)

## pulse
pulse.phenos <- batch_simulation_pheno(ffl_model, default.vars,
                                       seq(0, 5, 0.005),
                                       global.params.pulse,
                                       iffl1.and.params,
                                       c(1, 3), c(2, -2),
                                       'pulse',
                                       n=100, noise.sd=0.002,
                                       noise.thres=0.005)
hist(pulse.phenos)

pulse.ctrl.phenos <- batch_simulation_pheno(ffl_model, default.vars,
                                            seq(0, 5, 0.005),
                                            global.params.pulse.noy,
                                            iffl1.and.params,
                                            c(1, 3), c(2, -2),
                                            'pulse',
                                            n=100, noise.sd=0.002,
                                            noise.thres=0.005)
hist(pulse.ctrl.phenos)

t.test(pulse.phenos, pulse.ctrl.phenos)

pulse.phenos.mut <- batch_simulation_pheno(ffl_model, default.vars,
                                           seq(0, 5, 0.005),
                                           global.params.pulse.mut,
                                           iffl1.and.params,
                                           c(1, 3), c(2, -2),
                                           'pulse',
                                           n=100, noise.sd=0.002,
                                           noise.thres=0.005)
hist(pulse.phenos.mut)

t.test(pulse.phenos, pulse.phenos.mut)


## FCD
example.fcd <- batch_simulation_pheno(ffl_fcd_model, default.vars.fcd,
                                      seq(0, 20, 0.01),
                                      global.params.cffl1.and,
                                      cffl1.and.params,
                                      c(5, 10), c(4, 8),
                                      'FCD',
                                      n=100, noise.sd=0.002, noise.thres = 0)
hist(example.fcd)
t.test(example.fcd, mu=1)

