# flux 4D

# SIMULATE ECOLOGICAL PARAMETERS
# set simulation parameters : 
nsites <- 50
nalti <- 3
nday <- 10

# simulation of flux data following a random Poisson process
# w/ great variability around 150 bird/space-time unit
mu <- array(rnorm(nsites*nalti*nday, 150, sd= 30), 
            dim = c(nsites, nalti, nday))

# ordrer the simulate flux to create a site effect
for(i in 1:nalti){
  for(j in 1:nday){
  mu[,i,j] <- sort(mu[,i,j])
}}

# SIMUALTE DATA
# Costal radar data : 2 sites, all alti, all time
rmu <- mu
rmu[] <- rpois(nsites*nalti*nday, mu)
rmu[3:nsites,,] <- NA

# Sea radar data : all sites, all alti, 2 timesteps
mermu <- mu
mermu[] <- rpois(nsites*nalti*nday, mu)
mermu[,,1:8] <- NA

# Weather radar data : all sites, all timestep, 0 alti
metmu <- matrix(round(apply(mu,c(1,3),sum ),0), nrow = nsites, ncol = nday)

# simulate environmental covariate depending on site
env <- sort(rnorm(nsites,0,1))
env
# simulate weather covariate depending on site and time w/ no effect
meteo <- array(rnorm(nsites*nday,0,1), dim = c(nsites,nday))

# FIT THE MODEL

# NIMBLE code
f4 <- nimbleCode({
  
  # prior
  a0 ~ dnorm(0,1)
  a1 ~ dnorm(0,1)
  a2 ~ dnorm(0,1)
  
  # ecological process
  for(i in 1:nsites){
    log(mu[i,1:nalti,1:nday]) <-  a0 + a1 * env[i]  + a2 * meteo[i,1:nday]
  }
  
  #radar coast
  for(xy in 1:2){
    for(z in 1:nalti){
      for(t in 1:nday){
        rmu[xy,z,t] ~ dpois(mu[xy,z,t])
      }
    }
  }
  # sea radar
  for(xy in 1:nsites){
    for(z in 1:nalti){
      for(t in 9:nday){
        mermu[xy,z,t] ~ dpois(mu[xy,z,t])
      }
    }
  }
  
  
  # weather radar
  for(xy in 1:nsites){
      for(t in 1:nday){
        summu[xy,t] <- sum(mu[xy,1:nalti,t])
        metmu[xy,t] ~ dpois(summu[xy,t])
      }
    }
  
})

# BUNDLE MODEL 
data = list(rmu = rmu,
            mermu = mermu,
            metmu = metmu,
            env = env,
            meteo = meteo)

constants = list(nsites = nsites,
                 nalti = nalti,
                 nday = nday)

inits <- list(a2 = rnorm(1,0,1),
              a0 = rnorm(1,0,1),
              a1 = rnorm(1,0,1))


# Rmodel <- nimbleModel(f4, constants, data, inits)
# Rmodel$initializeInfo()
# Rmodel$calculate()

# RUN MODEL
 outint <- nimbleMCMC(code = f4,
                      data = data,
                      inits = inits,
                      constants = constants,
                      monitors = c("a2","a0","a1"),
                      niter = 1000,
                      nburnin = 20,
                      nchains = 2)

# OUTPUTS
samplesExti <- as_tibble(rbind(outint$chain1, outint$chain2))


mcmcplots::denplot(samplesExti)

summary(samplesExti)
# a1 significant --> postitive effect of env
# a2 non-significant --> no effect of weather

a0.post <- samplesExti$a0
a1.post <- samplesExti$a1
a2.post <- samplesExti$a2

mu.post <- exp(median(a0.post) + median(a1.post) * env + median(a2.post) * meteo)


# PLOT OUTPUTS
dim(mu.post) # nsites * nday

# true value of mu : sites * day
mu.true <- apply(mu,c(1,3),sum)

# plot flux per day
dmu.post <- apply(mu.post,2,sum)
dmu.true <- apply(mu.true,2,sum)
ggplot() + 
  geom_jitter(aes(x = 1:10, y = dmu.post, color ="Estimated")) +
  geom_jitter(aes(x = 1:10, y = dmu.true, color ="True")) +
  labs(title = "Flux per day",
       subtitle = "estimated vs. true",
       y = "Flux", x = "time steps",
       colour = "") + 
  theme(plot.title = element_text(face="bold"),
        plot.title.position = "plot",
        legend.position = "top")+
  theme_minimal()

# plot flux per sites
smu.post <- apply(mu.post,1,sum)
smu.true <- apply(mu.true,1,sum)
ggplot() + 
  geom_jitter(aes(x = 1:50, y = smu.post, color ="Estimated")) +
  geom_jitter(aes(x = 1:50, y = smu.true, color ="True")) +
  labs(title = "Flux per sites",
       subtitle = "estimated vs. true",
       y = "Flux", x = "time steps",
       colour = "") + 
  theme(plot.title = element_text(face="bold"),
        plot.title.position = "plot",
        legend.position = "top")+
  theme_minimal()


dim(mu.true)

# plot
ggplot() + 
  geom_jitter(aes(x = mu.post, y = mu.true)) + 
  geom_abline(intercept = 0, slope = 1) + 
  xlim(0 , max(mu.true)) + 
  ylim(0, max(mu.true)) +
  labs(title = "Scatterplot of flux w/ identity line",
       subtitle = "estimated vs. true",
       y = "True", x = "Estimated",
       colour = "") + 
  theme(plot.title = element_text(face="bold"),
        plot.title.position = "plot",
        legend.position = "top")+
  theme_minimal()
