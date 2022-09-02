# Spatial Distance Sampling model to fit terns counts from SAMM data 


# ----  Load required packages ----

library(ggplot2)
library(sf)
library(tidyverse)
library(here)
library(ggtext)
library(ggspatial)
library(lubridate)


# ---- load required data ----
load("Seabirds/gdlmap.rdata")
load("Seabirds/gdlhex.rdata")
load(here("Seabirds/Caugek/pelmed1721_occ.rdata"))
load(here("Seabirds/Caugek/samm1112_occ.rdata"))
load(here("Seabirds/Caugek/pnm1921_occ.rdata"))
load(here("Seabirds/Caugek/dataocc_telcaugek.Rdata"))

head(y_tel)
head(pelDS)
head(sam_obs)
head(pnm_obs)


# ----fit the GAM for the Caugek DS model ----

# load packages
library(nimble)
library(mgcv)

# make data for GAM and DS

# bind pelmed and telemetry data
# keep only 1/0 data
# remove duplicated cells

# Pelmed data for occupancy 
pelocc <- pelDS %>% 
  mutate( y  = case_when(obstot > 0 ~ 1,
                         obstot == 0 ~0) )
samocc <- sam_obs %>% 
  mutate( y  = case_when(obstot > 0 ~ 1,
                         obstot == 0 ~0) )
pnmocc <- pnm_obs %>% 
  mutate( y  = case_when(obstot > 0 ~ 1,
                         obstot == 0 ~0) )

yb <- bind_rows(pelocc%>% 
                  transmute(depth, depth.sc, id, y) %>% 
                  mutate(data = "pelmed"),
                samocc%>% 
                  transmute(depth, depth.sc, id, y) %>% 
                  mutate(data = "pelmed"),
                pnmocc%>% 
                  transmute(depth, depth.sc, id, y) %>% 
                  mutate(data = "pelmed"),
                y_tel %>% 
                  mutate(data = "telemetry"))%>%
  group_by(id) %>% 
  summarise(yn = max(y))  # long à tourener 4 min

head(yb)

yb <- yb %>% 
  mutate(idb = 1:nrow(yb))


# cherche quel site est echnatillonné par quel programme
dim(pelocc)
idpel <- pelocc %>% 
  st_centroid() %>% 
  st_intersection(yb) %>% 
  select(idb)

idpel$idb

idsam <- samocc %>% 
  st_centroid() %>% 
  st_intersection(yb) %>% 
  select(idb)

idsam$idb

idpnm <- pnmocc %>% 
  st_centroid() %>% 
  st_intersection(yb) %>% 
  select(idb)

idpnm$idb

# idtel <- y_tel %>% 
#   st_centroid() %>% 
#   st_intersection(yb) %>% 
#   select(idb)

dim(y_tel)
dim(yb)
head(idpel)

# n obs
sum(y_tel$y)
sum(pelocc$y)
sum(pnmocc$y)
sum(samocc$y)
sum(yb$yn) # doit être > autres


# set GAM ----
ygami <- yb$yn

stbathyi <- sea.gdl2 %>% 
  mutate(id = 1:nrow(sea.gdl2)) %>% 
  filter(id %in% yb$id) %>% 
  pull(depth.sc)

coordi <- yb %>% 
  st_centroid() %>% 
  st_coordinates() %>% 
  as_tibble() %>% 
  mutate(easting = (X - mean(X)) / sd(X), 
         northing = (Y - mean(Y)) / sd(Y)) %>%
  select(easting, northing) %>%
  as.matrix()

coordXi <- coordi[,1]
coordYi <- coordi[,2]
head(coordYi)
# use jagam 

terngami <- jagam( ygami ~ stbathyi + s(coordXi, coordYi, bs = "gp")  , 
                   family = binomial(link = "logit"),
                   file = "psi.txt") 
length(ygami)
length(terngami$jags.data$y)


# ---- INTEGRATED model w/ telemetry and PELMED ----

im_1 <-  nimbleCode({
  
  # Ecological density
  
  logit(mu[1:nsites]) <- X[1:nsites, 1:34] %*% b[1:34]
  
  for(i in 1:nsites){
    z[i] ~ dbern(mu[i])
  }
  
  # prior for intercept
  b[1] ~ dnorm(0,sd = 100)
  
  # prior for bethy linear effect
  b[2] ~ dnorm(0,sd = 10)
  
  ## prior for s(coordx,coordy) 
  K1[1:32,1:32] <- S1[1:32,1:32] * lambda[1]  + S1[1:32,33:64] * lambda[2]
  b[3:34] ~ dmnorm(zero[3:34], K1[1:32,1:32])
  
  ## smoothing parameter priors
  for (i in 1:2) {
    lambda[i] ~ dgamma(.05,.005)
    rho[i] <- log(lambda[i])
  }
  
  # Observation processes
  
  # PELMED detection process
  # prior
  b0 ~ dnorm(0, 1)
  b1 ~ dnorm(0, 1)
  
  # logit regression
  for(i in 1:246){
    logit(p_pel[i]) <- b0 + b1 * seffpel[i]
  }
  
  # Link Pelmed data and density surface
  for(i in 1:246){
    # ecological density layer
    y_pel[i] ~ dbern(p_pel[i] * z[idpel[i]])
  }
  
  # TELEMETRY observation process
  for(i in 1:ntel){
    y_tel[i] ~ dbern(mu[idtel[i]])
  }
  
})

# --- run the INTEGRATED model ----

data <- list(y_pel = pelocc$y, # pelmed data
             seffpel = c(scale(pelocc$efftot)),
             X = terngami$jags.data$X,
             S1 = terngami$jags.data$S1,
             zero = terngami$jags.data$zero,
             y_tel = y_tel$y)

constants <- list(npel = as.numeric(length(idpel)),
                  idpel =  idpel$idb,
                  idtel = yb$idb,
                  ntel = as.numeric(length(yb$idb)),
                  nsites = as.numeric(length(ygami)))

inits <- list(lambda = terngami$jags.ini$lambda,
              b = terngami$jags.ini$b,
              b0 = rnorm(1,0,1),
              b1 = rnorm(1,0,1),
              z = ygami)


outint <- nimbleMCMC(code = im_1,
                     data = data,
                     inits = inits,
                     constants = constants,
                     monitors = c("b","b0","b1"),
                     niter = 1000,
                     nburnin = 200,
                     nchains = 2)

# full nimble run 
#   Rmodel <- nimbleModel(im_1, constants, data, inits)
#   Rmodel$calculate()
#   
#   conf <- configureMCMC(Rmodel)
#   
#   conf$printMonitors() # see wht parameters are monitored
#   
#   # Custom samplers OPTIONNAL
#   # add Random Walk Block Samples for correlated parameters on sigma log-linear regression, and the IPP
#   conf$printSamplers(byType= TRUE)
#   
#   # Build and compile MCMC
#   Rmcmc <- buildMCMC(conf)
#   Cmodel <- compileNimble(Rmodel)
#   Cmcmc <- compileNimble(Rmcmc, project = Cmodel)
#   # Run 
#   # ART : 1000 it x 2 = 11 sec
#   t <- system.time(samples <- runMCMC(Cmcmc, niter = 10000, nburnin = 100, nchains = 3 , samplesAsCodaMCMC = TRUE))  ## DT: use runMCMC
#   
#   str(samples)
#   
#   mcmcplots::traplot(samples[,"b[13]"])
#   coda::effectiveSize(samples)
# ---- Plot outputs ----

# bind the 2 chains
samplesExti <- as_tibble(rbind(outint$chain1, outint$chain2))
#  samplesExti <- as_tibble(rbind(samples$chain1, samples$chain2))

summary(samplesExti)

grid_coord <- sea.gdl2 %>% 
  st_centroid() %>%
  st_coordinates() %>%
  as_tibble() %>%
  mutate(easting = (X - mean(X))/sd(X),
         northing = (Y - mean(Y))/sd(Y))

grid_bathy <-  c(sea.gdl2$depth.sc)

sm_xy <- smoothCon(s(coordx, coordy, bs = "gp"), 
                   data = data.frame(bathy = grid_bathy,
                                     coordx = grid_coord$easting, 
                                     coordy = grid_coord$northing), 
                   absorb.cons = TRUE) 


Xp_xy <- PredictMat(sm_xy[[1]], data.frame(stbathy = grid_bathy,
                                           coordx = grid_coord$easting, 
                                           coordy = grid_coord$northing))


# matrix nit x 33 
names(samplesExti)
b <- as.matrix(samplesExti[,1:34])
dim(b)


Xp <- cbind(1,median(b[2]), Xp_xy)
dim(Xp)


thetai <- matrix(NA, nrow = nrow(Xp), ncol = nrow(b))

for (i in 1:nrow(b)){
  thetai[1:nrow(Xp), i] <- Xp %*% b[i,]
}

mui <- apply(plogis(thetai),1,median)

range(mui)

# plot depth effect
ggplot()+
  geom_smooth(aes(sea.gdl2$depth,mui), method = "lm") +
  labs(y = "use of space probability",
       x = "seafloor depth")+
  theme_minimal()


# plot a density map
pi <- ggplot()+
  geom_sf(data = gdlmap) + 
  geom_sf(data = sea.gdl2, aes(fill = mui), lwd = 0.1) + 
  scale_fill_gradient(low = "white", high = "#1d3557") +
  theme_minimal() + 
  labs(title = "Utilisation de l'espace par les Sternes caugek",
       subtitle = "Modèle intégré comptages et télémétrie GPS",
       caption = "Source : PELMED2017, 22 GPS",
       fill = "Porbabilité d'utilisation de l'espace",
       x = "", y = "") + 
  theme(plot.title = element_text(face="bold"),
        plot.title.position = "plot",
        legend.position = "top") + 
  annotation_scale(location = "tl", width_hint = 0.2,  pad_y = unit(20, "pt"), pad_x = unit(40, "pt")) +
  guides(fill = guide_colourbar(title.position = 'top', title.hjust = .5,
                                barwidth = unit(10, 'lines'), barheight = unit(.5, 'lines')),
         colour = "none")

pi