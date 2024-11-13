## -----------------------------------------------------------------------------
##
## Script name: BowmanPerry_S-T_model.R
##
## Purpose of script: Recreate the Bowman and Perry (2017) model of savanna-forest transitions (for fun)
##
##        Bowman, D.M.J.S., Perry, G.L.W. (2017). Soil or fire: what causes treeless sedgelands in Tasmanian wet forests?. 
##        Plant Soil 420, 1–18. https://doi.org/10.1007/s11104-017-3386-7
##
## Author: Tom Keeble
##
## Date Created: 2024-10-31
##
## Copyright (c) Tom Keeble, 2024
## Email: tomkeeble0@gmail.com
##
## -----------------------------------------------------------------------------
##
## Notes: Should be self explanatory.
##   
## -----------------------------------------------------------------------------

# 1. PACKAGES

libs <- c(
  "terra"
)

installed_libraries <- libs %in% rownames(installed.packages())

if(any(installed_libraries == F)){
  install.packages(libs[!installed_libraries])
}

invisible(lapply(libs, library, character.only = T))

rm(libs, installed_libraries)

## -----------------------------------------------------------------------------

# 2. THE MODEL

# Step 1: Define the landscape grid
row_size <- 50; col_size = 200 # Define grid dimensions

# Step 2a: Fire ignition with fixed recurrence interval
fire_ignition <- function(landscape, burning) {
  # Randomly pick a cell in the savanna (0 or 3), or if no savanna exists, in forest (1)
  savanna_cells <- rbind(which(landscape == 0, arr.ind = TRUE), which(landscape == 3, arr.ind = TRUE))
  forest_cells <- which(landscape == 1, arr.ind = TRUE)
  
  if (length(savanna_cells) > 0) {
    ignition_cell <- savanna_cells[sample(1:nrow(savanna_cells), 1), ]
  } else if (length(forest_cells) > 0) {
    ignition_cell <- forest_cells[sample(1:nrow(forest_cells), 1), ]
  }
  
  burning[ignition_cell[1], ignition_cell[2]] <- 1  # Mark the selected cell as burning
  
  return(burning)
}

# Step 2b: Percolation fire spread
fire_spread <- function(landscape, burning, fire_probability_forest=0.035, fire_probability_savanna=0.3, plotting = T) {
  nrows <- nrow(landscape)
  ncols <- ncol(landscape)
  
  # Find all currently burning cells
  burning_cells <- which(burning == 1, arr.ind = TRUE)
  
  # List to store newly burning cells
  new_burning <- burning_cells
  
  # While there are new burning cells, continue to propagate fire
  while (length(new_burning) > 0) {
    new_burning_cells <- numeric() # Track newly ignited cells
    
    for (b in 1:nrow(new_burning)) {
      i <- new_burning[b, 1]
      j <- new_burning[b, 2]
      
      # Define neighbouring 8 cells
      neighbours <- list(c(i - 1, j - 1), c(i - 1, j), c(i - 1, j + 1),
                          c(i, j - 1), c(i, j + 1),
                          c(i + 1, j - 1), c(i + 1, j), c(i + 1, j + 1))
      
      for (n in neighbours) {
        ni <- n[1]
        nj <- n[2]
        
        if (ni > 0 && ni <= nrows && nj > 0 && nj <= ncols) {  # Ensure we don't assess cells out of bounds
          if (landscape[ni, nj] %in% c(0, 3) && burning[ni, nj] == 0 && runif(1) < fire_probability_savanna) {
            burning[ni, nj] <- 1  # Savanna/colonised savanna catches fire
            new_burning_cells <- rbind(new_burning_cells, c(ni, nj)) # Add to list of newly burning cells
          } else if (landscape[ni, nj] == 1 && burning[ni, nj] == 0 && runif(1) < fire_probability_forest) {
            burning[ni, nj] <- 1  # Forest catches fire
            new_burning_cells <- rbind(new_burning_cells, c(ni, nj)) # Add to list of newly burning cells
          }
        }
      }
    }
    
    # Update the burning cells for the next round of propagation
    new_burning <- new_burning_cells
  }
  
  if(plotting){
    plot(rast(burning), col = c("white", "red"))
    Sys.sleep(0.05)
  }
  return(burning)
}

# Step 3: Implement soil fertility feedbacks
update_soil_fertility <- function(soil_fertility, landscape, burning, fire_impact=0.2, recovery_rate=0.001) {
  new_soil_fertility <- soil_fertility
  for (i in 1:nrow(landscape)) {
    for (j in 1:ncol(landscape)) {
      if (burning[i, j] == 1 && landscape[i, j] %in% c(0, 3)) {
        new_soil_fertility[i, j] <- new_soil_fertility[i, j] + fire_impact # Decrease fertility due to fire in savanna
      } else if (landscape[i, j] == 1 && burning[i, j] == 0) {
        new_soil_fertility[i, j] <- max(1, new_soil_fertility[i, j] - recovery_rate) # Slow recovery in unburned forest
      }
    }
  }
  return(new_soil_fertility)
}


# Step 4: Update the landscape based on cells that burnt
update_landscape <- function(landscape, burning) {
  # Change any landscape forest (1) or colonised savanna (3) cells that were burnt to savanna
  for (i in 1:nrow(landscape)) {
    for (j in 1:ncol(landscape)) {
      if (burning[i, j] == 1 && landscape[i, j] %in% c(1, 3)) {
        landscape[i, j] <- 0  # Convert to savanna
      }
    }
  }
  return(landscape)
}

# Step 5: Forest expansion and regrowth after fire-free period
forest_expansion <- function(landscape, colonisation_time, soil_fertility, base_fire_recovery_time=15, dispersal_rate=1) {
  
  # Transition colonised savanna (3) to forest (1) after the colonisation period, adjusted by soil infertility
  colonised_cells <- which(landscape == 3, arr.ind = TRUE)
  
  if (length(colonised_cells) != 0) {
    for (c in 1:nrow(colonised_cells)) {
      ci <- colonised_cells[c, 1]
      cj <- colonised_cells[c, 2]
      
      # Increment colonisation time
      colonisation_time[ci, cj] <- colonisation_time[ci, cj] + 1
      
      # Calculate adjusted colonisation time based on soil fertility
      adjusted_colonisation_time <- base_fire_recovery_time * soil_fertility[ci, cj]
      
      # Convert to forest if the colonisation period has passed
      if (colonisation_time[ci, cj] >= adjusted_colonisation_time) {
        landscape[ci, cj] <- 1  # Transition to forest
      }
    }
  }
  
  
  nrows <- nrow(landscape)
  ncols <- ncol(landscape)
  
  # Propagule dispersal: each forest cell (1) disperses one propagule per year
  forest_cells <- which(landscape == 1, arr.ind = TRUE)
  
  if (length(forest_cells) != 0) {
    for (f in 1:nrow(forest_cells)) {
      i <- forest_cells[f, 1]
      j <- forest_cells[f, 2]
      
      # Dispersal: distance is a random deviate from a negative exponential distribution (mean=1)
      dispersal_distance <- rpois(1, dispersal_rate)
      angle <- runif(1, 0, 2 * pi)  # Random direction
      
      # Calculate the new cell based on dispersal distance and angle
      ni <- round(i + dispersal_distance * cos(angle))
      nj <- round(j + dispersal_distance * sin(angle))
      
      # Ensure the new cell is within bounds
      if (ni > 0 && ni <= nrows && nj > 0 && nj <= ncols) {
        # If the propagule lands in a savanna cell, colonise it (convert to state 3)
        if (landscape[ni, nj] == 0) {
          landscape[ni, nj] <- 3  # Convert to colonized savanna
          colonisation_time[ni, nj] <- 0  # Reset colonization time
        }
      }
    }
  }
  
  return(list(landscape = landscape, colonisation_time = colonisation_time))
}


# Step 6: Run the simulation
run_simulation <- function(n_steps=2500, recurrence_interval=15, base_fire_recovery_time=15, dispersal_rate=1, fire_soil_feedback = T, edaphic_boundary = T, plotting = T) {
  # Initial conditions: Create a sharp transition between forest and savanna
  landscape <- matrix(0, nrow = row_size, ncol = col_size) # 0: savanna, 1: forest, 3: colonised savanna
  landscape[ , 1:(col_size / 2)] <- 0  # Savanna on the left side
  landscape[ , (col_size / 2 + 1):col_size] <- 1  # Forest on the right side
  
  time_since_last_fire <- matrix(0, nrow = row_size, ncol = col_size) # Time since last fire
  colonisation_time <- matrix(0, nrow = row_size, ncol = col_size)  # Time since propagule colonisation
  
  if(edaphic_boundary){
    # Fertility of 5 in the savanna, 1 in forest
    soil_fertility <- landscape
    soil_fertility[ , 1:(col_size / 2)] <- 5  # Savanna cells set to 5
  } else {
    soil_fertility <- matrix(1, nrow = row_size, ncol = col_size)  # Fertility starts at 1 across the domain
  }
  
  burning <- matrix(0, nrow = row_size, ncol = col_size) # Binary matrix of burning cells (1: burning, 0: not burning)
  
  for (step in 1:n_steps) {
    
    # Fire ignition every `recurrence_interval` steps
    if (step %% recurrence_interval == 0) { # Make stochastic about the mean of 15 years <------------------------------------- TO DO!
      
      
      burning <- fire_ignition(landscape, burning)
      
      # Fire spread using percolation process
      burning <- fire_spread(landscape, burning, plotting = plotting)
    }
    
    if(fire_soil_feedback){
      # Update soil fertility based on fire activity
      soil_fertility <- update_soil_fertility(soil_fertility, landscape, burning)
    }
    
    # Update landscape based on cells that burnt
    landscape <- update_landscape(landscape, burning)
    
    # Forest expansion and regrowth after fire-free period, with soil infertility effects
    expansion_results <- forest_expansion(landscape, colonisation_time, soil_fertility)
    landscape <- expansion_results$landscape
    colonisation_time <- expansion_results$colonisation_time
    
    if (step %% recurrence_interval == 0) { # Make stochastic about the mean of 15 years <------------------------------------- TO DO!
      # Reset burning cells
      burning <- matrix(0, nrow = row_size, ncol = col_size)
    }
    
    # Optional: Visualise results at each step
    if(plotting){
      plot(rast(landscape), col = data.frame(value = c(0, 1, 3), color = c("brown", "green", "yellow"))$color)
      Sys.sleep(0.05)
    }
    print(step)
  }
  
  return(list(final_landscape=landscape, soil_fertility=soil_fertility))
}

# Run the simulation
results <- run_simulation(fire_soil_feedback = F, edaphic_boundary = F, plotting = T)

# Visualize the final landscape
plot(raster(results$final_landscape))
