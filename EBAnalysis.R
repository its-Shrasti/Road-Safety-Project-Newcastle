#--------------------STEP 1---------------

#Fit the nb model and extract coefficients

# Load necessary package

library(MASS)

ref_data <- read.csv("control.csv")

# Fit a Negative Binomial regression model

nb_model <- glm.nb(Slight+Serious+ Fatal~ Sp.Lim +Av.Sp+ ef.p+ ov.lim+ fifteen.ov.lim+  Flow+ Class +Type, data = ref_data)


# View summary of the model
summary(nb_model)

# Extract all coefficients
coeff <- coef(nb_model)

# Extract the dispersion parameter (theta)
dispersion_param <- nb_model$theta

shape_param<-(1/dispersion_param)

#--------------------STEP 2------------

#cal the prior means for each treated site j

treated_data<-read.csv("SC100.csv", stringsAsFactors = FALSE)
treated_data



# Convert the relevant columns to numeric explicitly
treated_data[, 2] <- as.numeric(treated_data[, 2])
treated_data[, 3] <- as.numeric(treated_data[, 3])
treated_data[, 4] <- as.numeric(treated_data[, 4])
treated_data[, 5] <- as.numeric(treated_data[, 5])
treated_data[, 6] <- as.numeric(treated_data[, 6])
treated_data[, 7] <- as.numeric(treated_data[, 7])

# Create new columns
treated_data$Total_Before <- treated_data[, 2] + treated_data[, 3] + treated_data[, 4]
treated_data$Total_After  <- treated_data[, 5] + treated_data[, 6] + treated_data[, 7]

# Select the required columns
treated_data <- data.frame(
  treated_data[ , 1],
  Total_Before = treated_data$Total_Before,
  Total_After = treated_data$Total_After,
  treated_data[, 8:13],
  treated_data[, 19:20]
)

# View the result
head(treated_data)



x_vars <- treated_data[, c(4:11)]
x_vars
# Swap columns by name
cols <- names(x_vars)
colA_index <- which(cols == "Column19")
colB_index <- which(cols == "Column20")

# Swap their positions
cols[c(colA_index, colB_index)] <- cols[c(colB_index, colA_index)]

# Reorder the data frame
x_vars <- x_vars[, cols]

# Calculate μ_j = exp(β_0 + Σ(β_p * x_pj))
x_vars[] <- lapply(x_vars, function(col) as.numeric(as.character(col)))
mu_j <- exp(coeff[1] + as.matrix(x_vars) %*% coeff[2:9])



# Calculate α_j = γ/(γ + μ_j)
alpha_j <- shape_param / (shape_param + mu_j)

y_j=treated_data$Total_Before

#--------------STEP 3-----------

posterior_mean=alpha_j*mu_j + (1-alpha_j)*y_j

# Calculate shape and rate
shape <- shape_param + y_j
rate <- shape_param/ mu_j + 1

# Standard deviation formula for Gamma(shape, rate)
posterior_std_dev <- sqrt(shape) / rate

#--------------------STEP 4-------------

#build the table

gamma <- shape_param 

treated_data$E_mj_yj <- posterior_mean
treated_data$SD_mj_yj <- posterior_std_dev

# Observed difference
treated_data$Observed <- treated_data$Total_After - treated_data$Total_Before

# After RTM difference
treated_data$After_RTM <- treated_data$Total_After - treated_data$E_mj_yj

treated_data$y_j=treated_data$Total_Before
treated_data$mu_j=mu_j
treated_data$alpha_j=alpha_j
summary(treated_data)

names(treated_data)[1]="site"

# Select and rename columns for the table
table_data <- treated_data[, c("site","y_j", "mu_j", "alpha_j", "E_mj_yj", "SD_mj_yj", "Total_After", "Observed", "After_RTM")]
names(table_data) <- c("site","y_j", "mu_j", "alpha_j", "E(m_j|y_j)", "SD(m_j|y_j)", "y_j,after", "Observed", "After RTM")

table_data





table_data<-table_data[-1,]
table_data<-table_data[-57,]
rownames(table_data) <- NULL

# Print the table
print(table_data)
write.csv(table_data, "table.csv", row.names = FALSE)
