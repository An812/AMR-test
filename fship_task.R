library(tidyverse)
library(ggplot2)
library(deSolve)

library(readr)
data <- read.csv("/Users/dodangan/Downloads/interviewtask/test_data_fship.csv")

natural_history <- function(times, states, parameters) {
  with(as.list(c(states, parameters)), {
    N = U + S + D + R + W + M
    
    foi_S = beta * (S + alpha * D + M) / N
    foi_R = beta * f * (R + (1-alpha) * D + W) / N
    
    dU = dur * N - U * (foi_S + foi_R) + gamma_s * (t * S + M) + gamma_r * (t * R + W + t * D) - dur * U
    dS = U * foi_S - gamma_s * t * S - theta * S - foi_R * S  - dur * S 
    dD = foi_R * S + foi_S * R - foi_R * D  - (gamma_s + gamma_r) * t * D - dur * D 
    dR = (U + D) * foi_R  - foi_S * R - gamma_r * t * R - theta * R + gamma_s * t * D  - dur * R 
    dW = theta * R - gamma_r * W - dur * W
    dM = theta * S - gamma_s * M  - dur * M
    
    list(c(dU, dS, dD, dR, dW, dM))
  })
} 

states <- c(U = 1240,
            S = 250,
            D = 45,
            R = 10, 
            W = 2,
            M = 3)

parameters <- c(beta = 0.6,  ## transmission parameter
                alpha = 0.5, ## reduced transmission from dual carriage 
                f = 0.8, ## relative fitness of resistant strain 
                
                gamma_s = 0.8, ## antibiotic exposure rate with susceptible effective drugs
                gamma_r = 0.7, ## antibiotic exposure rate with resistant effective drugs
                t = 0.5,  ## relative lower exposure of treatment for colonised 
                
                theta = 0.01, ## progression to infection rate from colonisation 
                dur = 1 / 8 ## average length of stay
)

run <- as.data.frame(ode(func = natural_history, 
                         y = states, 
                         parms = parameters, 
                         times = seq(0,1500,1))) %>% 
  pivot_longer(cols = "U":"M")

ggplot(run, aes(x = time, y = value, group = name)) + geom_line(aes(col = name), lwd = 2)
#ggplot(run, aes(x = time, y = value, group = name)) + geom_line(aes(col = name), lwd = 2) + scale_y_continuous(lim = c(0,100))


# Task 1

library(DiagrammeR)
grViz("
digraph {
U -> S [label=' '];
  U -> R [label=' '];
  S -> D [label='*'];
  S -> M [label='*'];
  S -> U [label=' '];
  D -> R [label=' '];
  D -> U [label='*'];
  R -> D [label='*'];
  R -> W [label='*'];
  R -> U [label=' '];
  W -> U [label=' '];
  M -> U [label=' '];
}
")

# I conduct stratified age by using S → M, R → W, S → D, R → D, and D → U because age affects progression to 
# infection (theta), driven by immune differences (e.g., weaker in elderly), and colonization 
# transitions (foi_S, foi_R), influenced by age-varying healthcare exposure and antibiotic use, 
# as seen in the data’s resistance trends.


# Task 2
# Task 2.1

# Calculate proportion resistant by age and bacteria
prop_data <- data %>%
  group_by(age, bacteria) %>%
  summarise(
    total = n(),
    resistant = sum(result == 1),
    prop_resistant = resistant / total,
    .groups = "drop"
  )

# Plot proportion resistant by age for each bacteria
plot <- ggplot(prop_data, aes(x = age, y = prop_resistant, color = bacteria)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = FALSE) + 
  labs(
    title = "Proportion of Resistant Infections by Age and Bacteria",
    x = "Age (years)",
    y = "Proportion Resistant",
    color = "Bacteria"
  ) +
  theme_minimal()

# Save 
ggsave("prop_resistant_by_age.png", plot, width = 8, height = 6)

# Print the plot
print(plot)

#The data shows infections from two bacteria, S. aureus and E. coli, in people of different ages. 
#For S. aureus, more older people (like 70+ years) have resistant bacteria, going up to half of them. 
#For E. coli, resistance changes a lot but is higher in some middle-aged (40–60 years) and older people (80+ years). 
#The lines on the chart show that S. aureus resistance goes up steadily, but E. coli resistance goes up and down.

#Task 2.2 
# Model 1

#Logistic Regression Model (Linear Age Effect):
#This model assumes a linear relationship between age and the log-odds of resistance, 
#with an intercept, an age effect, 
#and a bacteria-specific effect to account for differences between S. aureus and E. coli.

#Model 2
#Logistic Regression Model with Quadratic Age Effect:
# This extends the first model by adding a quadratic term (age^2) 
# to capture non-linear trends in resistance with age, and an 
# interaction term (age * bacteria) to allow the age effect to vary by bacteria type.

# Task 2.3. 
# Fit model

# Model 1: Linear Age Effect
model1 <- glm(result ~ age + bacteria, data = data, family = binomial)
summary(model1)

# Model 2: Quadratic Age Effect with Interaction
data$age_sq <- data$age^2  # Add squared term
model2 <- glm(result ~ age + age_sq + bacteria + age:bacteria, data = data, family = binomial)
summary(model2)

# Compare models (if time permits)
AIC(model1, model2)


#Both models confirm resistance increases with age, 
# but Model 2 provides a more nuanced view, showing a non-linear 
#(potentially U-shaped or accelerating) trend and a stronger age 
#effect for S. aureus compared to E. coli.
#The plot aligns with these findings: S. aureus resistance rises steadily with age, while E. coli shows more variability but higher baseline resistance.
