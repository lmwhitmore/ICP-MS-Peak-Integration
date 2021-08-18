## Laura M. Whitmore
## Peak Integration
## Original: August 17, 2021

## load useful functions
source(file = './R/source_peak_integration.R')

output = calculate.integrations()

write.xlsx(x = output, file = './Output/output.xlsx')
