## Laura M. Whitmore, T Bryce Kelly
## Peak Integration
## Original: August 17, 2021

## load useful toolboxes
library(data.table)
library(openxlsx)

## Set path
#If you open the R.Project file "ICP-MS Peak Integration" the path gets set to that folder

calculate.integrations = function(input.files = list.files(path = './Input/', pattern = '*.TXT', full.names = T),
                                  time_start = 27,
                                  time_end = 70,
                                  dt = 1e-4,
                                  verbose = T) {
  

  ## Start function here
  # Name manipulations and improvements
  file.names = strsplit(input.files, '/')
  for (i in 1:length(file.names)) {
    file.names[[i]] = file.names[[i]][length(file.names[[i]])]
  }
  
  ## Create empty output
  output = data.frame(filename = file.names)
  
  if (verbose) { message('Found ', length(input.files),' files for integration analysis.')}
  
  
  ## Setup archive
  if (verbose) { message('Setting up archive directory...', appendLF = F)}
  time = as.character(Sys.time())
  time = gsub(':', '', time)
  
  if (!dir.exists('./.archive')) {
    dir.create('./.archive')
  }
  archive.dir = paste0('./.archive/', time)
  dir.create(archive.dir)
  
  status = file.copy(input.files, to = archive.dir, overwrite = F)
  
  if (all(status)) {
    if (verbose) { message(' copied files.', appendLF = T) }
  } else {
    if (verbose) { message(' copy failed.', appendLF = T) }
  }
  
  
  ## Process Data
  # For each file
  for (index in 1:length(input.files)) {
    ## load raw file
    if (verbose) { message('Attempting to load file ', file.names[index], appendLF = F)}
    
    headers = fread(file = input.files[index], nrows = 1, header = FALSE)
    data = fread(file = input.files[index], skip = 11, header = FALSE)
    colnames(data) = c('Time', unlist(headers[1, 2:ncol(headers)]), 'Empty')
    
    headers = as.data.frame(headers)
    data = as.data.frame(data)
    
    if (verbose) { message(' (success)')}
    
    ## Fix column names
    for (i in 2:ncol(headers)) {
      headers[1, i] = gsub(pattern = '\\(', replacement = '.', x = headers[1, i])
      headers[1, i] = gsub(pattern = '\\)', replacement = '', x = headers[1, i])
        
      ## Add to output object if not there already
      if (!headers[1, i] %in% colnames(output)) {
        output[,headers[1, i]] = NA
        output[,paste0(headers[1, i], '.baseline')] = NA
      }
    }
    
    if (verbose) { message(' Column names added.')}
    
    if (verbose) { message(' Starting integrations...', appendLF = F)}
    
    ## integral evaluation times 
    xout = seq(time_start, time_end, by = dt)
    
    for (i in 2:ncol(headers)) {
      baseline = sum(approx(x = data$Time, y = data[,i], xout = c(time_start, time_end), rule = 2)$y) * (time_end - time_start) / 2
      output[index, (i-1)*2] = sum(approx(x = data$Time, y = data[,i], xout = xout, rule = 2)$y) * dt #- baseline
      output[index, (i-1)*2+1] = baseline
      
      if (T) {
        png(paste0(archive.dir, '/', gsub('.TXT', paste0(' - ', headers[i], '.png'), file.names[index])))
        plot(data$Time, data[,i],
             ylab = 'Counts', xlab = 'Time (s)', pch = 20,
             ylim = range(pretty(data[,i])), yaxs = 'i',
             xlim = c(0, max(pretty(data$Time))), xaxs = 'i')
        grid(); box()
        mtext(paste0('File: ', file.names[index]), side = 3, adj = 0, cex = 0.8)
        mtext(paste0('Sample: ', headers[i]), side = 3, adj = 1, cex = 0.8)
        
        baseline = approx(x = data$Time, y = data[,i], xout = c(time_start, time_end), rule = 2)
        curve = approx(x = data$Time, y = data[,i], xout = xout, rule = 2)$y
        polygon(x = c(xout, rev(baseline$x)), y= c(curve, rev(baseline$y)), col = '#00000020')
        polygon(x = c(baseline$x, time_end, time_start), y= c(baseline$y, 0, 0), col = '#eebbbb30')
        lines(baseline$x, baseline$y, col = 'red', lwd = 3)
        abline(v = time_start, lty = 2)
        abline(v = time_end, lty = 2)
        
        dev.off()
      }
    }
    
    ## Perform checks
    if (!any(data$Time > time_end)) {
      message('  Warning: file ends early! Check integration.')
    }
    if (data$Time[1] > time_start) {
      message('  Warning: file starts late! Check integration.')
    }
    
    ## Copy image files to plot directory
    if (verbose) { message(' Saving images...', appendLF = F) }
    file.copy(list.files(archive.dir, pattern = '*.png', full.names = T), to = 'Plots')
    if (verbose) { message(' Finished.') }
  }
  
  ## Save archive files
  save(list = ls(), file = paste0(archive.dir, '/workspace.rdata'))
  write.xlsx(output, paste0(archive.dir, '/output.xlsx'))
  if (verbose) { message(' output saved.\nArchiving complete.')}
  
  output
}

## Function: 
# message saying what file we're looking at 
# plot option (peak, peak integral, baseline integral)
# flag if time_end is greater than maximum time in data
# regularize column names


## load data function
## approx function
## plot function?
     
