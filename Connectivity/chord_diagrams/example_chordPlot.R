library(circlize)
circos.clear()
pathname = "/Users/stepeter/Documents/Neuroimage_revisions/event_timing_data/Connectivity_median_plots/"
png(paste(pathname,"example_chord.png",sep=""),bg="transparent",width = 720, height = 720)
color_vals = c('purple','orange','blue','magenta','gold','red','cyan','green')#'#FF00FF','#FF6600','#0000FF','#FF0000','#D4AA00','#800080','#00FF00','#00FFFF')
color_vals_links = adjustcolor(color_vals, alpha.f = 0.6)
factors = c("1","2","3","4","5","6","7","8")
circos.par(cell.padding = c(0, 0, 0, 0),start.degree = 24) 
circos.initialize(factors, xlim = cbind(c(0,0,0,0,0,0,0,0), c(12,12,12,12,12,12,12,12))) 


circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.05,
                       #panel.fun for each sector
                       panel.fun = function(x, y) {
                         offset = 24
                         #select details of current sector
                         name = get.cell.meta.data("sector.index")
                         i = get.cell.meta.data("sector.numeric.index")
                         xlim = get.cell.meta.data("xlim")
                         ylim = get.cell.meta.data("ylim")
                         
                         #text direction (dd) and adjusmtents (aa)
                         theta = circlize(mean(xlim), 1.3)[1, 1] %% 360
                         dd <- ifelse(theta < 90 || theta > 270, "clockwise", "reverse.clockwise")
                         aa = c(1, 0.5)
                         if(theta < 90 || theta > 270)  aa = c(0, 0.5)
                         
                         #plot cortical labels
                         #circos.text(x=mean(xlim), y=1.7, labels=name, facing = dd, cex=0.6,  adj = aa)
                         
                         #plot main sector
                         circos.rect(xleft=xlim[1], ybottom=ylim[1], xright=xlim[2], ytop=ylim[2], 
                                     col = color_vals[i], border=color_vals[i])
                         
                         #plot axis
                         circos.axis(labels.cex=2.8, direction = "outside", major.at=c(0,4,8,13), #seq(from=0,to=floor(12),by=4), 
                                     minor.ticks=1, labels.away.percentage = 0.15, lwd=2,labels.niceFacing=FALSE)
                       }) 

# #Add nonsignificant links first
# curr_starts_gray=matrix(0L, nrow = 1, ncol = 8)
# for(i in 1:8) {
#   for(j in 1:8){
#     if((A[i,j]!=0 | A[j,i]!=0) & (A[j,i]>A[i,j])){
#       if(sign_vals[i,j]==0){
#         if(plotBand_cort=='theta'){
#           data_temp = abs(netData_ave$netVals.ave[[1]])*10^4
#           circos.link(factors[i], c(curr_starts_gray[i],curr_starts_gray[i]+data_temp[i,j]), factors[j], c(curr_starts_gray[j],curr_starts_gray[j]+data_temp[j,i]), col = adjustcolor('gray', alpha.f = 0.4)) 
#         }else{
#           data_temp = abs(netData_ave$netVals.ave[[2]])*10^4
#           circos.link(factors[i], c(curr_starts_gray[i],curr_starts_gray[i]+data_temp[i,j]), factors[j], c(curr_starts_gray[j],curr_starts_gray[j]+data_temp[j,i]), col = adjustcolor('gray', alpha.f = 0.4)) 
#         }
#         curr_starts_gray[i] = curr_starts_gray[i]+data_temp[i,j]
#         curr_starts_gray[j] = curr_starts_gray[j]+data_temp[j,i]
#       }
#     }
#   }
# }

circos.link(factors[3], c(0,4), factors[7], c(0,8), col = adjustcolor('blue', alpha.f = 0.6))
circos.link(factors[1], c(0,6), factors[5], c(0,0), col = adjustcolor('red', alpha.f = 0.6))
circos.link(factors[1], c(6,8), factors[2], c(0,2), col = adjustcolor('magenta', alpha.f = 0.6))

# curr_starts=matrix(0L, nrow = 1, ncol = 8)
# for(i in 1:8) {
#   for(j in 1:8){
#     if((A[i,j]!=0 | A[j,i]!=0) & (A[j,i]>A[i,j])){
#       if(sign_vals[i,j]>0 & sign_vals[j,i]>0){
#         circos.link(factors[i], c(curr_starts[i],curr_starts[i]+A[i,j]), factors[j], c(curr_starts[j],curr_starts[j]+A[j,i]), col = adjustcolor('red', alpha.f = 0.6))#color_vals_links[8]) 
#         curr_starts[i] = curr_starts[i]+A[i,j]
#         curr_starts[j] = curr_starts[j]+A[j,i]
#       }else if(sign_vals[i,j]<0 & sign_vals[j,i]<0){
#         circos.link(factors[i], c(curr_starts[i],curr_starts[i]+A[i,j]), factors[j], c(curr_starts[j],curr_starts[j]+A[j,i]), col = adjustcolor('blue', alpha.f = 0.6)) 
#         curr_starts[i] = curr_starts[i]+A[i,j]
#         curr_starts[j] = curr_starts[j]+A[j,i]
#       }else if(sign_vals[i,j]!=0 & sign_vals[j,i]!=0){
#         circos.link(factors[i], c(curr_starts[i],curr_starts[i]+A[i,j]), factors[j], c(curr_starts[j],curr_starts[j]+A[j,i]), col = adjustcolor('magenta', alpha.f = 0.6)) 
#         curr_starts[i] = curr_starts[i]+A[i,j]
#         curr_starts[j] = curr_starts[j]+A[j,i]
#       }
#     }
#   }
# }
dev.off()
