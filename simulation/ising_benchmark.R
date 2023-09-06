source("ising_dos.R")

L = 6
h = runif(L)*.2 - .1

schemes = c("211", "221", "431", "441", "442", "443", "451", "452", "453", "454", "455", "456", "457", "458", "461", "471", "691", "6101", "6151", "6251", "671", "8171", "8501", "81251", "0171", "041", "0521", "0522", "0881", "0882", "03041", "03042", "01521", "01522")
cycles = c(1, 2, 3, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 6, 7, 9, 10, 15, 25, 7, 17, 50, 125, 3, .6, 8, 8, 13, 13, 47, 47, 30, 23)
scheme.names = c("Leap-Frog", "Omelyan2", "Forest-Ruth", "FR-Type", "small-B", "Non-Unitary1", "Suzuki4", "Algorithm-30", "Optimal-4th-order", "Non-Unitary2", "Omelyan-ST-4", "Non-Unitary3", "Non-Unitary Blanes", "Uniform non-unit", "Blanes-4", "Symplectic-6", "FR-Squared", "Blanes-6", "FR-Suzuki6", "Suzuk",    "i6", "Yoshida", "Morales", "Blanes6-Suzuki8", "Suzuki8", "Taylor-17-sum", "Taylor-4-sum", "Taylor-52-sum", "Taylor-52-prod", "Taylor-88-sum", "Taylor-88-prod", "Taylor-304-sum", "Taylor-304-prod", "Chebyshev-152-sum", "Chebyshev-152-prod")

#dummy = sapply(c(1:7, 9,10, 13, 16), function(s){
#				   cat(paste0("\\subsubsection{", scheme.names[s], "}\n"))
#				   noisy_trace(L, c(1,0,1), h, t=1, cutoff=1, time.series=FALSE, full.data=FALSE, operator.id=-1, first.all.x=TRUE, scheme=schemes[s], n_sources=1)
#})
#q()

## Look at different errors after constant time
t = 100i
effort = exp((0:40)/10)

## XY model split into 2 stages
errors = sapply(seq(schemes), function(i){
					sapply(effort, function(e){
							   Re(noisy_trace(L, c(1,0,1), h, t, t_step=cycles[i]/e, time.series=FALSE, full.data=FALSE, operator.id=-1, first.all.x=TRUE, scheme=schemes[i])$sd)
						})
			})
colnames(errors) = schemes

write.table(cbind(effort, errors), file="benchmark/2-stage_fix-t.txt", row.names=FALSE)

## XY model split into 2L stages
errors = sapply(seq(schemes), function(i){
					sapply(effort, function(e){
							   Re(noisy_trace(L, c(1,0,1), h, t, t_step=cycles[i]/e, time.series=FALSE, full.data=FALSE, operator.id=-1, first.all.x=FALSE, scheme=schemes[i])$sd)
						})
			})
colnames(errors) = schemes

write.table(cbind(effort, errors), file="benchmark/2L-stage_fix-t.txt", row.names=FALSE)

## XXZ model split into 3 stages
errors = sapply(seq(schemes), function(i){
					sapply(effort, function(e){
							   Re(noisy_trace(L, c(1,1,1), h, t, t_step=cycles[i]/e, time.series=FALSE, full.data=FALSE, operator.id=-1, first.all.x=TRUE, scheme=schemes[i])$sd)
						})
			})
colnames(errors) = schemes

write.table(cbind(effort, errors), file="benchmark/3-stage_fix-t.txt", row.names=FALSE)

## XXZ model split into 3L stages
errors = sapply(seq(schemes), function(i){
					sapply(effort, function(e){
							   Re(noisy_trace(L, c(1,1,1), h, t, t_step=cycles[i]/e, time.series=FALSE, full.data=FALSE, operator.id=-1, first.all.x=FALSE, scheme=schemes[i])$sd)
						})
			})
colnames(errors) = schemes

write.table(cbind(effort, errors), file="benchmark/3L-stage_fix-t.txt", row.names=FALSE)

## Look at different errors with constant time-step
t = 10i
dt = .01
times = seq(0, abs(t), by=dt)

## XY model split into 2 stages
errors = sapply(seq(schemes), function(i){
					Re(noisy_trace(L, c(1,0,1), h, t*cycles[i], t_step=dt*cycles[i], time.series=TRUE, full.data=FALSE, operator.id=-1, first.all.x=TRUE, scheme=schemes[i])$sd)
			})
colnames(errors) = schemes

write.table(cbind(times, errors), file="benchmark/2-stage_fix-dt.txt", row.names=FALSE)


## XY model split into 2L stages
errors = sapply(seq(schemes), function(i){
					Re(noisy_trace(L, c(1,0,1), h, t*cycles[i], t_step=dt*cycles[i], time.series=TRUE, full.data=FALSE, operator.id=-1, first.all.x=FALSE, scheme=schemes[i])$sd)
			})
colnames(errors) = schemes

write.table(cbind(times, errors), file="benchmark/2L-stage_fix-dt.txt", row.names=FALSE)

## XXZ model split into 3 stages
errors = sapply(seq(schemes), function(i){
					Re(noisy_trace(L, c(1,1,1), h, t*cycles[i], t_step=dt*cycles[i], time.series=TRUE, full.data=FALSE, operator.id=-1, first.all.x=TRUE, scheme=schemes[i])$sd)
			})
colnames(errors) = schemes

write.table(cbind(times, errors), file="benchmark/3-stage_fix-dt.txt", row.names=FALSE)


## XXZ model split into 3L stages
errors = sapply(seq(schemes), function(i){
					Re(noisy_trace(L, c(1,1,1), h, t*cycles[i], t_step=dt*cycles[i], time.series=TRUE, full.data=FALSE, operator.id=-1, first.all.x=FALSE, scheme=schemes[i])$sd)
			})
colnames(errors) = schemes

write.table(cbind(times, errors), file="benchmark/3L-stage_fix-dt.txt", row.names=FALSE)

### Look at different errors after constant time
#t = 10i
#effort = exp((0:60)/10)
#
### XY model split into 2 stages
#errors = sapply(seq(schemes), function(i){
#					sapply(effort, function(e){
#							   Re(noisy_trace(L, c(.5,0,1), h, t, t_step=cycles[i]/e, time.series=FALSE, full.data=FALSE, operator.id=-1, first.all.x=TRUE, scheme=schemes[i])$sd)
#						})
#			})
#colnames(errors) = schemes
#
#write.table(cbind(effort, errors), file="benchmark/2-stage_fix-t_small-x.txt", row.names=FALSE)
#
### XY model split into 2L stages
#errors = sapply(seq(schemes), function(i){
#					sapply(effort, function(e){
#							   Re(noisy_trace(L, c(.5,0,1), h, t, t_step=cycles[i]/e, time.series=FALSE, full.data=FALSE, operator.id=-1, first.all.x=FALSE, scheme=schemes[i])$sd)
#						})
#			})
#colnames(errors) = schemes
#
#write.table(cbind(effort, errors), file="benchmark/2L-stage_fix-t_small-x.txt", row.names=FALSE)
#
### XXZ model split into 3 stages
#errors = sapply(seq(schemes), function(i){
#					sapply(effort, function(e){
#							   Re(noisy_trace(L, c(.5,.5,1), h, t, t_step=cycles[i]/e, time.series=FALSE, full.data=FALSE, operator.id=-1, first.all.x=TRUE, scheme=schemes[i])$sd)
#						})
#			})
#colnames(errors) = schemes
#
#write.table(cbind(effort, errors), file="benchmark/3-stage_fix-t_small-x.txt", row.names=FALSE)
#
### XXZ model split into 3L stages
#errors = sapply(seq(schemes), function(i){
#					sapply(effort, function(e){
#							   Re(noisy_trace(L, c(.5,.5,1), h, t, t_step=cycles[i]/e, time.series=FALSE, full.data=FALSE, operator.id=-1, first.all.x=FALSE, scheme=schemes[i])$sd)
#						})
#			})
#colnames(errors) = schemes
#
#write.table(cbind(effort, errors), file="benchmark/3L-stage_fix-t_small-x.txt", row.names=FALSE)
