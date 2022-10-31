dyn.load("ising_dos.so")

noisy_trace = function(L, J, h, t, t_step, cutoff, relative.t_step=1, n_sources=2^L, scheme="443", time.series=FALSE, full.data=!time.series, operator.id=0, first.all.x=TRUE, seed=0){
	if(missing(t_step)){
		if(missing(cutoff)){
			volume = (sum(abs(J))*L + sum(abs(h))) * abs(t)
			order = as.numeric(substr(scheme, 1, 1))
			cutoff = ceiling(ifelse(order, abs(t) * volume^(1/order), volume))
		}
		t_step = t / ceiling(cutoff / relative.t_step)
	}

	trace = .Call("noisy_trace", L, J, h, as.complex(t), abs(t_step), n_sources, as.numeric(scheme), time.series, full.data, operator.id, first.all.x, seed)

	if(full.data){
		trace = matrix(trace, ncol=n_sources)
		res = c(2^L, apply(trace, 1, mean))
		dev = c(0, apply(Re(trace), 1, sd) + 1i*apply(Im(trace), 1, sd))
	}else{
		n = length(trace)/2
		res = c(2^L, trace[1:n])
		dev = c(0, trace[-(1:n)])
		trace = NULL
	}
	if(operator.id < 0) res[1] = 0 # results are exact at t=0, i.e. zero deviation from exact solution
	if(operator.id == 1) res[1] = 1 # perfect correlation at t=0

	times = seq(0, t, length.out=length(res))
	if(!time.series){
		times = times[[2]]
		res = res[[2]]
		dev = dev[[2]]
	}

	scale = sqrt((1 - (n_sources-1)/(2^L-1)) / n_sources)
	err = dev * scale
	sff = abs(res)
	d.sff = abs(err)
	
	return(list(samples=trace, times=times, mean=res, sd=dev, scale=scale, err=err, sff=sff, d.sff=d.sff, dt=t_step, L=L, J=J, h=h, scheme=scheme, op.id=operator.id))
}
