source("worker.common.R")

worker({
	print(session$session.key)
	data.frame(a=c(1,2,3),b=c("a", "b", "c"))
})
