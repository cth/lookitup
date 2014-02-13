source("worker.common.R")

worker({
	print(session$session.key)
	Sys.sleep(10)
	read.table("workers/test.annotation.txt",head=T)
	#data.frame(a=c(1,2,3),b=c("a", "b", "c"))
})
