#grafici Capitolo1, per quelli vm, basta cambiare la function
dev.new()
plot.function.circular(function(x) dwrappednormal(x, mu=circular(0), rho=0.5), xlim=c(-1, 2.2), col="red")

plot.function.circular(function(x) dwrappednormal(x, mu=circular(0), rho=0.9), xlim=c(-1, 2.2), col="blue", add=T)
plot.function.circular(function(x) dwrappednormal(x, mu=circular(0), rho=0.25), xlim=c(-1, 2.2), col="pink", add=T)

