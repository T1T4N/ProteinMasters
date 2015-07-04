library(lbfgs)
supervised.rw = function(graph, source.node, destinations, alpha = 0.3, 
                         lambda = 1, margin.loss = 0.4, max.iter = 100, quiet = FALSE) {

  features.names = list.edge.attributes(graph)
  features = matrix(nrow = length(features.names), ncol = ecount(graph))
  for(i in 1:length(features.names)) {
    features[i, ] = get.edge.attribute(graph, features.names[i])
  }
  w = rep(1, times=16)
  epsilon = 1e-12
  small.epsilon = 1e-18
  
  logistic.func = function(x, b) {
    1/(1+exp((-1)*x/b))
  }
  
  derivative.logistic.func = function(x, b) {
    (logistic.func(x, b) * (1 - logistic.func(x, b)))/b
  }
  
  strength.function = function(w) {
    logistic.func(w %*% features, 1)
  }
  
  dstrength.function = function(w) {
    as.vector(derivative.logistic.func(w %*% features, 1)) * t(features)
  }
  
  get.Q = function(strength) {
    g = set.edge.attribute(graph, name = "weight", value = strength)
    A = get.adjacency(g, attr = "weight")
    Q = (1-alpha) * (Matrix::Diagonal(x = (Matrix::rowSums(A))^(-1)) %*% A)
    Q[, source.node] = Q[, source.node] + alpha
    Q
  }
  
  get.p = function(Q) {
    p = rep(1/vcount(graph), vcount(graph))
    t = 0
    repeat{
      t = t + 1
      p.new = p %*% Q
      pre = p
      p = p.new
      if((max(abs(pre - p)) < epsilon) || t > max.iter){
        break
      }
    }
    p = as.vector(p)
  }
  
  get.dp = function(w, strength, Q, p) {
    
    get.dQk = function(k, dstrength.w) {
      g = set.edge.attribute(graph, "strength", value = strength)
      g = set.edge.attribute(g, "temp", value = dstrength.w)
      A = get.adjacency(g, attr = "strength")
      dA = get.adjacency(g, attr = "temp")
      A.rowsum = Matrix::Diagonal(x = Matrix::rowSums(A))
      dA.rowsum = Matrix::Diagonal(x = Matrix::rowSums(dA))
      # 1/(A.rowsum)^2
      rec = Matrix::Diagonal(x = Matrix::rowSums(A)^(-2))
      (1-alpha) * (rec %*% (A.rowsum %*% dA - dA.rowsum %*% A))
    }
    
    dp = matrix(0, nrow = vcount(graph), ncol = length(w))
    dstrength = dstrength.function(w)
    for(k in 1:length(w)) {
      t = 0
      dQk = get.dQk(k, dstrength[, k])
      repeat{
        t = t + 1
        dp.new = dp[, k] %*% Q + p %*% dQk
        pre = dp[, k]
        dp[, k] = as.vector(dp.new)
        if((max(abs(pre - dp[, k])) < epsilon) || t > max.iter){
          break
        }
      }
    }
    dp
  }
  
  get.differences = function(p) {
    pd = p[destinations]
    pl = p[-c(destinations, source.node)]
    pd.m = matrix(rep(pd, times = length(pl)), nrow = length(pl), byrow = T)
    pl.m = matrix(rep(pl, times = length(pd)), nrow = length(pl))
    differences = pl.m - pd.m
  }
  
  cost.function = function(w) {
    strength = as.vector(strength.function(w)) + small.epsilon
    Q = get.Q(strength)
    p = get.p(Q)
    diff = get.differences(p)
    h = logistic.func(diff, 0.4)
    sum(w^2) + lambda * sum(h)
  }
  
  gradient.function = function(w) {
    gr = vector(length = length(w))
    strength = as.vector(strength.function(w)) + small.epsilon
    Q = get.Q(strength)
    p = get.p(Q)
    dp = get.dp(w, strength, Q, p)
    diff = get.differences(p)
    dh = derivative.logistic.func(diff, 0.4)
    for(k in 1:length(w)) {
      gr[k] = 2*w[k] + lambda * sum(dh*get.differences(dp[, k]))
    }
    gr
  }
  
  t1 = Sys.time()
  x = lbfgs(cost.function, gradient.function, w, invisible = quiet)
  t2 = Sys.time()
  print(t2 - t1)

  new.p = get.p(get.Q(as.vector(strength.function(x$par)))) 
}
