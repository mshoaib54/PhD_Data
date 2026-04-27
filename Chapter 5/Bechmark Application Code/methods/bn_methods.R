bn_tabu <- function(train, test){
  bn <- bnlearn::tabu(train)
  bn.fit <- bnlearn::bn.fit(bn, train)
  predict(bn.fit, node = "answer", data = test)
}

bn_hc <- function(train, test) {
  bn <- bnlearn::hc(train)
  bn.fit <- bnlearn::bn.fit(bn, train)
  predict(bn.fit, node = "answer", data = test)
}