fun.smap <- function(area,  # área da bacia
                     param, # vetor c/ 6 parâmetros de entrada e 2 condiçoes iniciais
                     epq,   # tabela com dados observados
                     NSE = TRUE){  
  
  m <- length(param)
  
  # Pega precipitação, evapotranspiração e vazão observadas
  e <- epq[,1]    # coluna com a evapotranspiração
  p <- epq[,2]    # coluna com a precipitação
  qobs <- epq[,3] # coluna com a vazão observada
  ndias <- length(p)
  ad <- area
  
  # Parâmetros
  str <- param[1]
  k2t <- param[2]
  ai <- param[3]
  crec <- param[4]
  capc <- param[5]
  kkt <- param[6]
  tuin <- param[7] # condição inicial
  ebin <- param[8] # condição inicial
  
  qca <- rep(0, nrow(epq) + 1) # vetor nulo de vazões calculadas
  es <- rep(0, nrow(epq) + 1)  # vetor nulo de escoamento superficial
  er <- rep(0, nrow(epq) + 1)  # vetor nulo de evapotranspiração real 
  rec <- rep(0, nrow(epq) + 1) # vetor nulo de recarga
  ed <- rep(0, nrow(epq) + 1)  # vetor nulo do fluxo no reservatório rápido
  eb <- rep(0, nrow(epq) + 1)  # vetor nulo de escoamento básico
  rsolo <- rep(0, nrow(epq) + 1)
  rss <- rep(0, nrow(epq) + 1)
  res <- rep(0, nrow(epq) + 1)
  tu <- rep(0, nrow(epq) + 1)
  
  # Conversões de unidade
  capc <- capc/100
  crec <- crec/100
  kkt <- 0.5^(1/kkt);	
  k2t <- 0.5^(1/k2t);	
  
  # Inicialização dos reservatórios
  rsolo[1] <- tuin/100*str         # reservatório do solo
  rss[1] <- ebin/(1 - kkt)/ad*86.4 # reservatório subterrâneo
  res[1] <- 0                      # escoamento
  
  # Itera entre os dias e calcula vazão
  for(i in 2:ndias){
    
    tu[i-1] <- rsolo[i-1]/str # teor de umidade
    
    # Calcula o escoamento superficial
    if(p[i] > ai){
      es[i] <- (p[i] - ai)^2/(p[i] - ai + str - rsolo[i-1])
    }
    else{
      es[i] <- 0
    }
    
    # Calcula a evapotranspiração real
    if((p[i] - es[i]) > e[i]){
      er[i] <- e[i]
    }
    else{
      er[i] <- (p[i] - es[i]) + (e[i] - (p[i] - es[i]))*tu[i-1]
    }
    
    # Calcula a recarga
    if(rsolo[i-1] > capc*str){
      rec[i] <- crec*tu[i-1]*(rsolo[i-1] - (capc*str))
    }
    else{
      rec[i] <- 0
    }
    
    # Atualiza o reservatório do solo
    rsolo[i] <- rsolo[i-1] + p[i] - es[i] - er[i] - rec[i]
    
    # Testa a saturação do solo
    if(rsolo[i] > str){
      es[i] <- es[i] + rsolo[i] - str
      rsolo[i] <- str
    }
    
    ed[i] <- (res[i-1] + es[i])*(1 - k2t) # escoamento direto
    res[i] <- res[i-1] + es[i] - ed[i]    # atualiza reservatório da superficie
    eb[i] <- rss[i-1]*(1 - kkt)           # escoamento básico
    rss[i] <- rss[i-1] - eb[i] + rec[i]   # atualiza reservatório subterrâneo
    qca[i] <- (ed[i] + eb[i])*ad/86.4     # vazão calculada
    
  }
  
  # Dataframe resultante
  resultados <- data.frame(
    # data = as.Date(epq[,1]),
    evap = e,
    precip = p,
    tu = tu[-1],
    rsolo = rsolo[-1],
    rss = rss[-1],
    rec = rec[-1],
    es = es[-1],
    er = er[-1],
    ed = ed[-1],
    eb = eb[-1],
    qobs = qobs,
    qcalc = qca[-1])
  
  # Calcula o coeficiente de Nash-Sutcliff ou retorna o dataframe com os dados simulados
  if(NSE == TRUE){
    ns <- -(1 - sum((qca[-1] - qobs)^2)/sum((qobs - mean(qobs))^2))
    return(ns)
  }
  else{
    ns <- -(1 - sum((qca[-1] - qobs)^2)/sum((qobs - mean(qobs))^2))
    cat(paste("Coeficiente Nash-Sutcliff igual à", round(-ns, digits = 3)))
    return(resultados)
  }
  
}