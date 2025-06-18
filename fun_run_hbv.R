# PACOTES ####
pacman::p_load(pacman,
               tidyverse,
               lubridate,
               janitor)

# DADOS ####

# Temperatura e evapotranspiração médias mensais
aux.monthly.evap <- read.table(file = "Temperatura e Evapotranspiração.txt",
                               dec = ",",
                               sep = "\t",
                               header = TRUE) %>%
  rename("Tm" = "Monthly.Tave.", "PEd" = "Daily.PEm") %>%
  select(-PEm)

# Temperatura, precipitação e vazão observados
aux.obs.data <- read.table(file = "Dados Observados.txt",
                           dec = ",",
                           sep = "\t",
                           header = TRUE)

# HBV Concept p/ fazer comparação
# aux.hbv.concept <- read.table(file = "HBV Concept.txt",
#                               dec = ",",
#                               sep = "\t",
#                               header = TRUE) %>% clean_names()

# Dataframe c/ parâmetros
df.param <- {data.frame(ac = 410, # área de contribuição [km²]
                        dd = 3,
                        fc = 150, # capacidade de campo
                        pwp = 120, # ponto de murcha
                        c = 0.02,
                        beta = 2,
                        k0 = 0.5,
                        k1 = 0.15,
                        k2 = 0.01,
                        l1 = 21,
                        kperc = 0.01)}

# Dataframe c/ condições iniciais
df.cond.iniciais <- {data.frame(sn = 25 %>% as.numeric,   # snow [mm]
                                sm = 100 %>% as.numeric,  # soil moisture
                                s1 = 2 %>% as.numeric,    # reservatório superior
                                s2 = 200 %>% as.numeric)}

# FUNÇÕES ####

# Aplica as funções anteriores a partir de um dataframe com dados observados
fun.hbv.run <- function(df,            # dataframe c/ Date, Temp, Precip, Q observados
                        param,         # dataframe c/ parâmetros
                        cond.iniciais, # dataframe c/ condições iniciais
                        monthly.evap   # dataframe c/ evapotranspiração potencial média mensal
                        ){
  
  # Carregar pacotes necessários
  if(!require("pacman")) install.packages("pacman")
  pacman::p_load(tidyverse, lubridate)
  
  # Calcular PEa
  df <- df %>% 
    merge(monthly.evap, by = "Month") # os dados de evapotranspiração potencial mensal
    
  
  # Adicionar linha em branco no início do df
  row.na <- data.frame(matrix(data = NA, nrow = 1, ncol = ncol(df))) # cria um dataframe somente c/ linhas vazias e a mesma qtde de colunas que df
  colnames(row.na) <- colnames(df)                                   # atribui os mesmos nomes
  df <- rbind(row.na, df) %>% 
    mutate(Date = lubridate::dmy(Date)) %>% 
    arrange(Date)
  
  # Criar vetores vazios do tamanho dos dados, mas com uma linha a mais para armazenar os valores iniciais
  sn <- rep(0, nrow(df) + 1)   # snow
  lw <- rep(0, nrow(df) + 1)   # liquid water
  pe <- rep(0, nrow(df) + 1)   # potential evap
  dq <- rep(0, nrow(df) + 1)   # effective precip
  ea <- rep(0, nrow(df) + 1)   # evapotranspiração real
  sm <- rep(0, nrow(df) + 1)   # soil moisture
  s1 <- rep(0, nrow(df) + 1)   # reservatório superior
  s2 <- rep(0, nrow(df) + 1)   # reservatório inferiror
  qs <- rep(0, nrow(df) + 1)   # vazão total simulada [m³/s]
  dif1 <- rep(0, nrow(df) + 1) # (Qsim - Qobs)²
  dif2 <- rep(0, nrow(df) + 1) # (Qobs - Qm)²
  
  # Condições iniciais
  sn[1] <- cond.iniciais$sn
  sm[1] <- cond.iniciais$sm
  s1[1] <- cond.iniciais$s1
  s2[1] <- cond.iniciais$s2
  
  # Cálculo da média observada
  media.q.obs <- mean(df$Q, na.rm = TRUE)
  
  # Loop p/ calcular coisas
  for(i in 1:nrow(df) + 1){
    
    # Snow
    sn[i] <- ifelse(df$Temp[i] > 0, # checa se a temperatura no dia é maior que 0
                    max(sn[i-1] - param$dd*(df$Temp[i]), 0),
                    sn[i-1] + df$Precip[i])
    
    # Liquid Water
    lw[i] <- ifelse(df$Temp[i] > 0, # checa se a temperatura no dia é maior que 0
                    df$Precip[i] + min(sn[i-1], param$dd*(df$Temp[i])),
                    0)
    
    # Effective Precipitation
    dq[i] <- lw[i]*(sm[i-1]/param$fc)^param$beta
    
    # Potencial Evap
    pe[i] <- (1 + param$c*(df$Temp[i] - df$Tm[i]))*df$PEd[i]
    
    # Real Evap
    ea[i] <- ifelse(sm[i-1] > param$pwp, # checa se a umidade do solo é maior que o ponto de murcha
                    pe[i],
                    pe[i]*(sm[i-1]/param$pwp))
    
    # Balanço no solo
    sm[i] <- sm[i-1] + lw[i] - dq[i] - ea[i]
    
    # Balanço no reservatório superior
    s1[i] <- s1[i-1] + dq[i] - max(0, s1[i-1] - param$l1)*param$k0 - s1[i-1]*param$k1 - s1[i-1]*param$kperc
    
    # Balanço no inferior
    s2[i] <- s2[i-1] + s1[i-1]*param$kperc - s2[i-1]*param$k2
    
    # Vazão total simulada [m³/s]
    qs[i-1] <- (max(0, s1[i-1] - param$l1)*param$k0 + s1[i-1]*param$k1 + s2[i-1]*param$k2)*param$ac*1000/86400
    
  }
  
    # Resultados
  resultados <- df %>% 
    mutate(Snow = lag(sn[-1]),
           LiquidWater = lag(lw[-1]),
           DQ = lag(dq[-1]),
           SoilMoisture = lag(sm[-1]),
           PotentialEvap = lag(pe[-1]),
           Evap = lag(ea[-1]),
           S1 = lag(s1[-1]),
           S2 = lag(s2[-1]),
           QSim = lag(qs[-1]),
           Dif1 = ((Q - QSim)^2),
           Dif2 = (Q - media.q.obs)^2) %>% 
    select(-c(PEd, Month))
  
  return(resultados[-1,])
  
}

# SIMULAÇÃO ####

# Run
df.hbv <- fun.hbv.run(df = aux.obs.data,
                      param = df.param,
                      cond.iniciais = df.cond.iniciais,
                      monthly.evap = aux.monthly.evap)

# Totais
total.evap <- sum(df.hbv$Evap, na.rm = TRUE)                                   # mm
total.precip <- sum(df.hbv$Precip, na.rm = TRUE)                               # mm
total.discharge <- total.precip - total.evap                                   # m/h.km²
total.sim.discharge <- sum(df.hbv$QSim, na.rm = TRUE)*86400/(df.param$ac*1000) # m/h.km²
total.obs.discharge <- sum(df.hbv$Q, na.rm = TRUE)*86400/(df.param$ac*1000)    # m/h.km²
erro <- abs(total.obs.discharge - total.sim.discharge)*100/total.obs.discharge # %
square.diff <- sum(df.hbv$Dif1, na.rm = TRUE)
obs.deviation <- sum(df.hbv$Dif2, na.rm = TRUE)
corr <- cor(df.hbv$QSim, df.hbv$Q, use = "complete.obs")
NS <- 1 - square.diff/obs.deviation


# write.table(x = df.hbv$QSim,
#             file = "Vazão Simulada Teste.txt",
#             row.names = FALSE,
#             fileEncoding = "UTF-8")
