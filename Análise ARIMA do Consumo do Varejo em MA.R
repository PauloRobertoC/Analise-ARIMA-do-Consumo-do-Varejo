### Séries Temporais com R: Analise ARIMA do Consumo do Varejo em MA
### Elaboração: Paulo Roberto Carneiro de Sa

library(BETS)

# Pegando as series a partir do site do Banco Central do Brasil
# Índice de volume de vendas no varejo Total do maranhão
# mensal a partir de jan/2000 até fev/2020 
# 242 observações mensais
# pacotes necessarios:
library(BETS)
library(dygraphs)
library(ggplot2)
library(fpp2)
library(FitAR)
library(tseries)
library(forecast)

######################################################################################################
varejoma <- BETSget(1463) 
print(varejoma)
class(varejoma)
dput(varejoma)  # opção para ter os dados como na structure abaixo

#A rotina de dados obtidos pelo BETS já retorna a série em formato ts, ou seja, série temporal. 
#Farei então a criação de uma série em diferenças para observar o comportamento da série em nível e em diferenças.
#Inicialmente olharei as estatísticas descritivas da série. Em seguida farei um plot básico da série e o plot pelo pacote dygraphs,
#útil para ver os pontos de picos e momentos específicos.

dvarejo <- diff(varejoma)
# estatisticas basicas
summary(varejoma)

### plot basico lembrar que em class(), ele já indicou que era ts = serie
# temporal
plot(varejoma)

# pelo pacote dygraph dá mais opções
library(dygraphs)
help("dygraph")
dygraph(varejoma, main = "Índice de volume de vendas no varejo total do Maranhão <br> (Mensal)  (2011=100) BCB 1463") %>% 
  dyAxis("x", drawGrid = TRUE) %>% dyEvent("2005-1-01", "2005", labelLoc = "bottom") %>% 
  dyEvent("2015-1-01", "2015", labelLoc = "bottom") %>% dyEvent("2018-1-01", 
                                                                "2018", labelLoc = "bottom") %>% dyEvent("2019-2-20", "2020", labelLoc = "bottom") %>% 
  dyOptions(drawPoints = TRUE, pointSize = 2)

### Função de Autocorrelação (FAC) e Autocorrelação parcial (FACp) com defasagem 36
# Série em nível - Usarei a rotina do Hyndman e Athanasopoulos (2018).

library(ggplot2)
varejo <- varejoma
varejo %>% ggtsdisplay(main = "")

### Série em primeira diferença

#Como a série apresenta variações sazonais importantes, 
#assim como tendência importante, claramente não-estacionárias, 
#vou olhar também em primeira diferença.

varejo %>% diff() %>% ggtsdisplay(main = "Série varejo do MA em primeira diferença")

### Primeira diferença sazonal e ACF e PACF
#Farei agora a diferença sazonal.

varejo %>% diff(lag = 12) %>% ggtsdisplay(main = "Série varejo do MA em primeira diferença sazonal")

### Ajuste sazonal e ACF e PACF
#Aplicaremos o ajuste sazonal tipo STL (Seasonal Decomposition of Time Series by Loess) aos dados.

library(fpp2)
varejoadj <- varejo %>% stl(s.window = "periodic") %>% seasadj()
autoplot(varejoadj)

########################### analise ACF e PACF ############################################################# 
#E agora os plots de ACF e PACF na série ajustada sazonalmente.
varejoadj %>% diff() %>% ggtsdisplay(main = "Série varejo do MA em primeira diferença e ajuste sazonal")
#Parece indicar para algum AR de ordem 4 na PACF,
#alguma sazonalidade marcante na ACF (indicando um termo SMA), 
#e indicacao de um ARIMA (4,1,1)(0,1,1)[12].

#Uma estimação inicial desse modelo daria algo como: 
(fit <- Arima(varejoadj, order = c(4, 1, 1), seasonal = c(0, 1, 1)))
 
#############################################################################################################

autoplot(fit)

################### analise ARIMA ################################

checkresiduals(fit)

##################################################################

library(FitAR)  # pratico para plotar graficos de Q
# LBQPlot(res, lag.max = 30) res é a série que temos interesse de realizar
# os testes lag.max é o número de lags que se deseja imprimir o gráfico
LBQPlot(residuals(fit), 36)

# checar normalidade dos resíduos
library(tseries)
jarque.bera.test(residuals(fit))

########################## ARIMA #########################################
#Uma estimação automatizada pela função auto.arima {forecast} padrão indica um modelo diferente.
library(forecast)
fit2 <- auto.arima(varejo)
summary(fit2)

# p-values dos coeficientes
pvalues <- (1 - pnorm(abs(fit2$coef)/sqrt(diag(fit2$var.coef)))) * 2
pvalues

library(FitAR)  # pratico para plotar graficos de Q
LBQPlot(residuals(fit2), 36)

# checar normalidade dos resíduos
library(tseries)
jarque.bera.test(residuals(fit2))

######################## Modelo ARIMA automático mais geral ###########################################
require(fpp2)
fit3 <- auto.arima(varejoma, stepwise = FALSE, approximation = FALSE)
summary(fit3)

# p-values dos coeficientes
pvalues <- (1 - pnorm(abs(fit3$coef)/sqrt(diag(fit3$var.coef)))) * 2
pvalues

##################### Checando resíduos ##################################################
#Os modelos podem ser considerados white noise, mas não normais.

checkresiduals(fit3)

require(FitAR)  # pratico para plotar graficos de Q
LBQPlot(residuals(fit3), 36)

# checar normalidade dos resíduos
require(tseries)
jarque.bera.test(residuals(fit3))

#Pode-se fazer o gráfico dos forecasts deste modelo:
autoplot(forecast(fit3, h = 24), title = "Forecasts de volume de vendas no varejo total do Maranhão", 
         xlab = "Ano", ylab = "Índice (2011=100)")

######################################################################################################
# Teste de raiz unitária
# O teste de Dickey-Fuller aumentado (ADF) é nossa primeira escolha, 
# O padrão é que H0: série não estacionária ou série tem raiz unitária.
ADF.test <- adf.test(varejoma)
ADF.test

######################################################################################################
# Teste ADF invertendo a hipótese nula
#Neste caso, o padrão é que H0: a série é estacionária ou a série não tem raiz unitária, 
#Neste caso, a hipótese alternativa é que a série é "explosiva".

ADF.test12 <- adf.test(varejoma, k = 12, alternative = c("explosive"))
ADF.test12

######################################################################################################
# ARIMA - Codigo do Hyndman e Athanasopoulos 
x <- varejoma
ns <- nsdiffs(x)  #numero de diferencas sazonais
ns

if (ns > 0) {
  xstar <- diff(x, lag = frequency(x), differences = ns)
} else {
  xstar <- x
}

nd <- ndiffs(xstar)  # numero de diferencas
nd

if (nd > 0) {
  xstar <- diff(xstar, differences = nd)
}
# xstar será a série devidamente em diferenças sazonais e diferenças normais
autoplot(xstar)

auto.arima(xstar) # ver que o resultado é igual, mas agora nao tem diferencas

modelo <- auto.arima(xstar, stepwise = FALSE, approximation = FALSE)
summary(modelo)

autoplot(modelo)








































































































































