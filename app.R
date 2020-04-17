#-----------------------Code Description---------------------------------------#
# Notes:
# ver1.0, data: 20200414, by MaoYan
# ShinyApp, 交通运行与多源数据融合的安全风险特征数据平台原型系统。
#------------------------------------------------------------------------------#

if (!require(shiny)) {install.packages("shiny", repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")}
if (!require(leaflet)){install.packages("leaflet", repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")}
if (!require(readxl)){install.packages("readxl", repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")}
if (!require(magrittr)) {install.packages("magrittr", repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")}
if (!require(plyr)) {install.packages("plyr", repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")}
if (!require(ggplot2)) {install.packages("ggplot2", repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")}

# library(shiny)
# library(leaflet)
# library(magrittr)
# library(plyr)
# library(ggplot2)


CalcDistance <- function(Lat_A, Lng_A, Lat_B, Lng_B){
    # 考虑赤道与两极半径不同，地球是个椭球体
    ra  =  6378.140  # 赤道半径 (km)
    rb  =  6356.755  # 极半径 (km)
    flatten  =  (ra - rb) / ra  # 地球扁率
    rad_lat_A  =  Lat_A * pi/180
    rad_lng_A  =  Lng_A * pi/180
    rad_lat_B  =  Lat_B * pi/180
    rad_lng_B  =  Lng_B * pi/180
    pA  =  atan(rb / ra * tan(rad_lat_A))
    pB  =  atan(rb / ra * tan(rad_lat_B))
    xx  =  acos(sin(pA) * sin(pB) +
                    cos(pA) * cos(pB) * cos(rad_lng_A - rad_lng_B))
    c1  =  (sin(xx) - xx) * (sin(pA) + sin(pB)) ** 2 / cos(xx / 2) ** 2
    c2  =  (sin(xx) + xx) * (sin(pA) - sin(pB)) ** 2 / sin(xx / 2) ** 2
    dr  =  flatten / 8 * (c1 - c2)
    distance  =  ra * (xx + dr)
    return (distance)
}


Hageodist <- function(L1, phi1, L2, phi2){  # 只考虑地球半径，假设地球是个球体
    
    a = 6378.14
    f = 1/298.257
    F = (phi1+phi2)/2
    G = (phi1 - phi2)/2
    ramda  =  (L1 - L2)/2
    S = (sin(G*pi/180)^2)*(cos(ramda*pi/180)^2) +
        (cos(F*pi/180)^2)*(sin(ramda*pi/180)^2)
    C= (cos(G*pi/180)^2)*(cos(ramda*pi/180)^2) +
        (sin(F*pi/180)^2)*(sin(ramda*pi/180)^2)
    omega = atan(sqrt(S/C))
    R = sqrt(S*C)/omega
    D = 2*omega*a
    H1 = (3*R-1)/(2*C)
    H2 = (3*R+1)/(2*S)
    res = D*(1 + f*H1*(sin(F*pi/180)^2)*(cos(G*pi/180)^2) -
                 f*H2*(cos(F*pi/180)^2)*(sin(G*pi/180)^2))
    return(round(res,3))
}


SingleDataINI <- function(GPSData){  # 数据初始化处理
    
    GPSData <- GPSData[!duplicated(GPSData$gpsTime),]  # 删除重复的数据
    
    GPSData$gpsTime <- as.POSIXlt(GPSData$gpsTime)  # 将时间列转化为时间类型
    GPSData$gpsTime_UTC <- as.numeric(GPSData$gpsTime)
    
    GPSData <- GPSData[order(GPSData$gpsTime),]  # 按GPS时间重新排序
    GPSData$gpsTime_diff <- c(0, abs(diff(GPSData$gpsTime)))  # 相邻数据时间差
    
    # 计算相对于起点的距离坐标
    begin_x <- GPSData$locLat[1]
    begin_y <- GPSData$locLong[1]
    GPSData$coordsX <- CalcDistance(Lat_A = begin_x,
                                    Lng_A = begin_y,
                                    Lat_B = GPSData$locLat,
                                    Lng_B = begin_y) * 1000
    GPSData$coordsY <- CalcDistance(Lat_A = begin_x,
                                    Lng_A = begin_y,
                                    Lat_B = begin_x,
                                    Lng_B = GPSData$locLong) * 1000
    
    GPSData$speedChange <- c(0, diff(GPSData$gpsSpeed))  # 速度变化
    GPSData$calcAcc <- GPSData$speedChange / GPSData$gpsTime_diff  # 加速度
    GPSData$speedSplit <- GPSData$gpsSpeed %/% 10 + 1  # 速度分组
    
    # 计算间距，相邻数据行
    GPSData$calcSpacing <- 0
    for (i in 1:length(GPSData$vehID)) {
        if (i > 1) {
            GPSData$calcSpacing[i] <- CalcDistance(Lat_A = GPSData$locLat[i-1],
                                                   Lng_A = GPSData$locLong[i-1],
                                                   Lat_B = GPSData$locLat[i],
                                                   Lng_B = GPSData$locLong[i]) * 1000
            if (is.nan(GPSData$calcSpacing[i])) {
                GPSData$calcSpacing[i] <- 0
            }
        }
    }
    
    # 方位角变化
    GPSData$angleChange <- c(0, diff(GPSData$vehDir))
    # 方位角相对行驶距离的变化
    GPSData$angleChangeRate <- abs(GPSData$angleChange / GPSData$calcSpacing)
    
    # Nan值处理
    GPSData$coordsX[is.nan(GPSData$coordsX)] <- 0
    GPSData$coordsY[is.nan(GPSData$coordsY)] <- 0
    GPSData$speedChange[is.nan(GPSData$speedChange)] <- 0
    GPSData$calcAcc[is.nan(GPSData$calcAcc)] <- 0
    GPSData$calcSpacing[is.nan(GPSData$calcSpacing)] <- 0
    GPSData$angleChange[is.nan(GPSData$angleChange)] <- 0
    GPSData$angleChangeRate[is.nan(GPSData$angleChangeRate)] <- 0
    
    # 利用lubridate包函数提取日期和时间
    #  GPSData$tYear <- year(GPSData$gpsTime)
    #  GPSData$tMonth <- month(GPSData$gpsTime)
    #  GPSData$tDay <- day(GPSData$gpsTime)
    #  GPSData$tHour <- hour(GPSData$gpsTime)
    #  GPSData$weekDay <- wday(GPSData$gpsTime)
    
    return(GPSData)
}


fun_abnormalACC <- function(data, probs = 0.95){  # 异常加减速行为标准分析
    
    AAC <- subset(data, calcAcc > 0)
    DAC <- subset(data, calcAcc < 0)
    
    b1 <- subset(AAC, select = c("speedSplit", "calcAcc"))
    b1 <- ddply(b1, .(speedSplit),
                numcolwise(quantile), probs = c(probs), na.rm = TRUE)
    
    b2 <- subset(DAC, select = c("speedSplit", "calcAcc"))
    b2 <- ddply(abs(b2), .(speedSplit),
                numcolwise(quantile), probs = c(probs), na.rm = TRUE)
    
    names(b1) <- c("speedSplit", "abnlAAC")
    names(b2) <- c("speedSplit", "abnlDAC")
    
    c = merge(b1, b2, all = TRUE)
    c$speedBottom <- (c$speedSplit - 1) * 10
    c$speedTop <- c$speedSplit * 10
    
    return(c)
}


ShowProgress <- function(LoopIdx, kLoopLength){  # 可用于查看循环执行进度
    
    if (LoopIdx == kLoopLength) {
        print("--------->100% Complete")
    } else if (LoopIdx == kLoopLength %/% 10 * 9) {
        print("-------->-90% Complete")
    } else if (LoopIdx == kLoopLength %/% 10 * 8) {
        print("------->--80% Complete")
    } else if (LoopIdx == kLoopLength %/% 10 * 7) {
        print("------>---70% Complete")
    } else if (LoopIdx == kLoopLength %/% 10 * 6) {
        print("----->----60% Complete")
    } else if (LoopIdx == kLoopLength %/% 10 * 5) {
        print("---->-----50% Complete")
    } else if (LoopIdx == kLoopLength %/% 10 * 4) {
        print("--->------40% Complete")
    } else if (LoopIdx == kLoopLength %/% 10 * 3) {
        print("-->-------30% Complete")
    } else if (LoopIdx == kLoopLength %/% 10 * 2) {
        print("->--------20% Complete")
    } else if (LoopIdx == kLoopLength %/% 10 * 1) {
        print(">---------10% Complete")
    }
}


# 将GPS坐标转换为高德火星坐标，主函数是GPSToGaoDecoords( GPSData)----
transformLon <- function(x, y) {
    ret <- 300.0 + x + 2.0 * y + 0.1 * x * x + 0.1 * x * y + 0.1 * sqrt(abs(x))
    ret <- ret + (20.0 * sin(6.0 * x * pi) +
                      20.0 * sin(2.0 * x * pi)) * 2.0 / 3.0
    ret <- ret + (20.0 * sin(x * pi) +
                      40.0 * sin(x / 3.0 * pi)) * 2.0 / 3.0
    ret <- ret + (150.0 * sin(x / 12.0 * pi) +
                      300.0 * sin(x / 30.0 * pi))* 2.0 / 3.0
    return (ret)
}


transformLat <- function(x, y) {
    ret <- (-100.0) + 2.0 * x + 3.0 * y + 0.2 * y * y +
        0.1 * x * y + 0.2 *sqrt(abs(x))
    ret <- ret + (20.0 * sin(6.0 * x * pi) +
                      20.0 * sin(2.0 * x * pi)) * 2.0 / 3.0
    ret <- ret + (20.0 * sin(y * pi) +
                      40.0 * sin(y / 3.0 * pi)) * 2.0 / 3.0
    ret <- ret + (160.0 * sin(y / 12.0 * pi) +
                      320 * sin(y * pi / 30.0)) * 2.0 / 3.0
    return (ret)
}


GPSToGaoDecoords <- function(GPSData) {
    a <- 6378245.0
    ee <- 0.00669342162296594323
    # colnames(GPSData) <- c("long", "lat")
    GPSData$dLat <- transformLat(GPSData$long - 105.0, GPSData$lat - 35.0)
    GPSData$dLong <- transformLon(GPSData$long - 105.0, GPSData$lat - 35.0)
    GPSData$radLat <- GPSData$lat / 180.0 * pi
    GPSData$magic <- sin(GPSData$radLat)
    GPSData$magic <- 1 - ee * GPSData$magic * GPSData$magic
    GPSData$sqrtMagic <- sqrt(GPSData$magic)
    GPSData$dLat <- (GPSData$dLat * 180.0) /
        ((a * (1 - ee)) / (GPSData$magic * GPSData$sqrtMagic) * pi)
    GPSData$dLong <- (GPSData$dLong * 180.0) /
        (a / GPSData$sqrtMagic * cos(GPSData$radLat) * pi)
    GPSData$latitude <- GPSData$lat + GPSData$dLat
    GPSData$longitude <- GPSData$long + GPSData$dLong
    return(GPSData)
}


# 原始数据导入----
df.gpsdatasichuang <- read.csv(file = "Data/GPSDataSichuang.csv")


# 异常行为标记函数----
AddAbnlAccTag <- function(data){  # 给异常加减速行为加标签
    
    data$abnlAccTag <- ifelse(is.na(data$abnlAAC),
                              "Normal",
                              ifelse(data$calcAcc > data$abnlAAC,
                                     "Abnormal",
                                     "Normal"))  # 异常加速行为标签
    
    data$abnlDacTag <- ifelse(is.na(data$abnlDAC),
                              "Normal",
                              ifelse(data$calcAcc < data$abnlDAC,
                                     "Abnormal",
                                     "Normal"))  # 异常减速行为标签
    
    return(data)
}


# 异常加速行为标签
df.abnlaccsichuang <- subset(x = AddAbnlAccTag(data = df.gpsdatasichuang),
                             abnlAccTag == "Abnormal")
# 异常减速行为标签
df.abnldacsichuang <- subset(x = AddAbnlAccTag(data = df.gpsdatasichuang),
                             abnlDacTag == "Abnormal")
# 异常加速行为、异常减速行为数据集
df.abnlsichuang <- rbind(df.abnlaccsichuang, df.abnldacsichuang)
# 正常行为数据集
df.nldatasichuang <- subset(x = AddAbnlAccTag(data = df.gpsdatasichuang),
                            abnlAccTag == "Normal" & abnlDacTag == "Normal")


ui <- bootstrapPage(
    tags$style(type = "text/css", "html, body {width:100%;height:100%}"),
    leafletOutput("mymap", width = "100%", height = "100%"),
    absolutePanel(top = 10, right = 20, draggable = TRUE,
                  width = 500,
                  class = "panel panel-default",

                  h3("交通运行安全风险特征数据平台原型系统", align = "center"),
                  br(),
                  h4(textOutput("reactiveDataVol"), align = "left"),
                  h4(textOutput("reactiveVehNum"), align = "left"),
                  br(),
                  h4("异常行为数据分布", align = "left"),
                  plotOutput("reactivePlot",
                             height = "130px", width = "100%"),
                  br(),
                  h4("异常加速行为统计", align = "left"),
                  plotOutput("reactiveAccPlot",
                             height = "130px", width = "100%"),
                  br(),
                  h4("异常减速行为统计", align = "left"),
                  plotOutput("reactiveDacPlot",
                             height = "130px", width = "100%"),
                  br(),
                  selectInput("ptsTyp", label = "选择异常行为",
                              choices = list("全部异常行为点" = "abnlpts",
                                             "异常加速行为点" = "abnlaccpts",
                                             "异常减速行为点" = "abnldacpts"),
                              selected = "abnlpts",
                  ),
                  br(),
                  p("Copyright 2020 RIOH | All Rights Reserved", align = "center")
    )
)


server <- function(input, output, session) {
    
    output$mymap <- renderLeaflet({
        leaflet() %>%
            addTiles() %>%
            setView(lng = 105.3367, lat = 30.1124, zoom = 7)
    })
    
    filteredData <- reactive({
        
        if (input$ptsTyp == "abnlpts") {
            df.gpsdata <- df.abnlsichuang
        } else if (input$ptsTyp == "abnlaccpts") {
            df.gpsdata <- df.abnlaccsichuang
        } else if (input$ptsTyp == "abnldacpts") {
            df.gpsdata <- df.abnldacsichuang
        }
        
        df.gpsdata
    })
    
    colorPal <- reactive({
        
        if (input$ptsTyp == "abnlpts") {
            kcolorPal <- "black"
        } else if (input$ptsTyp == "abnlaccpts") {
            kcolorPal <- "red"
        } else if (input$ptsTyp == "abnldacpts") {
            kcolorPal <- "blue"
        }
        
        kcolorPal
    })
    
    observe({
        leafletProxy("mymap", data = filteredData()) %>%
            clearMarkers() %>%
            addCircleMarkers(lng = ~locLong, lat = ~locLat,
                             radius = 3, color = colorPal())
    })
    
    output$reactiveDataVol <- renderText({
        paste0("数据样本容量：",
               " ",
               nrow(df.gpsdatasichuang),
               "个")
    })
    
    output$reactiveVehNum <- renderText({
        paste0("驾驶人/车辆数：",
               " ",
               length(unique(df.gpsdatasichuang$vehID)),
               "人/辆")
    })
    
    output$reactivePlot <- renderPlot({
        
        ggplot(data = df.nldatasichuang, aes(x = speedBottom, y = calcAcc)) +
            geom_jitter(colour = "grey", size = 2, alpha = 0.5) +
            geom_jitter(data = df.abnlaccsichuang,
                        aes(x = speedBottom, y = calcAcc),
                        colour = "red", size = 2, alpha = 0.5) +
            geom_jitter(data = df.abnldacsichuang,
                        aes(x = speedBottom, y = calcAcc),
                        colour = "blue", size = 2, alpha = 0.5) +
            scale_x_continuous(name = "Speed Zone, km/h",
                               limits = c(0, 90),
                               breaks = seq(0, 90, 10),
                               labels = c("0-10", "10-20", "20-30", "30-40",
                                          "40-50", "50-60", "60-70", "70-80",
                                          "80-90", "90-100")) +
            scale_y_continuous(name = "Calc Acc, m/s2",
                               limits = c(-10, 10),
                               breaks = seq(-10, 10, 5)) +
            theme(axis.text.x = element_text(face = "bold", size = 10),
                  axis.text.y = element_text(face = "bold", size = 10),
                  axis.title.x = element_text(face = "bold", size = 12),
                  axis.title.y = element_text(face = "bold", size = 12))
    })
    
    output$reactiveAccPlot <- renderPlot({
        
        plotdf.abnlacc <- ddply(.data = df.abnlaccsichuang,
                                .variables = .(speedBottom),
                                .fun = summarise,
                                number = length(abnlAccTag))
        
        ggplot(data = plotdf.abnlacc, aes(x = speedBottom, y = number)) +
            geom_bar(stat = "identity", fill = "red") +
            scale_x_continuous(name = "Speed Zone, km/h",
                               breaks = seq(0, 90, 10),
                               labels = c("0-10", "10-20", "20-30", "30-40",
                                          "40-50", "50-60", "60-70", "70-80",
                                          "80-90", "90-100")) +
            scale_y_continuous(name = "Number of Abnormal Behavior",
                               limits = c(0, 100),
                               breaks = seq(0, 100, 20)) +
            theme(axis.text.x = element_text(face = "bold", size = 10),
                  axis.text.y = element_text(face = "bold", size = 10),
                  axis.title.x = element_text(face = "bold", size = 12),
                  axis.title.y = element_text(face = "bold", size = 12))
    })
    
    output$reactiveDacPlot <- renderPlot({
        
        plotdf.abnldac <- ddply(.data = df.abnldacsichuang,
                                .variables = .(speedBottom),
                                .fun = summarise,
                                number = length(abnlDacTag))
        
        ggplot(data = plotdf.abnldac, aes(x = speedBottom, y = number)) +
            geom_bar(stat = "identity", fill = "blue") +
            scale_x_continuous(name = "Speed Zone, km/h",
                               breaks = seq(0, 90, 10),
                               labels = c("0-10", "10-20", "20-30", "30-40",
                                          "40-50", "50-60", "60-70", "70-80",
                                          "80-90", "90-100")) +
            scale_y_continuous(name = "Number of Abnormal Behavior",
                               limits = c(0, 100),
                               breaks = seq(0, 100, 20)) +
            theme(axis.text.x = element_text(face = "bold", size = 10),
                  axis.text.y = element_text(face = "bold", size = 10),
                  axis.title.x = element_text(face = "bold", size = 12),
                  axis.title.y = element_text(face = "bold", size = 12))
    })
}


shinyApp(ui, server)

