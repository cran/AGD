# extractLMS, ageGrid

extractLMS <- function(fit, data, sex="M", grid="classic",
            decimals = c(4,4,4), 
            flatAge = NULL)
{
  # Extracts the LMS table after the 'Cole-transformation'
  # or any other transformation in t.age
  # fit  final gamlss object
  # data should contain both age and t.age
  
  check.names(df=data, needed=c("age","t.age"))
  if (!is.gamlss(fit)) stop("fit not a gamlss object.")
  tm <- data$t.age[which.min(data$age)]
  # if (abs(tm)>0.0001) stop("wrong offset of transformed age")
  # if (min(data$t.age) < 0) warning("Negative transformed age found. Results are unpredictable.")

  grd <- ageGrid(grid)
  grid.age <- grd$year

  minage <- min(data$age, na.rm=TRUE)
  maxage <- max(data$age, na.rm=TRUE)
  outside <- grid.age < minage | grid.age > maxage
  grid.age <- grid.age[!outside]

  t.grid.age <- approx(x=data$age, y=data$t.age, xout=grid.age, ties=mean)$y
  if (length(t.grid.age)==0) stop("No overlap between age grid and data.")
  newdata <- data.frame(t.age=t.grid.age)

  # lms <- predictAll(fit, newdata=newdata, data=data)
  # print(lms)
  lms <- predictAll(fit, newdata=newdata, data=data)
  lms$mu    <- round(lms$mu, decimals[2])
  lms$sigma <- round(lms$sigma/lms$mu, decimals[3])
  if (length(lms)>2) lms$nu <- round(lms$nu, decimals[1])
  else lms$nu <- 1
  
  lms <- as.data.frame(lms)
  result <- data.frame(sex=sex, x=grd$year, L=NA, M=NA, S=NA)
  result[!outside, "L"] <- lms$nu
  result[!outside, "M"] <- lms$mu
  result[!outside, "S"] <- lms$sigma
  
  if (!is.null(flatAge)) {
    if (!(flatAge %in% result$x)) stop("FlatAge value (', FlatAge,') not found in age grid")
    flatLMS <- result[flatAge==result$x, c("L","M","S")]
    result[!outside & flatAge<result$x, c("L","M","S")] <- flatLMS
  }
  
  return(result)
}

ageGrid <- function(grid="compact"){
  formats <- c("compact", "classic", "extensive", "0-104w", "0-24m", "0-21y", "0-21yd", "0-21yc")
  fmi <- pmatch(grid, formats)
  if (is.na(fmi)) stop("Grid format ",grid," unknown.")
  grid <- switch(fmi,
    compact =
      c((0:14)/365.25,
        (3:13)*7/365.25,
        seq(3,11.5,0.5)/12,
        (12:23)/12,
        seq(2, 21, 0.5)
      ),
    classic =
    {
      grid.weeks <- c(0,1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,22,24,26,
                      28,32,36,40,44,48,52,56,60,64)
      c(grid.weeks*7/365.25,seq(1.5,21,0.5))
    },
    extensive = (0:(365.25*21+1))/365.25,
    week = (0:104)*7/365.25,
    month = (0:24)/12,
    year = 0:21,
    dyear = seq(0, 21, 0.1),
    cyear = seq(0, 21, 0.01)
  )
  
  year    <- round(grid, 4)
  month   <- round(grid*12, 4)
  week    <- round(grid*365.25/7, 4)
  day     <- round(grid*365.25, 4)
  return(list(format = formats[fmi],
   year = year, month = month, week = week, day = day))
}


