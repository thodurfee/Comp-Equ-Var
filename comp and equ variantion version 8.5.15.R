 xname <-  "calculators"
 yname <-  "zebras"
 a1 <- .5
 a2 <- .5
 p1 <- 1
 p2 <- 2
 u <- 1000
 w <- 10
 b <- 1
 method <- "m"
#  leading <- NULL
#  lagging <- NULL
 w_color <- "#aa0000aa"
 u_color <- "#0000aaaa"
 w_type <- "solid"
 u_type <- "solid"
 u_thick <- 2
 w_thick <- 2
 x_star_type <- 8
 x_star_size <- 8
 dp1 <- 1
 dp2 <- 1.1
 variant_color <- "#00000077"
 variant_type <- "#dotted"
 variant_thick <- 1

# library(ggplot2)


Cobb_Douglas <- function(xname
                         , yname
                         , a1
                         , a2
                         , b
                         , p1
                         , p2
                         , u
                         , w
                         , method
                         , leading
                         , lagging
                         , w_color
                         , u_color
                         , w_type
                         , u_type
                         , u_thick
                         , w_thick
                         , x_star_type
                         , x_star_size
                         , x_star_color
                         , dp1
                         , dp2
                         , variant_color
                         , variant_type
                         , variant_thick){

  require(ggplot2)
  ###################################################
  #define functions derived from first order conditions

  m_lagrange <- function(x1,a1,x2,a2,lambda,w,p1,p2,b){
    (((x1)^a1)*((x2)^a2)*b) + lambda*(w - p1*x1 - p2*x2)
  }

  h_lagrange <- function(x1,a1,x2,a2,lambda,w,p1,p2,b){
    (p1*x1 + p2*x2) + lambda*(u - ((x1^a1)*(x2^a2))*b)
  }

  x1_con_d <- function(a1,a2,p1,p2,x2){
    (a1/a2)*(p2/p1)*x2
  }

  x2_con_d <- function(a1,a2,p1,p2,x1){
    (a2/a1)*(p1/p2)*x1
  }

  m_x1_unco <- function(a1,a2,w,p1,p2){
    (a1/a2)*(w/p1)*(1/(1+(a1/a2)))
  }

  m_x2_unco <- function(a1,a2,w,p1,p2){
    (a2/a1)*(w/p2)*(1/(1+(a2/a2)))
  }

  h_x1_unco <- function(a1,a2,p1,p2,u){
    ((a1/a2)*(p2/p1)*(u^(1/a2)))^(1/(1+(a1/a2)))
  }

  h_x2_unco <- function(a1,a2,p1,p2,u){
    ((a2/a1)*(p1/p2)*(u^(1/a1)))^(1/(1+(a2/a1)))
  }

  w_x1 <- function(p1,p2,x2,w){
    (-p2/p1)*x2 + (w/p1)
  }

  u_x1 <- function(u,p1,p2,a1,a2,b){
    (u/(b*(x2^a2)))^(1/a1)
  }

  w_x2 <- function(p1,p2,x1,w){
    (-p1/p2)*x1+(w/p2)
  }

  u_x2 <- function(u,x1,a1,a2,b){
    (u/(b*(x1^a1)))^(1/a2)
  }

  x_intercept <- function(w,p1){
    w/p1
  }

  y_intercept <- function(w,p2){
    w/p2
  }

  find_u <- function(x1,x2,a1,a2,b){
    (x1^a1)*(x2^a2)*b
  }

  find_w <- function(x1, x2, p1, p2){
    (p1*x1) + (p2*x2)
  }


  ###################################################
  ##define defaults


  if(  missing( method )   ){
    method <- "m"
  }

  if( missing(xname) ){
    xname <- "x.1"
  }

  if( missing(yname)  ){
    yname <- "x.2"
  }

  if( missing(w_color) ){
    w_color <- "#33aa33ff"
  }

  if( missing(w_type) ){
    w_type <- "solid"
  }

  if( missing(w_thick)  ){
    w_thick <- 2
  }

  if( missing(u_color)  ){
    u_color <- "#aa3333ff"
  }

  if( missing(u_type)  ){
    u_type <- "solid"
  }

  if( missing(u_thick) ){
    u_thick <- 2
  }

  if( missing(x_star_type)  ){
    x_star_type <- 8
  }

  if( missing(x_star_size)  ){
    x_star_size <- 6
  }

  if( missing(leading)  ){
    leading <- 0
  }

  if( missing(lagging)  ){
    lagging <- 0
  }

  if( missing(b)  ){
    b <- 1
  }

  if( missing(a1)  ){
    return("define a1 first")
  }

  if( missing(a2)  ){
    return("define a2 first")
  }

  if( missing(p1)  ){
    return("define prices")
  }

  if( missing(p2)  ){
    return("define prices")
  }

  if( missing(dp2)  ){
    dp2 <- 1
  }

  if( missing(dp1)  ){
    dp1 <- 1
  }

  if( missing(variant_color)  ){
    vairant_color <- "#aa333377"
  }

  if( missing(variant_type)  ){
    variant_type <- "dotted"
  }



  #################################
  ####################################
  ###################################
  ####################################
  ###################################################
  #find if I am using the hicks or the marshall method

  hicks_group <- c("Hicks", "hicks", "h", "H", 2, "hikcs", "hkisc", "hi", "hic")

  marshall_group <- c("marshall", "marshal", "Marshall", "Marshal", "m", "M",1, "masrlal", "mlsahr", "mar")


  verify_method <- function(in_method){
    for(i in 1:length(hicks_group)){
      if(hicks_group[i] == in_method){
        return("hicks")
      }
    }
    for(j in 1:length(marshall_group)){
      if(marshall_group[j] == in_method){
        return("marshall")
      }
    }
    return("marshall")
  }


  method <- verify_method(in_method = method)



  ###################################################
  #calculate x1 and x2 star with the hicks method
  find_x_star <- function(in_method){
    if(in_method == "hicks"){
      x1_star <- h_x1_unco(a1 = a1,
                           a2 = a2,
                           p1 = p1,
                           p2 = p2,
                           u = u)

      x2_star <- h_x2_unco(a1 = a1,
                           a2 = a2,
                           p1 = p1,
                           p2 = p2,
                           u = u)

      return(data.frame(x1_star, x2_star))
    }

    else if(in_method == "marshall"){

      x1_star <- m_x1_unco(a1 = a1,
                           a2 = a2,
                           p1 = p1,
                           p2 = p2,
                           w = w)

      x2_star <- m_x2_unco(a1 = a1,
                           a2 = a2,
                           p1 = p1,
                           p2 = p2,
                           w = w)

      return(data.frame(x1_star, x2_star))
    }
    else{
      return("improper definition of 'method'")
    }
  }


  xstar <- find_x_star(in_method = method)


  ###################################################
  #find w and u levels
  util <- find_u( x1 = as.numeric(xstar["x1_star"]),
                  x2 = as.numeric(xstar["x2_star"]),
                  a1 = a1,
                  a2 = a2,
                  b = b)

  wealth <- find_w( x1 = as.numeric(xstar["x1_star"]),
                    x2 = as.numeric(xstar["x2_star"]),
                    p1 = p1,
                    p2 = p2)



  yintercept <- y_intercept( w = wealth,
                             p2 = p2)

  xintercept <- x_intercept(w = wealth,
                            p1 = p1)



  xspace <- seq(from = 0, to = xintercept, by = .1)
  yspace <- seq(from = 0, to = yintercept, by = .1)




  ###################################################
  #build budget curve

  budget_curve <- sapply(xspace,
                         FUN = function(xspace) w_x2(x1 = xspace,
                                                     p1 = p1,
                                                     p2 = p2,
                                                     w = wealth)
  )

  bspace <- cbind.data.frame(budget_curve, xspace)



  names(bspace) <- c("yspace","xspace")



  ###################################################
  #build utility curve
  utility_curve <- sapply(xspace,
                          FUN = function(xspace) u_x2(x1 = xspace,
                                                      a1 = a1,
                                                      a2 = a2,
                                                      u = util,
                                                      b = b)
  )

  uspace <- cbind.data.frame(utility_curve, xspace)


  names(uspace) <- c("yspace","xspace")



  ######################################
  #####################################
  #####################################
  ###################################
  #now to add the vairants

  find_x_star2 <- function(in_method){
    if(in_method == "hicks"){
      x1_star <- h_x1_unco(a1 = a1,
                           a2 = a2,
                           p1 = p1*dp1,
                           p2 = p2*dp2,
                           u = u)

      x2_star <- h_x2_unco(a1 = a1,
                           a2 = a2,
                           p1 = p1*dp1,
                           p2 = p2*dp2,
                           u = u)

      return(data.frame(x1_star, x2_star))
    }

    else if(in_method == "marshall"){

      x1_star <- m_x1_unco(a1 = a1,
                           a2 = a2,
                           p1 = p1*dp1,
                           p2 = p2*dp2,
                           w = w)

      x2_star <- m_x2_unco(a1 = a1,
                           a2 = a2,
                           p1 = p1*dp1,
                           p2 = p2*dp2,
                           w = w)

      return(data.frame(x1_star, x2_star))
    }
    else{
      return("improper definition of 'method'")
    }
  }


  xstar2 <- find_x_star2(in_method = method)


  ###################################################
  #find w and u levels
  util2 <- find_u( x1 = as.numeric(xstar["x1_star"]),
                  x2 = as.numeric(xstar["x2_star"]),
                  a1 = a1,
                  a2 = a2,
                  b = b)

  wealth2 <- find_w( x1 = as.numeric(xstar["x1_star"]),
                    x2 = as.numeric(xstar["x2_star"]),
                    p1 = p1*dp1,
                    p2 = p2*dp2)



  yintercept2 <- y_intercept( w = wealth,
                             p2 = p2*dp2)

  xintercept2 <- x_intercept(w = wealth,
                            p1 = p1*dp1)



  xspace2 <- seq(from = 0, to = xintercept2, by = .1)
  yspace2 <- seq(from = 0, to = yintercept2, by = .1)




  ###################################################
  #build budget curve

  budget_curve2 <- sapply(xspace2,
                         FUN = function(xspace2) w_x2(x1 = xspace2,
                                                     p1 = p1*dp1,
                                                     p2 = p2*dp2,
                                                     w = wealth2)
                         )

  bspace2 <- cbind.data.frame(budget_curve2, xspace2)



  names(bspace2) <- c("yspace","xspace")



  ###################################################
  #build utility curve
  utility_curve2 <- sapply(xspace2,
                          FUN = function(xspace2) u_x2(x1 = xspace2,
                                                      a1 = a1,
                                                      a2 = a2,
                                                      u = util2,
                                                      b = b)
                          )

  uspace2 <- cbind.data.frame(utility_curve2, xspace2)


  names(uspace2) <- c("yspace","xspace")


  #################################
  #################################
  ################################
  ###################################################
  #graph all of the things

  wbase <- ggplot( data = bspace,
                   aes( x =  xspace,
                        y =  yspace
                   )

  )

  layer1 <- wbase + geom_line(data = bspace,
                              color = w_color,
                              linetype = w_type,
                              size = w_thick
  )


  layer2 <- layer1 + labs( title = paste(xname, " vs ", yname),
                           x = xname,
                           y = yname)

  layer3 <- layer2 + geom_hline() + geom_vline()

  layer4 <-  layer3 + geom_line(data = uspace,
                                color = u_color,
                                aes( x =  xspace,
                                     y =  yspace
                                ),
                                linetype = u_type,
                                size = u_thick
  )


  layer5 <- layer4# + xlim(0,xintercept*1.25) + ylim(0,yintercept*1.25)

  layer6 <- layer5+ geom_point(data = xstar,
                               aes(x = x1_star,
                                   y = x2_star),
                               size = x_star_size,
                               shape = x_star_type,
                               color = x_star_color)

  layer1.2 <- layer6 + geom_line(data = bspace2,
                                 color = variant_color,
                                 linetype = variant_type,
                                 size = variant_thick,
                                 aes(x = xspace2,
                                     y = yspace2))



  length(yspace)
  length(yspace2)

  layer4.2 <- layer1.2 + geom_line(data = uspace2,
                                   color = variant_color,
                                   aes(x = xspace2,
                                       y = yspace2),
                                   linetype = variant_type,
                                   size = variant_thick)

  layer5.2 <- layer4.2 + xlim(0, max(xintercept*1.25, xintercept2*1.25)) + ylim(0, max(yintercept*1.25, yintercept2*1.25))

  layer6.2 <- layer5.2 + geom_point(data = xstar2,
                                    aes(x = x1_star,
                                        y = x2_star),
                                    size = x_star_size,
                                    shape = x_star_type)


  layerf <- layer6.2



  print(paste(as.numeric(xstar$x1_star), "= x1 Star"))
  print(paste(as.numeric(xstar$x2_star), "= x2 Star"))
  print(paste(util, "= Utility"))
  print(paste(wealth, "= Wealth"))

  return(layerf)
}





Cobb_Douglas( w = 2000,
              p1 = 1,
              p2 = 1,
              u = 10,
              a1 = .5,
              a2 = .5,
              b = 2,
              method = "m",
              xname = "interns",
              yname = "printers",
              leading = NULL,
              lagging = NULL,
              w_color = "#00aa00aa",
              u_color = "#aa0000aa",
              u_type = "11",
              w_type = "solid",
              x_star_type = 8,
              x_star_size = 6,
              x_star_color = "black",
              u_thick = 2,
              w_thick = 2,
              dp1 = 1,
              dp2 = 2,
              variant_type = "solid",
              variant_color = "#00000088",
              variant_thick = 2)


