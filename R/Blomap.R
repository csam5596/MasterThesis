blomap <- function(aa){
  if(is.na(aa)) return(c(NA,NA,NA,NA,NA))
  switch(aa, 
         'A' = {return(c(-0.57, 0.39, -0.96, -0.61, -0.69))},
         'R' = {return(c(-0.4, -0.83, -0.61, 1.26, -0.28))},
         'N' = {return(c(-0.7, -0.63, -1.47, 1.02, 1.06))},
         'D' = {return(c(-1.62, -0.52, -0.67, 1.02, 1.47))},
         'C' = {return(c(0.07, 2.04, 0.65, -1.13, -0.39))},
         
         'Q' = {return(c(-0.05, -1.5, -0.67, 0.49, 0.21))},
         'E' = {return(c(-0.64, -1.59, -0.39, 0.69, 1.04))},
         'G' = {return(c(-0.90, 0.87, -0.36, 1.08, 1.95))},
         'H' = {return(c(0.73, -0.67, -0.42, 1.13, 0.99))},
         'I' = {return(c(0.59, 0.79, 1.44, -1.90, -0.93))},
         
         'L' = {return(c(0.65, 0.84, 1.25, -0.99, -1.90))},
         'K' = {return(c(-0.64, -1.19, -0.65, 0.68, -0.13))},
         'M' = {return(c(0.76, 0.05, 0.06, -0.62, -1.59))},
         'F' = {return(c(1.87, 1.04, 1.28, -0.61, -0.16))},
         'P' = {return(c(-1.82, -0.63, 0.32, 0.03, 0.68))},
         
         'S' = {return(c(-0.39, -0.27, -1.51, -0.25, 0.31))},
         'T' = {return(c(-0.04, -0.30, -0.82, -1.02, -0.04))},
         'W' = {return(c(1.38, 1.69, 1.91, 1.07, -0.05))},
         'Y' = {return(c(1.75, 0.11, 0.65, 0.21, -0.41))},
         'V' = {return(c(-0.02, 0.3, 0.97, -1.55, -1.16))},
         
         'X' = {return(c(0,0,0,0,0))}
  )
}
