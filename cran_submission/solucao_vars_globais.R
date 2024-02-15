x <- 2

func1 <- function(){
  #x <- 3
  x <- 2
  assign("x", 3, envir=parent.frame())
  func2 <- function(){
    x <- 4
    #assign("x", 4, envir=parent.frame())
    print(x)
  }
  func2()
  print(x)
}

func1()
print(x)

s[i] <- 0

s_s <- get("s")
s_s[i] <- 0
assign("s", s_s, envir=parent.frame())




