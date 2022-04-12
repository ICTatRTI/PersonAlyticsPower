

aNBI <- .arima.sim(12345,
                   model = list(ar=.1, ma=0),
                   n = 1000,
                   rand.gen = rNBI,
                   mu = 10, sigma = 0.5)

hist(floor(aNBI))
plot(floor(aNBI), aNBI)
cor(floor(aNBI), aNBI)
hist(scale(floor(aNBI)))
plot(scale(floor(aNBI)), aNBI)
cor(scale(floor(aNBI)), aNBI)


aNBII <- .arima.sim(12345,
                    model = list(ar=.1, ma=0),
                    n = 1000,
                    rand.gen = rNBII,
                    mu = 10, sigma = 0.5)

aLO <- .arima.sim(12345,
                  model = list(ar=.1, ma=0),
                  n = 1000,
                  rand.gen = rLO,
                  mu = 10, sigma = 0.5)

aBI <- .arima.sim(12345,
                  model = list(ar=.1, ma=0),
                  n = 1000,
                  rand.gen = rBI,
                  mu = .34)

# beta doesn't work
aBE <- .arima.sim(12345,
                   model = list(ar=.1, ma=0),
                   n = 1000,
                   rand.gen = BE,
                   mu = .5, sigma = 0.5)
head(table(aBE))




arima.sim <- function (model, n, rand.gen = rnorm, innov = rand.gen(n, ...),
          n.start = 10, start.innov = rand.gen(n.start, ...), ...)
{
  if (!is.list(model))
    stop("'model' must be list")
  if (n <= 0L)
    stop("'n' must be strictly positive")
  p <- length(model$ar)
  if (p) {
    minroots <- min(Mod(polyroot(c(1, -model$ar))))
    if (minroots <= 1)
      stop("'ar' part of model is not stationary")
  }
  q <- length(model$ma)
  if (is.na(n.start))
    n.start <- p + q + ifelse(p > 0, ceiling(6/log(minroots)),
                              0)
  if (n.start < p + q)
    stop("burn-in 'n.start' must be as long as 'ar + ma'")
  d <- 0
  if (!is.null(ord <- model$order)) {
    if (length(ord) != 3L)
      stop("'model$order' must be of length 3")
    if (p != ord[1L])
      stop("inconsistent specification of 'ar' order")
    if (q != ord[3L])
      stop("inconsistent specification of 'ma' order")
    d <- ord[2L]
    if (d != round(d) || d < 0)
      stop("number of differences must be a positive integer")
  }
  if (!missing(start.innov) && length(start.innov) < n.start)
    stop(sprintf(ngettext(n.start, "'start.innov' is too short: need %d point",
                          "'start.innov' is too short: need %d points"), n.start),
         domain = NA)
  x <- ts(c(start.innov[seq_len(n.start)], innov[1L:n]), start = 1 -
            n.start)
  if (length(model$ma)) {
    x <- filter(x, c(1, model$ma), sides = 1L)
    x[seq_along(model$ma)] <- 0
  }
  if (length(model$ar))
    x <- filter(x, model$ar, method = "recursive")
  if (n.start > 0)
    x <- x[-(seq_len(n.start))]
  if (d > 0)
    x <- diffinv(x, differences = d)
  as.ts(x)
}
