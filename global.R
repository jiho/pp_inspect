#
# Preprocess PlanktonPortal classification log to ease plotting
#
# (c) 2018 Jean-Olivier Irisson, GNU General Public License v3

library("tidyverse")
library("lubridate") # for manipulation of dates

library("ggrepel")   # for intelligent label plotting
library("dbscan")    # for HDBSCAN clustering of positions

library("shiny")     # for the web app
library("shinythemes")


message("Load data")

# list data files
# tar.gz files sent by zooniverse should be put in "classifications"
tar_files <- list.files("classifications", pattern=glob2rx("*plankton_classifications.csv.tar.gz"), full.names=TRUE)

# get latest one
latest_tar <- tail(tar_files, 1)
latest_rdata <- str_replace(latest_tar, fixed("tar.gz"), "Rdata")

# if it has been processed already, read the result
if (file.exists(latest_rdata)) {
  load(latest_rdata)

# if not, process it
} else {
  ## Read data
  message("Read csv file")
  untar(latest_tar)
  latest_csv <- str_replace(latest_tar, fixed(".tar.gz"), "") %>% basename()
  # read the columns of interest of the latest file
  dfull <- read_csv(latest_csv,
    col_types=cols_only(
      created_at=col_datetime(format="%Y-%m-%d %H:%M:%S %Z"),
      image_name=col_character(),
      user_name=col_character(),
      species=col_character(),
      p0_x=col_double(), p0_y=col_double()
    ))
  file.remove(latest_csv)

  ## Pre-process data ----
  message("Extract Mediterranean data")
  d <- dfull %>%
    # filter Mediterranean data
    filter(str_detect(image_name, "2013-07")) %>%
    # rename some columns
    rename(x=p0_x, y=p0_y, classif_datetime=created_at) %>%
    # add some columns
    mutate(
      # round classification date
      classif_date=as.Date(classif_datetime),
      # get capture date+time from name
      image_datetime=ymd_hms(str_sub(image_name, 1, 19))
    )

  message("Precompute some counts")
  # compute number of classifications per user
  dwho <- d %>%
    group_by(user_name) %>% summarise(n=n()) %>% ungroup() %>%
    arrange(desc(n)) %>% mutate(user_name=factor(user_name, levels=user_name))

  # add a score inversely proportional to the number of classif per user
  d <- left_join(d, mutate(dwho, score=sqrt(n)) %>% select(user_name, score), by="user_name")

  message("Compute concensus identifications")
  # compute classification agreement per frame

  # Get majority vote
  maj <- function(class) {
    table(class) %>% which.max() %>% names()
  }

  # Get majority vote weighted by score
  maj_score <- function(class, score) {
    by_score <- data.frame(class, score) %>% group_by(class) %>% summarise(score=sum(score))
    by_score$class[which.max(by_score$score)]
  }

  # prepare multicore computation
  library("doParallel")
  registerDoParallel(cores=detectCores()-1)
  # mutlicore
  system.time(dcons <- plyr::ddply(d, ~image_name, .parallel=TRUE, .fun=function(xf) {

  # # single core
  # system.time(dcons <- d %>% group_by(image_name) %>% do({
  #   xf <- .

    # remove points with no label
    x <- filter(xf, !is.na(species), !is.na(x), !is.na(y))

    if (nrow(x) > 3) {
      # cluster classifications to detect groups of clicks
      k <- hdbscan(select(x, x, y), minPts=3)
      x$k <- k$cluster

      # get cluster centers and majority classification for each cluster
      # = this is the retained identification
      centers <- x %>%
        # cluster 0 = isolated points => remove them
        filter(x, k > 0) %>%
        # for each cluster
        group_by(k) %>% summarise(
          # basic identification info
          image_name=image_name[1], image_datetime=image_datetime[1],
          # center
          x=mean(x), y=mean(y),
          # majority classification
          # species_maj=maj(species),
          # species_max=maj_score(species, score)
          # -> majority with score is better
          species=maj_score(species, score),
          # date of identification
          classif_date=max(classif_date, na.rm=T)
        )

      # # check with a plot
      # ggplot() +
      #   # individual votes
      #   geom_point(aes(x, -y, shape=factor(k), colour=str_trunc(user_name, 15)), size=1, data=x) +
      #   geom_text_repel(aes(x, -y, colour=str_trunc(user_name, 15), label=species), data=x, size=3, segment.size=0.25, segment.alpha=0.5, min.segment.length=0.1) +
      #   # concensus
      #   geom_point(aes(x=xc, y=-yc), data=centers) +
      #   geom_text(aes(x=xc, y=-yc, label=species_max), data=centers, size=3, hjust="left", vjust=0.2, nudge_x=5) +
      #   # beautify
      #   coord_fixed(xlim=c(0, 1024), ylim=c(0, -560)) +
      #   labs(shape="cluster", colour="user")
    } else {
      centers <- data.frame(image_name=xf$image_name[1])
    }
    centers
  }) %>% ungroup())

  message("Add transect data")
  # NB: Assumes a directory named "transects", with sub-directories named
  #     after each transect name, and, in each, a file named `isiis.csv`
  #     that contains columns:
  #     - dateTime: date and time (to the second) of the data recorded by ISIIS
  #     - one column per data field (lon, lat, depth, temperature, etc.)

  # read transect data
  transect_files <- list.files("transects", pattern="isiis.csv", rec=T, full=T)
  t <- map_df(transect_files, function(f) {
      x <- read_csv(f, col_types=cols())
      x$transect <- str_split(f, "/")[[1]][2]
      return(x)
    }) %>%
    select(transect, dateTime, Depth.m, dist_from_shore, lat, lon, Pressure.dbar:Density) %>%
    filter(dateTime >= min(dcons$image_datetime, na.rm=T), dateTime <= max(dcons$image_datetime, na.rm=T)) %>%
    group_by(transect, dateTime) %>% summarise_all(mean, na.rm=T) %>% ungroup()

  # associate transect data with classifications by finding the nearest time in transects for each classification
  idx <- round(approx(t$dateTime, 1:nrow(t), xo=dcons$image_datetime, rule=2)$y)
  dmeta <- cbind(dcons, t[idx,] %>% select(-dateTime))
  # remove cases with no species identified
  dmeta <- filter(dmeta, !is.na(species))

  message("Save result")
  # save result
  save(d, dmeta, file=latest_rdata)
  # cleanup
  rm(dfull, dwho)
}


## Prepare app elements ----

# plot theme
theme_set(theme_gray(14))

# precompute some lists of elements for the UI
transects <- sort(unique(dmeta$transect))
species <- sort(unique(dmeta$species))
dates <- range(d$classif_date)
