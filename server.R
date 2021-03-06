#
# Define the interactive functions that power the app and connect to the UI
#
# (c) 2018 Jean-Olivier Irisson, GNU General Public License v3

server <- function(input, output, session) {

  # Subset data based on the selected date range
  # full data
  filter_d <- reactive({
    filter(d, classif_date >= input$dates[1], classif_date <= input$dates[2])
  })
  # data with concensus identifications and location + metadata
  filter_dmeta <- reactive({
    filter(dmeta, classif_date >= input$dates[1], classif_date <= input$dates[2])
  })

  # Inpect classifications through time
  output$p_time <- renderPlot({
    # count classifs per day
    dd <- count(filter_d(), classif_date)
    # plot them
    p <- ggplot(dd) +
      geom_path(aes(classif_date, n)) +
      labs(y="Nb of classifications per day", x="Date")
    # possibly with log transformation of y
    if (input$log) {
      p <- p + scale_y_log10()
    }
    p
  })

  output$t_stats <- renderText({
    # subset data
    dd <- filter_d()
    dm <- filter_dmeta()
    # compute stats
    t <- paste0(
      "Number of classifications: ", nrow(dd), "\n",
      "Number of frames: ", length(unique(dd$image_name)), "\n",
      "Number of identified organisms: ", nrow(dm), "\n"
    )
    t
  })

  # Inpect users
  data_who <- reactive({
    filter_d() %>%
      # count classifs per user
      group_by(user_name) %>% summarise(n=n()) %>% ungroup() %>%
      # shorten user names
      mutate(user=str_trunc(user_name, 15)) %>%
      # sort in descending order
      arrange(desc(n)) %>% mutate(user=factor(user_name, levels=user_name, labels=user)) %>%
      # subset by minimum number of classif
      filter(n>input$minn)
  })
  output$p_who <- renderPlot({
    dw <- data_who()
    # plot result
    ggplot(dw) +
      geom_col(aes(x=user, y=n)) +
      theme(axis.text.x=element_text(angle=30, hjust=1)) +
      labs(y="Nb of classifications", x="User name")
  })
  output$t_who <- renderText({
    dd <- filter_d()
    dw <- data_who()
    tot <- nrow(dd)
    by_these <- sum(dw$n)
    paste0(by_these, " classifications made by these users, which is ", round(by_these/tot*100), "% of total")
  })

  # Inspect effort (number of users per image)
  output$p_effort <- renderPlot({
    dd <- filter_d() %>%
      # count number of users per frame, note if frame is empty or not
      group_by(image_name) %>% summarise(
        n_users=length(unique(user_name)),
        empty=all(is.na(species))
      ) %>% ungroup() %>%
      # pre-compute barchart
      mutate(n_users=factor(n_users)) %>%
      count(empty, n_users)

    # plot result
    ggplot(dd) + geom_col(aes(x=n_users, y=n)) +
      facet_wrap(~empty, labeller=label_both) +
      labs(x="Nb of users per image", y="Nb of images")
  })

  # Inspect results
  # compute height of plot dynamically, to accomodate several species
  p_transects_height <- reactive({
    length(input$species) * 200 + 50
  })
  # plot distribution of various species in space
  output$p_transects <- renderPlot({
    dd <- filter(filter_dmeta(), transect %in% input$transects, species %in% input$species)
    # localisation of classifications
    p <- ggplot(dd, aes(x=dist_from_shore, y=-Depth.m)) + facet_grid(species~transect)
    if (input$geom %in% c("Points", "Both")) {
      p <- p + geom_point(shape=16, size=1, alpha=0.5)
    }
    if (input$geom %in% c("Density", "Both")) {
      p <- p + geom_density2d()
    }
    p
  }, height=p_transects_height)
  # compute stats for classifications
  output$t_transects_t <- renderTable({
    dd <- filter_dmeta()
    dd %>% group_by(transect) %>% summarise(n_organisms=n(), n_images=length(unique(image_name))) %>% ungroup()
  }, striped=T, spacing="xs")
  output$t_transects_s <- renderTable({
    dd <- filter_dmeta()
    dd %>% group_by(species) %>% summarise(n_organisms=n(), n_images=length(unique(image_name))) %>% ungroup()
  }, striped=T, spacing="xs")
}
