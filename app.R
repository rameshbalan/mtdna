library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("mtDNA Introgression"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            sliderInput("max",
                        "Number of Generations",
                        min = 2,
                        max = 1000,
                        value = 100),
            sliderInput("N",
                        "Population Size",
                        min = 500,
                        max = 1000000,
                        value = 500000),
            sliderInput("sims",
                        "Number of Simulations",
                        min = 1,
                        max = 1000,
                        value = 10),
            sliderInput("loss1",
                        "Mutation rate of MEDEA(M1)",
                        min = 0.1,
                        max = 1,
                        value = 0.1),
            sliderInput("loss4",
                        "Mutation rate of MEDEA(M4)",
                        min = 0.1,
                        max = 1,
                        value = 0.1),
            sliderInput("s",
                        "Selection for MEDEA(M1)",
                        min = 0.1,
                        max = 1,
                        value = 0.12),
            sliderInput("t",
                        "Selection for MEDEA(M4)",
                        min = 0.1,
                        max = 1,
                        value = 0.12),
            sliderInput("q",
                        "Selection for mtDNA",
                        min = 0.1,
                        max = 1,
                        value = 0.99),
            sliderInput("bottleneck",
                        "Bottleneck strength applied every 7th and 10th generation",
                        min = 0.1,
                        max = 1,
                        step = 0.00001,
                        value = 0.99995),
            sliderInput("max.pop.inc",
                        "Max. Increase factor for Population",
                        min = 1,
                        max = 200,
                        value = 100)
            
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("distPlot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    source("functions.R")
    output$distPlot <- renderPlot({
        #result for all simulations
        major.result1 <- data.frame(matrix(, input$sims, input$max))
        major.result4 <- data.frame(matrix(, input$sims, input$max))
        major.resultM <- data.frame(matrix(, input$sims, input$max))
        
        for(i in 1:input$sims){
            ## SETUP FOR SIMULATION
            # we will assume that initially are population isn't in a bottleneck
            Ne <- input$N
            male_pop <- vector(length=32, mode="numeric")
            female_pop <- vector(length=32, mode="numeric")
            male_pop <- c(0, 0, 0, 0, 0, 62500, 0, 0,
                          0, 0, 62500, 0, 0, 0, 0, 0,
                          0 ,0 ,0 ,0 ,0 ,62500 ,0 ,0,
                          0 ,0 ,62500 ,0 ,0 ,0 ,0, 0)
            
            female_pop <- c(0, 0, 0, 0, 0, 62500, 0, 0,
                            0, 0, 62500, 0, 0, 0, 0, 0,
                            0 ,0 ,0 ,0 ,0 ,62500 ,0 ,0,
                            0 ,0 ,62500 ,0 ,0 ,0 ,0, 0)
            
            names(male_pop) <- c("ABABM", "ABAbM", "ABaBM", "ABabM", "AbABM", "AbAbM", "AbaBM", "AbabM",
                                 "aBABM", "aBAbM", "aBaBM", "aBabM", "abABM", "abAbM", "abaBM", "ababM", 
                                 "ABABm", "ABAbm", "ABaBm", "ABabm", "AbABm", "AbAbm", "AbaBm", "Ababm",
                                 "aBABm", "aBAbm", "aBaBm", "aBabm", "abABm", "abAbm", "abaBm", "ababm")
            
            names(female_pop) <- c("ABABM", "ABAbM", "ABaBM", "ABabM", "AbABM", "AbAbM", "AbaBM", "AbabM",
                                   "aBABM", "aBAbM", "aBaBM", "aBabM", "abABM", "abAbM", "abaBM", "ababM", "ABABm",
                                   "ABAbm", "ABaBm", "ABabm", "AbABm", "AbAbm", "AbaBm", "Ababm",
                                   "aBABm", "aBAbm", "aBaBm", "aBabm", "abABm", "abAbm", "abaBm", "ababm")
            
            # gamete table
            male_gametes <- vector(length=8, mode="numeric")
            female_gametes <- vector(length=18, mode="numeric")
            names(male_gametes) <- c("ABM", "AbM", "aBM", "abM",
                                     "ABm", "Abm", "aBm", "abm")
            names(female_gametes) <- c("ABM", "AbM", "aBM", "abM", "Ab'M", "a'BM", "ab'M", "a'bM", "a'b'M",
                                       "ABm", "Abm", "aBm", "abm", "Ab'm", "a'Bm", "ab'm", "a'bm", "a'b'm")
            
            # a matrix to store results of a single simulation
            results <- data.frame(matrix(, 32, input$max))
            row.names(results) <- names(female_pop)
            
            # increment variable
            counter <- 1
            
            # record starting conditions
            results[, counter] <- 2*male_pop
            
            ## Simulation
            while(counter < input$max){
                counter <- counter + 1
                male_pop[1:32] <- male_selection(male_pop,
                                                 N = Ne, input$s)
                male_gametes[1:8] <- makeGametes_male(male_pop,
                                                      r = 0.5,
                                                      N = Ne,
                                                      U.a = input$loss1,
                                                      U.b = input$loss4
                )
                female_pop[1:32] <- female_selection(female_pop,
                                                     N = Ne, input$s)
                female_gametes[1:18] <- makeGametes_female(female_pop,
                                                           r = 0.5,
                                                           N = Ne,
                                                           U.a = input$loss1,
                                                           U.b = input$loss4
                )
                # this is the bottleneck occuring every 7 and 10 generations
                if(counter %% 10 == 7 | counter %% 10 == 0 ){
                    Ne <- Ne - (input$bottleneck * Ne)
                }else{
                    if(Ne < input$N) Ne <- Ne * input$max.pop.inc
                    if(Ne > input$N) Ne <- input$N
                }
                # perform the bottleneck
                pop <- makePop(male_gametes, female_gametes, N = Ne, s = input$s, q = input$q, t = input$t)
                results[, counter] <- pop
                male_pop <- pop/2
                female_pop <- pop/2
            }
            # storing the results
            major.result1[i, 1:input$max] <- mon1.fnx(results=results, max=input$max)
            major.result4[i, 1:input$max] <- mon4.fnx(results=results, max=input$max)
            major.resultM[i, 1:input$max] <- mtDNA.fnx(results=results, max=input$max)
            
        }
        
        #Plotting frequencies of M1/M4 vs M
        layout(matrix(c(1,2),nrow=1), width=c(5,1)) 
        par(mar=c(5,4,4,0)) #No margin on the right side
        matplot(t(major.resultM), type="l", col = "#ff7f7f", lty = 1, ylim = range(0:1), xlab="Generations", ylab="Frequency")
        par(mar=c(5,1,4,0)) #No margin on the left side
        matlines(t(major.result1),lty = 1, col = "#7fff7f", alpha = 0.5)
        matlines(t(major.result4),lty = 1, col = "#7F7FFF")
        plot(c(0,1),type="n", axes =F,xlab="", ylab="")
        legend("center",legend=c("mTDNA","M1","M4"), col = c("#ff7f7f","#7fff7f","#7F7FFF"),lty=1)
        
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
