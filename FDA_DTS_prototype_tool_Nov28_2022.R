# Shiny App to Explore FDA Drug Trials Snapshots 2015-2021 Data Set
# Author: Ariel Carmeli (carmeli.b.ariel@gmail.com)

library(shiny)
library(tidyverse)
library(lubridate)
library(dplyr)
library(scales) 
library(RColorBrewer)
library(plotly)
library(ggridges)
#library(kableExtra)

########################
### READ IN THE DATA ###
########################

print( "Loading data ..." )
fda_approvals <- read.csv('FDA_Drug_Trials_Snapshots_2015-21.csv')
disease_burden <- read.csv('Disease_burden.csv')
print( "Data loaded!" )

############################
### PROCESS FDA DTS DATA ###
############################

# Change class type of select variables to aid in data processing and visualization
fda_approvals$Enrollment <- as.numeric(as.character(fda_approvals$Enrollment))
fda_approvals$Therapeutic_Area <- as.character(fda_approvals$Therapeutic_Area)
fda_approvals$Brand_Name <- as.character(fda_approvals$Brand_Name)
fda_approvals$United_States <- as.numeric(as.character(fda_approvals$United_States))

# Add columns for non-hispanic, Men, Age under 65 
fda_approvals <- fda_approvals %>% mutate(Non_Hispanic = 100 - Hispanic, .after = Hispanic)
fda_approvals <- fda_approvals %>% mutate(Male = 100 - Female, .after = Female)
fda_approvals <- fda_approvals %>% mutate(Age_under_65 = 100 - Age_65_or_older, .after = Age_65_or_older)

# Create longer version, for plotting
fda_approvals_long <- pivot_longer(fda_approvals, cols = Female:Age_80_or_older, names_to = "Demographic", values_to = "Percentage")

# Change class type of variable in long
fda_approvals_long$Percentage <- as.numeric(as.character(fda_approvals_long$Percentage))

#################################
### CREATE SUMMARY STATISTICS ###
#################################

# Calculate mean representation treating weighting each trial by enrollment size
averages <- fda_approvals_long %>% filter(Percentage != "NA") # create new df with non NA participation

averages <- averages %>% # Calculate number of people (% representation * Trial Enrollment)
    mutate(Population = round((Percentage / 100) * Enrollment))

averages_all_years <- averages %>% # Calculate weighted average (Demographic groups' enrollment / Total Enrollment) - all years
    group_by(Demographic) %>% 
    summarize(Demographic_enrollment = sum(Population), Total_enrollment = sum(Enrollment)) %>% 
    mutate(Weighted_average = round(100* Demographic_enrollment / Total_enrollment),1)

averages <- averages %>% # Same as above but by year
    group_by(Approval_Year, Demographic) %>% 
    summarize(Demographic_enrollment = sum(Population), Total_enrollment = sum(Enrollment)) %>% 
    mutate(Weighted_average = round(100* Demographic_enrollment / Total_enrollment),1)

# Calculate median representation
approval_median <- fda_approvals_long %>% filter(Percentage != "NA") # create new df with non NA participation

approval_median_all_years <- approval_median %>% # Calculate median enrollment percentage across all trials
    group_by(Demographic) %>% 
    summarize(Median = round(median(Percentage),1))

approval_median <- approval_median %>% # Same as above but by year
    group_by(Approval_Year, Demographic) %>% 
    summarize(Median = round(median(Percentage),1))

# Calculate mean representation treating each trial equally
approval_mean <- fda_approvals_long %>% filter(Percentage != "NA") # create new df with non NA participation

approval_mean_all_years <- approval_mean %>% # Calculate mean enrollment percentage across all trials (non weighted by enrollment)
    group_by(Demographic) %>% 
    summarize(Average = round(mean(Percentage),1))

approval_mean <- approval_mean %>% # Same as above but by year
    group_by(Approval_Year, Demographic) %>% 
    summarize(Average = round(mean(Percentage),1))

# Create combined dataset by Demographic
summary_statistics_all_years <- data.frame(averages_all_years$Demographic, approval_median_all_years$Median, approval_mean_all_years$Average, averages_all_years$Weighted_average)
names(summary_statistics_all_years) <- c("Demographic", "Median", "Average", "Weighted_average")

# Create combined dataset by Year and Demographic
summary_statistics <- data.frame(averages$Approval_Year, averages$Demographic, approval_median$Median, approval_mean$Average, averages$Weighted_average)
names(summary_statistics) <- c("Approval_Year", "Demographic", "Median", "Average", "Weighted_average")

#############################
### QUANTIFY MISSING DATA ###
#############################

# Determine demographics to include in data quality assessment
demographics <- c("Female", 
                  "Age_65_or_older", 
                  "Asian",
                  "Black",
                  "Hispanic",
                  "White")

# Calculate total number of approvals by Year
num_approvals <- fda_approvals %>% 
    group_by(Approval_Year) %>% 
    summarise(Count = n())
num_approvals <- num_approvals$Count

# Grab number of non missing values by Demographic and Year
not_missing_values <- fda_approvals_long %>% 
    filter(Percentage != "NA") %>% 
    filter(Demographic %in% demographics)

not_missing_values <- not_missing_values %>% 
    group_by(Approval_Year, Demographic) %>% 
    summarise(Count = n()) %>% 
    pivot_wider(names_from = Demographic, values_from = Count, values_fill = 0)

not_missing_values <- not_missing_values[, c("Approval_Year", "Female", "Age_65_or_older", "Asian", "Black", "Hispanic", "White")]

# Create missing values table. Subtract non missing from total to get number of missing values
missing_values <- num_approvals - not_missing_values
missing_values$NumApprovals <- num_approvals
missing_values$Approval_Year <- not_missing_values$Approval_Year
missing_values <- missing_values %>% 
    relocate(Approval_Year, NumApprovals, Female, Age_65_or_older, Asian, Black, Hispanic, White) 

# Calculate % of entries which are missing
missing_values_percentage <- round(missing_values / num_approvals * 100)
missing_values_percentage$Approval_Year <- not_missing_values$Approval_Year # Copy the original Years to replace unnecessary math done on the Year column

#########################################
### PROCESS DISEASE BURDEN COMPARISON ###
#########################################

# Read in data
approvals_trim <- fda_approvals %>% select(Brand_Name, Approval_Year, Therapeutic_Area, TA_subgroup, Indication, Asian, Black, Hispanic, White)
disease_burden_long <- disease_burden %>% filter(White != "NA")

# Turn both data frames into long
approvals_trim <- approvals_trim %>% pivot_longer(cols = Asian:White,
                                                    names_to = "Demographic",
                                                    values_to ="Enrollment")
disease_burden_long <- disease_burden_long %>% pivot_longer(cols = Asian:White, 
                                                  names_to = "Demographic", 
                                                  values_to = "Burden")
disease_burden_long$Burden <- disease_burden_long$Burden * 100 # Multiply disease burden from percentage to ones

# Merge the two df's along indication and demographic
comparison_df <- right_join(approvals_trim, disease_burden_long, 
                            by = c("Therapeutic_Area" = "Therapeutic_Area", 
                                   "Indication" = "Indication",
                                   "Demographic" = "Demographic"))

# Make comparison
comparison_df <- comparison_df %>% 
    mutate(Comparison = case_when(
        Enrollment > Burden ~ "Over represented",
        Enrollment == Burden ~ "Appropriately represented",
        Enrollment < Burden ~ "Under represented",
        is.na(Enrollment) ~ "Demographic not collected in trial"
    ))

# Count number of each comparison types
comparison_df_summary <- comparison_df %>% 
    group_by(Therapeutic_Area, Indication, Demographic, Comparison) %>% 
    summarise(Count = n())


saveRDS(fda_approvals, file = "fda_approvals_db.rds")
saveRDS(fda_approvals_long, file = "fda_approvals_long_db.rds")

#######################
### DEFINE SHINY UI ###
#######################

ui <- fluidPage(
    
    titlePanel("2015-2021 FDA Drug Trials Snapshots Data Visualization Explorer"),
    
    tabsetPanel( 
        
        tabPanel("Welcome + Instructions",
            mainPanel(
                h3("Welcome"),
                p("Welcome to a Data Visualization Explorer of the Food and Drug Administration (FDA) Drug Trials Snapshots (DTS) data from 2015-2021"),
                p("The FDA DTS program was launched in 2015 and reports clinical trial demographic data for new molecular 
                  entities and original biologics."),
                p("FDA data can be explored by race, ethnicity, sex, age group, therapeutic area, pharmaceutical sponsor, and approval year 
                  for clinical trials that supported each of the 339 FDA drug and biologic approvals between 2015-2021
                  We hope this data empowers you and your organization to make evidence-based decisions to improve trial representation
                  and advance health equity."),
                p("All data shown here are publicly available on various locations on FDA's website and 
                  our team at Harvard Medical School has scraped and processed these data to serve as input to this 
                  interactive visualization tool."),
                downloadButton("download_raw_data", "Download 2015-2021 FDA Drug Trials Snapshots Data"),
                br(),
                
                h3("Contents of this Data Explorer"),
                tags$b("Descriptive Statistics tab"),
                p("To get warmed up. Static high level summaries and trends of FDA approvals and clinical trial enrollment, 
                  and overview of missing data."),
                tags$b("Explore FDA Approvals"),
                p("To explore. Dynamic inputs allow you to filter clinical trial enrollment by age, sex, 
                  race, ethnicity and visualize the distributions of these demographics over time and 
                  across Therapeutic Areas."),
                tags$b("Detail by Demographic Participation"),
                p("To hone in on specific thresholds of participation. Dynamic inputs select race or ethnicity and 
                  demographic participation rate to explore and compare the characteristics of trials above or below the chosen rate."),
                tags$b("Detail by Therapeutic Area"),
                p("To hone in on specific Therapeutic Areas. Dynamic inputs to explore clinical trial enrollment
                  across demographics and pharma sponsor, and compare enrollment to disease burden."),
                br(),

                h3("Authors"),
                p("Ariel Carmeli, MBI: Graduate student, Harvard Business School 
                  Please reach out with feedback or questions to carmeli.b.ariel@gmail.com"),
                p("Laura Meloney, MS, MPH: Program Director at the Multi-Regional Clinical Trials 
                  Center of Brigham and Women’s Hospital and Harvard (MRCT Center)"),
                p("Barbara E. Bierer, MD: Faculty Director, MRCT Center, and Professor of Medicine, 
                  Harvard Medical School and Brigham and Women’s Hospital, Boston MA"),
                
                h3("Sources"),
                tags$a(href = "https://www.fda.gov/media/143592/download", 
                       "FDA's 2015-19 DTS Drug Trials Snapshots Summary Report"),
                br(),
                tags$a(href = "https://www.fda.gov/drugs/drug-approvals-and-databases/drug-trials-snapshots", 
                       "FDA website with 2015, 2016, 2017, 2018, 2019, 2020, 2021 annual reports"),
                br(),
                tags$a(href = "https://www.fda.gov/drugs/drug-approvals-and-databases/drug-trials-snapshots", 
                       "FDA website with Drug Trials Snapshot for each individual drug/FDA Approval"),
                br(),
                tags$a(href = "https://www.accessdata.fda.gov/scripts/cder/daf/", 
                       "FDA website with Drug Approval Packages"),
                                
                h3("Disclaimer"),
                p("The views and findings expressed herein are those of the individuals contributing 
                  to the work, serving in their individual capacity, and do not imply endorsement or 
                  reflect the views or policies of the U.S. FDA or any affiliated organization or entity. 
                  The MRCT Center is supported by voluntary contributions (www.MRCTCenter.org) and grants"),
                br(),
                
            )
        ),
        
        tabPanel("Descriptive Statistics",
            sidebarLayout(
                sidebarPanel(
                    
                    #radioButtons( 
                    #    "2020_inclusion", 
                    #    "Include 2020 data in all graphs?", 
                    #    choiceNames=list( "Yes", "No" ), 
                    #    choiceValues=list(TRUE,FALSE),
                    #    selected=FALSE
                    #),
                    
                    radioButtons( 
                        "DS_TA_stratify", 
                        "Stratify graphs 3 and 4 by Therapeutic Area?", 
                        choiceNames=list( "Yes", "No" ), 
                        choiceValues=list(TRUE,FALSE),
                        selected=FALSE
                    ),
                    width = 2
                ),
                         
                mainPanel(
                  
                    h2("1. Demographics of Trial Participation"),
                    h6("Data: All FDA approvals from 2015-19. Excludes 2020 and 2021 data in order to validate against FDA's 2015-2019 Drug Trial Snapshot report"),
                    tags$a(href = "https://www.fda.gov/media/143592/download", "This graph closely recreates page 9 in FDA's 2015-19 DTS Drug Trial Snapshot report (link)"),
                    plotOutput("Validation_Demographics", height=300, width = 1000),
                  
                    h2("2. Clinical Trial Participation by Therapeutic Area"),
                    h6("Data: All FDA approvals from 2015-19. Excludes 2020 and 2021 data in order to validate against FDA's 2015-2019 Drug Trial Snapshot report"),
                    tags$a(href = "https://www.fda.gov/media/143592/download", "This graph closely recreates page 30 in FDA's 2015-19 DTS Drug Trial Snapshot report (link)"),
                    h6("Small differences for a given TA exist below, and can be better understood if 
                       FDA provides enrollment size and TA label for each FDA approval"),
                    plotOutput("Validation_Enrollment_by_TA", height=550, width = 1000),
                    
                    h2("3. Total Approvals by Year"),
                    h6("Data: All FDA approvals from 2015-21"),
                    plotOutput("Approvals_DS", height=450, width = 1000),
                    
                    h2("4. Total Patients by Year"),
                    h6("Data: All FDA approvals from 2015-21"),
                    plotOutput("Patients_DS", height=450, width = 1000),
                    
                    h2("5. Distribution of Enrollment Size per FDA Approval by Therapeutic Area"),
                    h6("Data: All FDA approvals from 2015-21"),
                    plotOutput("Enrollment_TA_boxplot_DS", height=600, width = 1000),
                     
                    h2("6. Data Missingness: Number of FDA approvals with missing demographic data"),
                    h6("Data: All FDA approvals from 2015-21"),
                    DT::dataTableOutput("Demographics_reported", height=300, width = 1000),
                    
                )
            )
                
        ),
        
        tabPanel("Explore FDA Approvals",
            sidebarLayout(
                sidebarPanel(
                     selectInput(
                         "race",
                         "Race",
                         choices= c("Asian", "Black", "White", "Other"),
                         selected = c("Asian", "Black", "White"),
                         multiple=T
                     ),
                     
                     selectInput(
                         "ethnicity",
                         "Ethnicity",
                         choices= c("Hispanic", "Non_Hispanic"),
                         selected = c("Hispanic"),
                         multiple=T
                     ),
                     
                     selectInput(
                         "age",
                         "Age",
                         choices= c("Age_under_65", "Age_65_or_older"),
                         selected = FALSE,
                         multiple=T
                     ),
                     
                     selectInput(
                         "sex",
                         "Sex",
                         choices= c("Female", "Male"),
                         selected = FALSE,
                         multiple=T
                     ),
                     
                     radioButtons( 
                         "is_TA_Stratified", 
                         "Stratify by Therapeutic Area?", 
                         choiceNames=list( "Yes", "No" ), 
                         choiceValues=list(TRUE,FALSE),
                         selected=FALSE
                     ),
                     
                     selectInput(
                         "year",
                         "Year(s)",
                         choices=unique(fda_approvals$Approval_Year),
                         selected = c(2015, 2016, 2017, 2018, 2019, 2020, 2021),
                         multiple=T
                     ),
                     
                     radioButtons( 
                         "is_Year_labelled", 
                         "Label by Year?", 
                         choiceNames=list( "Yes", "No" ), 
                         choiceValues=list(TRUE,FALSE),
                         selected=FALSE
                     )#,
                     
                     , width = 2
                     
                     #textOutput("summaryText")
                 ),
                     
                 mainPanel(
                     
                     h2("Distribution of Participation by Demographic in Clinical Trials Across FDA approvals"),
                     h6("Each dot represents represents the proportion of patient enrollment that a certain race or ethnicity 
                        represents in each clinical trial used to inform regulatory approval of each new molecular entity or original biology"),
                     #plotOutput("individualPlotGriswold", heigh=700, width = 1200),
                     plotlyOutput("individualPlot", height=750, width = 1200),

                     h2("Data Table: Median and Average of Clinical Trial Participation"),
                     h6("Median and Average are calculated without accounting for trial size."),
                     h6("Weighted average is calculated accounting for trial size i.e., it equals total number 
                        of individuals per demographic divided by total number of individuals in all demographics. 
                        This weighted average is what FDA Drug Trial Snapshots describes as average in 
                        Table 1 of each of its annual reports"),
                     DT::dataTableOutput("stats_summary_Table", width = 700),
                                          
                     h2("Trend in Clinical Trial Participation by Demographic Over Time"),
                     h6("Boxplots represents 5 points in the distribution: Middle line is median. 
                        Ends of the box are 1st and 3rd quartile.
                        Ends of the whiskers are 1.5 * inter-quartile range from the closer of 1st or 3rd quartile"),
                     plotOutput("change_over_time", height = 600),

                     h2("Count of FDA Approvals Per Participation in Clinical Trials"),
                     h6("Here we zoom in to the tail of the distribution and allow you to count how many 
                        FDA approvals were based on trial data with 0, 1, 2, etc. percent participation by a demographic."),
                     plotOutput("participationCountPlot", height=500),
                     
                     h2("Approval Details"),
                     DT::dataTableOutput("approvalsTable"),
                 )
             )
        ), 

        tabPanel("Detail by Demographic Participation",
                 sidebarLayout(
                   sidebarPanel(
                     
                     selectInput(
                       "demographic",
                       "Demographic",
                       choices= c("Asian", "Black", "White", "Hispanic"),
                       selected = "Black",
                       multiple=FALSE
                     ),
                     
                     sliderInput("Participation_rate",
                                 "Demographic participation rate",
                                 min=0,
                                 max=99,
                                 value=1
                     ),
                     
                     width = 2
                     
                   ),
                   
                   mainPanel(
                     
                     h2("FDA approvals over/under chosen demographic participation"),
                     h6("Total does not equal 339 FDA approvals because of missing data. See Descriptive Statistics tab for number of FDA approvals with missing data by demographic"),
                     DT::dataTableOutput("tail_enrollment_table", width = 700),

                     h2("FDA approvals over/under chosen demographic participation by enrollment size"),
                     h3("Count"),
                     plotOutput("tail_enrollment_count", height = 400, width = 1000),

                     h3("Percent"),
                     plotOutput("tail_enrollment_percent", height = 400, width = 1000),

                     h2("FDA approvals over/under chosen demographic participation by Therapeutic Area"),
                     h3("Count"),
                     plotOutput("tail_TA_count", height = 600, width = 1000),
                     
                     h3("Percent"),
                     plotOutput("tail_TA_percent", height = 600, width = 1000),
                     
                                                               
                   )
                 )
        ), 
        
        
        
        tabPanel("Detail by Therapeutic Area",
            sidebarLayout(
                sidebarPanel(
                    
                    selectInput(
                        "therapeutic_area",
                        "Therapeutic Area",
                        choices=sort(unique(fda_approvals$Therapeutic_Area)),
                        selected = "Oncology",
                        multiple=FALSE
                    ),
                    
                    #selectInput(
                    #    "Stratify_by",
                    #    "Stratify by (Pharma Sponsor, Therapeutic Area subgroup):",
                    #    choices = list("None", "Sponsor", "TA_subgroup"),
                    #    selected = "None",
                    #    multiple = FALSE
                    #),
                    
                    radioButtons( 
                        "Stratify_by", 
                        "Stratify by", 
                        choiceNames=list( "None", "Indication", "Sponsor"), 
                        choiceValues=list("None", "Indication", "Sponsor"),
                        selected="None"
                    ),
                    
                    
                     selectInput(
                         "race_TA_page",
                         "Race",
                         #choices=str_to_title(unique(fda_approvals_long$Demographic)), 
                         choices= c("Asian", "Black", "White", "Other"),
                         selected = c("Asian", "Black", "White"),
                         multiple=T
                     ),
                     
                     selectInput(
                         "ethnicity_TA_page",
                         "Ethnicity",
                         choices= c("Hispanic", "Non_Hispanic"),
                         selected = c("Hispanic"),
                         multiple=T
                     ),
                     
                     selectInput(
                         "age_TA_page",
                         "Age",
                         choices= c("Age_under_65", "Age_65_or_older"),
                         selected = FALSE,
                         multiple=T
                     ),
                     
                     selectInput(
                         "sex_TA_page",
                         "Sex",
                         choices= c("Female", "Male"),
                         selected = FALSE,
                         multiple=T
                     ),
                     
                     selectInput(
                         "year_TA_page",
                         "Year(s)",
                         choices=unique(fda_approvals$Approval_Year),
                         selected = c(2015, 2016, 2017, 2018, 2019, 2020, 2021),
                         multiple=T
                     ),
                     
                     sliderInput("Sponsor_size",
                                 "If too many sponsors, filter >= X approvals",
                                 min=1,
                                 max=3,
                                 value=1
                     ),
                     
                     radioButtons( 
                         "is_Enrollment_Stratified", 
                         "View Enrollment size?", 
                         choiceNames=list( "Yes", "No" ), 
                         choiceValues=list(TRUE,FALSE),
                         selected=FALSE
                     ),
                     
                     width = 2
                     
                 ), 
                 mainPanel(
                     
                     h2("Distribution of Participation by Demographic in Clinical Trials Across FDA approvals"),
                     h6("Each dot represents represents the proportion of patient enrollment that a certain race or ethnicity 
                        represents in each clinical trial used to inform regulatory approval of each new molecular entity or original biology"),
                     plotlyOutput("TA_individualPlot", height = 800),
                     
                     h2("Trend in Clinical Trial Participation by Demographic Over Time"),
                     h6("Boxplots represents 5 points in the distribution: Middle line is median. 
                        Ends of the box are 1st and 3rd quartile.
                        Ends of the whiskers are 1.5 * inter-quartile range from the closer of 1st or 3rd quartile"),
                     plotOutput("TA_change_over_time", height=600),
                     
                     h2("Comparison of clinical trial representation and disease incidence"),
                     h6("Bar chart represents all 2015-2021 FDA approvals"),
                     h6("Note: This analysis is only available with oncology and infectious disease. 
                        You wil see a red error message below if you have selected a different TA"),
                     h6("Table below the graph reflects incidence data collected from NCI SEER and the CDC"),
                     plotOutput("TA_Disease_Burden_Comparison", height = 700), #height = 325, width - 800
                     DT::dataTableOutput("Disease_Burden_table"),
                     
                     h2("Approval Details"),
                     DT::dataTableOutput("TA_approvalsTable"),
                 )
             )
        )
    ),
    tags$head(tags$style(HTML('* {font-family: "Arial"};')))
)

################################
### DYNAMIC DATA INTERACTION ###
################################

server <- function(input, output) {

    # reactive for disease burden
    disease_burden_table <- reactive({
        
        selection <- disease_burden %>% 
            filter(White != "NA") %>% # filter out NAs
            filter( Therapeutic_Area %in% input$therapeutic_area ) # Filter per user input on TA
        
        selection
        
    })
    
    # reactive for descriptive statistics
    approvals_DS <- reactive({
        selection <- fda_approvals
        selection
    })
    
    missing <- reactive({
        selection <- missing_values #missing_values_percentage
        selection
    })
    
    # default no adjustment
    approvals_15_19 <- reactive({
        selection <- fda_approvals %>% filter(Approval_Year == 2015 |
                                                Approval_Year == 2016 |
                                                Approval_Year == 2017 |
                                                Approval_Year == 2018 |
                                                Approval_Year == 2019)
        selection
    })
    
    # reactive expression to filter selected approvals
    approvals <- reactive({
        selection <- fda_approvals_long %>%
            select(Brand_Name, Therapeutic_Area, TA_subgroup, Indication, Indication_long, Enrollment, Demographic, Percentage, Approval_Year, Sponsor) %>%
            filter(Percentage != "NA") %>% 
            unique()

        if ( !is.null( input$race ) | !is.null( input$ethnicity ) | !is.null( input$age ) | !is.null( input$sex ) ) {
            selection <- selection %>% filter( Demographic %in% input$race | Demographic %in% input$ethnicity | Demographic %in% input$age | Demographic %in% input$sex)
        }
        
        if ( !is.null( input$year ) ) {
            selection <- selection %>% filter( Approval_Year %in% input$year )
        }
        
        # Find unique set of drugs that have participation values in the bounds from input
        drugs <- selection %>% 
            group_by(Brand_Name) %>% 
            unique()
        
        # Filter selection to include this set of drugs
        selection <- selection %>% 
            filter( Brand_Name %in% drugs$Brand_Name )
        
        selection
        
    })

    # reactive for tail tab: enrollment
    approvals_tail_enrollment <- reactive({
      
      selection <- fda_approvals_long %>% 
        mutate(Enrollment_bucket = case_when(
          Enrollment <= 100 ~ "<=100",
          Enrollment > 100 & Enrollment <=250 ~ "101-250",
          Enrollment > 250 & Enrollment <=500 ~ "251-500",
          Enrollment > 500 & Enrollment <=1000 ~ "501-1000",
          Enrollment > 1000 & Enrollment <=3000 ~ "1001-3000",
          Enrollment > 3000 ~ ">3000",
          TRUE ~ "unknown"
        ))
      selection <- selection %>% relocate(Enrollment_bucket, .after = Enrollment)
      
      selection$Enrollment_bucket <- factor(selection$Enrollment_bucket,
                                                levels = c(">3000",
                                                           "1001-3000", 
                                                           "501-1000", 
                                                           "251-500", 
                                                           "101-250",
                                                           "<=100"))
      
      categories <- selection %>% select(Enrollment_bucket) %>% unique()
      
      equal_or_less <- selection %>% filter(Demographic %in% input$demographic) %>% filter(Percentage <= input$Participation_rate) %>% group_by(Enrollment_bucket) %>% summarize(count = n())
      more <- selection %>% filter(Demographic %in% input$demographic) %>% filter(Percentage > input$Participation_rate) %>% group_by(Enrollment_bucket) %>% summarize(count = n())
      
      equal_or_less_all_categories <- categories
      equal_or_less_all_categories$count <- equal_or_less$count[match(equal_or_less_all_categories$Enrollment_bucket,
                                                              equal_or_less$Enrollment_bucket)]
      equal_or_less_all_categories <- equal_or_less_all_categories %>% replace(is.na(.),0)
      equal_or_less_all_categories <- equal_or_less_all_categories %>% mutate(percent = count / sum(equal_or_less_all_categories$count))
      equal_or_less_all_categories$over_under_label <- "Equal_or_under"
      
      more_all_categories <- categories
      more_all_categories$count <- more$count[match(more_all_categories$Enrollment_bucket,
                                           more$Enrollment_bucket)]
      more_all_categories <- more_all_categories %>% replace(is.na(.),0)
      more_all_categories <- more_all_categories %>% mutate(percent = count / sum(more_all_categories$count))
      more_all_categories$over_under_label <- "Over"
      
      combine <- rbind(equal_or_less_all_categories, more_all_categories)
      combine$percent <- round(combine$percent, 2) * 100
      
      combine
      
    })

    # reactive for tail tab: therapeutic area
    approvals_tail_TA <- reactive({
      
      selection <- fda_approvals_long
      
      categories <- selection %>% select(Therapeutic_Area) %>% unique()
      
      equal_or_less <- selection %>% filter(Demographic %in% input$demographic) %>% filter(Percentage <= input$Participation_rate) %>% group_by(Therapeutic_Area) %>% summarize(count = n())
      more <- selection %>% filter(Demographic %in% input$demographic) %>% filter(Percentage > input$Participation_rate) %>% group_by(Therapeutic_Area) %>% summarize(count = n())
      
      equal_or_less_all_categories <- categories
      equal_or_less_all_categories$count <- equal_or_less$count[match(equal_or_less_all_categories$Therapeutic_Area,
                                                                      equal_or_less$Therapeutic_Area)]
      equal_or_less_all_categories <- equal_or_less_all_categories %>% replace(is.na(.),0)
      equal_or_less_all_categories <- equal_or_less_all_categories %>% mutate(percent = count / sum(equal_or_less_all_categories$count))
      equal_or_less_all_categories$over_under_label <- "Equal_or_under"
      
      more_all_categories <- categories
      more_all_categories$count <- more$count[match(more_all_categories$Therapeutic_Area,
                                                    more$Therapeutic_Area)]
      more_all_categories <- more_all_categories %>% replace(is.na(.),0)
      more_all_categories <- more_all_categories %>% mutate(percent = count / sum(more_all_categories$count))
      more_all_categories$over_under_label <- "Over"
      
      combine <- rbind(equal_or_less_all_categories, more_all_categories)
      combine$percent <- round(combine$percent, 2) * 100
      
      combine
      
    })
            
    # reactive expression for the count chart
    approvals_count <- reactive({
        selection <- fda_approvals_long %>% 
            filter(Percentage != "NA")
        
        if ( !is.null( input$race ) | !is.null( input$ethnicity ) | !is.null( input$age ) | !is.null( input$sex ) ) {
            selection <- selection %>% filter( Demographic %in% input$race | Demographic %in% input$ethnicity | Demographic %in% input$age | Demographic %in% input$sex )
        }
        
        if ( !is.null( input$year ) ) {
            selection <- selection %>% filter( Approval_Year %in% input$year )
        }
        
        selection <- selection %>% 
            group_by(Brand_Name, Demographic, Percentage, Approval_Year) %>% 
            unique() %>% 
            summarize(Count = n()) %>% 
            group_by(Demographic, Percentage, Approval_Year) %>% 
            summarize(Count = sum(Count))
        
        selection
    })
    
    # reactive expression for the cumulative chart
    approvals_cum_count <- reactive({
        selection <- fda_approvals_long %>% 
            filter(Percentage != "NA")
        
        if ( !is.null( input$race ) | !is.null( input$ethnicity ) | !is.null( input$age ) | !is.null( input$sex ) ) {
            selection <- selection %>% filter( Demographic %in% input$race | Demographic %in% input$ethnicity | Demographic %in% input$age | Demographic %in% input$sex)
        }
        
        if ( !is.null( input$year ) ) {
            selection <- selection %>% filter( Approval_Year %in% input$year )
        }
        
        selection
    })
    
    
    # Reactive expression to filter selected approvals on the TA/disease tab 
    approvals_TA <- reactive({
        
        selection <- fda_approvals_long %>%
            select(Brand_Name, Therapeutic_Area, TA_subgroup, Indication, Indication_long, Enrollment, Demographic, Percentage, Approval_Year, Sponsor) %>%
            filter(Percentage != "NA") %>% 
            unique()
        
        if ( !is.null( input$race_TA_page ) | !is.null( input$ethnicity_TA_page ) | !is.null( input$age_TA_page ) | !is.null( input$sex_TA_page ) ) {
            selection <- selection %>% filter( Demographic %in% input$race_TA_page | Demographic %in% input$ethnicity_TA_page | Demographic %in% input$age_TA_page | Demographic %in% input$sex_TA_page)
        }
        
        if ( !is.null( input$year_TA_page ) ) {
            selection <- selection %>% filter( Approval_Year %in% input$year_TA_page )
        }
        
        if ( !is.null( input$therapeutic_area ) ) {
            selection <- selection %>% filter( Therapeutic_Area %in% input$therapeutic_area )
        }
        
        if ( "Sponsor" %in% input$Stratify_by ) {
            
            # Identify sponsors with > X drugs
            sponsors <- selection %>% 
                pivot_wider(names_from = Demographic, values_from = Percentage, values_fill = -1) %>% 
                select(Sponsor) %>% 
                group_by(Sponsor) %>% 
                summarize(Count = n()) %>% 
                filter(Count >= input$Sponsor_size) %>% 
                select(Sponsor)
            
            selection <- selection %>% filter( Sponsor %in% sponsors$Sponsor )
            
        }
        
        selection
        
    })
    
    # Reactive expression to make comparison of trial enrollment and disease burden
    disease_burden_comparison <- reactive({
        
        selection <- comparison_df_summary %>% filter( Therapeutic_Area %in% input$therapeutic_area )
        
        selection
        
    })
    
    output$download_raw_data <- downloadHandler(
        filename = function() {
            paste("FDA_Drug_Trials_Snapshots_2015_to_2021.csv", sep = "")
        },
        content = function(file) {
            write.csv(approvals_DS(), file, row.names = FALSE)
        }
    )
    
    output$Approvals_DS <- renderPlot({
        
        plot <- approvals_DS() %>% 
            ggplot(aes(x=Approval_Year)) + 
            geom_bar(width = 0.7) +
            geom_text(stat='count', aes(label=..count..), vjust=-1, size=5) +
            xlab("Year") + 
            ylab("Number of Approvals") +
            theme(strip.text.x = element_text(size = 15),
                  axis.text.x=element_text(size=12),
                  axis.text.y=element_text(size=12),
                  axis.title.x=element_text(size=15, face="bold"), 
                  axis.title.y=element_text(size=15, face="bold")) +
            #theme(axis.title.x=element_text(size=12, face="bold"), axis.title.y=element_text(size=12, face="bold")) + 
            scale_x_continuous(breaks = round(seq(min(2015), max(2021), by = 1),1)) #+ 
            #expand_limits(y=c(0, 70))
        
        if ( input$DS_TA_stratify ){
            plot <- plot + 
                facet_wrap(~Therapeutic_Area) +
                expand_limits(y=c(0, 30))

        }
        
        else{
            plot <- plot + expand_limits(y=c(0, 70))
        }
        
        plot
        
    })
    
    output$Patients_DS <- renderPlot({
        
        plot <- approvals_DS() %>% 
            ggplot(aes(x=Approval_Year, y=Enrollment)) + 
            geom_col(width = 0.7) +
            #geom_text(aes(label=Enrollment), position = position_stack(vjust = 0.5), size = 4) +
            #geom_text(aes(label = Enrollment)) +
            #geom_text(aes(label= sum(Enrollment)), vjust=-1) +
            #geom_text(size = 3, position = position_stack(vjust = 0.5)) +
            xlab("Year") + 
            ylab("Number of patients enrolled in pivotal clinical trial") +
            theme(axis.text.x=element_text(size=12),
                  axis.text.y=element_text(size=12),
                  axis.title.x=element_text(size=15, face="bold"), 
                  axis.title.y=element_text(size=15, face="bold")) +
            #theme(axis.title.x=element_text(size=12, face="bold"), axis.title.y=element_text(size=12, face="bold")) + 
            scale_x_continuous(breaks = round(seq(min(2015), max(2021), by = 1),1))
        
        if ( input$DS_TA_stratify ){
            plot <- plot + 
                facet_wrap(~Therapeutic_Area) + 
                expand_limits(y=c(0, 75000))
        }
        else{
            plot <- plot + expand_limits(y=c(0, 120000))
        }
        
        plot
        
    })
    
    output$Enrollment_TA_boxplot_DS <- renderPlot({
        plot <- approvals_DS() %>% 
            select(Therapeutic_Area, Enrollment) %>% 
            filter(Enrollment != "NA") %>% 
            mutate(Therapeutic_Area = fct_reorder(Therapeutic_Area, Enrollment, .fun='median')) %>%
            ggplot(aes(Therapeutic_Area, Enrollment)) +
            geom_boxplot() +
            geom_jitter() +
            theme(axis.text.x=element_text(size=12),
                  axis.text.y=element_text(size=12),
                  axis.title.x=element_text(size=15, face="bold"), 
                  axis.title.y=element_text(size=15, face="bold")) +
            #theme(axis.title.x = element_text(size=12, face="bold"), axis.title.y=element_text(size=12, face="bold")) +
            #theme(text = element_text(size=12)) +
            coord_flip() +
            xlab("Therapeutic Area") +
            ylab("Number of patients")
        
        plot
        
    })
    
    
    output$Validation_Enrollment_by_TA <- renderPlot({
        
        # Create data frame with patient count by TA - using our scraped data
        df <- approvals_15_19() %>% 
            #filter(Approval_Year != 2020) %>% 
            select(Therapeutic_Area, Enrollment) %>% 
            mutate(Therapeutic_Area = case_when(
                Therapeutic_Area == "Oncology" ~ "Oncology and Hematology", 
                Therapeutic_Area == "Hematology" ~ "Oncology and Hematology",
                Therapeutic_Area == "Endocrinology and Metabolism" ~ "Endocrinology and Metabolism",
                Therapeutic_Area == "Infectious Disease" ~ "Infectious Disease",
                Therapeutic_Area == "Neurology" ~ "Neurology",
                Therapeutic_Area == "Gynecology" ~ "Gynecology",
                Therapeutic_Area == "Dermatology" ~ "Dermatology",
                Therapeutic_Area == "Pulmonology and Rheumatology" ~ "Pulmonology and Rheumatology",
                Therapeutic_Area == "Gastroenterology" ~ "Gastroenterology",
                Therapeutic_Area == "Psychiatry" ~ "Psychiatry",
                Therapeutic_Area == "Sleep Disorders" ~ "Psychiatry",
                Therapeutic_Area == "Ophthalmology" ~ "Ophthalmology",
                Therapeutic_Area == "Anesthesia and Analgesia" ~ "Anesthesia and Analgesia",
                Therapeutic_Area == "Medical Imaging" ~ "Medical Imaging",
                Therapeutic_Area == "Cardiovascular Diseases" ~ "Cardiovascular Diseases"
            )) %>% 
            group_by(Therapeutic_Area) %>% 
            summarize(Our_database = sum(Enrollment)) %>% 
            arrange(desc(Our_database))
        
        validation_df <- data.frame(df$Therapeutic_Area, df$Our_database) #, fda_snapshots)
        names(validation_df) <- c("Therapeutic_Area", "Our_Database") #, "FDA_Snapshots")
        
        # Create data frame with patient count by TA - using estimate from FDA DTS 5 year report 
        # Cardio, Onc, Endocrin, ID, Neurology, Derm, Gyn, Pulm/Rheum, GI, Psychiatry, Anaesthesia, Ophthal, Medical Imaging
        fda_snapshots <- c(59000, 
                           41000, 
                           35000, 
                           32500, 
                           26000, 
                           21000, 
                           20500, 
                           20000, 
                           15000, 
                           10500, 
                           6500, 
                           4500, 
                           1000)
        fda_snapshots_TAname <- c("Cardiovascular Diseases",
                                  "Endocrinology and Metabolism",
                                  "Oncology and Hematology",
                                  "Infectious Disease",
                                  "Neurology",
                                  "Gynecology",
                                  "Dermatology",
                                  "Pulmonology and Rheumatology",
                                  "Gastroenterology",
                                  "Psychiatry",
                                  "Ophthalmology",
                                  "Anesthesia and Analgesia",
                                  "Medical Imaging")
        
        fda_dts_df <- data.frame(fda_snapshots_TAname, fda_snapshots)
        names(fda_dts_df) <- c("Therapeutic_Area", "FDA_Snapshots")
        
        validation_df$fda_snapshots <- fda_dts_df$FDA_Snapshots[match(validation_df$Therapeutic_Area, 
                                                                      fda_dts_df$Therapeutic_Area)]
        names(validation_df) <- c("Therapeutic_Area", "Our_Database", "FDA_Snapshots")
        
        fda_approvals_check_long <- pivot_longer(validation_df, cols = Our_Database:FDA_Snapshots, 
                                                 names_to = "Database", values_to = "Participants")
        
        #fda_approvals_check_long <- pivot_longer(validation_df, cols = Our_Database, 
        #                                         names_to = "Database", values_to = "Participants")
        
        Patient_participation_graph <- fda_approvals_check_long %>% 
            ggplot(aes(reorder(Therapeutic_Area, Participants), Participants, fill=Database)) +
            geom_bar(stat = "identity", position = 'dodge') +
            geom_text(aes(label=Participants), position = position_dodge(0.9), hjust=-0.2, size=4) +
            coord_flip() + 
            expand_limits(y=c(0, 70000)) +
            xlab("Therapeutic Area") +
            ylab("Number of Participants") + 
            scale_y_continuous(breaks = round(seq(min(0), max(70000), by = 5000),1)) +
            #theme(text = element_text(size=20)) +
            theme(legend.text = element_text(size=12),
                  legend.title = element_text(size=12),
                  axis.text.x=element_text(size=12),
                  axis.text.y=element_text(size=12),
                  axis.title.x=element_text(size=15, face="bold"), 
                  axis.title.y=element_text(size=15, face="bold"))
        
        Patient_participation_graph
        
    })
    
    output$Validation_Demographics <- renderPlot({
        ## To account for missing data, we capture % representation of each demographic only in trials where that demographic is captured
        
        # Race
        
        global_black_participation <- approvals_15_19() %>% 
            select(Brand_Name, Enrollment, Black) %>% 
            filter(Black != "NA") %>% 
            mutate(Black_participants = (Black/100) * Enrollment)
        
        global_white_participation <- approvals_15_19() %>% 
            select(Brand_Name, Enrollment, White) %>% 
            filter(White != "NA") %>% 
            mutate(White_participants = (White/100) * Enrollment)
        
        global_asian_participation <- approvals_15_19() %>% 
            select(Brand_Name, Enrollment, Asian) %>% 
            filter(Asian != "NA") %>% 
            mutate(Asian_participants = (Asian/100) * Enrollment)
        
        race_participation <- c(round(sum(100 * global_asian_participation$Asian_participants) /
                                          sum(global_asian_participation$Enrollment), 0), 
                                round(sum(100 * global_black_participation$Black_participants) / 
                                          sum(global_black_participation$Enrollment), 0),
                                round(sum(100 * global_white_participation$White_participants) / 
                                          sum(global_white_participation$Enrollment), 0))
        
        races <- c("Asian", "Black", "White")
        race_participation <- data.frame("Race", races, race_participation)
        names(race_participation) <- c("Demographic", "Category", "Value")
        
        # Sex
        
        global_female_participation <- approvals_15_19() %>% 
            select(Brand_Name, Enrollment, Female) %>% 
            filter(Female != "NA") %>% 
            mutate(Female_participants = (Female/100) * Enrollment)
        
        global_male_participation <- approvals_15_19() %>%
            select(Brand_Name, Enrollment, Male) %>% 
            filter(Male != "NA") %>% 
            mutate(Male_participants = (Male/100) * Enrollment)
        
        sex_participation <- c(round(sum(100 * global_female_participation$Female_participants) /
                                         sum(global_female_participation$Enrollment), 0), 
                               round(sum(100 * global_male_participation$Male_participants) / 
                                         sum(global_male_participation$Enrollment), 0))
        
        sexes <- c("Female", "Male")
        sex_participation <- data.frame("Sex", sexes, sex_participation)
        names(sex_participation) <- c("Demographic", "Category", "Value")
        
        # Age
        
        global_Age_under_65_participation <- approvals_15_19() %>% 
            select(Brand_Name, Enrollment, Age_under_65) %>% 
            filter(Age_under_65 != "NA") %>% 
            mutate(Age_under_65_participants = (Age_under_65/100) * Enrollment)
        
        global_Age_65_or_older_participation <- approvals_15_19() %>%
            select(Brand_Name, Enrollment, Age_65_or_older) %>% 
            filter(Age_65_or_older != "NA") %>% 
            mutate(Age_65_or_older_participants = (Age_65_or_older/100) * Enrollment)
        
        age_participation <- c(round(sum(100 * global_Age_under_65_participation$Age_under_65_participants) /
                                         sum(global_Age_under_65_participation$Enrollment), 0), 
                               round(sum(100 * global_Age_65_or_older_participation$Age_65_or_older_participants) / 
                                         sum(global_Age_65_or_older_participation$Enrollment), 0))
        
        ages <- c("Age_under_65", "Age_65_or_older")
        age_participation <- data.frame("Age", ages, age_participation)
        names(age_participation) <- c("Demographic", "Category", "Value")
        
        # Combine the tables
        demographic_participation <- rbind(age_participation, sex_participation, race_participation)
        
        demographic_participation_graph <- demographic_participation %>% 
            ggplot(aes(x = Category, y = Value)) +
            geom_col(width = 0.5) + 
            geom_text(aes(label = Value), position = position_dodge(0.9), vjust=-0.4, size=5) +
            expand_limits(y=c(0, 100)) +
            theme(axis.title.x=element_text(size=15, face="bold"), 
                  axis.title.y=element_text(size=15, face="bold"),
                  axis.text.x=element_text(size=15),
                  axis.text.y=element_text(size=15),
                  strip.text.x = element_text(size = 15)) + 
            facet_wrap(~Demographic, scales = "free") +
            xlab("") +
            ylab("% Participation") #+
            #ggtitle("Demographics of Trial Participation [FDA approvals 2015-19; excludes 2020 for purposes of dataset validation]")
        
        demographic_participation_graph
        
    })

    output$Demographics_reported <- DT::renderDataTable({
        missing()
    })
    
    output$individualPlotGriswold <- renderPlot({
        
        plot <- approvals() %>% 
            ggplot(aes(x=Percentage,
                       y=fct_rev(Demographic), 
                       group=Demographic, 
                       text = paste("Total trial enrollment (all demographics):", Enrollment,
                                    "<br>Drug:", Brand_Name,
                                    "<br>Sponsor:", Sponsor,
                                    "<br>Therapeutic Area:", Therapeutic_Area,
                                    "<br>Indication:", Indication,
                                    "<br>Approval Year :", Approval_Year))) +
            geom_density_ridges2(aes(#point_color = Demographic,
                                     #point_fill = Demographic,
                                     point_shape = 'circle'),
                                 alpha = .2, 
                                 point_alpha = 1, 
                                 jittered_points = TRUE,
                                 point_size = 3,
                                 scale=0.9,
                                 trim=TRUE) +
            scale_discrete_manual(aesthetics = "point_shape", values = c(21, 22, 23, 24)) +
            scale_point_color_hue(l = 40) +
            scale_x_continuous(breaks = round(seq(min(0), max(100), by = 10),1)) +
            theme(legend.position = "none",
                  axis.text.x = element_text(size=20),
                  axis.text.y = element_text(size=20),
                  axis.title.x=element_text(size=20, face="bold"),
                  axis.title.y=element_text(size=20, face="bold")) +
            ylab("Demographic")
        
        if ( input$is_TA_Stratified ) {
            plot <- plot + 
                facet_wrap(~Therapeutic_Area) #+ 
                #scale_y_continuous(breaks = round(seq(min(0), max(100), by = 10),1))
        }
        
        if ( input$is_Year_labelled ) {
            plot <- plot #+ 
                #geom_jitter(width = 0.2, aes(colour=factor(Approval_Year))) +
                #scale_color_brewer(palette="PuRd")
        }
        else{
            plot <- plot #+ 
                #geom_jitter(width = 0.2, aes(colour=Demographic))
        }
        
        plot
        
    })
    
    output$individualPlot <- renderPlotly({
        
        plot <- approvals() %>% 
            ggplot(aes(Demographic, 
                       Percentage, 
                       text = paste("Total trial enrollment (all demographics):", Enrollment,
                                    "<br>Drug:", Brand_Name,
                                    "<br>Sponsor:", Sponsor,
                                    "<br>Therapeutic Area:", Therapeutic_Area,
                                    "<br>Indication:", Indication,
                                    "<br>Approval Year :", Approval_Year))) +
            scale_y_continuous(breaks = round(seq(min(0), max(100), by = 5),1)) +
            theme(axis.text.x = element_text(size=12),
                  axis.text.y = element_text(size=8),
                  axis.title.x=element_text(size=12, face="bold"),
                  axis.title.y=element_text(size=12, face="bold"))
            #geom_jitter(width = 0.2, aes(colour=Demographic)) + 
            #theme(legend.position = "top")+#, legend.title = element_blank()) +
             #+
            #theme(axis.title.x=element_text(size=12, face="bold"), axis.title.y=element_text(size=12, face="bold"))
            #ggtitle("Distribution of clinical trial participation")
        
        if ( input$is_TA_Stratified ) {
            plot <- plot + 
                facet_wrap(~Therapeutic_Area) + 
                scale_y_continuous(breaks = round(seq(min(0), max(100), by = 10),1))
        }
        
        if ( input$is_Year_labelled ) {
            plot <- plot + 
                geom_jitter(width = 0.2, aes(colour=factor(Approval_Year))) +
                scale_color_brewer(palette="PuRd")
        }
        else{
            plot <- plot + 
                geom_jitter(width = 0.2, aes(colour=Demographic))
        }
        
        #plot <- plot + theme(legend.position = "top")
    
        ggplotly(plot) %>% 
            layout(legend = list(orientation = "h", y = 1.1, x = 0.03))
        
    })
    
    output$stats_summary_Table <- DT::renderDataTable(
        
        if ( input$is_TA_Stratified ) { 
            if (input$is_Year_labelled) { # Year and TA
                demographics_selected <- approvals() %>% select(Demographic) %>% unique()
                
                summary_statistics_all_years %>% filter(Demographic %in% demographics_selected$Demographic)
                
                #approvals() %>% 
                #    group_by(Demographic, Therapeutic_Area, Approval_Year) %>% 
                #    summarize(Median = round(median(Percentage), 1), 
                #              Average = round(mean(Percentage), 1))
            }
            else{ # TA only
                approvals() %>% 
                    group_by(Demographic, Therapeutic_Area) %>% 
                    summarize(Median = round(median(Percentage), 1), 
                              Average = round(mean(Percentage), 1))
            }
        }
        else if ( input$is_Year_labelled) { # Year only
            demographics_selected <- approvals() %>% select(Demographic) %>% unique()
            
            summary_statistics %>% filter(Demographic %in% demographics_selected$Demographic)
            
            #approvals() %>% 
            #    group_by(Demographic, Approval_Year) %>% 
            #    summarize(Median = round(median(Percentage), 1), 
            #              Average = round(mean(Percentage), 1))
        }
        else{
            
            demographics_selected <- approvals() %>% select(Demographic) %>% unique()
            
            summary_statistics_all_years %>% filter(Demographic %in% demographics_selected$Demographic)
            
            #approvals() %>% 
            #    group_by(Demographic) %>% 
            #    summarize(Median = round(median(Percentage), 1), 
            #              Average = round(mean(Percentage), 1))
            
        }
    )
    
    output$participationCountPlot <- renderPlot({
        
        # bar chart - for all years
        plot <- approvals_count() %>% 
            group_by(Demographic, Percentage) %>% 
            summarize(Count = sum(Count)) %>% 
            filter(Percentage <= 20) %>%
            ggplot(aes(Percentage, y=Count)) +
            geom_bar(stat = "identity", position = 'dodge', width = 0.8) +
            geom_text(aes(label=Count), position = position_dodge(0.9), vjust=-0.3, size=6) +
            #ylim(0,65) +
            expand_limits(y=c(0, 65)) +
            scale_x_continuous(breaks = round(seq(0, 20, by = 1),1)) +
            theme(axis.title.x=element_text(size=15, face="bold"), 
                  axis.title.y=element_text(size=15, face="bold"),
                  axis.text.x=element_text(size=15),
                  strip.text.x = element_text(size = 15),
                  axis.ticks.y = element_blank(),
                  axis.text.y = element_blank()) + 
            xlab("Percent participation") +
            ylab("Count of FDA approvals") +
            #ggtitle("Number of FDA approvals with under 25 percent participation by demographic in its trials") +
            facet_wrap(~Demographic, ncol=2)
        
        plot
    })
    
    output$change_over_time <- renderPlot({
        
        plot <- approvals() %>% 
            #filter(Demographic %in% dems) %>% 
            filter(Percentage != "NA") %>% 
            ggplot(aes(x = Demographic, y = Percentage, fill = factor(Approval_Year))) +
            geom_boxplot(outlier.shape = NA, lwd=0.8) +
            scale_fill_brewer(palette="PuRd")+
            #geom_point(position=position_jitterdodge(),alpha=0.1) +
            scale_y_continuous(breaks = round(seq(0, 100, by = 10),1)) +
            theme(legend.position = "top",
                  legend.title = element_blank(),
                  legend.text = element_text(size=20),
                  axis.text.x = element_text(size=16),
                  axis.text.y = element_text(size=16),
                  axis.title.x=element_text(size=16, face="bold"), 
                  axis.title.y=element_text(size=16, face="bold")) #+ 
            #coord_flip()
        
        if ( input$is_TA_Stratified ) {
            plot <- plot + 
                facet_wrap(~Therapeutic_Area, scales = "free")
        }
        
        plot
        
    })
    
    output$approvalsTable <- DT::renderDataTable(
        approvals() %>% 
            pivot_wider(names_from = Demographic, values_from = Percentage), options=list(pageLength=10)
    )

    output$tail_enrollment_table <- DT::renderDataTable(

      approvals_tail_enrollment() %>% group_by(over_under_label) %>% summarize(sum = sum(count))
            
    )
    
    output$tail_enrollment_count <- renderPlot({

      plot <- approvals_tail_enrollment() %>% 
        ggplot(aes(Enrollment_bucket, count, fill=over_under_label)) +
        geom_bar(stat = "identity", position = 'dodge') +
        geom_text(aes(label=count), position = position_dodge(0.9), hjust=-0.2, size=4) +
        coord_flip() + 
        xlab("Enrollment bucket") +
        ylab("Number of FDA approvals") + 
        theme(legend.text = element_text(size=12),
              legend.title = element_blank(),
              legend.position = "top",
              axis.text.x=element_text(size=12),
              axis.text.y=element_text(size=12),
              axis.title.x=element_text(size=15, face="bold"), 
              axis.title.y=element_text(size=15, face="bold"))
      
      plot
      
    })

    output$tail_enrollment_percent <- renderPlot({
      
      plot <- approvals_tail_enrollment() %>% 
        ggplot(aes(Enrollment_bucket, percent, fill=over_under_label)) +
        geom_bar(stat = "identity", position = 'dodge') +
        geom_text(aes(label=percent), position = position_dodge(0.9), hjust=-0.2, size=4) +
        coord_flip() + 
        expand_limits(y=c(0, 100)) +
        xlab("Enrollment bucket") +
        ylab("Percent") + 
        scale_y_continuous(breaks = round(seq(min(0), max(100), by = 20),1)) +
        theme(legend.text = element_text(size=12),
              legend.title = element_blank(),
              legend.position = "top",
              axis.text.x=element_text(size=12),
              axis.text.y=element_text(size=12),
              axis.title.x=element_text(size=15, face="bold"), 
              axis.title.y=element_text(size=15, face="bold"))
      
      plot
      
    })
    
    output$tail_TA_count <- renderPlot({
      
      plot <- approvals_tail_TA() %>% 
        ggplot(aes(reorder(Therapeutic_Area, count), count, fill=over_under_label)) +
        geom_bar(stat = "identity", position = 'dodge') +
        geom_text(aes(label=count), position = position_dodge(0.9), hjust=-0.2, size=4) +
        coord_flip() + 
        xlab("Enrollment bucket") +
        ylab("Number of FDA approvals") + 
        theme(legend.text = element_text(size=12),
              legend.title = element_blank(),
              legend.position = "top",
              axis.text.x=element_text(size=12),
              axis.text.y=element_text(size=12),
              axis.title.x=element_text(size=15, face="bold"), 
              axis.title.y=element_text(size=15, face="bold"))
      
      plot
      
    })
    
    output$tail_TA_percent <- renderPlot({
      
      plot <- approvals_tail_TA() %>% 
        ggplot(aes(reorder(Therapeutic_Area, count), percent, fill=over_under_label)) +
        geom_bar(stat = "identity", position = 'dodge') +
        geom_text(aes(label=percent), position = position_dodge(0.9), hjust=-0.2, size=4) +
        coord_flip() + 
        expand_limits(y=c(0, 100)) +
        xlab("Enrollment bucket") +
        ylab("Percent") + 
        scale_y_continuous(breaks = round(seq(min(0), max(100), by = 20),1)) +
        theme(legend.text = element_text(size=12),
              legend.title = element_blank(),
              legend.position = "top",
              axis.text.x=element_text(size=12),
              axis.text.y=element_text(size=12),
              axis.title.x=element_text(size=15, face="bold"), 
              axis.title.y=element_text(size=15, face="bold"))
      
      plot
      
    })

    
                
    output$TA_individualPlot <- renderPlotly({
        
        plot <- approvals_TA() %>% 
            ggplot(aes(Demographic, Percentage, text = paste("Total trial enrollment (all demographics):", Enrollment,
                                                             "<br>Drug:", Brand_Name,
                                                             "<br>Sponsor:", Sponsor,
                                                             "<br>Therapeutic Area:", Therapeutic_Area,
                                                             "<br>Indication:", Indication,
                                                             "<br>Approval Year :", Approval_Year)))
        
        if ( input$is_Enrollment_Stratified ) {
            plot <- plot + 
                geom_jitter(width = 0.2, aes(colour=Demographic, size = Enrollment)) + 
                theme(legend.position = "top", legend.title = element_blank()) +
                scale_y_continuous(breaks = round(seq(min(0), max(100), by = 5),1))
            
        }   
        else{
            plot <- plot + 
                geom_jitter(width = 0.2, aes(colour=Demographic)) + 
                theme(legend.position = "top", legend.title = element_blank()) +
                scale_y_continuous(breaks = round(seq(min(0), max(100), by = 5),1)) 
        }
        
        if ("None" %in% input$Stratify_by){
            plot <- plot
        }
        
        if ("Sponsor" %in% input$Stratify_by) {
            plot <- plot + 
                facet_wrap(~Sponsor) + 
                scale_y_continuous(breaks = round(seq(min(0), max(100), by = 10),1))
        }
        
        if ("Indication" %in% input$Stratify_by) {
            plot <- plot + 
                facet_wrap(~TA_subgroup) + 
                scale_y_continuous(breaks = round(seq(min(0), max(100), by = 10),1))
        }
        
        plot <- plot + 
            theme(axis.text.x = element_text(size=8),
                  axis.text.y = element_text(size=8),
                  axis.title.x=element_text(size=12, face="bold"),
                  axis.title.y=element_text(size=12, face="bold"))
        
        #plot
        ggplotly(plot) %>% 
            layout(legend = list(orientation = "h", y = 1.1, x = 0.03))
        
    })
    
    output$TA_Disease_Burden_Comparison <- renderPlot({
        
        # Plot stacked bar
        # https://www.datanovia.com/en/blog/awesome-list-of-hexadecimal-colors-you-should-have/
        group.colors <- c("Over represented" = "#6699FF", # Darker Blue
                          "Appropriately represented" = "#99CCFF", # Medium Blue
                          "Under represented" = "#CC3333", #Red
                          "Demographic not collected in trial" = "#999999") #Grey
        
        plot <- disease_burden_comparison() %>% 
            ggplot(aes(fill=factor(Comparison, levels = c("Demographic not collected in trial",
                                                          "Under represented",
                                                          "Appropriately represented",
                                                          "Over represented")), 
                                   y=Count, 
                                   x=Demographic, 
                                   label = Count)) + 
            geom_bar(position="stack", stat="identity") +
            facet_wrap(~Indication) + 
            geom_text(size = 7, position = position_stack(vjust = 0.5)) + 
            scale_fill_manual(values=group.colors) + 
            theme(axis.text = element_text(size=16),
                  axis.title=element_text(size=16, face="bold"),
                  strip.text.x = element_text(size = 16),
                  legend.text = element_text(size=14),
                  legend.position = "top",
                  legend.title = element_blank())
        
        plot
        
    })
    
    output$Disease_Burden_table <- DT::renderDataTable(
        
        disease_burden_table()

    )
    
    output$TA_change_over_time <- renderPlot({

        plot <- approvals_TA() %>% 
            #filter(Demographic %in% dems) %>% 
            filter(Percentage != "NA") %>% 
            ggplot(aes(x = Demographic, y = Percentage, fill = factor(Approval_Year))) +
            geom_boxplot(outlier.shape = NA) +
            scale_fill_brewer(palette="PuRd")+
            geom_point(position=position_jitterdodge(),alpha=0.1) +
            scale_y_continuous(breaks = round(seq(0, 100, by = 5),1)) +
            theme(legend.position = "top",
                  legend.title = element_blank(),
                  legend.text = element_text(size=12),
                  axis.text.x = element_text(size=15),
                  axis.text.y = element_text(size=15),
                  axis.title.x=element_text(size=15, face="bold"), 
                  axis.title.y=element_text(size=15, face="bold"))
                  
        plot
        
    })

    output$TA_approvalsTable <- DT::renderDataTable(
        approvals_TA() %>% 
            #select(-Enrollment_bucket) %>%
            pivot_wider(names_from = Demographic, values_from = Percentage), options=list(pageLength=10)
    )
    
}

# Run the application 
shinyApp(ui = ui, server = server)
