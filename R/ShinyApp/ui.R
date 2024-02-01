
ui <- fluidPage(
  # Application title
  titlePanel("Multiomics Atlas of Biomarkers for cancer diagnostics"),

  sidebarLayout(
    # Sidebar with a slider and selection inputs
    sidebarPanel(
      selectInput("disease", "Disease",
                  choices = c("AllCancers", "Bladder", "Breast", "Colorectal",
                            "Leukaemia", "Lung", "Lymphoma", "Malignant_melanoma", "Ovarian",
                            "Prostate", "Thyroid", "Uterine", "Colon", "HeadNeck",
                            "Liver", "Pancreas") ),
      helpText("Please note that some cancers are available either only in UKBB or only in CPTAC"),
      
      selectInput("omics", "Omics",
                  choices = c('proteomics', 'metabolomics', 'genomics',
                              'clinical', 'CPTAC' = 'tissue')),
      
      
      radioButtons("HC_type", "Type of controls",
                   choices = c('HC' = "HC",
                               'Non-cancer' = "NonCancer"),
                   selected = "HC"),
      
      sliderInput("N",
                  "# Features",
                  min = 1,  max = 15, value = 5, step = 1),
      helpText("Slide to change the number of molecules")
      
    ),
    
    
    # Show Word Cloud
    mainPanel(fluidRow(
      splitLayout(cellWidths = c("50%", "50%"), plotOutput("plot1"), plotOutput("plot2"))
    ),
    fluidRow(column(12,
                    helpText(HTML('The figure on the left shows the coefficients for each biomarker. Biomarkers with positive coefficients may have disease-inducing roles, while negative coefficient indicate protective roles. 
                                        <br/>
                                        <br/>
                                       The figure on the right shows AUCs for different numbers of genomic, proteomic or metabolic potential biomarkers'))
                    
    ))
    )
    
  )
)