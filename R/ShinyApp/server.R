# Rely on the 'WorldPhones' dataset in the datasets
# package (which generally comes preloaded).
library(shiny)
library(ggplot2)
library(pROC)
library(plotROC)



ShinyData = read.csv('ShinyData_cancer_diagnosis.csv', row.names = 1)


ShinyData_gen = read.csv('ShinyData_cancer_diagnosis_genomics.csv', row.names = 1)
ShinyData_CPTAC = read.csv('ShinyData_cancer_diagnosis_CPTAC.csv', row.names = 1)

ShinyData_CPTAC$Disease = unlist(lapply(ShinyData_CPTAC$Disease, function(x) strsplit(x, split = ' ')[[1]][1]))
ShinyData_CPTAC[ShinyData_CPTAC=="UCEC"]<- 'Uterine'

ShinyData[ShinyData=="met"]<- 'metabolomics'
ShinyData[ShinyData=="prot"]<- 'proteomics'
ShinyData_CPTAC[ShinyData_CPTAC=="Ovary"]<- 'Ovarian'


ShinyData_genomics = data.frame()
for(disease in unique(ShinyData$Disease)){
  disease_dataset = read.csv(paste('use_', disease, '_plink.txt', sep = ''), sep = '\t')
  colnames(disease_dataset) = c('SNP', 'effect_allele', 'effect_weight')
  disease_dataset$Disease = disease
  ShinyData_genomics = rbind(ShinyData_genomics, disease_dataset)
}

server <- function(input, output) {
  # Define a reactive expression for the document term matrix
  
  
  output$plot1 <- renderPlot({
    
    if(input$omics == 'genomics'){
      
      ShinyData_genomics_subset = ShinyData_genomics[ShinyData_genomics['Disease'] == input$disease,]
      print(dim(ShinyData_genomics_subset))
      
      ShinyData_genomics_subset = ShinyData_genomics_subset[sort(abs(ShinyData_genomics_subset$effect_weight),decreasing=T,index.return=T)[[2]],]
      ShinyData_genomics_subset = ShinyData_genomics_subset[1:min(input$N, dim(ShinyData_genomics_subset)[1]) ,]
      print(dim(ShinyData_genomics_subset))
      
      PlotData = data.frame('Features' = ShinyData_genomics_subset$SNP,
                            'Coef' = ShinyData_genomics_subset$effect_weight)
    }
    else if(input$omics == 'tissue'){
      SelectedRow = ShinyData_CPTAC['Disease'] == input$disease &
        ShinyData_CPTAC['N'] == input$N 

      Features = ShinyData_CPTAC[SelectedRow, 'Features']
      Coef = ShinyData_CPTAC[SelectedRow, 'Coef']
      
      # Remove brackets and single quotes and split the string
      Features <- unlist(strsplit(Features, "', '"))
      Features[1] <- gsub("\\['", "",Features[1])
      Features[length(Features)] <- gsub("'\\]", "",Features[length(Features)])
      
      Coef <- unlist(strsplit(Coef, ", "))
      Coef[1] <- gsub("\\[", "",Coef[1])
      Coef[length(Coef)] <- gsub("\\]", "",Coef[length(Coef)])
      Coef <- as.numeric(Coef)
      
      PlotData = data.frame('Features' = Features, 'Coef' = Coef)
    } 
    
    else{
      SelectedRow = ShinyData['Disease'] == input$disease &
        ShinyData['N'] == input$N &
        ShinyData['Data'] == input$omics &
        ShinyData['control_data'] == input$HC_type
        
      Features = ShinyData[SelectedRow, 'Features']
      Coef = ShinyData[SelectedRow, 'Coef']
        
      # Remove brackets and single quotes and split the string
      Features <- unlist(strsplit(Features, "', '"))
      Features[1] <- gsub("\\['", "",Features[1])
      Features[length(Features)] <- gsub("'\\]", "",Features[length(Features)])
        
      Coef <- unlist(strsplit(Coef, ", "))
      Coef[1] <- gsub("\\[", "",Coef[1])
      Coef[length(Coef)] <- gsub("\\]", "",Coef[length(Coef)])
      Coef <- as.numeric(Coef)
        
      PlotData = data.frame('Features' = Features, 'Coef' = Coef)
    }
    
    ggplot(PlotData, aes(x=Coef, y=Features)) +
      geom_segment( aes(y=Features, yend=Features, x=0, xend=Coef), color="grey") +
      geom_point( color="orange", size=4) +
      theme_light() +
      theme(
        panel.grid.major.y = element_blank(),
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 12)
      ) +
      ylab("") +
      xlab("Coefficient")
  })
  
  
  output$plot2 <- renderPlot({

    if(input$omics == 'genomics'){
      SelectedRow = ShinyData_gen['Disease'] == input$disease &
        ShinyData_gen['control_data'] == input$HC_type
      
      Prob = ShinyData_gen[SelectedRow, 'prob']
      Reference = ShinyData_gen[SelectedRow, 'reference']
    }
    else if(input$omics == 'tissue'){
      SelectedRow = ShinyData_CPTAC['Disease'] == input$disease &
        ShinyData_CPTAC['N'] == input$N 

      Prob = ShinyData_CPTAC[SelectedRow, 'prob']
      Reference = ShinyData_CPTAC[SelectedRow, 'reference']
    }
    else{
      SelectedRow = ShinyData['Disease'] == input$disease &
        ShinyData['N'] == input$N &
        ShinyData['Data'] == input$omics &
        ShinyData['control_data'] == input$HC_type
      
      Prob = ShinyData[SelectedRow, 'prob']
      Reference = ShinyData[SelectedRow, 'reference']
    }
    

    # Remove brackets and single quotes and split the string
    Reference <- unlist(strsplit(Reference, "', '"))
    Reference[1] <- gsub("\\['", "",Reference[1])
    Reference[length(Reference)] <- gsub("'\\]", "",Reference[length(Reference)])
    
    Prob <- unlist(strsplit(Prob, ", "))
    Prob[1] <- gsub("\\[", "",Prob[1])
    Prob[length(Prob)] <- gsub("\\]", "",Prob[length(Prob)])
    Prob <- as.numeric(Prob)
    

    
    PlotData = data.frame('Reference' = Reference, 'Prob' = Prob)
    
    print(dim(PlotData))
    basicplot <- ggplot(PlotData, aes(d = Reference, m = Prob)) +
      geom_roc(n.cuts = 0)+
      style_roc(theme = theme_gray())
    basicplot
  })
}

server
