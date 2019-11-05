\name{DAPAR-NEWS}
\title{DAPAR News}



\section{CHANGES IN VERSION 1.11.15}{

\itemize{
   \item See changes in version 1.11.13 of the package Prostar
   }
}



\section{CHANGES IN VERSION 1.9.15}{

\subsection{BUG FIXES}{
\itemize{
   \item When the aggregation step has been performed, the interface switches to the first tab of the 'Descriptive Statistics' in order to view informations aout the new dataset (the protein one).
         
      \item  A new package (readxl) is used to read xls or xlxs files. In certain circumstances, the functions of the previsous package openxlsx is not able to decode properly Excel files.
      \item When converting a new (text or Excel) file in Prostar : the missing values were not registered as expected. Especially, they did not appear in blue in the table above the volcanoplot. Bug fixed
}
  }
  \subsection{NEW FEATURES}{
  
  \itemize{
      \item Implementation of a parallel version of the function which saves the (new) protein dataset after the aggregation step.
            \item  Enhancement of the string-based filtering UI
      \item The automatic generation of an analysis report has been integrated in the Dataset Manager (menu 'Export'). It allows the user to download plots and parameters used in Prostar ont their dataset. 
      \item Implementation of a parallel version of the function which saves the (new) protein dataset after the aggregation step.
            \item  Added a Gene Ontology (GO) analysis module in Data Processing. This module allows to perform GO grouping and GO Enrichment
      \item  Several plots are now based on the package highcharter which is a wrapper to the highcharts graphical library. It provides interactivity with the user.
      }
       }
       }
       