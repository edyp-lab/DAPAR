setClass("PipelineTemplate", contains='MultiAssayExperiment')
setClass("Pipelinepeptide", contains='PipelineTemplate')
setClass("PipelineProtein", contains='PipelineTemplate')
