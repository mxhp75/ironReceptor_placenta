## Required packages
req_packages <- c("REDCapR","Hmisc","progress") #required packages
if(!all(req_packages%in%installed.packages())){install.packages(req_packages[!req_packages%in%installed.packages()],repos="http://cran.us.r-project.org")} #install if required
sapply(req_packages,require,character.only=T) #load packages


#remove HTML tags
#https://stackoverflow.com/questions/17227294/removing-html-tags-from-a-string-in-r
rmHTMLtags <- function(htmlString) {
  return(gsub("<.*?>", "", htmlString))
}

# get_choices
#select_choices_or_calculations
get_choices <- function(lab_code){
	lab_str <- strsplit(lab_code,"\\|")
	lab_str <- lapply(lab_str,trimws)

	label <- unlist(lapply(lab_str,function(x) trimws(gsub("^(.*),","",x))))
	code <- unlist(lapply(lab_str,function(x) trimws(gsub("^(.*?),.*","\\1",x))))
	return(data.frame(code,label))
}


## apply_data_dict
apply_data_dict <- function(rc_data,rc_metadata){
	rc_metadata2 <- rc_metadata
	if(all(!names(rc_metadata)%in%c("field_name","field_type","field_label","select_choices_or_calculations"))){
		names(rc_metadata2)[c(1,4,5,6)]<-c("field_name","field_type","field_label","select_choices_or_calculations")
		cat("\nData dictionary CSV import...\n")}
	rc_metadata2$field_label <- rmHTMLtags(rc_metadata2$field_label)
	rc_metadata2$field_label <- unlist(lapply(rc_metadata2$field_label,function(x) gsub("\"","",x))) #remove extra "s
	
	# set factors (radio/dropdown)
	choices_vars <- rc_metadata2$field_name[rc_metadata2$field_type%in%c("dropdown","radio")]
	if(length(choices_vars)>0){
	cat("\nSetting factors (radio/dropdown)\n")
	pb2 <- txtProgressBar(min=0,max=length(choices_vars),initial=0,style=3,width=60)
	for(k in 1:length(choices_vars)){
		setTxtProgressBar(pb2,k)
		choices_dat <- get_choices(rc_metadata2$select_choices_or_calculations[rc_metadata2$field_name%in%choices_vars[k]])
		if(choices_vars[k] %in% names(rc_data)){
		eval(parse(text=paste0("rc_data$",choices_vars[k]," <- factor(rc_data$",choices_vars[k],",levels=choices_dat$code)")))
		eval(parse(text=paste0("levels(rc_data$",choices_vars[k],") <- choices_dat$label")))
	}}}

	# set factors (checkbox)
	checkbox_vars <- rc_metadata2$field_name[rc_metadata2$field_type%in%"checkbox"]
	if(length(checkbox_vars)>0){
	cat("\nSetting factors (checkbox)\n")
	pb3 <- txtProgressBar(min=0,max=length(checkbox_vars),initial=0,style=3,width=60)
	for(j in 1:length(checkbox_vars)){
		setTxtProgressBar(pb3,j)
		x <- checkbox_vars[j]
		if(x %in% names(rc_data)){ #only if data has variable
		choices_dat <- get_choices(rc_metadata2$select_choices_or_calculations[rc_metadata2$field_name%in%checkbox_vars[j]])
		for(j2 in 1:nrow(choices_dat)){
		  eval(parse(text=paste0("rc_data$",x,"___",choices_dat$code[j2]," <- factor(rc_data$",x,"___",choices_dat$code[j2],",levels=c(0,1))")))
		eval(parse(text=paste0("levels(rc_data$",x,"___",choices_dat$code[j2],") <- c(\"Unchecked\",\"Checked\")")))
	}}}}

	# set factor (yesno)
	yesno_vars <- rc_metadata2$field_name[rc_metadata2$field_type%in%"yesno"]
	if(length(yesno_vars)>0){
	cat("\nSetting factors (yes/no)\n")
	pb4 <- txtProgressBar(min=0,max=length(yesno_vars),initial=0,style=3,width=60)
	for(k in 1:length(yesno_vars)){
		setTxtProgressBar(pb4,k)
		x <- yesno_vars[k]
		if(x %in% names(rc_data)){ #only if data has variable
		eval(parse(text=paste0("rc_data$",x," <- factor(rc_data$",x,",levels=c(1,0))")))
		eval(parse(text=paste0("levels(rc_data$",x,") <- c(\"Yes\",\"No\")")))
	}}}
	
	# set labels
	cat("\nSetting variable labels\n")
	pb1 <- txtProgressBar(min=0,max=nrow(rc_metadata2),initial=0,style=3,width=60)
	for(i in 1:nrow(rc_metadata2)){
		setTxtProgressBar(pb1,i)
		x <- rc_metadata2$field_name[i]
		if(x %in% names(rc_data)){ #only if data has variable
		if(!rc_metadata2$field_type[i]%in%c("descriptive","checkbox")){
  		  eval(parse(text=paste0("label(rc_data$",x,") <- \"", rc_metadata2$field_label[rc_metadata2$field_name%in%x],"\"")))
		}
		if(rc_metadata2$field_type[i]%in%"checkbox"){
		  choices_dat <- get_choices(rc_metadata2$select_choices_or_calculations[rc_metadata2$field_name%in%x])
		  for(i2 in 1:nrow(choices_dat)){
			eval(parse(text=paste0("label(rc_data$",x,"___",choices_dat$code[i2],") <- \"",rc_metadata2$field_label[rc_metadata2$field_name%in%x]," (choice=",choices_dat$label[i2],")\"")))}
		}
	}}
	#label(rc_data)
	
	cat("\n\n\t\t>>>>>> DATA FORMAT COMPLETE <<<<<<\n")
	return(rc_data) #return formatted data
}


## Wrapper function with redcap_read
redcap_api_import <- function(redcap_url,token){
	cat("\n\t>>>>>>>>>>>>>> REDCap Import <<<<<<<<<<<<<<\n")
	rc_data <- redcap_read(redcap_uri=redcap_url,token=token)$data
	rc_metadata <- redcap_metadata_read(redcap_uri=redcap_url,token=token)$data
	cat("\n\t>>>>>>>>>>>>>> Imported data <<<<<<<<<<<<<<\n")
	cat("\tNumber of rows = ",nrow(rc_data),"\n\tNumber of columns = ",ncol(rc_data),"\n")
	cat("\n\t>>>>>>>>>>>>>> Apply Data Dictionary <<<<<<<<<<<<<<\n")
	output_dat <- apply_data_dict(rc_data,rc_metadata)
	return(output_dat)

}
                         
## From data dictionary file
#names(datadict)[c(1,4,5,6)] <- c("field_name","field_type","field_label","select_choices_or_calculations")
#apply_data_dict(main_dat,datadict)
