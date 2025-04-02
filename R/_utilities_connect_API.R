library("httr2")
library("jsonlite")

api_url <- "http://biodiv-app.biol.lu.se:8088/v1/"

req <- request(paste0(api_url, "api-info")) |> 
  req_headers( "Accept" = "application/json")

resp <- req_perform(req)

resp_body_json(resp)

Hmsc::computePredictedValues()