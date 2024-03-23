#!/bin/bash

for FILE in *.pdf; do 
  #apply pdfcrop and use the same filename for output
  pdfcrop $FILE $FILE
done;
