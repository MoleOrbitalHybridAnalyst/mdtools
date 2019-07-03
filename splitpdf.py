from PyPDF2 import PdfFileWriter, PdfFileReader
from sys import argv

inputpdf = PdfFileReader(open(argv[1], "rb"))

for i in range(inputpdf.numPages):
    output = PdfFileWriter()
    output.addPage(inputpdf.getPage(i))
    with open("page%.2d.pdf"%i, "wb") as outputStream:
        output.write(outputStream)
