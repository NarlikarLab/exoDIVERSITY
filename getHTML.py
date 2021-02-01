import sys
import os,commands
import re
def makehtml(outdir,modelno,hmflag):
    f = outdir+"/"+str(modelno)+"modes/"+str(modelno)+"modes.html"
    of = open(f,'w')
    header = "<!DOCTYPE html>\n<html>\n<meta charset=\"UTF-8\">\n\t<body>\n"
    footer = "\t</body>\n</html>"
    of.write(header)
    of.write("\t\t<div style\"margin:10px\">\n\t\t\t<h3> Output for model with "+str(modelno)+" modes</h3>\n\t\t</div>\n")
    if (hmflag):
        captions = ["Sequence alignments","Positive strand reads","Negative strand reads"]
        filenames = ["alignedMotifs.png","posreadsHeatmap.png","negreadsHeatmap.png"]
        altnames = ["dna","posreads","negreads"]
        of.write("\t\t<div class=\"row\">\n")
        for i in range(3):
            of.write("\t\t\t<div class=\"column\">\n\t\t\t\t<figure>\n\t\t\t\t<figcaption>"+captions[i]+"</figcaption>\n\t\t\t\t<img src=\""+filenames[i]+"\" alt=\""+altnames[i]+"\" style=\"width: 90%\">\n\t\t\t\t</figure>\n\t\t\t</div>\n")
        of.write("\t\t</div>\n")

    of.write("\t\t<div style=\"clear: both;\">\n")
    of.write("\t\t\t<p>For details about the sequences and the modes please refer to <a href=\"info.txt\">info.txt</a></p>\n")
    of.write("\t\t\t<h3>Motif logos and read profiles:</h3>\n")
    of.write("\t\t</div>\n\t\t<div style=\"margin: 10px\">\n\t\t\t<table style=\"width: 100%\">\n")
    of.write("\t\t\t\t<tr>\n\t\t\t\t\t<th>Motif</th>\n\t\t\t\t\t<th>Motif reverse complement</th>\n\t\t\t\t\t<th>Reads</th>\n\t\t\t\t</tr>\n")
    for i in range(modelno):
        if not(os.path.isfile(outdir+"/"+str(modelno)+"modes/logo_"+str(i)+".png")): continue
        status,out = commands.getstatusoutput("file "+outdir+"/"+str(modelno)+"modes/logo_"+str(i)+".png")
        if status:
            print status,out
            exit()
        w,h = map(int,re.findall('(\d+)\s*x\s*(\d+)',out)[-1])
        logowidth = w*0.75
        logoheight = h*0.75

        of.write("\t\t\t\t<tr>\n")
        of.write("\t\t\t\t\t<td><img src=\"logo_"+str(i)+".png\" width=\""+str(logowidth)+"\"px height=\""+str(logoheight)+"\"px></td>\n")
        of.write("\t\t\t\t\t<td><img src=\"logo_"+str(i)+"_rc.png\" width=\""+str(logowidth)+"\"px height=\""+str(logoheight)+"\"px></td>\n")
        of.write("\t\t\t\t\t<td><img src=\"reads_"+str(i)+".png\" width=\"400\"px height=\"300\"px></td>\n")
        of.write("\t\t\t\t</tr>\n")
    of.write("\t\t\t</table>\n\t\t</div>\n")


    if hmflag:
        style="\t\t<style type=\"text/css\">\n\t\t.column {\n\t\t\tfloat: left;\n\t\t\twidth: 25%;\n\t\t\tpadding: 1px;\n\t\t\ttext-align: center;\n\t\t}\n\t\t.row::after{\n\t\t\tcontent: \"\";\n\t\t\tclear: both;\n\t\t\tdisplay: inline-block;\n\t\t}\n\t\tfigcaption{\n\t\t\tfont-weight: bold;\n\t\t\ttext-align:center;\n\t\t}\n\t\ttd{ \n\t\t\ttext-align: center;\n\t\t}\n\t\t</style>\n"
    else:
        style="\t\t<style type=\"text/css\">\n\t\ttd{\n\t\t\ttext-align: center;\n\t\t}\n\t\t</style>\n"
    of.write(style)
    of.write(footer)
    of.close()

def makeBestmodelHTML(outdir,bm,bestmodel,hmflag):
    modelno = int(re.findall(r"(\d+)",bm)[0])
    f = outdir+"/exoDiversity.html"
    of = open(f,'w')
    header = "<!DOCTYPE html>\n<html>\n<meta charset=\"UTF-8\">\n\t<body>\n"
    footer = "\t</body>\n</html>"
    of.write(header)
    of.write("\t\t<div style\"margin:10px\">\n\t\t\t<h3> Output for best model with "+str(modelno)+" modes</h3>\n\t\t</div>\n")

    if (hmflag):
        captions = ["Sequence alignments","Positive strand reads","Negative strand reads"]
        filenames = [bm+"/alignedMotifs.png",bm+"/posreadsHeatmap.png",bm+"/negreadsHeatmap.png"]
        altnames = ["dna","posreads","negreads"]
        of.write("\t\t<div class=\"row\">\n")
        for i in range(3):
            of.write("\t\t\t<div class=\"column\">\n\t\t\t\t<figure>\n\t\t\t\t<figcaption>"+captions[i]+"</figcaption>\n\t\t\t\t<img src=\""+filenames[i]+"\" alt=\""+altnames[i]+"\" style=\"width: 90%\">\n\t\t\t\t</figure>\n\t\t\t</div>\n")
        of.write("\t\t</div>\n")

    of.write("\t\t<div style=\"clear: both;\">\n")
    of.write("\t\t\t<p>For details about the sequences and the modes please refer to <a href=\""+bm+"/info.txt\">info.txt</a></p>\n")
    of.write("\t\t\t<p>The output for other models can be viewed from <a href=\"links.html\">Other models</a></p>\n")


    of.write("\t\t\t<h3>Motif logos and read profiles:</h3>\n")
    of.write("\t\t</div>\n\t\t<div style=\"margin: 10px\">\n\t\t\t<table style=\"width: 100%\">\n")
    of.write("\t\t\t\t<tr>\n\t\t\t\t\t<th>Motif</th>\n\t\t\t\t\t<th>Motif reverse complement</th>\n\t\t\t\t\t<th>Reads</th>\n\t\t\t\t</tr>\n")
    for i in range(modelno):
        if not(os.path.isfile(outdir+"/"+str(modelno)+"modes/logo_"+str(i)+".png")): continue
        status,out = commands.getstatusoutput("file "+outdir+"/"+str(modelno)+"modes/logo_"+str(i)+".png")
        if status:
            print status,out
            exit()
        w,h = map(int,re.findall('(\d+)\s*x\s*(\d+)',out)[-1])
        logowidth = w*0.75
        logoheight = h*0.75

        of.write("\t\t\t\t<tr>\n")
        of.write("\t\t\t\t\t<td><img src=\""+bm+"/logo_"+str(i)+".png\" width=\""+str(logowidth)+"\"px height=\""+str(logoheight)+"\"px></td>\n")
        of.write("\t\t\t\t\t<td><img src=\""+bm+"/logo_"+str(i)+"_rc.png\" width=\""+str(logowidth)+"\"px height=\""+str(logoheight)+"\"px></td>\n")
        of.write("\t\t\t\t\t<td><img src=\""+bm+"/reads_"+str(i)+".png\" width=\"400\"px height=\"300\"px></td>\n")
        of.write("\t\t\t\t</tr>\n")
    of.write("\t\t\t</table>\n\t\t</div>\n")


    if hmflag:
        style="\t\t<style type=\"text/css\">\n\t\t.column {\n\t\t\tfloat: left;\n\t\t\twidth: 25%;\n\t\t\tpadding: 1px;\n\t\t\ttext-align: center;\n\t\t}\n\t\t.row::after{\n\t\t\tcontent: \"\";\n\t\t\tclear: both;\n\t\t\tdisplay: inline-block;\n\t\t}\n\t\tfigcaption{\n\t\t\tfont-weight: bold;\n\t\t\ttext-align:center;\n\t\t}\n\t\ttd{ \n\t\t\ttext-align: center;\n\t\t}\n\t\t</style>\n"
    else:
        style="\t\t<style type=\"text/css\">\n\t\ttd{\n\t\t\ttext-align: center;\n\t\t}\n\t\t</style>\n"
    of.write(style)
    of.write(footer)
    of.close()


def makelinksHTML(outdir,minMode,maxMode):
    of = open(outdir+"/links.html","w")
    header = "<!DOCTYPE html>\n<html>\n<meta charset=\"UTF-8\">\n\t<body>\n"
    footer = "\t</body>\n</html>"
    of.write(header)
    of.write("\t<p>Links to HTML files for different models</p>\n")
    of.write("\t<div style=\"clear: both;\">\n")
    for m in range(minMode,maxMode+1):
        of.write("\t<a href="+str(m)+"modes/"+str(m)+"modes.html"+"> Model with "+str(m)+"modes</a><br>\n")
    of.write("\t</div>\n")
    of.write(footer)
    of.close()

if __name__=='__main__':
    outdir = sys.argv[1]
    modelno = int(sys.argv[2])
    hmflag = 1
    #makehtml(outdir,modelno,hmflag)
    makeBestmodelHTML(outdir,"2modes","bestModel_2modes",1)
    makelinksHTML(outdir,2,5)
