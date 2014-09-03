#!/usr/bin/python
import sys
import os
import re
import tempfile

HIDE_FROM_SPAG = '!!!SPAGIGNORE'
HIDE_CONTINUE  = '!!!&&&'
def encodeCpp(lines):
    continuesPattern = re.compile('(\\\\\s*\n)')
    newLines = []

    newLine = ''
    continues = False
    for line in lines:
        line = line.replace('#',HIDE_FROM_SPAG+'#')
        m = re.search(continuesPattern, line)
        if (m != None):
            line = re.sub(continuesPattern,' ', line)
            newLine = newLine + line
            continues = True
        else:
            continues = False
            newLine = newLine + line
            newLines.append(newLine)
            newLine = ''
    return newLines

def unhide(lines):
    newLines = []
    for line in lines:
        newLines.append(line.replace(HIDE_FROM_SPAG,''))
    return newLines


# SPAG does not support MODULE nor CONTAINS statements.
# We skirt the issue by replacing with subroutine/end subroutine
# statements.
def encodeModule(lines):
    beginRegexp = re.compile('^(\s*)\s+(module|subroutine)\s+([^(]+)(.*)\n',re.IGNORECASE)
    endRegexp = re.compile('^(\s*end)\s+(module|subroutine)\s+([^(]+)',re.IGNORECASE)
    containsRegexp = re.compile('^\s*contains',re.IGNORECASE)
    newLines = [];
    programUnits = []
    for line in lines:
        if (beginRegexp.match(line)):
            m = re.search(beginRegexp,line)
            name = m.group(2)+"_"+m.group(3)
            programUnits.append(name)
            newLines.append("      subroutine "+name+m.group(4)+"\n")
        elif containsRegexp.match(line):
            name = programUnits.pop()
            newLines.append("      end subroutine "+name+" !contains \n")
            programUnits.append('TEMP_'+name)
        elif endRegexp.match(line):
            try:
                name = programUnits.pop()
            except:
                print "*******"
                print "uh oh: ",line
                print " "
            if re.match('^TEMP_',name) != None:
                newLines.append("      subroutine "+name+"\n")
            newLines.append("      end subroutine "+name+'\n')
        else:
            newLines.append(line)
    return newLines

def decodeModule(lines):
    beginRegexp = re.compile('^(\s*)subroutine\s+(subroutine|function|module)_(\S+.*)$',re.IGNORECASE)
    endRegexp = re.compile('^(\s*end)\s+subroutine\s+(temp_)?(subroutine|function|module)_?(\S+)',re.IGNORECASE)
    containsRegexp = re.compile('^(\s*end)\s+subroutine\s+(\S+)\s*!contains',re.IGNORECASE)
    ignoreRegexp   = re.compile('^\s*subroutine\s+temp_(\S*)',re.IGNORECASE)
    newLines = [];
    for line in lines:
        if (beginRegexp.match(line)):
            m = re.search(beginRegexp,line)
            newLines.append('      '+m.group(2)+" "+m.group(3) + "\n")
        elif containsRegexp.match(line):
            newLines.append("      contains"+"\n")
        elif ignoreRegexp.match(line):
            pass
        elif endRegexp.match(line):
            m = re.search(endRegexp,line)
            newLines.append("      end "+m.group(3)+" "+m.group(4)+"\n")
        else:
            newLines.append(line)
    return newLines

# Comment out F90 syntax unknown to SPAG, including
# proper handling of continuation lines and intervening
# comment and CPP lines.
def encodeUnknownF90(lines):
    unknownPattern = re.compile('\s*(use|procedure|(end)?\s*interface|public|(end)?\s*type|allocate|deallocate)',
                                re.IGNORECASE)
    continuePattern = re.compile('^     \S')
    commentPattern = re.compile('^\s*[!cC#]')
    newLines = []
    previousLineWasUnknown = False
    for line in lines:
        if previousLineWasUnknown and (continuePattern.match(line)!=None):
            line = HIDE_CONTINUE + line[6:]
        elif (commentPattern.match(line) != None):
            pass
        else:
            previousLineWasUnknown = (unknownPattern.match(line) != None)
        newLines.append(line)
    return newLines
    

def decodeUnknownF90(lines):
    continuePattern = re.compile('^'+HIDE_CONTINUE)
    commentPattern = re.compile('[!cC#]')
    trailingCommentPattern = re.compile('\S*\s*!')
    newLines = []
    for line in lines:
        if continuePattern.match(line):
            for j in range(len(newLines)-1,-1,-1):
                if commentPattern.match(newLines[j]) == None:
                    if trailingCommentPattern.search(newLines[j]):
                        newLines[j] = newLines[j].replace('!','& !')
                    else:
                        newLines[j] = newLines[j][:-1]+' &' + newLines[j][-1]
                    break
            line = '      ' + line[6:]
        newLines.append(line)
    return newLines

def runSpag(lines):
    inTmp = 'spagin.f'
    outTmp = 'spagout.F90'
    file = open(inTmp,mode='w')
    for line in lines:
        file.write(line)
    file.close()
    os.system('spag ' + inTmp + ' TO=' + outTmp)
    file = open(outTmp,"r")
    newLines = file.readlines()
    file.close()
    return newLines

def stripSpagComments(lines):
    newLines = []
    for line in lines:
        if re.match('^!\*==',line):
            continue
        newLines.append(line)
    return newLines
    
        


    
filters = []
filters.append(encodeCpp)
filters.append(encodeUnknownF90)
filters.append(encodeModule)
filters.append(runSpag)
filters.append(decodeModule)
filters.append(decodeUnknownF90)
filters.append(unhide)
filters.append(stripSpagComments)


inFileName = sys.argv[1]
outFileName = inFileName.replace('.f','.F90')

inFile = open(inFileName,"r")
lines = inFile.readlines()
inFile.close()
for filter in filters:
    lines = filter(lines)

outFile = open(outFileName,"w")
for line in lines:
    outFile.write(line)
outFile.close()    

os.system('my-indent '+outFileName)
