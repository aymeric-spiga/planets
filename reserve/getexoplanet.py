#! /usr/bin/env python
import urllib,re
import numpy

###########################################
name="hd_189733_b"
#name="gl_581_c"
###########################################


######################
url = "http://exoplanet.eu/catalog/"+name+"/"
html = urllib.urlopen(url).read()

######################
find = re.compile(r'<td class="label">') ; html = find.sub("",html)
find = re.compile(r"</td>\n") ; html = find.sub(" ___ ",html)
find = re.compile(r"<td>") ; html = find.sub(" ___ ",html)
find = re.compile(r"&mdash;") ; html = find.sub(" NA ",html)
find = re.compile(r"<a href") ; html = find.sub(" \n YORGL ",html)
ff = open("temp.txt",'w') ; ff.write(html) ; ff.close()

######################
ff = open("temp.txt",'r')
ff2 = open("temp2.txt",'w')
for line in ff:
    if re.search(r'___', line) is not None:
      if re.search(r'YORGL', line) is None:
         line = " ".join(line.split())
         find = re.compile(r"___ ___") ; line = find.sub("YYYY",line)
         find = re.compile(r" ___ </tr>") ; line = find.sub("",line)
         line = line + '\n'
         ########
         try: l1,l2,l3 = line.strip().split('YYYY')
         except: 
           try: l1,l2 = line.strip().split('YYYY') ; l3 = ''
           except: l1 = line.strip().split('YYYY') ; l2 = '' ; l3 = ''
         ########
         if l2 == '': 
           pass
         else: 
           line = 'ZZZZ %s / %s \n' % (l1,l2)
           find = re.compile(r" M<sub>") ; line = find.sub("\n",line)
           find = re.compile(r" \(") ; line = find.sub("\n",line)
           find = re.compile(r" AU") ; line = find.sub("\n",line)
           ff2.write(line)        
ff.close()
ff2.close()

######################
ff2 = open("temp2.txt",'r')
ff3 = open(name+".txt",'w')
for line in ff2:
   if 'ZZZZ' in line:
       find = re.compile(r"ZZZZ") ; line = find.sub("",line)
       ff3.write(line) 
ff2.close()
ff3.close()

               


exit()


