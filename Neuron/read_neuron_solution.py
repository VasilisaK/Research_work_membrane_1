import re

f_input = open('input.txt', 'r')
text = f_input.read().strip()
f_input.close()

m1 = re.search('GraphicsComplex', text)
text1 = text[m1.end():]
m2 = re.search('RGBColor', text1)
text2 = text1[:m2.start()].strip()
text3 = text2[2:-6]
text4 = text3.replace('*^', '*10^')
text5 = re.findall(r"\{(.*?)\}", text4) 

f_output = open('output.txt', 'w')
for elem in text5:
    f_output.write(elem.replace(',','') + '\n')
f_output.close()
