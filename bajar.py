import requests

regionPoints = [4.500833, 4.466944, -75.404720, -75.370000, "CERRO MACHÍN"]



def downloadData(regionPoints):
    variable = []
    for j in range(2):
        i = regionPoints[j]
        if i>0:
            if int(i)<10:
               variable.append('N0'+str(int(i)))
            else:
               variable.append('N'+str(int(i)))
        if i<0:
            if int(i)<10:
               variable.append('S0'+str(abs(int(i))))
            else:
               variable.append('S'+str(abs(int(i))))
    for j in range(2,4):
        i = regionPoints[j]
        print(len(str(int(i))))
        print(str(int(i)))
        if i>0:
            if len(str(abs(int(i))))<3:
               variable.append('E0'+str(int(i)))
            else:
               variable.append('E'+str(int(i)))
        if i<0:
            if len(str(abs(int(i))))<3:
               variable.append('W0'+str(abs(int(i))))
            else:
               variable.append('W'+str(abs(int(i))))

    return list(set(variable))

def download_file(downloadUrl):
    req = requests.get(downloadUrl)
    filename = req.url[downloadUrl.rfind('/')+1:]

    with open(filename, 'wb') as f:
        for chunk in req.iter_content(chunk_size=8192):
            if chunk:
                f.write(chunk)



coor = downloadData(regionPoints)
print(coor[0])
print(coor[1])
downloadLink = 'http://step.esa.int/auxdata/dem/SRTMGL1/'+coor[0]+coor[1]+'.SRTMGL1.hgt.zip'
downloadLink = 'http://step.esa.int/auxdata/dem/SRTMGL1/N45E003.SRTMGL1.hgt.zip'
print(downloadLink)
download_file(downloadLink)   
