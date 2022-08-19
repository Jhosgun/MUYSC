import requests

downloadUrl = 'http://step.esa.int/auxdata/dem/SRTMGL1/N04W072.SRTMGL1.hgt.zip'

req = requests.get(downloadUrl)
filename = req.url[downloadUrl.rfind('/')+1:]

with open(filename, 'wb') as f:
    for chunk in req.iter_content(chunk_size=8192):
        if chunk:
            f.write(chunk)
