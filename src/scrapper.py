from urllib.request import urlopen

url = 'https://celestrak.org/NORAD/elements/gp.php?GROUP=last-30-days&FORMAT=tle'

page = urlopen(url)
html = page.read().decode("utf-8")
print(html)