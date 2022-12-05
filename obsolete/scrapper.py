from urllib.request import urlopen

url_candidates = {
    'last30days': 'https://celestrak.org/NORAD/elements/gp.php?GROUP=last-30-days&FORMAT=tle'
}


def getRawTLEPage(url):
    page = urlopen(url)
    html = page.read().decode("utf-8")
    return html
# print(html)
