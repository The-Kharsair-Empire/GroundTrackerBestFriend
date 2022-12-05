from obsolete import getRawTLEPage, url_candidates
from obsolete.tle import TLE


def main():
    # print(html)
    html = getRawTLEPage(url_candidates['last30days'])
    lines = html.split('\r\n')
    lines = list(map(lambda x:x.strip(), lines))
    # print(lines)
    tles = []
    for i in range(0, len(lines) - 2, 3):
        # print(lines[i])
        # print(lines[i+1])
        # print(lines[i+2])
        tles.append(TLE(lines[i], lines[i+1], lines[i+2], True))

    # print(tles[0].get_semi_major_axis())
        
if __name__ == '__main__':
    main()