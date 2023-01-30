import os


def groundtrack():
    pass


def get_world_city_coor():
    data_file = os.path.join(os.path.dirname(__file__), '..', 'file', 'csv', 'world_cities.csv')
    with open(data_file, encoding="utf8") as f:
        lines = f.readlines()

    header = lines[0]
    print(header)
    cities = {}
    for line in lines[1:]:
        line = line.split(',')

        try:
            cities[line[1]] = (float(line[2]), float(line[3]))
        except Exception as e:
            print(e)

    return cities


# test
if __name__ == '__main__':
    get_world_city_coor()
