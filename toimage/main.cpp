#include <QImage>
#include <fstream>
#include <iostream>
#include <QImageWriter>
#include <QCoreApplication>

const int WH = 512;
const int size = 64;
const QString name = "data_train";
const int cntsq = WH/size;

int read_int(std::ifstream& in) {
    int val;
    in.read(reinterpret_cast<char *>(&val), sizeof(val));
    return val;
}

QColor Image[WH][WH];
int main(int argc, char *argv[]) {
    QCoreApplication a(argc, argv);

    QString img_read_path = R"(C:\Users\Main\Base\huawei\parsed\)" + name + "\\" +
            QString::number(size) + "\\1205.txt";
    QString ans_read_path = R"(C:\Users\Main\Base\huawei\parsed\)" + name + "\\" +
            QString::number(size) + "\\answers.txt";
    QString img_write_path = R"(C:\Users\Main\Base\huawei\RESULT.png)";

    std::cout << img_read_path.toStdString() << std::endl;

    std::ifstream img_in(img_read_path.toStdString(), std::ios::binary);
    for (int i = 0; i < WH; ++i)
        for (int j = 0; j < WH; ++j) {
            int cur = read_int(img_in);

            int r = cur>>16,
                g = (cur>>8)&255,
                b = cur&255;

            Image[i][j] = QColor(r, g, b);
        }

    std::cout << "read image" << std::endl;
    QImage res = QImage(WH, WH, QImage::Format_RGB32);
    img_in.close();

    std::ifstream ans_in(ans_read_path.toStdString());
    std::string trash;
    std::getline(ans_in, trash);

    std::cout << "reading permutation" << std::endl;
    for (int j = 0; j < cntsq; ++j, std::cout << std::endl)
        for (int i = 0; i < cntsq; ++i){
            int pos; ans_in >> pos;
            std::cout << pos << " ";

            int ni = pos/cntsq, nj = pos%cntsq;

            for (int x = 0; x < size; ++x)
                for (int y = 0; y < size; ++y) {
                    res.setPixelColor(i*size+y, j*size+x, Image[nj*size+y][ni*size+x]);
                }
        }

    QImageWriter writer(img_write_path);
    std::cout << "came here" << std::endl;
    writer.write(res);

}
