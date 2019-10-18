#include <QImage>
#include <fstream>
#include <iostream>
#include <QImageWriter>
#include <QCoreApplication>

//  image properties
//  wh - square's side of picture in pix.
//  size - square's side of puzzle's fragment in pix.
//  name - folder with parsed images
//  nl, nr - segment of what to solve
//  cntsq - square's side of picture in fragments
const int WH = 512;
const int size = 64;
const QString name = "data_train";
const int cntsq = WH/size;
const int NL = 1200, NR = 1230;

//  read int from bytes
int read_int(std::ifstream& in) {
    int val;
    in.read(reinterpret_cast<char*>(&val), sizeof(val));
    return val;
}

QColor Image[WH][WH];
int main(int argc, char *argv[]) {
    QCoreApplication a(argc, argv);

    //  path - directory with parsed images
    //          and file for results
    QString path = R"(C:\Users\Main\Base\huawei\parsed\)";
    for (int NUMBER = NL; NUMBER <= NR; ++NUMBER) {
        QString img_write_path = path + "img_results\\"+QString::number(NUMBER)+".png";
        QString img_read_path = path + name + "\\" +
                QString::number(size) + "\\" + QString::number(NUMBER) + ".txt";
        QString ans_read_path = path + name + "\\" +
                QString::number(size) + "\\answers.txt";

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
        while (std::getline(ans_in, trash))
            if (trash == QString::number(NUMBER).toStdString()+".png")
                break;

        //  gathering the picture
        std::cout << "reading permutation" << std::endl;
        for (int j = 0; j < cntsq; ++j)
            for (int i = 0; i < cntsq; ++i){
                int pos; ans_in >> pos;
                if (!j) std::cout << pos << " ";

                int ni = pos/cntsq, nj = pos%cntsq;

                for (int x = 0; x < size; ++x)
                    for (int y = 0; y < size; ++y) {
                        res.setPixelColor(i*size+y, j*size+x, Image[nj*size+y][ni*size+x]);
                    }
            }
        std::cout << std::endl;

        //  writing to png
        QImageWriter writer(img_write_path);
        writer.write(res);
    }
}
