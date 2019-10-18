#include <QImage>
#include <fstream>
#include <iostream>
#include <QCoreApplication>

//  read int from bytes
int read_int(std::ifstream& in) {
    int val;
    in.read(reinterpret_cast<char *>(&val), sizeof(val));
    return val;
}

//  write int in bytes
inline void write_int(std::ofstream& out, int val) {
    out.write(reinterpret_cast<char*>(&val), sizeof(val));
}

//  image properties
//  wh - square's side of picture in pix.
//  size - square's side of puzzle's fragment in pix.
//  name - folder with parsed images
//  nl, nr - segment of what to solve
//  cntsq - square's side of picture in fragments
const int WH = 512;
const int size = 32;
const int NL = 600, NR = 1200;
const QString name = "data_train";
const int cntsq = WH/size;
int main(int argc, char *argv[]) {
    QCoreApplication a(argc, argv);

    int mx = 0;

    //  path - directory with parsed images
    //          and file for results
    QString path = R"(C:\Users\Main\Base\huawei\)";
    std::ios::sync_with_stdio(0);
    QString read_path = path+name+"\\"+QString::number(size);
    QString write_path = path+"parsed\\"+name+"\\"+QString::number(size);
    for (int n = NL; n < NR; ++n) {
        QString lol = QString::number(n);
        while (lol.size() < 4)
            lol = "0"+lol;
        lol = read_path+"\\"+lol+".png";
        QImage cur(lol);


        //  writing parsed image
        std::ofstream out((write_path+"\\"+QString::number(n)+".txt").toStdString(), std::ios::binary);
        for (int i = 0; i < 512; ++i)
            for (int j = 0; j < 512; ++j) {
                int r, g, b;
                QColor pix = cur.pixel(i, j);
                pix.getRgb(&r, &g, &b);

                write_int(out, r<<16|g<<8|b);
                mx = std::max(mx, r);
                mx = std::max(mx, g);
                mx = std::max(mx, b);
            }
    }

    return 0;
}
