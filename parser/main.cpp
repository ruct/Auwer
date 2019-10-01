#include <QImage>
#include <fstream>
#include <iostream>
#include <QCoreApplication>

int read_int(std::ifstream& in) {
    int val;
    in.read(reinterpret_cast<char *>(&val), sizeof(val));
    return val;
}

inline void write_int(std::ofstream& out, int val) {
    out.write(reinterpret_cast<char*>(&val), sizeof(val));
}

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    int mx = 0;

    std::ios::sync_with_stdio(0);
    QString read_path = R"(C:\Users\Main\Base\huawei\data_test1_blank\64)";
    QString write_path = R"(C:\Users\Main\Base\huawei\parsed\data_test1_blank\64)";
    for (int n = 2400; n < 2700; ++n) {
        QString lol = QString::number(n);
        while (lol.size() < 4)
            lol = "0"+lol;
        lol = read_path+"\\"+lol+".png";
        QImage cur(lol);

        std::ofstream out((write_path+"\\"+QString::number(n)+".txt").toStdString(), std::ios::binary);
        for (int i = 0; i < 512; ++i)
            for (int j = 0; j < 512; ++j) {
                int r, g, b;
                QColor pix = cur.pixel(i, j);
                pix.getRgb(&r, &g, &b);

                write_int(out, r<<16|g<<8|b);
//                out << (r<<16|g<<8|b) << " ";
                mx = std::max(mx, r);
                mx = std::max(mx, g);
                mx = std::max(mx, b);
                /*write_int(out, r);
                write_int(out, g);
                write_int(out, b);*/
            }
    }

    return 0;
}
