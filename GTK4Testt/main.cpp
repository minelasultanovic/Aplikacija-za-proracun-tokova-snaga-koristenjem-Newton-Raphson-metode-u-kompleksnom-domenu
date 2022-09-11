// Include gtk
#include <gtk/gtk.h>
#include <string>
#include <complex>
#include <iostream>
#include <vector>
#include <tuple>
#include <functional>
#include <algorithm>

using namespace std;

struct Node {
    string type;
    complex<double> values;
    vector<int> adjacentNodes;
    string dependence;
};

struct Admittance {
    int node1;
    int node2;
    complex<double> impedance;
};


vector<vector<complex<double>>> getCofactor(vector<vector<complex<double>>> A, int p, int q, int n) {
    int i = 0, j = 0;
    vector<vector<complex<double>>> temp(n, vector<complex<double>>(n, complex<double>(0, 0)));
    // Looping for each element of the matrix
    for (int row = 0; row < n; row++) {
        for (int col = 0; col < n; col++) {
            //  Copying into temporary matrix only those
            //  element which are not in given row and
            //  column
            if (row != p && col != q) {
                temp[i][j++] = A[row][col];

                // Row is filled, so increase row index and
                // reset col index
                if (j == n - 1) {
                    j = 0;
                    i++;
                }
            }
        }
    }

    return temp;
}



/* Recursive function for finding determinant of matrix.
  n is current dimension of A[][]. */
complex<double> determinant(vector<vector<complex<double>>> A, int n)
{
    complex<double> D = complex<double>(0); // Initialize result

    //  Base case : if matrix contains single element
    if (n == 1)
        return A[0][0];

    //vector<vector<complex<double>>> temp(n, vector<complex<double>>(n, complex<double>(0,0)));

    double sign = double(1); // To store sign multiplier

    // Iterate for each element of first row
    for (int f = 0; f < n; f++) {
        // Getting Cofactor of A[0][f]
        auto temp = getCofactor(A, 0, f, n);
        D += sign * A[0][f] * determinant(temp, n - 1);

        // terms are to be added with alternate sign
        sign = -sign;
    }

    return D;
}

vector<vector<complex<double>>> adjoint(vector<vector<complex<double>>> A)
{
    auto N = A.size();
    vector<vector<complex<double>>> adj(N, vector<complex<double>>(N, complex<double>(0, 0)));

    if (A.size() == 1) {
        adj[0][0] = 1;
        return adj;
    }

    // temp is used to store cofactors of A[][]
    double sign = double(1);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            // Get cofactor of A[i][j]
            auto temp = getCofactor(A, i, j, N);

            // sign of adj[j][i] positive if sum of row
            // and column indexes is even.
            sign = ((i + j) % 2 == 0) ? double(1) : double(-1);

            // Interchanging rows and columns to get the
            // transpose of the cofactor matrix
            adj[j][i] = (sign) * (determinant(temp, N - 1));
        }
    }

    return adj;
}

vector<vector<complex<double>>> inverse(vector<vector<complex<double>>> A) {
    // Find determinant of A[][]
    /*
        for(int i = 0; i<A.size(); i++) {
        for(int j=0; j<A[i].size(); j++) {
            cout<<A[i][j]<<" ";
        }
    cout<<endl;
    }
    */

    auto N = A.size();
    auto det = determinant(A, N);

    vector<vector<complex<double>>> inverse(N, vector<complex<double>>(N, complex<double>(0, 0)));

    if (det == complex<double>(0, 0)) {
        cout << "Singular matrix, can't find its inverse";
        return inverse;
    }

    // Find adjoint
    auto adj = adjoint(A);

    // Find Inverse using formula "inverse(A) =
    // adj(A)/det(A)"
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            inverse[i][j] = adj[i][j] / det;


    /*
    for (int i = 0; i < N; i++) {
        cout<<endl;
        for (int j = 0; j < N; j++) {
                cout<<inverse[i][j];
        }
    }
    */


    return inverse;
}

int brojIteracija;
double greska;

vector<vector<complex<double>>> computeHandleFunctions(vector<function<complex<double>(vector<complex<double>>)>> jednacine, vector<complex<double>> x_0) {
    auto n = jednacine.size();
    vector<vector<complex<double>>> y(n, vector<complex<double>>(n, complex<double>(0, 0)));

    for (int i = 0; i < n; i++) {
        y[i][0] = jednacine.at(i)(x_0);
    }

    return y;
}

function<complex<double>(vector<complex<double>>)> partialDerivativeEq(function<complex<double>(vector<complex<double>>)> f, int n, int i) {
    complex<double> dx = complex<double>(0.000001, 0);
    vector<complex<double>> dX(n, complex<double>(0, 0));
    dX.at(i) = dx;

    auto f1 = [dx, dX, f](vector<complex<double>> x) {
        vector<complex<double>> sum;

        for (int k = 0; k < x.size(); k++) {
            sum.push_back(x.at(k) + dX.at(k));
        }

        return (f(sum) - f(x)) / dx;
    };

    return f1;
}

vector<vector<complex<double>>> jacobian(vector<function<complex<double>(vector<complex<double>>)>> jednacine, vector<complex<double>> x_0) {
    vector<vector<complex<double>>> y(jednacine.size(), vector<complex<double>>(x_0.size(), complex<double>(0, 0)));

    for (int i = 0; i < jednacine.size(); i++) {
        for (int j = 0; j < x_0.size(); j++) {
            y[i][j] = partialDerivativeEq(jednacine[i], x_0.size(), j)(x_0);
        }
    }

    return y;
}

vector<complex<double>> getFirstColumn(vector<vector<complex<double>>> A) {
    vector<complex<double>> column(A.size());

    for (int i = 0; i < A.size(); i++) {
        column[i] = A[i][0];
    }

    return column;
}

vector<vector<complex<double>>> mnozenjeMatrica(vector<vector<complex<double>>> A, vector<vector<complex<double>>> B) {
    auto m = A.size();
    auto c = A.at(0).size();
    vector<vector<complex<double>>> mul(m, vector<complex<double>>(c, complex<double>(0, 0)));

    for (int p = 0; p < m; p++) {
        for (int j = 0; j < c; j++) {
            mul[p][j] = 0;
            for (int k = 0; k < c; k++) {
                mul[p][j] += A[p][k] * B[k][j];
            }
        }
    }

    return mul;
}

vector<vector<complex<double>>> sabiranjeMatrica(vector<vector<complex<double>>> A, vector<vector<complex<double>>> B) {
    auto m = A.size();
    auto c = A.at(0).size();
    vector<vector<complex<double>>> mul(m, vector<complex<double>>(c, complex<double>(0, 0)));

    for (int p = 0; p < m; p++) {
        for (int j = 0; j < c; j++) {
            mul[p][j] = A[p][j] + B[p][j];
        }
    }

    return mul;
}

vector<vector<complex<double>>> NewtonRaphson(vector<function<complex<double>(vector<complex<double>>)>> jednacine, vector<complex<double>> x_0, int maxIter, double eps) {
    brojIteracija = maxIter;
    vector<vector<complex<double>>> x_min(x_0.size(), vector<complex<double>>(x_0.size(), complex<double>(0, 0)));
    vector<vector<complex<double>>> x(x_0.size(), vector<complex<double>>(x_0.size(), complex<double>(0, 0)));

    for (int i = 0; i < x_0.size(); i++) {
        x[i][0] = x_0[i];
    }

    for (int i = 0; i < maxIter; i++) {
        vector<double> absValues = vector<double>();
        auto kolona = getFirstColumn(x);
        auto results = computeHandleFunctions(jednacine, kolona);

        for (int j = 0; j < results.size(); j++) {
            absValues.push_back(abs(results[j][0]));
        }

        auto maxValue = *max_element(begin(absValues), end(absValues));

        if (maxValue < eps) {
            //error
            brojIteracija = i;
            greska = maxValue;
            x_min = x;
            break;
        }

        ////Racunaj Newton-ov korak
        auto dx = mnozenjeMatrica(inverse(jacobian(jednacine, kolona)), results);
        for (int k = 0; k < dx.size(); k++) {
            //cout<<endl;
            for (int g = 0; g < dx[0].size(); g++) {
                dx[k][g] *= double(-1);
                // cout<<dx[k][g];
            }
        }

        x = sabiranjeMatrica(x, dx);
    }

    return x_min;
}

complex<double> getAdm(vector<Admittance> admittance, int n1, int n2) {
    complex<double> y = complex<double>(0, 0);
    if (n1 == n2) {
        for (int i = 0; i < admittance.size(); i++) {
            if (admittance.at(i).node1 == n1 || admittance.at(i).node2 == n1 || (admittance.at(i).node1 == n1 && admittance.at(i).node2 == -1)) {
                y = y + double(1) / admittance.at(i).impedance;
            }
        }
    }
    else {
        for (int i = 0; i < admittance.size(); i++) {
            if ((admittance.at(i).node1 == n1 && admittance.at(i).node2 == n2) || (admittance.at(i).node1 == n2 && admittance.at(i).node2 == n1)) {
                y = double(1) / admittance.at(i).impedance;
                break;
            }
        }
    }

    return y;
}

tuple<function<complex<double>(vector<complex<double>>)>, function<complex<double>(vector<complex<double>>)>> SLACKcvor(int i, int n, complex<double> V) {
    auto f1 = [i, V](vector<complex<double>> x) {
        return x.at(i) - V;
    };

    auto f2 = [i, n, V](vector<complex<double>> x) {
        return x.at(n + i) - conj(V);
    };

    return  { f1, f2 };
}

tuple<function<complex<double>(vector<complex<double>>)>, function<complex<double>(vector<complex<double>>)>>  PQ(vector<vector<complex<double>>> Y, int i, int n, complex<double> s_pq, string dependence) {
    function<complex<double>(vector<complex<double>>)> f1 = [](vector<complex<double>> x) {
        return complex<double>(0, 0);
    };

    for (int j = 0; j < n; j++) {
        f1 = [i, j, Y, f1, n](vector<complex<double>> x) {
            return conj(Y[i][j]) * x[n + j] + f1(x);
        };
    }

    function<complex<double>(vector<complex<double>>)> f2 = [](vector<complex<double>> x) {
        return complex<double>(0, 0);
    };

    for (int j = 0; j < n; j++) {
        f2 = [f2, j, i, Y](vector<complex<double>> x) {
            return Y[i][j] * x[j] + f2(x);
        };
    }
    f1 = [i, f1, s_pq](vector<complex<double>> x) {
        return x[i] * f1(x) - s_pq;
    };

    f2 = [f2, s_pq, i, n](vector<complex<double>> x) {
        return x[n + i] * f2(x) - conj(s_pq);
    };


    return { f1, f2 };
}

tuple<function<complex<double>(vector<complex<double>>)>, function<complex<double>(vector<complex<double>>)>>  PV(vector<vector<complex<double>>> Y, int i, int n, complex<double> pv) {
    function<complex<double>(vector<complex<double>>)> f1 = [](vector<complex<double>> x) {
        return complex<double>(0, 0);
    };

    for (int j = 0; j < n; j++) {
        f1 = [i, j, Y, f1, n](vector<complex<double>> x) {
            return conj(Y[i][j]) * x[n + j] + f1(x);
        };
    }

    function<complex<double>(vector<complex<double>>)> f2 = [](vector<complex<double>> x) {
        return complex<double>(0, 0);
    };
    for (int j = 0; j < n; j++) {
        f2 = [f2, j, i, Y](vector<complex<double>> x) {
            return Y[i][j] * x[j] + f2(x);
        };
    }

    f1 = [i, f1](vector<complex<double>> x) {
        return x[i] * f1(x);
    };

    f2 = [f2, i, n](vector<complex<double>> x) {
        return x[n + i] * f2(x);
    };

    auto f3 = [f1, f2, i, Y, pv](vector<complex<double>> x) {
        return (f1(x) + f2(x)) / double(2) - pv.real();
    };

    auto f4 = [i, n, pv](vector<complex<double>> x) {
        return x[i] * x[n + i] - (pv.imag() * pv.imag());
    };

    return { f3, f4 };
}

vector<vector<complex<double>>> admSistem(vector<Node> nodes, vector<Admittance> admittance) {
    vector<vector<complex<double>>> Y(nodes.size(), vector<complex<double>>(nodes.size(), complex<double>(0, 0)));

    for (int i = 0; i < nodes.size(); i++) {
        if (nodes.at(i).type == "SLACK") {
            Y[i][i] = complex<double>(1, 0);
        }
        else {
            for (int j = 0; j < nodes.at(i).adjacentNodes.size(); j++) {
                auto n = nodes.at(i).adjacentNodes.at(j);
                Y[i][n - 1] = double(-1) * getAdm(admittance, i + 1, n);
            }
            Y[i][i] = getAdm(admittance, i + 1, i + 1);
        }
    }

    return Y;
}

vector<complex<double>> getFinal(vector<Node> nodes, vector<Admittance> admittancees, vector<complex<double>> x_0, int maxIter, double eps) {
    auto n = nodes.size();
    vector<function<complex<double>(vector<complex<double>>)>> jednacine;
    auto Y = admSistem(nodes, admittancees);

    for (int i = 0; i < n; i++) {
        if (nodes.at(i).type == "SLACK") {
            function<complex<double>(vector<complex<double>>)> f1, f2;
            tie(f1, f2) = SLACKcvor(i, int(n), nodes.at(i).values);
            jednacine.push_back(f1);
            jednacine.push_back(f2);
        }
        else if (nodes.at(i).type == "PQ") {
            function<complex<double>(vector<complex<double>>)> f1, f2;
            tie(f1, f2) = PQ(Y, i, int(n), nodes.at(i).values, nodes.at(i).dependence);
            jednacine.push_back(f1);
            jednacine.push_back(f2);
        }
        else if (nodes.at(i).type == "PV") {
            function<complex<double>(vector<complex<double>>)> f1, f2;
            tie(f1, f2) = PV(Y, i, int(n), nodes.at(i).values);
            jednacine.push_back(f1);
            jednacine.push_back(f2);
        }
    }

    auto v = NewtonRaphson(jednacine, x_0, maxIter, eps);
    auto rjesenje = getFirstColumn(v);
    return  rjesenje;
}

vector<Node> nodes;
vector<Admittance> admittances;

/// 
/// GUI -------------------------------------------------------
/// 
/// 


GtkWidget* grid;
GtkWidget* t_0, * t_1, * nodeTypeSelector, * prvaVrijednostRealniDioEntryBox, * prvaVrijednostImaginarniDioEntryBox,
* susjedniCvoroviEntryBox, * dodajNodeButton;
GtkWidget* mojaLista, * admitansaLista, * vrstaCvoraLabel, * prvaVrijednostRealniDioLabel, * prvaVrijednostImaginarniDioLabel, * susjedniCvoroviLabel,
* drugaVrijednostRealniDioLabel, * drugaVrijednostImaginarniDioLabel, * drugaVrijednostRealniDioEntryBox, * drugaVrijednostImaginarniDioEntryBox;
GtkWidget* node1Label, * node2Label, * node1EntryBox, * node2EntryBox, * impedansaRealniDioLabel, * impedansaImaginarniDioLabel, * impedansaRealniDioEntryBox,
* impedansaImaginarniDioEntryBox, * dodajAdmitansuButton, * pocetneVrijednostiLabel, * pocetneVrijednostiEntryBox, * izracunajButton;
GtkWidget* rezultatLabel, * iteracijeGreskaLabel;


void setEntryText(GtkWidget* pEntry, const char* text, int len)
{
    GtkEntryBuffer* buffer = gtk_entry_get_buffer(GTK_ENTRY(pEntry));

    size_t nch = len;
    if (nch <= 0)
        nch = strlen(text);

    gtk_entry_buffer_set_text(buffer, text, (int)nch);
}

const char* getEntryText(GtkWidget* pEntry)
{
    GtkEntryBuffer* buffer = gtk_entry_get_buffer(GTK_ENTRY(pEntry));
    const char* text = gtk_entry_buffer_get_text(buffer);
    return text;
}

enum {

    LIST_ITEM = 0,
    N_COLUMNS
};

void init_list(GtkWidget* list) {

    GtkCellRenderer* renderer;
    GtkTreeViewColumn* column;
    GtkListStore* store;

    renderer = gtk_cell_renderer_text_new();
    column = gtk_tree_view_column_new_with_attributes("Lista:",
        renderer, "text", LIST_ITEM, NULL);
    gtk_tree_view_append_column(GTK_TREE_VIEW(list), column);

    store = gtk_list_store_new(N_COLUMNS, G_TYPE_STRING);

    gtk_tree_view_set_model(GTK_TREE_VIEW(list),
        GTK_TREE_MODEL(store));

    g_object_unref(store);

}

void add_to_list(GtkWidget* list, const gchar* str) {

    GtkListStore* store;
    GtkTreeIter iter;

    store = GTK_LIST_STORE(gtk_tree_view_get_model
    (GTK_TREE_VIEW(list)));

    gtk_list_store_append(store, &iter);
    gtk_list_store_set(store, &iter, LIST_ITEM, str, -1);
}

static void izracunajClicked() {
    string vrijednostiString(getEntryText(pocetneVrijednostiEntryBox));
    string delimiter = ",";

    vector<complex<double>> vrijednosti;

    size_t pos = 0;
    std::string token;
    while ((pos = vrijednostiString.find(delimiter)) != string::npos) {
        token = vrijednostiString.substr(0, pos);
        vrijednosti.push_back(complex<double>(stod(token), 0));
        vrijednostiString.erase(0, pos + delimiter.length());
    }
    vrijednosti.push_back(complex<double>(stod(vrijednostiString), 0));

    auto x_0 = vector<complex<double>>();

    for (int i = 0; i < vrijednosti.size(); i++) {
        x_0.push_back(vrijednosti.at(i));
    }

    for (int i = 0; i < vrijednosti.size(); i++) {
        x_0.push_back(conj(vrijednosti.at(i)));
    }

    auto rez = getFinal(nodes, admittances, x_0, 1000, 0.000001);
    string rezultatString = "Rezultat:\n";

    for (int i = 0; i < rez.size(); i++) {
        rezultatString += to_string(rez.at(i).real()) + " + i" + to_string(rez.at(i).imag()) + "\n";
    }

    gtk_label_set_text(GTK_LABEL(rezultatLabel), rezultatString.c_str());

    string iteracijeGreskaString = "Iteracije: " + to_string(brojIteracija) + "\nGreska: " + to_string(1000000 * greska) + "e-6";

    gtk_label_set_text(GTK_LABEL(iteracijeGreskaLabel), iteracijeGreskaString.c_str());
}

static void buttonClicked() {
    string nodeType(gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(nodeTypeSelector)));
    string prvaVrijednostRealni(getEntryText(prvaVrijednostRealniDioEntryBox));
    string prvaVrijednostImaginarni(getEntryText(prvaVrijednostImaginarniDioEntryBox));

    string cvorovi(getEntryText(susjedniCvoroviEntryBox));
    auto cvorovi2 = cvorovi;
    string delimiter = ",";

    vector<string> cvoroviVektor;
    size_t pos = 0;
    std::string token;
    while ((pos = cvorovi.find(delimiter)) != string::npos) {
        token = cvorovi.substr(0, pos);
        cvoroviVektor.push_back(token);
        cvorovi.erase(0, pos + delimiter.length());
    }

    cvoroviVektor.push_back(cvorovi);

    //string type;
    //vector<complex<double>> values;
    //vector<int> adjacentNodes;
    //string dependence;

    auto value = complex<double>(stod(prvaVrijednostRealni), stod(prvaVrijednostImaginarni));

    auto cvoroviInt = vector<int>();
    for (int i = 0; i < cvoroviVektor.size(); i++) {
        cvoroviInt.push_back(stoi(cvoroviVektor.at(i)));
    }

    Node node;
    node.type = nodeType;
    node.values = value;
    node.adjacentNodes = cvoroviInt;
    node.dependence = "no";

    string drugaVrijednostRealni, drugaVrijednostImaginarni;

    if (nodeType != "PQ") {
        drugaVrijednostRealni = getEntryText(drugaVrijednostRealniDioEntryBox);
        drugaVrijednostImaginarni = getEntryText(drugaVrijednostImaginarniDioEntryBox);
        if (nodeType == "PV")
            node.values = complex<double>(stod(prvaVrijednostRealni), stod(drugaVrijednostRealni));
        else
            node.values = complex<double>(stod(drugaVrijednostRealni), stod(drugaVrijednostImaginarni));
    }

    nodes.push_back(node);

    string konacni;
    if (nodeType == "PQ")
        konacni = node.type + "    " + "Vrijednosti: " + prvaVrijednostRealni + " + i" + prvaVrijednostImaginarni + "    " + "Susjedni cvorovi: " + cvorovi2;
    else if (nodeType == "PV")
        konacni = node.type + "    " + "Vrijednosti: " + prvaVrijednostRealni + ", " + drugaVrijednostRealni + "    " + "Susjedni cvorovi: " + cvorovi2;
    else
        konacni = node.type + "    " + "Vrijednosti: " + drugaVrijednostRealni + ", + i" + drugaVrijednostImaginarni + "    " + "Susjedni cvorovi: " + cvorovi2;
    add_to_list(mojaLista, konacni.c_str());
}

static void admitansaClicked() {
    string cvor1(getEntryText(node1EntryBox));
    string cvor2(getEntryText(node2EntryBox));
    string realniDio(getEntryText(impedansaRealniDioEntryBox));
    string imaginarniDio(getEntryText(impedansaImaginarniDioEntryBox));

    Admittance admit;
    admit.node1 = stoi(cvor1);
    admit.node2 = stoi(cvor2);
    admit.impedance = complex<double>(stod(realniDio), stod(imaginarniDio));
    admittances.push_back(admit);

    string konacni = "Cvor1: " + to_string(admit.node1) + " Cvor 2: " + cvor2 + "    " + "Impedansa: " + realniDio + " + i" + imaginarniDio;
    add_to_list(admitansaLista, konacni.c_str());
}

static void nodeTypeChanged() {
    string active(gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(nodeTypeSelector)));

    if (active != "PQ") {
        gtk_widget_show(drugaVrijednostRealniDioLabel);
        gtk_widget_show(drugaVrijednostImaginarniDioLabel);
        gtk_widget_show(drugaVrijednostRealniDioEntryBox);
        gtk_widget_show(drugaVrijednostImaginarniDioEntryBox);
    }
    else {
        gtk_widget_hide(drugaVrijednostRealniDioLabel);
        gtk_widget_hide(drugaVrijednostImaginarniDioLabel);
        gtk_widget_hide(drugaVrijednostRealniDioEntryBox);
        gtk_widget_hide(drugaVrijednostImaginarniDioEntryBox);
    }

}

static void on_activate(GtkApplication* app) {

    grid = gtk_grid_new();
    gtk_grid_set_column_homogeneous(GTK_GRID(grid), TRUE);
    gtk_grid_set_row_spacing(GTK_GRID(grid), 10);
    gtk_grid_set_column_spacing(GTK_GRID(grid), 10);

    //gtk_window_set_child(GTK_WINDOW(win), mainVertBox);
    //gtk_box_append((GtkBox*)mainVertBox, grid);

    //setEntryText(unost_1, "0.1", 3);

    //gtk_grid_attach(GTK_GRID(grid), t_0, 0, 1, 1, 1);
    //gtk_grid_attach(GTK_GRID(grid), unost_0, 1, i, 1, 1);
    //(GTK_GRID(grid), t_1, 0, 2, 1, 1);
    //gtk_grid_attach(GTK_GRID(grid), unost_1, 3, i++, 1, 1);
    int redIndex = 0;

    vrstaCvoraLabel = gtk_label_new("Vrsta cvora");
    gtk_grid_attach(GTK_GRID(grid), vrstaCvoraLabel, 0, redIndex, 1, 1);

    prvaVrijednostRealniDioLabel = gtk_label_new("S ili P realni dio");
    gtk_grid_attach(GTK_GRID(grid), prvaVrijednostRealniDioLabel, 1, redIndex, 1, 1);

    prvaVrijednostImaginarniDioLabel = gtk_label_new("S ili P imaginarni dio");
    gtk_grid_attach(GTK_GRID(grid), prvaVrijednostImaginarniDioLabel, 2, redIndex, 1, 1);

    susjedniCvoroviLabel = gtk_label_new("Susjedni cv. format: 1,2,3");
    gtk_grid_attach(GTK_GRID(grid), susjedniCvoroviLabel, 3, redIndex, 1, 1);

    redIndex++;
    nodeTypeSelector = gtk_combo_box_text_new();
    gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(nodeTypeSelector), "0", "PQ");
    gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(nodeTypeSelector), "1", "PV");
    gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(nodeTypeSelector), "2", "SLACK");
    gtk_grid_attach(GTK_GRID(grid), nodeTypeSelector, 0, redIndex, 1, 1);
    gtk_combo_box_set_active(GTK_COMBO_BOX(nodeTypeSelector), 0);

    prvaVrijednostRealniDioEntryBox = gtk_entry_new();
    gtk_grid_attach(GTK_GRID(grid), prvaVrijednostRealniDioEntryBox, 1, redIndex, 1, 1);

    prvaVrijednostImaginarniDioEntryBox = gtk_entry_new();
    gtk_grid_attach(GTK_GRID(grid), prvaVrijednostImaginarniDioEntryBox, 2, redIndex, 1, 1);

    susjedniCvoroviEntryBox = gtk_entry_new();
    gtk_grid_attach(GTK_GRID(grid), susjedniCvoroviEntryBox, 3, redIndex, 1, 1);

    dodajNodeButton = gtk_button_new_with_label("Dodaj cvor");
    gtk_grid_attach(GTK_GRID(grid), dodajNodeButton, 4, redIndex, 1, 1);


    redIndex++;

    drugaVrijednostRealniDioLabel = gtk_label_new("V realni dio");
    gtk_grid_attach(GTK_GRID(grid), drugaVrijednostRealniDioLabel, 1, redIndex, 1, 1);
    gtk_widget_hide(drugaVrijednostRealniDioLabel);

    drugaVrijednostImaginarniDioLabel = gtk_label_new("V imaginarni dio");
    gtk_grid_attach(GTK_GRID(grid), drugaVrijednostImaginarniDioLabel, 2, redIndex, 1, 1);
    gtk_widget_hide(drugaVrijednostImaginarniDioLabel);

    redIndex++;
    drugaVrijednostRealniDioEntryBox = gtk_entry_new();
    gtk_grid_attach(GTK_GRID(grid), drugaVrijednostRealniDioEntryBox, 1, redIndex, 1, 1);

    drugaVrijednostImaginarniDioEntryBox = gtk_entry_new();
    gtk_grid_attach(GTK_GRID(grid), drugaVrijednostImaginarniDioEntryBox, 2, redIndex, 1, 1);
    gtk_widget_hide(drugaVrijednostRealniDioEntryBox);
    gtk_widget_hide(drugaVrijednostImaginarniDioEntryBox);

    redIndex++;
    mojaLista = gtk_tree_view_new();
    init_list(mojaLista);
    gtk_grid_attach(GTK_GRID(grid), mojaLista, 0, redIndex, 4, 3);


    redIndex += 3;

    node1Label = gtk_label_new("Cvor 1");
    node2Label = gtk_label_new("Cvor 2");
    impedansaRealniDioLabel = gtk_label_new("Impedansa realni dio");
    impedansaImaginarniDioLabel = gtk_label_new("Impedansa imaginarni dio");

    gtk_grid_attach(GTK_GRID(grid), node1Label, 0, redIndex, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), node2Label, 1, redIndex, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), impedansaRealniDioLabel, 2, redIndex, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), impedansaImaginarniDioLabel, 3, redIndex, 1, 1);

    redIndex++;
    impedansaImaginarniDioEntryBox = gtk_entry_new();
    dodajAdmitansuButton = gtk_button_new_with_label("Dodaj admitansu");
    node1EntryBox = gtk_entry_new();
    node2EntryBox = gtk_entry_new();
    impedansaRealniDioEntryBox = gtk_entry_new();



    gtk_grid_attach(GTK_GRID(grid), node1EntryBox, 0, redIndex, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), node2EntryBox, 1, redIndex, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), impedansaRealniDioEntryBox, 2, redIndex, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), impedansaImaginarniDioEntryBox, 3, redIndex, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), dodajAdmitansuButton, 4, redIndex, 1, 1);

    redIndex++;
    admitansaLista = gtk_tree_view_new();
    init_list(admitansaLista);
    gtk_grid_attach(GTK_GRID(grid), admitansaLista, 0, redIndex, 4, 3);

    redIndex += 3;
    pocetneVrijednostiLabel = gtk_label_new("Pocetne vrijednosti: primjer: 1.2,2.3");
    gtk_grid_attach(GTK_GRID(grid), pocetneVrijednostiLabel, 0, redIndex, 1, 1);
    pocetneVrijednostiEntryBox = gtk_entry_new();
    gtk_grid_attach(GTK_GRID(grid), pocetneVrijednostiEntryBox, 1, redIndex, 1, 1);
    izracunajButton = gtk_button_new_with_label("Izracunaj");
    gtk_grid_attach(GTK_GRID(grid), izracunajButton, 2, redIndex, 1, 1);


    redIndex++;

    rezultatLabel = gtk_label_new("");
    gtk_grid_attach(GTK_GRID(grid), rezultatLabel, 0, redIndex, 1, 1);

    iteracijeGreskaLabel = gtk_label_new("");
    gtk_grid_attach(GTK_GRID(grid), iteracijeGreskaLabel, 1, redIndex, 1, 1);

    GtkWidget* window = gtk_application_window_new(app);
    g_signal_connect_swapped(GTK_BUTTON(dodajNodeButton), "clicked", G_CALLBACK(buttonClicked), window);
    g_signal_connect_swapped(GTK_BUTTON(dodajAdmitansuButton), "clicked", G_CALLBACK(admitansaClicked), window);
    g_signal_connect_swapped(GTK_BUTTON(izracunajButton), "clicked", G_CALLBACK(izracunajClicked), window);
    g_signal_connect(GTK_COMBO_BOX_TEXT(nodeTypeSelector), "changed", G_CALLBACK(nodeTypeChanged), window);
    gtk_window_set_child(GTK_WINDOW(window), grid);
    //gtk_window_set_child(GTK_WINDOW(window), GTK_WIDGET(store));
    gtk_window_present(GTK_WINDOW(window));

}

int main(int argc, char* argv[]) {
    // Create a new application
    GtkApplication* app = gtk_application_new("com.example.GtkApplication",
        G_APPLICATION_FLAGS_NONE);
    g_signal_connect(app, "activate", G_CALLBACK(on_activate), NULL);
    return g_application_run(G_APPLICATION(app), argc, argv);
}