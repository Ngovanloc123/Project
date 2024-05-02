#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stdbool.h>

#define MAX_SIZE 100
#define MAX_REPE 100
#define esp 0.001

struct Node {
    float data;
    struct Node *next;
    struct Node *prev;
};
typedef struct Node *List;
typedef struct Node *Position;

struct ArrayList {
    int max_size;
    List *PointToHeaderNode;
};
typedef struct ArrayList *Matrix;

Matrix createMatrix(int max_size) { // Tạo ma trận
    Matrix mt = malloc(sizeof(struct ArrayList));
    mt->max_size = max_size;
    mt->PointToHeaderNode = malloc((max_size + 1)*sizeof(struct Node));
    return mt;
}

List createNode(float data) { // Tạo Node
    List lt = malloc(sizeof(struct Node));
    lt->data = data;
    lt->next = NULL;
    lt->prev = NULL;
    return lt;
}

List createHeaderNode() { // Tạo HeaderNode
    List headerNode = malloc(sizeof(struct Node));
    headerNode->next = NULL;
    headerNode->prev = NULL;
    return headerNode;
}

List addNodetoList(int n) { // Thêm Node vào List
    List headerNode = createHeaderNode();
    Position p = headerNode;
    for(int i = 1; i <= n + 1; i++) {
        float x; scanf("%f", &x);
        List newNode = createNode(x);
        newNode->next = p->next;
        newNode->prev = p;
        p->next = newNode;
        if(newNode->next != NULL){
            newNode->next->prev = newNode;
        }
        p = p->next;
    }
    return headerNode;
}

void addListtoMatrix(Matrix A, int n) {// Thêm list vào Matrix
    for(int i = 1; i <= n; i++) {
        List headerNode = createHeaderNode();
        headerNode = addNodetoList(n);
        A->PointToHeaderNode[i] = headerNode;
    }
}

List createList(int n) {
    return addNodetoList(n - 1);
}

List createListZero(int n) {
    List headerNode = createHeaderNode();
    Position p = headerNode;
    for(int i = 1; i <= n; i++) {
        List newNode = createNode(0);
        newNode->next = p->next;
        newNode->prev = p;
        p->next = newNode;
        if(newNode->next != NULL) {
            newNode->next->prev = newNode;
        }
        p = p->next;
    }
    return headerNode;
}

void createMatrixZero(Matrix Temp, int n) {
    for(int i = 1; i <= n; i++) {
        List headerNode = createHeaderNode();
        headerNode = createListZero(n);
        Temp->PointToHeaderNode[i] = headerNode;
    }
}

void copyMatrix(Matrix Temp, Matrix A, int n) {
    for(int i = 1; i <= n; i++) {
        Position p = A->PointToHeaderNode[i]->next;
        Position q = Temp->PointToHeaderNode[i]->next;
        while(p != NULL) {
            q->data = p->data;
            p = p->next;
            q = q->next;
        }
    }
}

void copyList(List Temp, List B, int n) {
    Position p = B->next;
    Position q = Temp->next;
    while(p != NULL) {
        q->data = p->data;
        p = p->next;
        q = q->next;
    }
}

void deleteMiddle(Matrix A, int n) {
    for(int i = 1; i <= n; i++) {
        Position p = A->PointToHeaderNode[i];
        for(int j = 1; j <= (n + 1)/2 + 1; j++) {
            p = p->next;
        }
        p->prev->next = p->next;
        p->next->prev = p->prev;
        free(p);
    }
}

List getNodeFromMarix(Matrix A, int x, int y) {
    Position p = A->PointToHeaderNode[x];
    for(int i = 1; i <= y; i++){
        p = p->next;
    }
    return p;
}

List getNodeFromList(List headerNode, int x) {
    for(int i = 1; i <= x; i++){
        headerNode = headerNode->next;
    }
    return headerNode;
}

void displayMatrix(Matrix A, int n) {  
    for(int i = 1; i <= n; i++) {
        Position p = A->PointToHeaderNode[i]->next;
        while(p != NULL) {
            printf("%.3f\t", p->data);
            p = p->next;
        }
        printf("\n");
    }
}

void displayList(List headerNode) {
    headerNode = headerNode->next;
    while(headerNode != NULL) {
        printf("%.3f\t", headerNode->data);
        headerNode = headerNode->next;
    }
    printf("\n");
}

void displayN0(List N0) {// In ra nghiệm của phương trình
    N0 = N0->next;
    int i = 1;
    while(N0 != NULL) {
        printf("x%d = %.3f\t", i, N0->data);
        N0 = N0->next;
        i++;
    }
    printf("\n");
}

void Swap(float *a, float *b);

void Swap_row(Matrix A, List B, int n, int row1, int row2) {
    List t = A->PointToHeaderNode[row1];
    A->PointToHeaderNode[row1] = A->PointToHeaderNode[row2];
    A->PointToHeaderNode[row2] = t;

    Swap(&getNodeFromList(B, row1)->data, &getNodeFromList(B, row2)->data);
}

bool check_Matrix(Matrix A, List B, int n) {
    List idx = createListZero(n);
    List maxidx = createListZero(n);
    List pointList_idx = idx->next;
    List pointList_maxidx = maxidx->next;
    for(int i = 1; i <= n; i++) {
        pointList_maxidx->data = fabs(getNodeFromMarix(A, i, 1)->data);
        pointList_idx->data = 1;
        for(int j = 2; j <= n; j++) {
            if(pointList_maxidx->data < fabs(getNodeFromMarix(A, i, j)->data)) {
                pointList_maxidx->data = fabs(getNodeFromMarix(A, i, j)->data);
                pointList_idx->data = j;
            }
        }
        pointList_maxidx = pointList_maxidx->next;
        pointList_idx = pointList_idx->next;
    }
    pointList_idx = idx->next;
    for(int i = 1; i <= n; i++) {
        if(pointList_idx->data != i) {
            Swap_row(A, B, n, i, pointList_idx->data);
            if(pointList_idx->data == getNodeFromList(idx, pointList_idx->data)->data) return 0;
            Swap(&pointList_idx->data, &getNodeFromList(idx, pointList_idx->data)->data);
        }
        pointList_idx = pointList_idx->next;
    }
    return 1;
}

bool Gauss_Siedel(int n, Matrix A, List B, List N0) {
    List N1 = createListZero(n);
    bool check;
    int res = 0;
    do{
        check = false;
        res++;
        List pointList_N1 = N1->next;
        List pointList_B = B->next;
        List pointList_N0 = N0->next;
        for(int i = 1; i <= n; i++) {
            float s = 0;
            List pointListTemp_N0 = N0->next;
            for(int j = 1; j <= n; j++) {
                if(i != j){
                    s += getNodeFromMarix(A, i, j)->data * pointListTemp_N0->data;
                }
                pointListTemp_N0 = pointListTemp_N0->next;
            }
            pointList_N1->data = (pointList_B->data - s) / getNodeFromMarix(A, i, i)->data;
            if(fabs(pointList_N1->data - pointList_N0->data) >= esp) check = true;
            if(res == MAX_REPE) return 0;
            pointList_N1 = pointList_N1->next;
            pointList_B = pointList_B->next;
            pointList_N0 = pointList_N0->next;
        }
        for(int i = 1; i <= n; i++) getNodeFromList(N0, i)->data = getNodeFromList(N1, i)->data;
    }while (check);
    return true;
}

int Gauss(int n, Matrix A, List B, List N0) {
    for(int i = 1; i <= n - 1; i++){
        if(getNodeFromMarix(A, i, i)->data  == 0){
            int check = 0;
            for(int j = i + 1; j <= n; j++){
                if(getNodeFromMarix(A, j, i)->data != 0){
                    for(int k = 1; k <= n; k++){
                        Swap(&getNodeFromMarix(A, i, j)->data, &getNodeFromMarix(A, j, k)->data);
                    }
                    Swap(&getNodeFromList(B, i)->data, &getNodeFromList(B, j)->data);
                    
                    check++;
                    break;
                }
            }
            if(check == 0) return 0;
        }
        for(int j = i + 1; j <= n; j++){
            float h = -getNodeFromMarix(A, j, i)->data / getNodeFromMarix(A, i, i)->data;
            for(int k = i; k <= n; k++) getNodeFromMarix(A, j, k)->data += h * getNodeFromMarix(A, i, k)->data;
            getNodeFromList(B, j)->data += h * getNodeFromList(B, i)->data;
        }
    }
    
	for(int i = n; i > 0; i--){
		float s = 0;
        if(getNodeFromMarix(A, i, i)->data == 0) {
			if(getNodeFromList(B, i)->data == 0) return 1;
			else return 0;
		}
		for(int j = i; j <= n; j++) s += getNodeFromMarix(A, i, j)->data * getNodeFromList(N0, j)->data;
		getNodeFromList(N0, i)->data = (getNodeFromList(B, i)->data - s) / getNodeFromMarix(A, i, i)->data;
	}
	return 2;
}

float det(int n, Matrix A) {
    float res = 1;
    for(int i = 1; i <= n - 1; i++){
        if(getNodeFromMarix(A, i, i)->data  == 0){
            int check = 0;
            for(int j = i + 1; j <= n; j++){
                if(getNodeFromMarix(A, j, i)->data != 0){
                    for(int k = 1; k <= n; k++){
                        Swap(&getNodeFromMarix(A, i, j)->data, &getNodeFromMarix(A, j, k)->data);
                    }                   
                    check++;
                    break;
                }
            }
            if(check == 0) return 0;
        }
        for(int j = i + 1; j <= n; j++){
            float h = -getNodeFromMarix(A, j, i)->data / getNodeFromMarix(A, i, i)->data;
            for(int k = i; k <= n; k++) getNodeFromMarix(A, j, k)->data += h * getNodeFromMarix(A, i, k)->data;
        }
        res *= getNodeFromMarix(A, i, i)->data;
    }
    return res*getNodeFromMarix(A, n, n)->data;
}

void Krame(int n, Matrix A, Matrix Temp, List B, List N0) {
    copyMatrix(Temp, A, n);
    float d = det(n, Temp);
    if(d == 0) {
        printf("Phuong trinh khong the giai bang phuong phap Krame\n");
        return;
    }
    printf("Nghiem cua he phuong trinh:\n");
    List pointList_N0 = N0->next;
    for(int i = 1; i <= n; i++) {
        copyMatrix(Temp, A, n);
        List pointList_B = B->next;
        for(int j = 1; j <= n; j++) {
            getNodeFromMarix(Temp, j, i)->data = pointList_B->data;
            pointList_B = pointList_B->next;
        } 
        pointList_N0->data = det(n, Temp)/d;
        printf("x%d = %.3f\t", i, pointList_N0->data);
        pointList_N0 = pointList_N0->next;
    }
}

void change_Matrix(Matrix A, List B, int n) {
    List idx = createListZero(n);
    List maxidx = createListZero(n);
    List pointList_idx = idx->next;
    List pointList_maxidx = maxidx->next;
    for(int i = 1; i <= n; i++) {
        pointList_maxidx->data = fabs(getNodeFromMarix(A, i, 1)->data);
        pointList_idx->data = 1;
        for(int j = 2; j <= n; j++) {
            if(pointList_maxidx->data < fabs(getNodeFromMarix(A, i, j)->data)) {
                pointList_maxidx->data = fabs(getNodeFromMarix(A, i, j)->data);
                pointList_idx->data = j;
            }
        }
        pointList_maxidx = pointList_maxidx->next;
        pointList_idx = pointList_idx->next;
    }
    pointList_idx = idx->next;
    for(int i = 1; i <= n; i++) {
        if(pointList_idx->data != i) {
            Swap_row(A, B, n, i, pointList_idx->data);
            Swap(&pointList_idx->data, &getNodeFromList(idx, pointList_idx->data)->data);
        }
        pointList_idx = pointList_idx->next;
    }
}

bool Giam_du(int n, Matrix A, List B, List N0) {
    List pointList_B;
    List pointList_r;
    List pointList_N0;
    pointList_B = B->next;
    for(int i = 1; i <= n; i++) {
        float t = getNodeFromMarix(A, i, i)->data;
        if(t == 0) return 0;
        for(int j = 1; j <= n; j++) {
            getNodeFromMarix(A, i, j)->data /= t;
        }
        pointList_B->data /= t;
        pointList_B = pointList_B->next;
    }

    List r = createListZero(n);
    pointList_r = r->next;
    pointList_B = B->next;
    // Tinh so du ban dau
    for(int i = 1; i <= n; i++) {
        pointList_r->data = pointList_B->data;
        pointList_N0 = N0->next;
        for(int j = 1; j <= n; j++) {
            pointList_r->data -= getNodeFromMarix(A, i, j)->data*pointList_N0->data;
            pointList_N0 = pointList_N0->next;
        }
        pointList_B = pointList_B->next;
        pointList_r = pointList_r->next;
    }

    int res = 0;
    bool check;
    do {
        pointList_r = r->next;
        res++;
        check = false;// Dieu kien lap
        displayList(r);
        // Tim max {|r[i]|} 
        float max = fabs(getNodeFromList(r, 1)->data);
        pointList_r = pointList_r->next;
        int idx = 1;
        for(int i = 2; i <= n; i++) {
            if(max < fabs(pointList_r->data)) {
                max = fabs(pointList_r->data);
                idx = i;
            }
            pointList_r = pointList_r->next;
        }
        getNodeFromList(N0, idx)->data += getNodeFromList(r, idx)->data;
        // Tinh lai r[i] kiem tra kha nang lap tiep theo
        pointList_r = r->next;
        float d = getNodeFromList(r, idx)->data;
        for(int i = 1; i <= n; i++) {
            pointList_r->data -= getNodeFromMarix(A, i, idx)->data*d;
            if(fabs(pointList_r->data) >= esp) check = true;
            pointList_r = pointList_r->next;
        }
        if(res == MAX_REPE) return 0;
    }  while(check);
    return 1;
}

void addListtoMatrixfromfile(Matrix A, int n, char *filename) {
    FILE *file;
    file = fopen(filename, "r");
    if(file == NULL) return;
    fscanf(file, "%d", &n);
    for(int i = 1; i <= n; i++) {
        List headerNode = createHeaderNode();
        Position p = headerNode;
        for(int j = 1; j <= n + 1; j++) {
            float x;
            fscanf(file, "%f", &x);
            List newNode = createNode(x);
            newNode->next = p->next;
            newNode->prev = p;
            p->next = newNode;
            if(newNode->next != NULL){
                newNode->next->prev = newNode;
            }
            p = p->next;
        }
        A->PointToHeaderNode[i] = headerNode;
    }
    fclose(file);
}

bool is_square_matrix(Matrix A, int n) {
    Position p = A->PointToHeaderNode[1]->next;
    int res = 0;
    while(p != NULL) {
        res++;
        p = p->next;
    }
    return res == n;
}

int main() {
    Matrix A = createMatrix(MAX_SIZE);
    List B;
    bool have_matrix_B = false;
    int n;// So an cua phuong trinh

    printf("----------------------------------------------------------\n");
    printf("----------------GIAI HE PHUONG TRINH Ax = B---------------\n");
    printf("----------------------------------------------------------\n\n");
    Repeat:// Nhap lai ma tran A
    have_matrix_B = false;
    printf("---------------------------MENU---------------------------\n");
    printf("1. Nhap ma tran A thu cong\n");
    printf("2. Nhap ma tran A bang file\n");
    printf("----------------------------------------------------------\n");
    printf("Nhap lua chon cua ban: ");

    int select; scanf("%d", &select);
    if(select == 1) {
        printf("Nhap so an cua cua phuong trinh: "); scanf("%d", &n);
        printf("Nhap Ma Tran %dx%d:\n", n, n + 1);
        addListtoMatrix(A, n);
        printf("Ma tran A ban nhap da duoc luu!!\n");
    }
    else if(select == 2) {
        char *filename;
        printf("-----------------FILE-----------------\n");
        printf("1. Matrix1.inp\n");
        printf("2. Matrix2.inp\n");
        printf("--------------------------------------\n");
        printf("Chon file ban muon: ");
        scanf("%d", &select);
        if(select == 1) {
            filename = "Matrix1.inp";
        } else {
            filename = "Matrix2.inp";
        }
        FILE *file;
        file = fopen(filename, "r");
        if(file != NULL) {
            fscanf(file, "%d", &n);

            addListtoMatrixfromfile(A, n, filename);
            printf("Du lieu file %s:\n", filename);
            displayMatrix(A, n);
            printf("Ma tran da duoc luu vao A\n");

        }
    }
    
    while(1) {
        printf("\n--------------------------MENU---------------------------\n");
        printf("1. Nhap lai ma tran A\n");
        printf("2. In ra ma tran A\n");
        printf("3. Xoa 1 cot o vi tri giua ma tran A\n");
        printf("4. Nhap/nhap lai (neu ban da co ma tran B) ma tran B \n");
        if(is_square_matrix(A, n) && have_matrix_B) {
            printf("5. Giai he phuong trinh Ax = B\n");
        } else {
            printf("Luu y:\n");
            printf("Ban chua the giai phuong trinh Ax = B !!\n");
            if(!is_square_matrix(A, n)) printf("!! Ma tran A chua phai ma tran vuong (nhan phim 3 de xoa cot)\n");
            if(!have_matrix_B) printf("!! Ban chua nhap ma tran B (nhan phim 4 de nhap ma tran B)\n");
        }
        printf("6. Thoat chuong trinh!!\n");
        printf("----------------------------------------------------------\n");

        printf("Nhap lua chon cua ban: ");
        scanf("%d", &select);
        if(select == 1) {
            goto Repeat;
        } else if(select == 2) {
            displayMatrix(A, n);
        } else if(select == 3) {
            deleteMiddle(A, n);
            printf("Ma tran A sau khi xoa cot o vi tri gua ma tran!!\n");
            displayMatrix(A, n);
        } else if(select == 4) {
            have_matrix_B = true;
            printf("Nhap ma tran B 1x%d:\n", n);
            B = createList(n); 
        } else if(select == 5) {
            Matrix TempMatrix = createMatrix(MAX_SIZE);// Matrix temporary
            createMatrixZero(TempMatrix, n); // Tạo ma trận chứa các phần tử bằng không
            List TempList = createListZero(n);// // Tạo list tepmorary chứa các phần tử bằng không
 
            printf("\n---------------------------MENU---------------------------\n");
            printf("Cac phuong phap giai he phuong trinh Ax = B\n");
            printf("1. Phuong phap Krame\n");
            printf("2. Phuong phap Gauss\n");
            printf("3. Phuong phap Gauss_Siedel\n");
            printf("4. Phuong phap Giam du\n");
            printf("----------------------------------------------------------\n");
            printf("Nhap lua chon cua ban: ");
            scanf("%d", &select);
            if(select == 1) {
                List N0 = createListZero(n);
                Krame(n, A, TempMatrix, B, N0);
            } else if(select == 2) {
                List N0 = createListZero(n);
                copyMatrix(TempMatrix, A, n);
                copyList(TempList, B, n);
                int check = Gauss(n, TempMatrix, TempList, N0);
                if(check == 2) {
                    displayN0(N0);
                } else if(check == 0) {
                    printf("He phuong trinh vo nghiem!!\n");
                } else {
                    printf("He phuong trinh vo so nghiem!!\n");
                }
            } else if(select == 3) {
                copyMatrix(TempMatrix, A, n);// Sao chép ma trận A cho ma trận TempMatrix
                copyList(TempList, B, n);// Sao chép ma trận B cho ma trận TempList
                printf("Nhap nghiem ban dau:\n");
                List N0 = createList(n);
                bool check = check_Matrix(TempMatrix, TempList, n);
                if(check) {
                    if(Gauss_Siedel(n, TempMatrix, TempList, N0)) {
                        printf("Nghiem cua he phuong trinh:\n");
                        displayN0(N0);
                    } else printf("1He phuong trinh khong the giai bang phuong phap Gauss_Siedel\n");
                } else printf("2He phuong trinh khong the giai bang phuong phap Gauss_Siedel\n");
            } else if(select == 4) {
                printf("Nhap nghiem ban dau:\n");
                List N0 = createList(n);
                copyMatrix(TempMatrix, A, n);
                copyList(TempList, B, n);
                change_Matrix(TempMatrix, TempList, n);
                if(Giam_du(n, TempMatrix, TempList, N0)) {
                    printf("Nghiem cua he phuong trinh:\n");
                    displayN0(N0);
                }
            }
        } else if(select == 6) {
            printf("Da thoat chuong trinh!!\n");
            break;
        } else {
            printf("Lua chon cua ban khong hop le!!\n");
        }
    }    
    return 0;
}

void Swap(float *a, float *b) {
    float t = *a;
    *a = *b;
    *b = t;
}
