#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define OK 0
#define ERR -1

typedef struct point
{
    double x;
    double y;
} point_t;

double func(double x, double y)
{
    return x * x + y * y;
}

void matrix_free(double **matrix, int n)
{
    for (int i = 0; i < n; i++)
    {
        free(matrix[i]);
    }

    free(matrix);
}

int matrix_allocate(double ***matrix, int n)
{
    *matrix = calloc(n, sizeof(double *));
    if (!*matrix)
    {
        return ERR;
    }

    for (int i = 0; i < n; i++)
    {
        (*matrix)[i] = malloc(n * sizeof(double));
        if (!(*matrix)[i])
        {
            matrix_free(*matrix, i);
            return ERR;
        }
    }

    return OK;
}

// Initial values table
int initial_table_create(FILE *f, point_t **table, double ***z, int n)
{
    *table = malloc(n * sizeof(point_t));
    if (!*table)
    {
        return ERR;
    }

    if (matrix_allocate(z, n) == ERR)
    {
        free(*table);
        return ERR;
    }

    for (int i = 0; i < n; i++)
    {
        //printf("Input %d point data (x y): ", i + 1);
        fscanf(f, "%lf", &((*table)[i].x));
    }
    fscanf(f, "\n");

    for (int i = 0; i < n; i++)
    {
        //printf("Input %d point data (x y): ", i + 1);
        fscanf(f, "%lf",  &((*table)[i].y));

        for (int j = 0; j < n; j++)
        {
            fscanf(f, "%lf", &((*z)[i][j]));
        }
        fscanf(f, "\n");
    }

    return OK;
}

void initial_table_print(FILE *f, point_t *table, double **z, int n)
{
   fprintf(stdout, "%9s", "y\\x");

    for (int i = 0; i < n; i++)
    {
        //printf("Input %d point data (x y): ", i + 1);
        fprintf(f, "%9.2lf", table[i].x);
    }
    fprintf(f, "\n");

    for (int i = 0; i < n; i++)
    {
        //printf("Input %d point data (x y): ", i + 1);
        fprintf(f, "%9.2lf",  table[i].y);

        for (int j = 0; j < n; j++)
        {
            fprintf(f, "%9.2lf", z[i][j]);
        }
        fprintf(f, "\n");
    }
}

int min_x_i(point_t *points, int n_init, double x)
{
    int min_i = 0, min;
    int flag = 0;

    for (int i = 0; i < n_init; i++)
    {
        if (!flag)
        {
            min = fabs(points[i].x - x);
            flag = 1;
        }
        else if (fabs(points[i].x - x) < min)
        {
            min = fabs(points[i].x - x);
            min_i = i;
        }
    }

    return min_i;
}

int min_point_i(point_t *points, int n_init, double x)
{
    int min_i = 0, min;
    int flag = 0;

    for (int i = 0; i < n_init; i++)
    {
        if (!flag)
        {
            min = abs(points[i].x - x);
            flag = 1;
        }
        else if (abs(points[i].x - x) < min)
        {
            min = abs(points[i].x - x);
            min_i = i;
        }
    }

    return min_i;
}

// Points selection
point_t *choose_points(point_t *points, double *z, int n_init, double x, int n, int flag)
{
    point_t *selected_points = malloc(n * sizeof(point_t));
    if (!selected_points)
    {
        return NULL;
    }

    int i_beg, i_flag = 0;
    int n_required = ceil(n / 2);
    int i_near = min_point_i(points, n_init, x);

    if (i_near + n_required + 1 > n_init)
    {
        i_beg = n_init - n;
    }
    else if (i_near < n_required)
    {
        i_beg = 0;
    }
    else
    {
        i_beg = i_near - n_required + 1;
    }

    if (flag)
    {
        i_flag = i_beg;
    }
    for (int i = 0; i < n; i++)
    {
        selected_points[i].x = points[i_beg + i].x;
        selected_points[i].y = z[i_beg - i_flag + i];
    }

    printf("\nSELECTED:\n");
    for (int i = 0; i < n; i++)
    {
        printf("%lf %lf\n", selected_points[i].x, selected_points[i].y);
    }

    return selected_points;
}

void arr_print(double *arr, char *s, int n)
{
    printf("\n%s\n", s);

    for (int i = 0; i < n; i++)
    {
        printf(" %lf ", arr[i]);
    }
}

// Divided difference matrix
void dif_matrix_free(double **dif_matrix, int n)
{
    for (int i = 0; i < n; i++)
    {
        free(dif_matrix[i]);
    }

    free(dif_matrix);
}

int dif_matrix_allocate(double ***dif_matrix, int n)
{
    *dif_matrix = calloc(n, sizeof(double *));
    if (!*dif_matrix)
    {
        return ERR;
    }

    for (int i = 0; i < n; i++)
    {
        (*dif_matrix)[i] = malloc((n - i) * sizeof(double));
        if (!(*dif_matrix)[i])
        {
            dif_matrix_free(*dif_matrix, i);
            return ERR;
        }
    }

    return OK;
}

void arr_null(double *arr, int n)
{
    for (int  i = 0; i < n; i++)
    {
        arr[i] = 0;
    }
}

void arr_copy(double *destination, double *source, int n)
{
    for (int i = 0; i < n; i++)
    {
        destination[i] = source[i];
    }
}

int dif_matrix_fill(double **dif_matrix, point_t *points, int n)
{
    double *tmp_arr = malloc(n *sizeof(double));
    if (!tmp_arr)
    {
        return ERR;
    }

    int cur_count;

    for (int i = 0; i < n; i++)
    {
        arr_null(tmp_arr, n);

        for (int j = 0; j < n - i; j++)
        {
            if (i == 0)
            {
                tmp_arr[j] = (points[j].y - points[j + 1].y) / (points[j].x - points[j + i + 1].x);
                
            }
            else
            {

                tmp_arr[j] = (dif_matrix[i - 1][j] - dif_matrix[i - 1][j + 1]) / (points[j].x - points[j + i + 1].x);

            }  
        }
        cur_count = n - i - 1;

        arr_copy(dif_matrix[i], tmp_arr, cur_count);
        arr_print(dif_matrix[i], "DIF", cur_count);
    }

    free(tmp_arr);

    return OK;
}

int newton_interpolation(point_t *points, double *z, int n_init, \
                        double x, int n, double *result, int flag)
{
    point_t *selected_points = calloc(n + 1, sizeof(point_t));
    selected_points = choose_points(points, z, n_init, x, n + 1, flag);
    if (!selected_points)
    {
        return ERR;
    }

    double **dif_matrix;
    dif_matrix_allocate(&dif_matrix, n);
    dif_matrix_fill(dif_matrix, selected_points, n + 1);

    double cur_k = 1;

    for (int i = 0; i < n + 1; i++)
    {
        if (i == 0)
        {
            *result += selected_points[i].y;
        }
        else
        {
            *result += cur_k * dif_matrix[i - 1][0];
        }

        cur_k *= x - selected_points[i].x;

        //arr_print(dif_matrix[i], "dif_matrix cur_i", n - i);
    }

    free(selected_points);
    dif_matrix_free(dif_matrix, n);

    return OK;
}

void change_y_x(point_t *points, int n_init)
{
    double tmp;

    for (int i = 0; i < n_init; i++)
    {
        tmp = points[i].x;
        points[i].x = points[i].y;
        points[i].y = tmp;
    }
}

int interpolation(point_t *points, double **z, int n, \
                double x, double y, int n_x, int n_y, double *result)
{
    double *cur_res = calloc((n_y + 1), sizeof(double));

    for (int i = 0; i < n_y + 1; i++)
    {
        if (newton_interpolation(points, z[i], n, x, n_y, &(cur_res[i]), 0) != OK)
        {
            free(points);
            return ERR;
        }
        //printf("\nCUR RES[%d]: %lf\n", i + 1, cur_res[i]);
    }

    change_y_x(points, n);
    if (newton_interpolation(points, cur_res, n, y, n_x, result, 1) != OK)
    {
        free(points);
        return ERR;
    }
    //printf("\nCUR RES[end]: %lf\n", *result);

    free(cur_res);

    return OK;
}

int cmp_table_print(point_t *initial_table, double **z, int n_init, \
                    double x, double y, int n_max)
{
    double cur_res = 0.0;

    printf("\n|%10s |%19s |\n", "n_x = n_y", "Intorpolated value");

    for (int i = 1; i < n_max + 1; i++)
    {
        cur_res = 0.0;

        if (interpolation(initial_table, z, n_init, x, y, i, i, &cur_res) != OK)
        {
            free(initial_table);
            return ERR;
        }

        printf("|%10d |", i);
        printf("%20lf|\n", cur_res);
    }

    return OK;
}

int main(int argc, char **argv)
{
    int n_init = 0;
    setbuf(stdout, NULL);

    if (argc != 2)
    {
        return ERR;
    }

    printf("\nINITIAL TABLE DATA\n");

    FILE *f = fopen(argv[1], "r");
    if (!f)
    {
        return ERR;
    }

    if (fscanf(f, "%d", &n_init) != 1 || n_init < 1)
    {
        return ERR;
    }
    fscanf(f, "\n");

    point_t *initial_table;
    double **z = NULL;

    if (initial_table_create(f, &initial_table, &z, n_init) != OK)
    {
        return ERR;
    }

    fclose(f);

    printf("\n");
    initial_table_print(stdout, initial_table, z, n_init);

    int n_x, n_y;
    double x, y;

    printf("\nInput n_x n_y: ");
    if (scanf("%d%d", &n_x, &n_y) != 2 || n_x < 0 || n_y < 0)
    {
        free(initial_table);
        matrix_free(z, n_init);
        return ERR;
    }

    printf("Input x y: ");
    if (scanf("%lf%lf", &x, &y) != 2)
    {
        free(initial_table);
        matrix_free(z, n_init);
        return ERR;
    }

    cmp_table_print(initial_table, z, n_init, x, y, n_y);

    printf("\nz(x,y)      : %lf", func(x, y));

    matrix_free(z, n_init);
    free(initial_table);

    return OK;
}