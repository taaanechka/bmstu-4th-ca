#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.1415
#define OK 0
#define ERR -1

typedef struct point
{
    double x;
    double y;
    double y_der;
} point_t;

double f(double x)
{
    return cos(x) - x;
}

/* double f_to_rad(double (*f)(double), double x)
{
    //return f(x * PI / 180);
    return f(x);
} */

// Initial values table
int initial_table_create(point_t **table, int n)
{
    *table = malloc(n * sizeof(point_t));
    if (!*table)
    {
        return ERR;
    }

    for (int i = 0; i < n; i++)
    {
        printf("Input %d point data (x y y'): ", i + 1);
        scanf("%lf%lf%lf", &((*table)[i].x), &((*table)[i].y), &((*table)[i].y_der));
    }

    return OK;
}

void initial_table_print(point_t *table, int n)
{
    printf("\n%10s | %10s | %10s\n", "X", "Y", "Y'");

    for (int i = 0; i < n; i++)
    {
        printf("%10.6lf | %10.6lf | %10.6lf\n", table[i].x, table[i].y, table[i].y_der);
    }
}

void points_sort_x(point_t *table, int n)
{
    int idx = 0;
    double val = 0.0;
    point_t cur_dot = {0, 0 , 0};

    for (int i = 1; i < n; i++)
    {
        val = table[i].x;
        idx = i - 1;
        cur_dot = table[i];

        while (idx >= 0 && table[idx].x > val)
        {
            table[idx + 1] = table[idx];
            idx--;
        }

        table[idx + 1] = cur_dot;
    }
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
point_t *choose_points(point_t *points, int n_init, double x, int n)
{
    point_t *selected_points = malloc(n * sizeof(point_t));
    if (!selected_points)
    {
        return NULL;
    }

    int i_beg;
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

    for (int i = 0; i < n; i++)
    {
        selected_points[i] = points[i_beg + i];
    }

    return selected_points;
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

int dif_matrix_fill(double **dif_matrix, point_t *points, int n, int hermit_flag)
{
    double *tmp_arr = malloc(n *sizeof(double));
    if (!tmp_arr)
    {
        return ERR;
    }

    int cur_j = 0, cur_count;

    for (int i = 0; i < n; i++)
    {
        arr_null(tmp_arr, n);

        for (int j = 0; j < n - i; j++)
        {
            if (i == 0)
            {
                if (!hermit_flag)
                    tmp_arr[j] = (points[j].y - points[j + 1].y) / (points[j].x - points[j + i + 1].x);
                else
                {
                    cur_j = j / 2;

                    if (j % 2 == 0)
                    {
                        tmp_arr[j] = points[cur_j].y_der;
                    }
                    else
                    {
                        tmp_arr[j] = (points[cur_j].y - points[cur_j + 1].y) / (points[cur_j].x - points[cur_j + i + 1].x);
                    } 
                }
                
            }
            else
            {
                if (!hermit_flag)
                {
                    tmp_arr[j] = (dif_matrix[i - 1][j] - dif_matrix[i - 1][j + 1]) / (points[j].x - points[j + i + 1].x);
                }
                else
                {
                    tmp_arr[j] = (dif_matrix[i - 1][j] - dif_matrix[i - 1][j + 1]) / (points[j / 2].x - points[(j + i - 1) / 2 + 1].x);
                }
            }  
        }

        if (!hermit_flag)
            cur_count = n - i - 1;
        else
            cur_count = n - i;
        
        arr_copy(dif_matrix[i], tmp_arr, cur_count);
    }

    free(tmp_arr);

    return OK;
}

int newton_interpolation(point_t *points, int n_init, double x, int n, double *result)
{
    point_t *selected_points = calloc(n + 1, sizeof(point_t));
    selected_points = choose_points(points, n_init, x, n + 1);
    if (!selected_points)
    {
        return ERR;
    }

    double **dif_matrix;
    dif_matrix_allocate(&dif_matrix, n);
    dif_matrix_fill(dif_matrix, selected_points, n + 1, 0);

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
    }

    free(selected_points);
    dif_matrix_free(dif_matrix, n);

    return OK;
}

int hermit_interpolation(point_t *points, int n_init, double x, int n, double *result)
{
    point_t *selected_points = calloc(n + 1, sizeof(point_t));
    selected_points = choose_points(points, n_init, x, n + 1);
    if (!selected_points)
    {
        return ERR;
    }

    double **dif_matrix;
    dif_matrix_allocate(&dif_matrix, 2 * n - 1);
    dif_matrix_fill(dif_matrix, selected_points, 2 * n - 1, 1);

    double cur_k = 1;

    for (int i = 0; i < 2 * n - 1; i++)
    {
        if (i == 0)
        {
            *result += selected_points[i].y;
        }
        else
        {
            *result += cur_k * dif_matrix[i - 1][0];
        }

        cur_k *= x - selected_points[i / 2].x;
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

int cmp_table_print(point_t *initial_table, int n_init, double x, int n_max)
{
    double cur_res_1 = 0.0, cur_res_2 = 0.0;
    int n_cur = 1;

    if (n_max == 0)
    {
        n_cur = 0;
        n_max = 1;
    }

    printf("\n|%3s |%27s |%27s |\n", "N", "Newton intorpolated value", "Hermit intorpolated value");

    for (int i = 0; i < n_max; i++)
    {
        cur_res_1 = 0.0;
        cur_res_2 = 0.0;

        if (n_cur != 0)
        {
            n_cur = i + 1;
        }

        if (newton_interpolation(initial_table, n_init, x, n_cur, &cur_res_1) != OK)
        {
            free(initial_table);
            return ERR;
        }

        if (hermit_interpolation(initial_table, n_init, x, n_cur, &cur_res_2) != OK)
        {
            free(initial_table);
            return ERR;
        }

        printf("|%3d |", i +1);
        printf("%27lf|", cur_res_1);
        printf("%27lf|\n", cur_res_2);
    }

    return OK;
}

int main(void)
{
    int n_init = 0;
    setbuf(stdout, NULL);

    printf("\nINITIAL TABLE DATA\n");

    printf("Input dots amount: ");
    if (scanf("%d", &n_init) != 1 || n_init < 1)
    {
        return ERR;
    }

    point_t *initial_table;
    if (initial_table_create(&initial_table, n_init) != OK)
    {
        return ERR;
    }

    //initial_table_print(initial_table, n_init);
    points_sort_x(initial_table, n_init);
    //initial_table_print(initial_table, n_init);

    // Task conditions
    int n;
    double x;

    printf("\nInput n: ");
    if (scanf("%d", &n) != 1 || n < 0)
    {
        free(initial_table);
        return ERR;
    }

    printf("Input x: ");
    if (scanf("%lf", &x) != 1)
    {
        free(initial_table);
        return ERR;
    }

    // Task solution
    double root = 0;

    cmp_table_print(initial_table, n_init, x, n);

    change_y_x(initial_table, n_init);
    if (newton_interpolation(initial_table, n_init, 0, n, &root) != OK)
    {
        free(initial_table);
        return ERR;
    }

    printf("\ny(x)              : %lf", f(x));
    printf("\nThis function root: %lf", root);

    free(initial_table);

    return OK;
}