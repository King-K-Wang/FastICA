#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

/* 信号数量 */
#define M 4
/* 信号长度 */
#define T 1024
/* 最大循环长度 */
#define MAX_CYCLE 1000000
/* 精度 */
#define ACCURATE 0.00001

int main(int argc, char *argv[]) {
    if (argc != 2 * M + 1) {
        printf("format error, format: %s \"input files\" \"output files\"\n", argv[0]);
        return 1;
    }

    /* 输入信号 */
    double X[M][T]; 
    for (uint16_t i = 0; i < M; i++) {
        FILE *sensor = fopen(argv[i + 1], "r");
        if (sensor == NULL) {
            printf("Failed to open file: %s\n", argv[i + 1]);
            return 1;
        }
        for (uint32_t j = 0; j < T; j++) {
            fscanf(sensor, "%lf", &X[i][j]);
        }
        fclose(sensor);
    }

    /* 计算时间 */
    clock_t tstart, tstop;
    tstart = clock();
    /* 去均值 */
    double average[M];
    for (uint16_t i = 0; i < M; i++) {
        double sum = 0;
        for (uint32_t j = 0; j < T; j++) {
            sum += X[i][j];
        }
        average[i] = sum / T;
        for (uint32_t j = 0; j < T; j++) {
            X[i][j] -= average[i];
        }
    }

    /* 计算协方差矩阵 */
    double covariance[M][M] = {0};
    for (uint16_t i = 0; i < M; i++) {
        for (uint16_t j = 0; j < M; j++) {
            for (uint32_t k = 0; k < T; k++) {
                covariance[i][j] += X[i][k] * X[j][k];
            }
            covariance[i][j] /= T;
        }
    }
    printf("Covariance matrix = \n");
    for (uint16_t i = 0; i < M; i++) {
        for (uint16_t j = 0; j < M; j++) {
            printf("%20.5lf,", covariance[i][j]);
        }
        printf("\n");
    }

    /* 初始化特征值 */
    double eigenvalue[M][M];
    double eigenvector[M][M];
    for (uint16_t i = 0; i < M; i++) {
        for (uint16_t j = 0; j < M; j++) {
            eigenvalue[i][j] = covariance[i][j];
            eigenvector[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }

    for (uint32_t i = 0; i < MAX_CYCLE; i++) {
        /* 查找最大非对角线元素 */
        double max_number = 0.0;
        int max_row = 0, max_column = 0;
        for (uint16_t j = 1; j < M; j++) {
            for (uint16_t k = 0; k < j; k++) {
                if ((fabs(eigenvalue[j][k])) > max_number) {
                    max_number = fabs(eigenvalue[j][k]);
                    max_row = j;
                    max_column = k;
                }
            }
        }
        if(max_number < ACCURATE && max_number > -ACCURATE) {
            // printf("Jacobi uses cycles = %d\n", i);
            break;
        }

        /* 旋转角计算 */
        double theta = atan(2.0 * max_number / (eigenvalue[max_row][max_row] - eigenvalue[max_column][max_column])) / 2.0;

        /* Givens矩阵 */
        double givens[M][M];
        for (uint16_t j = 0; j < M; j++) {
            for (uint16_t k = 0; k < M; k++) {
                givens[j][k] = (j == k) ? 1.0 : 0.0;
            }
        }
        givens[max_row][max_row] = cos(theta);
        givens[max_row][max_column] = -sin(theta);
        givens[max_column][max_row] = sin(theta);
        givens[max_column][max_column] = cos(theta);

        /* 旋转变换 */
        double temp[M][M];
        for (uint16_t j = 0; j < M; j++) {
            for (uint16_t k = 0; k < M; k++) {
                temp[j][k] = 0.0;
                for (uint16_t l = 0; l < M; l++) {
                    temp[j][k] += givens[j][l] * eigenvalue[l][k];
                }
            }
        }
        givens[max_row][max_column] = -givens[max_row][max_column];
        givens[max_column][max_row] = -givens[max_column][max_row];
        for (uint16_t j = 0; j < M; j++) {
            for (uint16_t k = 0; k < M; k++) {
                eigenvalue[j][k] = 0.0;
                for (uint16_t l = 0; l < M; l++) {
                    eigenvalue[j][k] += temp[j][l] * givens[l][k];
                }
            }
        }

        /* 更新特征矩阵 */
        for (uint16_t j = 0; j < M; j++) {
            for (uint16_t k = 0; k < M; k++) {
                temp[j][k] = eigenvector[j][k];
            }
        }
        for (uint16_t j = 0; j < M; j++) {
            for (uint16_t k = 0; k < M; k++) {
                eigenvector[j][k] = 0.0;
                for (uint16_t l = 0; l < M; l++) {
                    eigenvector[j][k] += temp[j][l] * givens[l][k];
                }
            }
        }
    }

    /* 特征值和特征矩阵 */
    printf("EigenValue = \n");
    for (uint16_t i = 0; i < M; i++) {
        for (uint16_t j = 0; j < M; j++) {
            printf("%20.5lf,", eigenvalue[i][j]);
        }
        printf("\n");
    }
    printf("EigenVector = \n");
    for (uint16_t i = 0; i < M; i++) {
        for (uint16_t j = 0; j < M; j++) {
            printf("%20.5lf,", eigenvector[i][j]);
        }
        printf("\n");
    }
    
    /* 白化 */
    double sqrt_eigenvalue[M][M];
    for (uint16_t i = 0; i < M; i++) {
        for (uint16_t j = 0; j < M; j++) {
            sqrt_eigenvalue[i][j] = (i == j) ? pow(eigenvalue[i][j], -0.5) : 0.0;
        }
    }
    double whiten_matrix[M][M];
    for (uint16_t i = 0; i < M; i++) {
        for (uint16_t j = 0; j < M; j++) {
            whiten_matrix[i][j] = 0;
            for (uint16_t k = 0; k < M; k++) {
                whiten_matrix[i][j] += eigenvector[j][k] * sqrt_eigenvalue[i][k];
            }
        }
    }
    double Z[M][T];
    for (uint16_t i = 0; i < M; i++) {
        for (uint32_t j = 0; j < T; j++) {
            Z[i][j] = 0.0;
            for (int k = 0; k < M; k++) {
                Z[i][j] += whiten_matrix[i][k] * X[k][j];
            }
        }
    }

    /* 验证白化后信号相关性 */
    for (uint16_t i = 0; i < M; i++) {
        for (uint16_t j = 0; j < M; j++) {
            covariance[i][j] = 0;
            for (uint32_t k = 0; k < T; k++) {
                covariance[i][j] += Z[i][k] * Z[j][k];
            }
            covariance[i][j] /= T;
        }
    }
    printf("Whiten Verification = \n");
    for (uint16_t i = 0; i < M; i++) {
        for (uint16_t j = 0; j < M; j++) {
            printf("%20.5lf,", covariance[i][j]);
        }
        printf("\n");
    }

    /* 恢复矩阵 */
    double W[M][M];
    for (uint16_t i = 0; i < M; i++) {
        double w[M] = {0};
        double w_new[M];
        for (uint16_t j = 0; j < M; j++) {
            w_new[j] = 0.5;
        }
        uint32_t cycle = 0;
        do {
            /* 更新w */
            for (uint16_t j = 0; j < M; j++) {
                w[j] = w_new[j];
            }
            for (uint16_t j = 0; j < M; j++) {
                double mean1 = 0.0;
                double mean2 = 0.0;
                for (int k = 0; k < T; k++) {
                    mean1 += Z[j][k] * tanh(w[0] * Z[0][k] + w[1] * Z[1][k] + w[2] * Z[2][k] + w[3] * Z[3][k]);
                    mean2 += (1 - pow(tanh(w[0]) * Z[0][k] + tanh(w[1]) * Z[1][k] + tanh(w[2]) * Z[2][k] + tanh(w[3]) * Z[3][k], 2));
                }
                w_new[j] = mean1 / T - mean2 / T * w[j];
            }

            /* 施密特正交化 */
            for (uint16_t j = 0; j < i; j++) {
                double dot_product = 0.0;
                for (uint16_t k = 0; k < M; k++) {
                    dot_product += w_new[k] * W[k][j];
                }
                for (uint16_t k = 0; k < M; k++) {
                    w_new[k] -= dot_product * W[k][j];
                }
            }
            double norm = 0.0;
            for (uint16_t j = 0; j < M; j++) {
                norm += w_new[j] * w_new[j];
            }
            norm = sqrt(norm);
            for (uint16_t j = 0; j < M; j++) {
                w_new[j] /= norm;
            }
            
            cycle++;
        } while (cycle < MAX_CYCLE && fabs(w_new[0] - w[0]) > ACCURATE && fabs(w_new[1] - w[1]) > ACCURATE && fabs(w_new[2] - w[2]) > ACCURATE && fabs(w_new[3] - w[3]) > ACCURATE);
        
        // printf("Column %d uses cycles %d\n", i, cycle);
        for (uint16_t j = 0; j < M; j++) {
            W[j][i] = w_new[j];
        }
    }

    /* 恢复矩阵 */
    printf("W = \n");
    for (uint16_t i = 0; i < M; i++) {
        for (uint16_t j = 0; j < M; j++) {
            printf("%20.5f,", W[i][j]);
        }
        printf("\n");
    }
    
    /* 恢复信号 */
    double Y[M][T];
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < T; j++) {
            Y[i][j] = 0;
            for (int k = 0; k < M; k++) {
                Y[i][j] += W[k][i] * Z[k][j];
            }
        }
    }

    tstop = clock();

    for (uint16_t i = 0; i < M; i++) {
        FILE *recovery = fopen(argv[i + M + 1], "wb");
        if (recovery == NULL) {
            printf("Failed to open file: %s\n", argv[i + 1]);
            return 1;
        }
        for (uint32_t j = 0; j < T; j++) {
            fprintf(recovery, "%.5lf,\n", Y[i][j]);
        }
        fclose(recovery);
    }
    printf("Done!\n");
    printf("Time usage = %lf\n", (double)(tstop - tstart) / CLOCKS_PER_SEC);
    return 0;
}
