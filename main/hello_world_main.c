#include <stdio.h>  
#include "esp_log.h"
#include "driver/i2c.h"
#include "sdkconfig.h"
#include "freertos/FreeRTOS.h"
#include "freertos/task.h"
#include <math.h>
#include "freertos/semphr.h"  // Pre mutex
#include "fft.h"
#include "butterworth.h"
#include <complex.h>
// I2C adresa pre MPU6050 (akcelerometer a gyroskop)
#define MPU6050_ADDR 0x68

// I2C nastavenia
#define I2C_MASTER_SCL_IO 7              /*!< GPIO number for I2C master clock */
#define I2C_MASTER_SDA_IO 6              /*!< GPIO number for I2C master data  */
#define I2C_MASTER_NUM I2C_NUM_0         /*!< I2C port number for master dev */
#define I2C_MASTER_FREQ_HZ 100000        /*!< I2C master clock frequency */
#define I2C_MASTER_TX_BUF_DISABLE 0      /*!< I2C master doesn't need buffer */
#define I2C_MASTER_RX_BUF_DISABLE 0      /*!< I2C master doesn't need buffer */

// Registro MPU6050 pre akcelerometer (X, Y, Z dáta)
#define MPU6050_REG_ACCEL_XOUT_H 0x3B
#define MPU6050_REG_ACCEL_XOUT_L 0x3C
#define MPU6050_REG_ACCEL_YOUT_H 0x3D
#define MPU6050_REG_ACCEL_YOUT_L 0x3E
#define MPU6050_REG_ACCEL_ZOUT_H 0x3F
#define MPU6050_REG_ACCEL_ZOUT_L 0x40

// Registro MPU6050 pre gyroskop (X, Y, Z dáta)
#define MPU6050_REG_GYRO_XOUT_H 0x43
#define MPU6050_REG_GYRO_XOUT_L 0x44
#define MPU6050_REG_GYRO_YOUT_H 0x45
#define MPU6050_REG_GYRO_YOUT_L 0x46
#define MPU6050_REG_GYRO_ZOUT_H 0x47
#define MPU6050_REG_GYRO_ZOUT_L 0x48

// Prevod akcelerometrických hodnôt (rozsah ±2g)
#define ACCEL_SCALE 2.0f / 32768.0f  // Konverzia na g pre ±2g rozsah

// Prevod gyroskopických hodnôt (rozsah ±250 °/s)
#define GYRO_SCALE 250.0f / 32768.0f  // Konverzia na °/s pre ±250 °/s rozsah

#define pocet_dat 10
// Struktúra na ukladanie údajov
typedef struct {
    float accel_x_g;
    float accel_y_g;
    float accel_z_g;
   // float gyro_x_dps;
    //float gyro_y_dps;
    //float gyro_z_dps;
} sensor_data_t;

//potrebné parametre

int vzorkovacia_freq = 10;
int cut_off_freq = 4;
int order = 20;
double senzor_data_array[3][pocet_dat];


// Globálna premenná pre zdieľanie dát medzi taskmi
sensor_data_t sensor_data;
bool write_flag = false;

// Globálny mutex
SemaphoreHandle_t xMutex = NULL;

// Funkcia na nastavenie I2C
esp_err_t i2c_master_init() {
    int i2c_master_port = I2C_MASTER_NUM;
    i2c_config_t conf = {
        .mode = I2C_MODE_MASTER,
        .sda_io_num = I2C_MASTER_SDA_IO,
        .scl_io_num = I2C_MASTER_SCL_IO,
        .sda_pullup_en = GPIO_PULLUP_ENABLE,
        .scl_pullup_en = GPIO_PULLUP_ENABLE,
        .master.clk_speed = I2C_MASTER_FREQ_HZ,
    };
    esp_err_t ret = i2c_param_config(i2c_master_port, &conf);
    if (ret != ESP_OK) {
        return ret;
    }
    ret = i2c_driver_install(i2c_master_port, conf.mode, I2C_MASTER_RX_BUF_DISABLE, I2C_MASTER_TX_BUF_DISABLE, 0);
    return ret;
}

// Funkcia na čítanie 16-bitových hodnôt (ako je akcelerometer alebo gyroskop)
int16_t read_word(i2c_port_t i2c_num, uint8_t reg) {
    uint8_t data[2];
    i2c_cmd_handle_t cmd = i2c_cmd_link_create();
    i2c_master_start(cmd);
    i2c_master_write_byte(cmd, (MPU6050_ADDR << 1) | I2C_MASTER_WRITE, true);
    i2c_master_write_byte(cmd, reg, true);
    i2c_master_start(cmd);
    i2c_master_write_byte(cmd, (MPU6050_ADDR << 1) | I2C_MASTER_READ, true);
    i2c_master_read_byte(cmd, &data[0], I2C_MASTER_ACK);
    i2c_master_read_byte(cmd, &data[1], I2C_MASTER_NACK);
    i2c_master_stop(cmd);
    esp_err_t ret = i2c_master_cmd_begin(i2c_num, cmd, 1000 / portTICK_PERIOD_MS);
    i2c_cmd_link_delete(cmd);

    if (ret == ESP_OK) {
        return (data[0] << 8) | data[1];
    }
    return 0;
}

void lowpass_filter_and_fft(double **acc_data, int n, double *frequencies, double *P1) {
    double *b, *a;
    butterworth_lowpass(order, cut_off_freq, vzorkovacia_freq, &b, &a); // Funkcia na generovanie Butterworth koeficientov
    
    double **acc_filtered = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++) {
        acc_filtered[i] = (double *)malloc(3 * sizeof(double));
    }
    
    for (int i = 0; i < 3; i++) {
        filtfilt(b, a, order, acc_data[i], acc_filtered[i], n); // Aplikácia filtra
    }
    
    double *magnitude_filtered = (double *)malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        magnitude_filtered[i] = sqrt(acc_filtered[i][0] * acc_filtered[i][0] +
                                     acc_filtered[i][1] * acc_filtered[i][1] +
                                     acc_filtered[i][2] * acc_filtered[i][2]);
    }
    
    // FFT výpočet
    int N = n;
    double *Y_real = (double *)malloc(N * sizeof(double));
    double *Y_imag = (double *)malloc(N * sizeof(double));
    
    fft(magnitude_filtered, Y_real, Y_imag, N); // FFT funkcia
    
    for (int i = 0; i <= N / 2; i++) {
        double magnitude = sqrt(Y_real[i] * Y_real[i] + Y_imag[i] * Y_imag[i]) / N;
        P1[i] = (i == 0 || i == N / 2) ? magnitude : 2 * magnitude;
        frequencies[i] = vzorkovacia_freq * i / (double)N;
    }
    
    // Free allocated memory
    free(b);
    free(a);
    free(magnitude_filtered);
    free(Y_real);
    free(Y_imag);
    for (int i = 0; i < n; i++) free(acc_filtered[i]);
    free(acc_filtered);
}

// Task na čítanie dát z MPU6050
void read_sensor_data_task(void *pvParameter) {
    int data_index = 0;

    while (data_index < pocet_dat) {
        // Získanie mutexu pred prístupom k zdieľanej štruktúre
        if (xSemaphoreTake(xMutex, portMAX_DELAY)) {
            // Čítanie hodnôt z akcelerometra a gyroskopu
            int16_t accel_x = read_word(I2C_MASTER_NUM, MPU6050_REG_ACCEL_XOUT_H);
            int16_t accel_y = read_word(I2C_MASTER_NUM, MPU6050_REG_ACCEL_YOUT_H);
            int16_t accel_z = read_word(I2C_MASTER_NUM, MPU6050_REG_ACCEL_ZOUT_H);

            /*int16_t gyro_x = read_word(I2C_MASTER_NUM, MPU6050_REG_GYRO_XOUT_H);
            int16_t gyro_y = read_word(I2C_MASTER_NUM, MPU6050_REG_GYRO_YOUT_H);
            int16_t gyro_z = read_word(I2C_MASTER_NUM, MPU6050_REG_GYRO_ZOUT_H);*/

            // Prevod na fyzikálne jednotky (g pre akcelerometer, °/s pre gyroskop)
            sensor_data.accel_x_g = accel_x * ACCEL_SCALE;
            sensor_data.accel_y_g = accel_y * ACCEL_SCALE;
            sensor_data.accel_z_g = accel_z * ACCEL_SCALE;

            //sensor_data.gyro_x_dps = gyro_x * GYRO_SCALE;
            //sensor_data.gyro_y_dps = gyro_y * GYRO_SCALE;
            //sensor_data.gyro_z_dps = gyro_z * GYRO_SCALE;

            write_flag = true;

            // Uvoľnenie mutexu po prístupe
            xSemaphoreGive(xMutex);
            senzor_data_array[0][data_index] = sensor_data.accel_x_g;
            senzor_data_array[1][data_index] = sensor_data.accel_y_g;
            senzor_data_array[2][data_index] = sensor_data.accel_z_g;
            data_index++;
        }
        vTaskDelay(100 / portTICK_PERIOD_MS); // Zrýchlenie / spomalenie čítania
    }
    vTaskDelete(NULL);
    printf("Koniec merania.\n");
    double * frequencies = (double*)malloc((pocet_dat / 2 + 1) * sizeof(double));
    double * P1 = (double*)malloc((pocet_dat / 2 + 1) * sizeof(double));
    lowpass_filter_and_fft(senzor_data_array, pocet_dat, frequencies, P1);
    for (int i = 0; i <= pocet_dat / 2; i++) {
        printf("%f Hz\n", frequencies[i]);
    }

    free(frequencies);
    free(P1);

}

// Task na vypisovanie dát
/*void print_sensor_data_task(void *pvParameter) {
    int data_index = 0;

    while (data_index < pocet_dat) {
        // Získanie mutexu pred prístupom k zdieľanej štruktúre
        if (xSemaphoreTake(xMutex, portMAX_DELAY)) {
            // Vypísanie hodnot z čítania
            if (write_flag) {
                printf("%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",
                       sensor_data.accel_x_g, sensor_data.accel_y_g, sensor_data.accel_z_g,
                       sensor_data.gyro_x_dps, sensor_data.gyro_y_dps, sensor_data.gyro_z_dps);
                data_index++;
                write_flag = false;
            }

            // Uvoľnenie mutexu po prístupe
            xSemaphoreGive(xMutex);
        }
    }
    printf("Koniec merania.\n");
    vTaskDelete(NULL);
}*/



void process_data_task(void *pvParameter){

}
// Hlavná funkcia
void app_main(void) {
    esp_err_t ret = i2c_master_init();
    if (ret != ESP_OK) {
        ESP_LOGE("I2C", "I2C initialization failed: %s", esp_err_to_name(ret));
        return;
    }

    // Inicializácia MPU6050 (vypnutie spánkového režimu)
    i2c_cmd_handle_t cmd = i2c_cmd_link_create();
    i2c_master_start(cmd);
    i2c_master_write_byte(cmd, (MPU6050_ADDR << 1) | I2C_MASTER_WRITE, true);
    i2c_master_write_byte(cmd, 0x6B, true);  // Register pre Sleep
    i2c_master_write_byte(cmd, 0x00, true);  // Vypneme Sleep režim
    i2c_master_stop(cmd);
    ret = i2c_master_cmd_begin(I2C_MASTER_NUM, cmd, 1000 / portTICK_PERIOD_MS);
    i2c_cmd_link_delete(cmd);

    if (ret != ESP_OK) {
        ESP_LOGE("I2C", "MPU6050 initialization failed: %s", esp_err_to_name(ret));
        return;
    }

    // Vytvorenie mutexu
    xMutex = xSemaphoreCreateMutex();
    if (xMutex == NULL) {
        ESP_LOGE("Mutex", "Mutex creation failed!");
        return;
    }

    // Vytvorenie taskov pre čítanie a vypisovanie
    xTaskCreate(read_sensor_data_task, "ReadSensorData", 2048, NULL, 5, NULL);
    //xTaskCreate(process_data_task, "ProcessData", 2048, NULL, 5, NULL);
    //xTaskCreate(print_sensor_data_task, "PrintSensorData", 2048, NULL, 5, NULL);
}
