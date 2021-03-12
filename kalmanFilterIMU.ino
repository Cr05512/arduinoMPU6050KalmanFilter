#include "math.h"
#include "I2Cdev.h"
#include "MPU6050.h"
#include "Wire.h"

MPU6050 imu;

#define degToRad(angleInDegrees) ((angleInDegrees) * M_PI / 180.0)
#define radToDeg(angleInRadians) ((angleInRadians) * 180.0 / M_PI)

int16_t ax, ay, az;
int16_t gx, gy, gz;

float accX, accY, accZ, gyroX, gyroY, gyroZ;

float accRes = 2.0; //2g default
float divider = 32767.0;
float gyroRes = 250.0; // 250°/s default
float phiXa=0.0, phiYa=0.0, phiX=0.0, phiY=0.0, phiZ=0.0, bx = 0.0, by = 0.0, bz = 0.0;
uint16_t t0 = 0, t1 = 0;
float dt = 0.0;

// kalman filter params
float Qa = 0.001;
float Qb = 0.003;
float Racc = 0.03;

float P00 = 1.0, P01 = 0.0, P02 = 0.0, P03 = 0.0, P04 = 0.0, P05 = 0.0;
float P10 = 0.0, P11 = 1.0, P12 = 0.0, P13 = 0.0, P14 = 0.0, P15 = 0.0;
float P20 = 0.0, P21 = 0.0, P22 = 1.0, P23 = 0.0, P24 = 0.0, P25 = 0.0;
float P30 = 0.0, P31 = 0.0, P32 = 0.0, P33 = 1.0, P34 = 0.0, P35 = 0.0;
float P40 = 0.0, P41 = 0.0, P42 = 0.0, P43 = 0.0, P44 = 1.0, P45 = 0.0;
float P50 = 0.0, P51 = 0.0, P52 = 0.0, P53 = 0.0, P54 = 0.0, P55 = 1.0;

float S00 = 1.0, S01 = 0.0;
float S10 = 1.0, S11 = 0.0;
float denS = 1.0;

float S00inv = 1.0, S01inv = 0.0, S10inv = 0.0, S11inv = 1.0;

float K00 = 0.0, K01 = 0.0;
float K10 = 0.0, K11 = 0.0;
float K20 = 0.0, K21 = 0.0;
float K30 = 0.0, K31 = 0.0;
float K40 = 0.0, K41 = 0.0;
float K50 = 0.0, K51 = 0.0;

float errX = 0.0, errY = 0.0;


void readAngles();
void kalmanFilter();

void setup() {
  Wire.begin();
  Serial.begin(115200);
  imu.initialize();

}

void loop() {
  readAngles(); // We get the measure
    
  kalmanFilter(); // We filter it and update our angles

  Serial.print("Roll: "); Serial.print(radToDeg(phiX)); Serial.print(", ");
  Serial.print("Pitch: "); Serial.print(radToDeg(phiY)); Serial.print(", ");
  Serial.print("Yaw: "); Serial.println(radToDeg(phiZ)); 

}

void readAngles()
{
  imu.getMotion6(&ax, &ay, &az, &gx, &gy, &gz);

  t1 = millis();
  dt = (t1 - t0)/1000.0;
  t0 = t1;
  
  accX = ax*accRes/divider;
  accY = ay*accRes/divider;
  accZ = az*accRes/divider;
  phiXa = atan2(accY,accZ);
  phiYa = atan2(-accX,sqrt(accY*accY+accZ*accZ));
    
  gyroX = degToRad(gx*gyroRes/divider);
  gyroY = degToRad(gy*gyroRes/divider);
  gyroZ = degToRad(gz*gyroRes/divider);

}

void kalmanFilter()
{
  //Measurement update

  //Innovavtion covariance S
  S00 = Racc + P00;
  S01 = P02;
  S10 = P20;
  S11 = Racc + P22;

  denS = S00*S11 - S01*S10;
  //Inverse Innovation covariance Sinv
  S00inv = S11/denS;
  S01inv = -S01/denS;
  S10inv = -S10/denS;
  S11inv = S00/denS;

  S00 = S00inv;
  S01 = S01inv;
  S10 = S10inv;
  S11 = S11inv;



  //Kalman Gain
  K00 = P00*S00 + P02*S10;
  K01 = P00*S01 + P02*S11;
  K10 = P10*S00 + P12*S10;
  K11 = P10*S01 + P12*S11;
  K20 = P20*S00 + P22*S10;
  K21 = P20*S01 + P22*S11;
  K30 = P30*S00 + P32*S10;
  K31 = P30*S01 + P32*S11;
  K40 = P40*S00 + P42*S10;
  K41 = P40*S01 + P42*S11;
  K50 = P50*S00 + P52*S10;
  K51 = P50*S01 + P52*S11;

  //State estimation
  errX = phiXa - phiX;
  errY = phiYa - phiY;

  
  phiX = phiX + K00*errX + K01*errY;
  bx = bx + K10*errX + K11*errY;
  phiY = phiY + K20*errX + K21*errY;
  by = by + K30*errX + K31*errY;
  phiZ = phiZ + K40*errX + K41*errY;
  bz = bz + K50*errX + K51*errY;
  //Covariance update
  

  P00 = - P00*(K00 - 1) - K01*P20;
  P01 = - P01*(K00 - 1) - K01*P21;
  P02 = - P02*(K00 - 1) - K01*P22;
  P03 = - P03*(K00 - 1) - K01*P23;
  P04 = - P04*(K00 - 1) - K01*P24;
  P05 = - P05*(K00 - 1) - K01*P25;

  P10 = P10 - K10*P00 - K11*P20;
  P11 = P11 - K10*P01 - K11*P21;
  P12 = P12 - K10*P02 - K11*P22;
  P13 = P13 - K10*P03 - K11*P23;
  P14 = P14 - K10*P04 - K11*P24;
  P15 = P15 - K10*P05 - K11*P25;

  P20 = P20*(1-K21) - K20*P00;
  P21 = P21*(1-K21) - K20*P01;
  P22 = P22*(1-K21) - K20*P02;
  P23 = P23*(1-K21) - K20*P03;
  P24 = P24*(1-K21) - K20*P04;
  P25 = P25*(1-K21) - K20*P05;

  P30 = P30 - K30*P00 - K31*P20;
  P31 = P31 - K30*P01 - K31*P21;
  P32 = P32 - K30*P02 - K31*P22;
  P33 = P33 - K30*P03 - K31*P23;
  P34 = P34 - K30*P04 - K31*P24;
  P35 = P35 - K30*P05 - K31*P25;

  P40 = P40 - K40*P00 - K41*P20;
  P41 = P41 - K40*P01 - K41*P21;
  P42 = P42 - K40*P02 - K41*P22;
  P43 = P43 - K40*P03 - K41*P23;
  P44 = P44 - K40*P04 - K41*P24;
  P45 = P45 - K40*P05 - K41*P25;

  P50 = P50 - K50*P00 - K51*P20;
  P51 = P51 - K50*P01 - K51*P21;
  P52 = P52 - K50*P02 - K51*P22;
  P53 = P53 - K50*P03 - K51*P23;
  P54 = P54 - K50*P04 - K51*P24;
  P55 = P55 - K50*P05 - K51*P25;

  //Model Prediction

  //State prediction
  phiX = phiX - bx*dt + gyroX*dt;
  phiY = phiY - by*dt + gyroY*dt;
  phiZ = phiZ - bz*dt + gyroZ*dt;

  //Covariance prediction

  P00 = P00 + Qa - P10*dt - dt*(P01 - P11*dt);
  P01 = P01 - P11*dt;
  P02 = P02 - P12*dt - dt*(P03 - P13*dt);
  P03 = P03 - P13*dt;
  P04 = P04 - P14*dt - dt*(P05 - P15*dt);
  P05 = P05 - P15*dt;

  P10 = P10 - P11*dt;
  P11 = P11 + Qb;
  P12 = P12 - P13*dt;
  P13 = P13;
  P14 = P14 - P15*dt;
  P15 = P15;

  P20 = P20 - P30*dt - dt*(P21 - P31*dt);
  P21 = P21 - P31*dt;
  P22 = P22 - P32*dt - dt*(P23 - P33*dt) + Qa;
  P23 = P23 - P33*dt;
  P24 = P24 - P34*dt - dt*(P25 - P35*dt);
  P25 = P25 - P35*dt;

  P30 = P30 - P31*dt;
  P31 = P31;
  P32 = P32 - P33*dt;
  P33 = P33 + Qb;
  P34 = P34 - P35*dt;
  P35 = P35;

  P40 = P40 - P50*dt - dt*(P41 - P51*dt);
  P41 = P41 - P51*dt;
  P42 = P42 - P52*dt - dt*(P43 - P53*dt);
  P43 = P43 - P53*dt;
  P44 = P44 - P54*dt - dt*(P45 - P55*dt) + Qa;
  P45 = P45 - P55*dt;

  P50 = P50 - P51*dt;
  P51 = P51;
  P52 = P52 - P53*dt;
  P53 = P53;
  P54 = P54 - P55*dt;
  P55 = P55 + Qb;
  
  
}
