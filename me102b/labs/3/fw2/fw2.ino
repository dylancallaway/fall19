#include <Encoder.h>

int led0 = 13;

int pwm_val = 0;
int pwm_pin = 6;
int dir_pin = 7;

int enc_pin1 = 9;
int enc_pin2 = 10;
int enc_val = 0;
int last_enc_val = 0;
int enc_diff = 0;
int enc_speed = 0;
Encoder myEnc(9, 10);
int t0 = 0;
int t1 = 0;
int t_diff;

void setup()
{
  pinMode(led0, OUTPUT);
  pinMode(pwm_pin, OUTPUT);
  pinMode(dir_pin, OUTPUT);
}


void loop()
{
  t0 = millis();
  enc_val = myEnc.read();
  
  enc_diff = last_enc_val - enc_val; // Turning positive gives positive diff
  t_diff = t1 - t0;

  pwm_val = -100 * enc_diff / t_diff;
  
  digitalWrite(dir_pin, (pwm_val / abs(pwm_val)) == 1);
  analogWrite(pwm_pin, abs(pwm_val));


  t1 = millis();
  last_enc_val = enc_val;
  
  delay(100);
}

void toggle(int pin)
{
  digitalWrite(pin, !digitalRead(pin));
}
