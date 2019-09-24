#include <Encoder.h>

int led0 = 13;

int pot = 14;
int pot_val = 0;

int pwm_val = 0;
int pwm_pin = 6;
int dir_pin = 7;

int enc_pin1 = 9;
int enc_pin2 = 10;
int enc_val;
Encoder myEnc(9, 10);


void setup()
{
  pinMode(led0, OUTPUT);
  pinMode(pot, INPUT);
  pinMode(pwm_pin, OUTPUT);
  pinMode(dir_pin, OUTPUT);
}


void loop()
{
  pot_val = analogRead(pot); // 0 - 1024
  pwm_val = map(pot_val, 0, 1024, -256, 256);
  digitalWrite(dir_pin, (pwm_val / abs(pwm_val)) == 1);
  analogWrite(pwm_pin, abs(pwm_val));

}

void toggle(int pin)
{
  digitalWrite(pin, !digitalRead(pin));
}
