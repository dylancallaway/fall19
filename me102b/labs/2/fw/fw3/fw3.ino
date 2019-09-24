#include <Encoder.h>

int led0 = 13;
int led1 = 3;
int led2 = 4;
int led3 = 5;

int pot = 14;
int pot_val = 0;
long pos  = 0;


void setup()
{
  pinMode(led0, OUTPUT);
  pinMode(led1, OUTPUT);
  pinMode(led2, OUTPUT);
  pinMode(led3, OUTPUT);

  pinMode(pot, INPUT);


}


void loop()
{

  pot_val = analogRead(pot);
  analogWrite(led1, pot_val / 4);
  analogWrite(led2, pot_val / 4);
  analogWrite(led3, pot_val / 4);


}
