#include <Encoder.h>

int led0 = 13;
int led1 = 3;
int led2 = 4;
int led3 = 5;
int enc_pin1 = 9;
int enc_pin2 = 10;
int enc_val;
Encoder myEnc(9, 10);
long pos  = 0;


void setup()
{
  pinMode(led0, OUTPUT);
  pinMode(led1, OUTPUT);
  pinMode(led2, OUTPUT);
  pinMode(led3, OUTPUT);

}


void loop()
{
  if (pos < 0)
  {
    digitalWrite(led1, LOW);
    digitalWrite(led2, LOW);
    digitalWrite(led3, LOW);
  }
  if (pos > 20)
  {
    digitalWrite(led1, HIGH);
    digitalWrite(led2, LOW);
    digitalWrite(led3, LOW);
  }
  if (pos > 40)
  {
    digitalWrite(led1, HIGH);
    digitalWrite(led2, HIGH);
    digitalWrite(led3, LOW);
  }
  if (pos > 60)
  {
    digitalWrite(led1, HIGH);
    digitalWrite(led2, HIGH);
    digitalWrite(led3, HIGH);
  }



  pos = myEnc.read();

}
