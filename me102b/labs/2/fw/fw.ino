int led0 = 13;
int led1 = 3;
int led2 = 4;
int led3 = 5;
volatile bool flag = 0;
int int_pin = 12;
int t0 = 0;
int t1 = 0;

void interrupt_cb() {
  t1 = millis();
  if ((t1 - t0) > 200)
  {
    flag = !flag;
  }
  t0 = t1;
}

void setup()
{
  pinMode(led0, OUTPUT);
  pinMode(led1, OUTPUT);
  pinMode(led2, OUTPUT);
  pinMode(led3, OUTPUT);
  pinMode(int_pin, INPUT_PULLUP);


  attachInterrupt(digitalPinToInterrupt(int_pin), interrupt_cb, FALLING);
}

void loop()
{
  if (flag)
  {
    toggle(led0);
    toggle(led1);
    toggle(led2);
    toggle(led3);

    delay(500);
  }
}

void toggle(int pin)
{
  digitalWrite(pin, !digitalRead(pin));
}
