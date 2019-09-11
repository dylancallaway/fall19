int led = 13;

void setup()
{
  pinMode(led, OUTPUT);
}

void loop()
{
  toggle(led);
  delay(500);
}

void toggle(int pin)
{
  digitalWrite(pin, !digitalRead(pin));
}
