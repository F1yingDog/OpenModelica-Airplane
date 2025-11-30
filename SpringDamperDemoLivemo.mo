package SpringDemper
  model Spring
    ForcePort forcePort annotation(
      Placement(visible = true, transformation(origin = {0, 100}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 100}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    parameter Real k = 1 "Stiffness factor";
    Real forceFactor(start=1);
    Modelica.Units.SI.Distance s(start = 3);
    Boolean e;
  equation
    s = forcePort.x;
    forcePort.f = s*forceFactor*k;
    
  
  algorithm
    when s<-10 then
        e:=true;
        forceFactor:=1;
    end when;
    
  
    annotation(
      Icon(graphics = {Line(origin = {0, 9}, points = {{0, 91}, {0, 71}, {40, 51}, {-40, 31}, {40, 9}, {-38, -9}, {40, -31}, {-40, -49}, {0, -69}, {0, -89}}), Line(origin = {-0.146447, -85.2071}, points = {{-19.8536, 5.20711}, {20.1464, 5.20711}, {10.1464, -4.79289}, {10.1464, -4.79289}}), Line(origin = {-5, -85}, points = {{5, 5}, {-5, -5}}), Line(origin = {-25, -85}, points = {{5, 5}, {-5, -5}})}, coordinateSystem(extent = {{-100, -100}, {100, 100}})));
  annotation(
      Icon(graphics = {Line(origin = {0.0963055, 34.7929}, points = {{-0.0963055, 63.2071}, {-0.0963055, 45.2071}, {19.9037, 35.2071}, {-20.0963, 25.2071}, {19.9037, 5.20711}, {-20.0963, -14.7929}, {-0.0963055, -26.7929}, {-0.0963055, -54.7929}, {-20.0963, -54.7929}, {19.9037, -54.7929}, {11.9037, -62.7929}, {11.9037, -62.7929}}), Line(origin = {5, -23}, points = {{3, 3}, {-3, -3}}), Line(origin = {-9, -23}, points = {{3, 3}, {-3, -3}}), Line(origin = {-21, -23}, points = {{3, 3}, {-3, -3}})}));
  end Spring;

  model Mass
    Modelica.Units.SI.Position p;
    Modelica.Units.SI.Velocity v;
    Modelica.Units.SI.Acceleration a;
    parameter Modelica.Units.SI.Mass m = 1 "Mass";
    constant Modelica.Units.SI.Acceleration g = Modelica.Constants.g_n;
    ForcePort forcePort annotation(
      Placement(visible = true, transformation(origin = {0, -62}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, -62}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    p=forcePort.x;
    der(v) = a;
    der(forcePort.x) = v;
    m*a= -m*g +forcePort.f;
    annotation(
      Diagram,
      Icon(graphics = {Rectangle(lineColor = {255, 255, 255}, fillPattern = FillPattern.CrossDiag, lineThickness = 0.5, extent = {{-64, 46}, {64, -46}}), Line(origin = {0, -53}, points = {{0, 7}, {0, -7}})}));
  end Mass;

  connector ForcePort
    Modelica.Units.SI.Position x;
    flow Modelica.Units.SI.Force f;
    annotation(
      Icon(graphics = {Rectangle(fillPattern = FillPattern.Solid, extent = {{-20, 20}, {20, -20}})}));
  end ForcePort;

  model Example
  SpringDemper.Mass mass(m = 1) annotation(
      Placement(visible = true, transformation(origin = {-12, 16}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    SpringDemper.Demper demper(d = 0.2)  annotation(
      Placement(visible = true, transformation(origin = {8, -16}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  SpringDemper.Spring spring(s(start = 0))  annotation(
      Placement(visible = true, transformation(origin = {-27, -15}, extent = {{-15, -15}, {15, 15}}, rotation = 0)));
  equation
    connect(demper.forcePort, mass.forcePort) annotation(
      Line(points = {{8, -10}, {8, 10}, {-12, 10}}));
  connect(spring.forcePort, mass.forcePort) annotation(
      Line(points = {{-27, 0}, {-13.5, 0}, {-13.5, 10}, {-12, 10}}));
  protected    annotation(
      Diagram);
  end Example;

  model Demper
    parameter Real d = 1 "Dampening factor";
    ForcePort forcePort annotation(
      Placement(visible = true, transformation(origin = {0, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    forcePort.f = d*der(forcePort.x);
    annotation(
      Icon(graphics = {Rectangle(fillPattern = FillPattern.Solid, extent = {{-40, 20}, {40, -20}}), Line(origin = {-40, 30}, points = {{0, 10}, {0, -10}}), Line(origin = {40, 30}, points = {{0, 10}, {0, -10}}), Line(origin = {0, 40}, points = {{0, 20}, {0, -20}}), Line(origin = {0.384111, -46.2042}, points = {{-0.384111, 26.2042}, {-0.384111, -13.7958}, {-40.3841, -13.7958}, {39.6159, -13.7958}, {29.6159, -25.7958}}), Line(origin = {15, -65}, points = {{5, 5}, {-5, -5}}), Line(origin = {-5, -64}, points = {{5, 4}, {-5, -4}}), Line(origin = {-25, -64}, points = {{5, 4}, {-5, -4}}), Line(origin = {-46, -64}, points = {{6, 4}, {-6, -4}})}));
  end Demper;

  model MassSpringDemper
  
    parameter Real k=1;
    parameter Real b=0.2;
    parameter Real g=-9.81;
    parameter Real m=1;
    
    Real Fk;
    Real Fb;
    Real Fg;
    Real FRes;
    Real x(start=0.3);
    Real v;
  equation
  Fg=m*g;
    Fk=k*x;
    Fb=b*der(x);
    FRes=Fg-Fk-Fb;
    v=der(x);
    FRes=m*der(v);
  
  
  annotation(
      Icon);
end MassSpringDemper;
  annotation(
    uses(Modelica(version = "4.0.0")));
end SpringDemper;

