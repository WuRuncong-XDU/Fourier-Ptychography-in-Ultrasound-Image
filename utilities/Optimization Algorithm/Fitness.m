function result = Fitness(params, des)
try
    % [out, ~] = Ptychography30MHzSimulation1(params);
    % [out, ~] = Ptychography30MHzPigEye(params);
    % [out, ~] = Ptychography30MHzTungsten(params);
    [out, ~] = Ptychography30MHzSimulation2(params);
catch message
    result = 1000;
    fprintf(['Fitness.m - ' message.message, ' : ', num2str(params), '\n']);
    return;
end

% ����Ȩֵ
a = 1/3;
b = 1/3;
c = 1/3;
% a = 3/20;
% b = 3/20;
% c = 7/10;

% ������������У��
result = a * ((out.resH-des.resH)*25)^2 + b * ((out.resR-des.resR)*25)^2 + c * ((out.contrast-des.contrast) / 10)^2;

% result = 0.5 * ((out.resH-des.resH)*25)^2 + 0.5 * ((out.resR-des.resR)*25)^2;

% �����¼�ļ�
writePg(params, result);
end
