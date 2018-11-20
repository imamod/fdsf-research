clc; clear all; close all;
hold on

% ������
bottom = -3*pi; top = 3*pi; span = 2*pi;
poles = bottom:span:top;

% �������
left_border = -6; right_border = 8; bottom_border = bottom - 1; top_border = top + 1;

% ���
line([left_border;0],[0; 0],'linewidth', 1.5, 'color', 'black');
line([0;right_border],[0; 0],'linewidth', 3, 'color', 'black');
line([0;0],[bottom_border; top_border],'linewidth', 3, 'color', 'black');
xlabel('Re\tau','Position', [right_border + 1 bottom_border]);
ylabel('Im\tau','Rotation',0,'Position', [left_border top_border]);

% ����� �����
left = -3;
line([left; left],[bottom_border; top_border], 'color', 'black');
plot(left, poles, 'k*','linewidth', 4);
text(-2.95, 2*pi, 'x<0');
% ������ �����
right = 5;
line([right; right],[bottom_border; top_border], 'color', 'black');
plot(right, poles, 'k*','linewidth', 4);
text(5.05, 2*pi, 'x>0');

% ��������� �������� ticks
set(gca, 'XTick', left_border:right_border)
% set(gca, 'XTickLabel', left_border:1:right_border)
% 
set(gca, 'YTick', poles)
% set(gca, 'YTickLabel', ['-3\pi' '-\pi' '\pi' '3\pi'])

% �������
axis([left_border right_border bottom_border top_border]);

%title('������ ��������������� ������� ��')
