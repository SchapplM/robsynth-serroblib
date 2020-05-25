% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4RRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% T_c_mdh [4x4x(4+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   5:  mdh base (link 0) -> mdh frame (5-1), link (5-1)
%   ...
%   4+1:  mdh base (link 0) -> mdh frame (4)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-01-31 13:16
% Revision: 9ef80adae39e3cd5824e7abdb6e4e1e7895c437e (2019-01-31)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S4RRPR1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-01-31 13:16:31
% EndTime: 2019-01-31 13:16:31
% DurationCPUTime: 0.06s
% Computational Cost: add. (58->17), mult. (12->6), div. (0->0), fcn. (28->8), ass. (0->19)
t13 = qJ(1) + qJ(2);
t22 = pkin(4) + 0;
t14 = sin(qJ(1));
t21 = t14 * pkin(1) + 0;
t15 = cos(qJ(1));
t20 = t15 * pkin(1) + 0;
t19 = pkin(5) + t22;
t9 = sin(t13);
t18 = pkin(2) * t9 + t21;
t10 = cos(t13);
t17 = pkin(2) * t10 + t20;
t8 = pkin(7) + t13;
t16 = qJ(3) + t19;
t7 = qJ(4) + t8;
t4 = cos(t8);
t3 = sin(t8);
t2 = cos(t7);
t1 = sin(t7);
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t15, -t14, 0, 0; t14, t15, 0, 0; 0, 0, 1, t22; 0, 0, 0, 1; t10, -t9, 0, t20; t9, t10, 0, t21; 0, 0, 1, t19; 0, 0, 0, 1; t4, -t3, 0, t17; t3, t4, 0, t18; 0, 0, 1, t16; 0, 0, 0, 1; t2, -t1, 0, pkin(3) * t4 + t17; t1, t2, 0, pkin(3) * t3 + t18; 0, 0, 1, pkin(6) + t16; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,4+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,4+1]); end % symbolisch
for i = 1:4+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
Tc_stack = NaN(3*size(T_c_mdh,3),4);
% Zusätzliche Ausgabe: Als 2D-array gestapelt, ohne Zeile mit 0001
for i = 1:size(T_c_mdh,3), Tc_stack((i-1)*3+1:3*i,1:4) = T_c_mdh(1:3,1:4,i); end
