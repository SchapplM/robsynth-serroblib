% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% T_c_mdh [4x4x(5+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   6:  mdh base (link 0) -> mdh frame (6-1), link (6-1)
%   ...
%   5+1:  mdh base (link 0) -> mdh frame (5)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:49
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S5PRRRP5_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP5_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP5_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:47:37
% EndTime: 2019-12-05 16:47:37
% DurationCPUTime: 0.15s
% Computational Cost: add. (110->50), mult. (117->56), div. (0->0), fcn. (168->8), ass. (0->36)
t26 = -pkin(7) - pkin(6);
t22 = sin(qJ(3));
t39 = t22 * pkin(3);
t24 = cos(qJ(3));
t10 = t24 * pkin(3) + pkin(2);
t20 = sin(pkin(8));
t38 = t20 * t22;
t25 = cos(qJ(2));
t37 = t20 * t25;
t21 = cos(pkin(8));
t36 = t21 * t25;
t35 = t22 * t25;
t19 = qJ(3) + qJ(4);
t11 = sin(t19);
t23 = sin(qJ(2));
t34 = t23 * t11;
t33 = t24 * t25;
t32 = t20 * pkin(1) + 0;
t17 = qJ(1) + 0;
t31 = t21 * pkin(1) + t20 * pkin(5) + 0;
t30 = pkin(2) * t25 + pkin(6) * t23;
t18 = -qJ(5) + t26;
t12 = cos(t19);
t5 = pkin(4) * t12 + t10;
t29 = -t18 * t23 + t25 * t5;
t28 = t10 * t25 - t23 * t26;
t27 = -t21 * pkin(5) + t32;
t9 = t21 * t23;
t8 = t20 * t23;
t7 = t23 * t12;
t6 = pkin(4) * t11 + t39;
t4 = t20 * t11 + t12 * t36;
t3 = -t11 * t36 + t20 * t12;
t2 = -t21 * t11 + t12 * t37;
t1 = -t11 * t37 - t21 * t12;
t13 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t21, -t20, 0, 0; t20, t21, 0, 0; 0, 0, 1, t17; 0, 0, 0, 1; t36, -t9, t20, t31; t37, -t8, -t21, t27; t23, t25, 0, t17; 0, 0, 0, 1; t21 * t33 + t38, t20 * t24 - t21 * t35, t9, t21 * t30 + t31; t20 * t33 - t21 * t22, -t20 * t35 - t21 * t24, t8, t20 * t30 + t27; t23 * t24, -t23 * t22, -t25, t23 * pkin(2) - t25 * pkin(6) + t17; 0, 0, 0, 1; t4, t3, t9, pkin(3) * t38 + t21 * t28 + t31; t2, t1, t8, (-pkin(5) - t39) * t21 + t28 * t20 + t32; t7, -t34, -t25, t23 * t10 + t25 * t26 + t17; 0, 0, 0, 1; t4, t3, t9, t20 * t6 + t21 * t29 + t31; t2, t1, t8, (-pkin(5) - t6) * t21 + t29 * t20 + t32; t7, -t34, -t25, t25 * t18 + t23 * t5 + t17; 0, 0, 0, 1;];
T_ges = t13;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
Tc_stack = NaN(3*size(T_c_mdh,3),4);
% Zusätzliche Ausgabe: Als 2D-array gestapelt, ohne Zeile mit 0001
for i = 1:size(T_c_mdh,3), Tc_stack((i-1)*3+1:3*i,1:4) = T_c_mdh(1:3,1:4,i); end
