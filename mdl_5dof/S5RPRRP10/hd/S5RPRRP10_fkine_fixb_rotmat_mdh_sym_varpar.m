% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S5RPRRP10_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP10_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP10_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:50:48
% EndTime: 2019-12-31 18:50:48
% DurationCPUTime: 0.12s
% Computational Cost: add. (111->43), mult. (92->44), div. (0->0), fcn. (138->8), ass. (0->33)
t16 = pkin(8) + qJ(3);
t13 = sin(t16);
t22 = sin(qJ(4));
t36 = t13 * t22;
t23 = sin(qJ(1));
t35 = t23 * t22;
t24 = cos(qJ(4));
t34 = t23 * t24;
t25 = cos(qJ(1));
t33 = t25 * t22;
t32 = t25 * t24;
t17 = pkin(5) + 0;
t19 = cos(pkin(8));
t10 = t19 * pkin(2) + pkin(1);
t31 = t25 * t10 + 0;
t18 = sin(pkin(8));
t30 = t18 * pkin(2) + t17;
t21 = -pkin(6) - qJ(2);
t29 = t23 * t10 + t25 * t21 + 0;
t14 = cos(t16);
t28 = pkin(3) * t14 + pkin(7) * t13;
t12 = t24 * pkin(4) + pkin(3);
t20 = -qJ(5) - pkin(7);
t27 = t12 * t14 - t13 * t20;
t26 = -t23 * t21 + t31;
t9 = t25 * t13;
t8 = t23 * t13;
t7 = t13 * t24;
t4 = t14 * t32 + t35;
t3 = -t14 * t33 + t34;
t2 = t14 * t34 - t33;
t1 = -t14 * t35 - t32;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t25, -t23, 0, 0; t23, t25, 0, 0; 0, 0, 1, t17; 0, 0, 0, 1; t25 * t19, -t25 * t18, t23, t25 * pkin(1) + t23 * qJ(2) + 0; t23 * t19, -t23 * t18, -t25, t23 * pkin(1) - t25 * qJ(2) + 0; t18, t19, 0, t17; 0, 0, 0, 1; t25 * t14, -t9, t23, t26; t23 * t14, -t8, -t25, t29; t13, t14, 0, t30; 0, 0, 0, 1; t4, t3, t9, t25 * t28 + t26; t2, t1, t8, t23 * t28 + t29; t7, -t36, -t14, t13 * pkin(3) - t14 * pkin(7) + t30; 0, 0, 0, 1; t4, t3, t9, t27 * t25 + (pkin(4) * t22 - t21) * t23 + t31; t2, t1, t8, -pkin(4) * t33 + t23 * t27 + t29; t7, -t36, -t14, t13 * t12 + t14 * t20 + t30; 0, 0, 0, 1;];
T_ges = t5;
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
