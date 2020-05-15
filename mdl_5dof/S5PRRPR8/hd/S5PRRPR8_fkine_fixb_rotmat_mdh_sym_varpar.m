% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-31 17:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S5PRRPR8_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR8_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR8_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:42:12
% EndTime: 2019-12-31 17:42:12
% DurationCPUTime: 0.10s
% Computational Cost: add. (122->42), mult. (69->44), div. (0->0), fcn. (110->10), ass. (0->29)
t23 = -pkin(6) - pkin(5);
t17 = sin(pkin(8));
t16 = qJ(2) + qJ(3);
t9 = pkin(9) + t16;
t5 = sin(t9);
t34 = t17 * t5;
t18 = cos(pkin(8));
t33 = t18 * t5;
t22 = cos(qJ(2));
t8 = t22 * pkin(2) + pkin(1);
t19 = sin(qJ(5));
t32 = t17 * t19;
t21 = cos(qJ(5));
t31 = t17 * t21;
t30 = t18 * t19;
t29 = t18 * t21;
t14 = qJ(1) + 0;
t15 = -qJ(4) + t23;
t11 = cos(t16);
t3 = pkin(3) * t11 + t8;
t28 = t18 * t15 + t17 * t3 + 0;
t20 = sin(qJ(2));
t27 = t20 * pkin(2) + t14;
t6 = cos(t9);
t26 = pkin(4) * t6 + pkin(7) * t5;
t10 = sin(t16);
t25 = pkin(3) * t10 + t27;
t24 = -t17 * t15 + t18 * t3 + 0;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t18, -t17, 0, 0; t17, t18, 0, 0; 0, 0, 1, t14; 0, 0, 0, 1; t18 * t22, -t18 * t20, t17, t18 * pkin(1) + t17 * pkin(5) + 0; t17 * t22, -t17 * t20, -t18, t17 * pkin(1) - t18 * pkin(5) + 0; t20, t22, 0, t14; 0, 0, 0, 1; t18 * t11, -t18 * t10, t17, -t17 * t23 + t18 * t8 + 0; t17 * t11, -t17 * t10, -t18, t17 * t8 + t18 * t23 + 0; t10, t11, 0, t27; 0, 0, 0, 1; t18 * t6, -t33, t17, t24; t17 * t6, -t34, -t18, t28; t5, t6, 0, t25; 0, 0, 0, 1; t6 * t29 + t32, -t6 * t30 + t31, t33, t26 * t18 + t24; t6 * t31 - t30, -t6 * t32 - t29, t34, t26 * t17 + t28; t5 * t21, -t5 * t19, -t6, t5 * pkin(4) - t6 * pkin(7) + t25; 0, 0, 0, 1;];
T_ges = t1;
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
