% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
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
% Datum: 2019-12-31 19:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S5RRPPR8_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR8_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR8_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:37:59
% EndTime: 2019-12-31 19:37:59
% DurationCPUTime: 0.13s
% Computational Cost: add. (89->44), mult. (119->46), div. (0->0), fcn. (167->8), ass. (0->32)
t20 = sin(qJ(1));
t19 = sin(qJ(2));
t35 = qJ(3) * t19;
t21 = cos(qJ(2));
t6 = t20 * t21;
t39 = pkin(2) * t6 + t20 * t35;
t16 = sin(pkin(8));
t38 = t19 * t16;
t37 = t20 * t19;
t22 = cos(qJ(1));
t36 = t22 * t19;
t7 = t22 * t21;
t15 = pkin(5) + 0;
t34 = t20 * pkin(1) + 0;
t33 = t19 * pkin(2) + t15;
t32 = t22 * pkin(1) + t20 * pkin(6) + 0;
t31 = t34 + t39;
t14 = pkin(8) + qJ(5);
t8 = sin(t14);
t9 = cos(t14);
t30 = t19 * t9 - t21 * t8;
t29 = t19 * t8 + t21 * t9;
t17 = cos(pkin(8));
t28 = -t21 * t16 + t19 * t17;
t27 = t21 * t17 + t38;
t26 = -t22 * pkin(6) + t34;
t25 = pkin(2) * t7 + t22 * t35 + t32;
t5 = t17 * pkin(4) + pkin(3);
t24 = pkin(4) * t38 + t21 * t5;
t23 = -t21 * qJ(3) + t33;
t18 = -pkin(7) - qJ(4);
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t22, -t20, 0, 0; t20, t22, 0, 0; 0, 0, 1, t15; 0, 0, 0, 1; t7, -t36, t20, t32; t6, -t37, -t22, t26; t19, t21, 0, t15; 0, 0, 0, 1; t7, t20, t36, t25; t6, -t22, t37, t26 + t39; t19, 0, -t21, t23; 0, 0, 0, 1; t27 * t22, t28 * t22, -t20, pkin(3) * t7 - t20 * qJ(4) + t25; t27 * t20, t28 * t20, t22, pkin(3) * t6 + (-pkin(6) + qJ(4)) * t22 + t31; t28, -t27, 0, t19 * pkin(3) + t23; 0, 0, 0, 1; t29 * t22, t30 * t22, -t20, t20 * t18 + t22 * t24 + t25; t29 * t20, t30 * t20, t22, (-pkin(6) - t18) * t22 + t24 * t20 + t31; t30, -t29, 0, t19 * t5 + (-pkin(4) * t16 - qJ(3)) * t21 + t33; 0, 0, 0, 1;];
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
