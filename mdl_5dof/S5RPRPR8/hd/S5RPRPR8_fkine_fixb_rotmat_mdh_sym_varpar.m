% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-31 18:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S5RPRPR8_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR8_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR8_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:21:11
% EndTime: 2019-12-31 18:21:11
% DurationCPUTime: 0.10s
% Computational Cost: add. (128->42), mult. (83->48), div. (0->0), fcn. (125->10), ass. (0->31)
t16 = sin(pkin(9));
t15 = qJ(1) + pkin(8);
t8 = sin(t15);
t35 = t8 * t16;
t21 = cos(qJ(3));
t34 = t8 * t21;
t10 = cos(t15);
t33 = t10 * t21;
t32 = t16 * t21;
t17 = cos(pkin(9));
t31 = t17 * t21;
t30 = pkin(5) + 0;
t20 = sin(qJ(1));
t29 = t20 * pkin(1) + 0;
t22 = cos(qJ(1));
t28 = t22 * pkin(1) + 0;
t27 = t8 * pkin(2) + t29;
t11 = qJ(2) + t30;
t26 = t10 * pkin(2) + t8 * pkin(6) + t28;
t18 = -pkin(7) - qJ(4);
t19 = sin(qJ(3));
t6 = t17 * pkin(4) + pkin(3);
t25 = -t18 * t19 + t21 * t6;
t24 = pkin(3) * t21 + qJ(4) * t19;
t23 = -t10 * pkin(6) + t27;
t14 = pkin(9) + qJ(5);
t9 = cos(t14);
t7 = sin(t14);
t2 = t10 * t19;
t1 = t8 * t19;
t3 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t22, -t20, 0, 0; t20, t22, 0, 0; 0, 0, 1, t30; 0, 0, 0, 1; t10, -t8, 0, t28; t8, t10, 0, t29; 0, 0, 1, t11; 0, 0, 0, 1; t33, -t2, t8, t26; t34, -t1, -t10, t23; t19, t21, 0, t11; 0, 0, 0, 1; t10 * t31 + t35, -t10 * t32 + t8 * t17, t2, t24 * t10 + t26; -t10 * t16 + t8 * t31, -t10 * t17 - t8 * t32, t1, t24 * t8 + t23; t19 * t17, -t19 * t16, -t21, t19 * pkin(3) - t21 * qJ(4) + t11; 0, 0, 0, 1; t9 * t33 + t8 * t7, -t7 * t33 + t8 * t9, t2, pkin(4) * t35 + t25 * t10 + t26; -t10 * t7 + t9 * t34, -t10 * t9 - t7 * t34, t1, t25 * t8 + (-pkin(4) * t16 - pkin(6)) * t10 + t27; t19 * t9, -t19 * t7, -t21, t21 * t18 + t19 * t6 + t11; 0, 0, 0, 1;];
T_ges = t3;
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
