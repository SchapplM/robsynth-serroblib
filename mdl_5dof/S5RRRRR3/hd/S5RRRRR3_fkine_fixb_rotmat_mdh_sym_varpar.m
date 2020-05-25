% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
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
% Datum: 2019-12-05 18:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S5RRRRR3_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR3_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:54:59
% EndTime: 2019-12-05 18:54:59
% DurationCPUTime: 0.10s
% Computational Cost: add. (96->37), mult. (82->43), div. (0->0), fcn. (128->10), ass. (0->31)
t18 = sin(qJ(1));
t15 = qJ(2) + qJ(3);
t9 = sin(t15);
t3 = t18 * t9;
t21 = cos(qJ(1));
t4 = t21 * t9;
t11 = cos(t15);
t35 = t18 * t11;
t16 = sin(qJ(4));
t34 = t18 * t16;
t19 = cos(qJ(4));
t33 = t18 * t19;
t20 = cos(qJ(2));
t32 = t18 * t20;
t31 = t21 * t11;
t30 = t21 * t16;
t29 = t21 * t19;
t28 = t21 * t20;
t13 = pkin(4) + 0;
t27 = pkin(1) * t32 + 0;
t26 = pkin(1) * t28 + 0;
t25 = pkin(5) * t3 + t27;
t24 = pkin(5) * t4 + t26;
t17 = sin(qJ(2));
t23 = t17 * pkin(1) + t13;
t22 = -t11 * pkin(5) + t23;
t14 = qJ(4) + qJ(5);
t10 = cos(t14);
t8 = sin(t14);
t7 = t19 * pkin(3) + pkin(2);
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t21, -t18, 0, 0; t18, t21, 0, 0; 0, 0, 1, t13; 0, 0, 0, 1; t28, -t21 * t17, t18, 0; t32, -t18 * t17, -t21, 0; t17, t20, 0, t13; 0, 0, 0, 1; t31, -t4, t18, t26; t35, -t3, -t21, t27; t9, t11, 0, t23; 0, 0, 0, 1; t11 * t29 + t34, -t11 * t30 + t33, t4, pkin(2) * t31 + t24; t11 * t33 - t30, -t11 * t34 - t29, t3, pkin(2) * t35 + t25; t9 * t19, -t9 * t16, -t11, t9 * pkin(2) + t22; 0, 0, 0, 1; t10 * t31 + t18 * t8, t18 * t10 - t8 * t31, t4, pkin(3) * t34 + t7 * t31 + t24; t10 * t35 - t21 * t8, -t21 * t10 - t8 * t35, t3, -pkin(3) * t30 + t7 * t35 + t25; t9 * t10, -t9 * t8, -t11, t9 * t7 + t22; 0, 0, 0, 1;];
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
