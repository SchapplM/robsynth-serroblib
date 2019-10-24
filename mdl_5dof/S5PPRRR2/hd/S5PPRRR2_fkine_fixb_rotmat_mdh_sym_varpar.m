% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-10-24 10:20
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5PPRRR2_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR2_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR2_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-10-24 10:20:52
% EndTime: 2019-10-24 10:20:52
% DurationCPUTime: 0.11s
% Computational Cost: add. (121->48), mult. (92->54), div. (0->0), fcn. (138->10), ass. (0->34)
t15 = qJ(4) + qJ(5);
t10 = sin(t15);
t17 = sin(pkin(8));
t37 = t17 * t10;
t11 = cos(t15);
t36 = t17 * t11;
t21 = sin(qJ(4));
t35 = t17 * t21;
t22 = cos(qJ(4));
t34 = t17 * t22;
t19 = cos(pkin(8));
t33 = t19 * t10;
t32 = t19 * t11;
t31 = t19 * t21;
t30 = t19 * t22;
t18 = cos(pkin(9));
t6 = t18 * pkin(2) + pkin(1);
t29 = t19 * t6 + 0;
t13 = qJ(1) + 0;
t20 = -pkin(5) - qJ(2);
t28 = t17 * t6 + t19 * t20 + 0;
t16 = sin(pkin(9));
t27 = t16 * pkin(2) + t13;
t14 = pkin(9) + qJ(3);
t8 = sin(t14);
t9 = cos(t14);
t26 = pkin(3) * t9 + pkin(6) * t8;
t23 = -pkin(7) - pkin(6);
t7 = t22 * pkin(4) + pkin(3);
t25 = -t23 * t8 + t7 * t9;
t24 = -t17 * t20 + t29;
t4 = t19 * t8;
t3 = t17 * t8;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t19, -t17, 0, 0; t17, t19, 0, 0; 0, 0, 1, t13; 0, 0, 0, 1; t19 * t18, -t19 * t16, t17, t19 * pkin(1) + t17 * qJ(2) + 0; t17 * t18, -t17 * t16, -t19, t17 * pkin(1) - t19 * qJ(2) + 0; t16, t18, 0, t13; 0, 0, 0, 1; t19 * t9, -t4, t17, t24; t17 * t9, -t3, -t19, t28; t8, t9, 0, t27; 0, 0, 0, 1; t9 * t30 + t35, -t9 * t31 + t34, t4, t26 * t19 + t24; t9 * t34 - t31, -t9 * t35 - t30, t3, t26 * t17 + t28; t8 * t22, -t8 * t21, -t9, t8 * pkin(3) - t9 * pkin(6) + t27; 0, 0, 0, 1; t9 * t32 + t37, -t9 * t33 + t36, t4, t25 * t19 + (pkin(4) * t21 - t20) * t17 + t29; t9 * t36 - t33, -t9 * t37 - t32, t3, -pkin(4) * t31 + t25 * t17 + t28; t8 * t11, -t8 * t10, -t9, t9 * t23 + t8 * t7 + t27; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
