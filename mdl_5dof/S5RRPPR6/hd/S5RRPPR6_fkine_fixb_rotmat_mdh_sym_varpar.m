% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2019-12-29 18:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RRPPR6_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR6_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR6_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-29 18:18:35
% EndTime: 2019-12-29 18:18:35
% DurationCPUTime: 0.23s
% Computational Cost: add. (121->48), mult. (92->54), div. (0->0), fcn. (138->10), ass. (0->32)
t14 = qJ(2) + pkin(8);
t11 = cos(t14);
t21 = sin(qJ(1));
t35 = t21 * t11;
t16 = sin(pkin(9));
t34 = t21 * t16;
t17 = cos(pkin(9));
t33 = t21 * t17;
t23 = cos(qJ(1));
t32 = t23 * t11;
t31 = t23 * t16;
t30 = t23 * t17;
t15 = pkin(5) + 0;
t22 = cos(qJ(2));
t7 = t22 * pkin(2) + pkin(1);
t29 = t23 * t7 + 0;
t18 = -qJ(3) - pkin(6);
t28 = t23 * t18 + t21 * t7 + 0;
t20 = sin(qJ(2));
t27 = t20 * pkin(2) + t15;
t19 = -pkin(7) - qJ(4);
t5 = t17 * pkin(4) + pkin(3);
t9 = sin(t14);
t26 = t11 * t5 - t19 * t9;
t25 = pkin(3) * t11 + qJ(4) * t9;
t24 = -t21 * t18 + t29;
t13 = pkin(9) + qJ(5);
t10 = cos(t13);
t8 = sin(t13);
t4 = t23 * t9;
t3 = t21 * t9;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t23, -t21, 0, 0; t21, t23, 0, 0; 0, 0, 1, t15; 0, 0, 0, 1; t23 * t22, -t23 * t20, t21, t23 * pkin(1) + t21 * pkin(6) + 0; t21 * t22, -t21 * t20, -t23, t21 * pkin(1) - t23 * pkin(6) + 0; t20, t22, 0, t15; 0, 0, 0, 1; t32, -t4, t21, t24; t35, -t3, -t23, t28; t9, t11, 0, t27; 0, 0, 0, 1; t11 * t30 + t34, -t11 * t31 + t33, t4, t25 * t23 + t24; t11 * t33 - t31, -t11 * t34 - t30, t3, t25 * t21 + t28; t9 * t17, -t9 * t16, -t11, t9 * pkin(3) - t11 * qJ(4) + t27; 0, 0, 0, 1; t10 * t32 + t21 * t8, t21 * t10 - t8 * t32, t4, t26 * t23 + (pkin(4) * t16 - t18) * t21 + t29; t10 * t35 - t23 * t8, -t23 * t10 - t8 * t35, t3, -pkin(4) * t31 + t26 * t21 + t28; t9 * t10, -t9 * t8, -t11, t11 * t19 + t9 * t5 + t27; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
