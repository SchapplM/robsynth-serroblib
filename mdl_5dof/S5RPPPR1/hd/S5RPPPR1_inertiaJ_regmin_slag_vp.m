% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x19]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPPPR1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t29 = cos(pkin(8));
t42 = -0.2e1 * t29;
t30 = cos(pkin(7));
t20 = -t30 * pkin(1) - pkin(2);
t26 = sin(pkin(8));
t12 = -t29 * pkin(3) - t26 * qJ(4) + t20;
t27 = sin(pkin(7));
t19 = t27 * pkin(1) + qJ(3);
t25 = sin(pkin(9));
t28 = cos(pkin(9));
t37 = t28 * t29;
t6 = t25 * t12 + t19 * t37;
t22 = t26 ^ 2;
t41 = t22 * t19;
t40 = t25 * t26;
t39 = t25 * t29;
t38 = t28 * t26;
t36 = t25 ^ 2 + t28 ^ 2;
t24 = t29 ^ 2;
t35 = t22 + t24;
t10 = t28 * t12;
t5 = -t19 * t39 + t10;
t34 = t6 * t25 + t5 * t28;
t31 = sin(qJ(5));
t32 = cos(qJ(5));
t14 = t32 * t25 + t31 * t28;
t13 = -t31 * t25 + t32 * t28;
t18 = t19 ^ 2;
t17 = t26 * t19;
t16 = t22 * t18;
t11 = pkin(4) * t40 + t17;
t8 = t13 * t26;
t7 = t14 * t26;
t4 = -pkin(6) * t40 + t6;
t3 = -pkin(6) * t38 + t10 + (-t19 * t25 - pkin(4)) * t29;
t2 = t31 * t3 + t32 * t4;
t1 = t32 * t3 - t31 * t4;
t9 = [1, 0, 0, (t27 ^ 2 + t30 ^ 2) * pkin(1) ^ 2, t20 * t42, 0.2e1 * t20 * t26, 0.2e1 * t35 * t19, t24 * t18 + t20 ^ 2 + t16, 0.2e1 * t25 * t41 - 0.2e1 * t5 * t29, 0.2e1 * t28 * t41 + 0.2e1 * t6 * t29, -0.2e1 * t34 * t26, t5 ^ 2 + t6 ^ 2 + t16, t8 ^ 2, -0.2e1 * t8 * t7, t8 * t42, 0.2e1 * t29 * t7, t24, -0.2e1 * t1 * t29 + 0.2e1 * t11 * t7, 0.2e1 * t11 * t8 + 0.2e1 * t2 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t19 * t29 - t25 * t5 + t28 * t6) * t26, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 1, 0, 0, 0, t35, 0, 0, 0, t36 * t22 + t24, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -t29, t26, 0, t20, -t37, t39, -t36 * t26, t34, 0, 0, 0, 0, 0, -t13 * t29, t14 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, t36, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t40, t38, 0, t17, 0, 0, 0, 0, 0, t7, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t7, -t29, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t9;
