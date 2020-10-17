% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x24]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRPRRP3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 23:44:44
% EndTime: 2019-05-04 23:44:46
% DurationCPUTime: 0.52s
% Computational Cost: add. (458->81), mult. (990->154), div. (0->0), fcn. (1207->10), ass. (0->57)
t31 = sin(pkin(11));
t33 = cos(pkin(11));
t36 = sin(qJ(4));
t59 = cos(qJ(4));
t19 = t59 * t31 + t36 * t33;
t62 = -0.2e1 * t19;
t25 = -t33 * pkin(3) - pkin(2);
t61 = 0.2e1 * t25;
t38 = cos(qJ(5));
t60 = t38 * pkin(5);
t32 = sin(pkin(6));
t58 = t32 * sin(qJ(2));
t39 = cos(qJ(2));
t57 = t32 * t39;
t18 = t36 * t31 - t59 * t33;
t35 = sin(qJ(5));
t56 = t35 * t18;
t55 = t35 * t19;
t54 = t35 * t38;
t51 = pkin(8) + qJ(3);
t20 = t51 * t31;
t21 = t51 * t33;
t13 = -t36 * t20 + t59 * t21;
t53 = t38 * t13;
t52 = t38 * t19;
t50 = -qJ(6) - pkin(9);
t49 = t31 ^ 2 + t33 ^ 2;
t29 = t35 ^ 2;
t30 = t38 ^ 2;
t48 = t29 + t30;
t47 = qJ(6) * t19;
t46 = t18 * t62;
t11 = t18 * pkin(4) - t19 * pkin(9) + t25;
t3 = t38 * t11 - t35 * t13;
t45 = -pkin(4) * t19 - pkin(9) * t18;
t1 = t18 * pkin(5) - t38 * t47 + t3;
t2 = t53 + (t11 - t47) * t35;
t44 = t1 * t38 + t2 * t35;
t34 = cos(pkin(6));
t14 = -t31 * t58 + t34 * t33;
t15 = t34 * t31 + t33 * t58;
t9 = t36 * t14 + t59 * t15;
t5 = -t35 * t9 - t38 * t57;
t6 = -t35 * t57 + t38 * t9;
t43 = t6 * t35 + t5 * t38;
t42 = -t14 * t31 + t15 * t33;
t22 = t50 * t35;
t23 = t50 * t38;
t41 = t38 * t22 - t35 * t23;
t12 = t59 * t20 + t36 * t21;
t26 = -pkin(4) - t60;
t17 = t19 ^ 2;
t16 = t38 * t18;
t8 = -t59 * t14 + t36 * t15;
t7 = pkin(5) * t55 + t12;
t4 = t35 * t11 + t53;
t10 = [1, 0, 0, 0, 0, 0, 0, t32 ^ 2 * t39 ^ 2 + t14 ^ 2 + t15 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5 ^ 2 + t6 ^ 2 + t8 ^ 2; 0, 0, t57, -t58, t33 * t57, -t31 * t57, t42, pkin(2) * t57 + t42 * qJ(3), 0, 0, 0, 0, 0, -t18 * t57, -t19 * t57, 0, 0, 0, 0, 0, t5 * t18 + t8 * t55, -t6 * t18 + t8 * t52, -t43 * t19, t5 * t1 + t6 * t2 + t8 * t7; 0, 1, 0, 0, 0.2e1 * pkin(2) * t33, -0.2e1 * pkin(2) * t31, 0.2e1 * t49 * qJ(3), t49 * qJ(3) ^ 2 + pkin(2) ^ 2, t17, t46, 0, 0, 0, t18 * t61, t19 * t61, t30 * t17, -0.2e1 * t17 * t54, 0.2e1 * t18 * t52, t35 * t46, t18 ^ 2, 0.2e1 * t12 * t55 + 0.2e1 * t3 * t18, 0.2e1 * t12 * t52 - 0.2e1 * t4 * t18, t44 * t62, t1 ^ 2 + t2 ^ 2 + t7 ^ 2; 0, 0, 0, 0, 0, 0, 0, -t57, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43; 0, 0, 0, 0, -t33, t31, 0, -pkin(2), 0, 0, 0, 0, 0, t18, t19, 0, 0, 0, 0, 0, t16, -t56, -t48 * t19, t44; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, -t9, 0, 0, 0, 0, 0, -t8 * t38, t8 * t35, -t5 * t35 + t6 * t38, t5 * t22 - t6 * t23 + t8 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t18, 0, -t12, -t13, t35 * t52 (-t29 + t30) * t19, t56, t16, 0, -t12 * t38 + t45 * t35, t12 * t35 + t45 * t38, -t1 * t35 - t41 * t19 + t2 * t38, t1 * t22 - t2 * t23 + t7 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t29, 0.2e1 * t54, 0, 0, 0, 0.2e1 * pkin(4) * t38, -0.2e1 * pkin(4) * t35, -0.2e1 * t22 * t35 - 0.2e1 * t23 * t38, t22 ^ 2 + t23 ^ 2 + t26 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, -t6, 0, t5 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, -t55, t18, t3, -t4, -pkin(5) * t52, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t35, 0, t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, t38, 0, -t35 * pkin(9), -t38 * pkin(9), -t35 * pkin(5), t22 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t10;
