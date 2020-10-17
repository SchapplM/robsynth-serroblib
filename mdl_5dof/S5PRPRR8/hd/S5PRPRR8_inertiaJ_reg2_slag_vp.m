% Calculate inertial parameters regressor of joint inertia matrix for
% S5PRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRPRR8_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR8_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR8_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:04:27
% EndTime: 2019-12-05 16:04:30
% DurationCPUTime: 0.59s
% Computational Cost: add. (202->69), mult. (469->125), div. (0->0), fcn. (508->8), ass. (0->55)
t24 = cos(pkin(5));
t26 = sin(qJ(4));
t29 = cos(qJ(4));
t23 = sin(pkin(5));
t30 = cos(qJ(2));
t49 = t23 * t30;
t7 = t24 * t26 + t29 * t49;
t56 = t7 ^ 2;
t28 = cos(qJ(5));
t55 = 0.2e1 * t28;
t54 = 2 * qJ(3);
t53 = t29 * pkin(4);
t52 = t7 * t29;
t21 = t29 ^ 2;
t31 = -pkin(2) - pkin(7);
t51 = t21 * t31;
t27 = sin(qJ(2));
t50 = t23 * t27;
t25 = sin(qJ(5));
t48 = t25 * t26;
t47 = t25 * t28;
t46 = t25 * t29;
t45 = t26 * t31;
t44 = t28 * t26;
t15 = t28 * t29;
t43 = t28 * t31;
t42 = t29 * t26;
t41 = t29 * t31;
t18 = t25 ^ 2;
t20 = t28 ^ 2;
t40 = t18 + t20;
t19 = t26 ^ 2;
t13 = t19 + t21;
t39 = -0.2e1 * t42;
t38 = t25 * t15;
t37 = t40 * t26;
t36 = -pkin(8) * t26 - t53;
t9 = t24 * t29 - t26 * t49;
t2 = -t9 * t25 + t28 * t50;
t3 = t25 * t50 + t9 * t28;
t35 = -t2 * t25 + t3 * t28;
t11 = t26 * pkin(4) - t29 * pkin(8) + qJ(3);
t4 = t28 * t11 - t25 * t45;
t5 = t25 * t11 + t26 * t43;
t34 = -t4 * t25 + t5 * t28;
t1 = t9 * t26 - t52;
t32 = qJ(3) ^ 2;
t22 = t31 ^ 2;
t17 = t23 ^ 2;
t16 = t21 * t22;
t14 = t17 * t27 ^ 2;
t12 = qJ(3) * t50;
t10 = t13 * t31;
t6 = t17 * t30 ^ 2 + t24 ^ 2 + t14;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9 ^ 2 + t14 + t56, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2 ^ 2 + t3 ^ 2 + t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, -t50, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, t50, pkin(2) * t49 + t12, 0, 0, 0, 0, 0, 0, t26 * t50, t29 * t50, -t1, t1 * t31 + t12, 0, 0, 0, 0, 0, 0, t2 * t26 + t7 * t46, t7 * t15 - t3 * t26, (-t2 * t28 - t25 * t3) * t29, t2 * t4 + t3 * t5 - t7 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -0.2e1 * pkin(2), t54, pkin(2) ^ 2 + t32, t21, t39, 0, t19, 0, 0, t26 * t54, t29 * t54, -0.2e1 * t10, t19 * t22 + t16 + t32, t20 * t21, -0.2e1 * t21 * t47, t42 * t55, t18 * t21, t25 * t39, t19, -0.2e1 * t25 * t51 + 0.2e1 * t4 * t26, -0.2e1 * t21 * t43 - 0.2e1 * t5 * t26, 0.2e1 * (-t25 * t5 - t28 * t4) * t29, t4 ^ 2 + t5 ^ 2 + t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26 * t35 - t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(2), 0, 0, 0, 0, 0, 0, 0, 0, -t13, t10, 0, 0, 0, 0, 0, 0, -t13 * t25, -t13 * t28, 0, t26 * t34 + t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40 * t19 + t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, -t9, 0, 0, 0, 0, 0, 0, 0, 0, -t7 * t28, t7 * t25, t35, -t7 * pkin(4) + pkin(8) * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, -t26, 0, t41, -t45, 0, 0, t38, (-t18 + t20) * t29, t48, -t38, t44, 0, t25 * t36 + t28 * t41, -t25 * t41 + t36 * t28, t34, pkin(4) * t41 + pkin(8) * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t26, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t46, t37, pkin(8) * t37 + t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t18, 0.2e1 * t47, 0, t20, 0, 0, pkin(4) * t55, -0.2e1 * pkin(4) * t25, 0.2e1 * t40 * pkin(8), t40 * pkin(8) ^ 2 + pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, -t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, 0, -t46, t26, t4, -t5, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, -t44, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, t28, 0, -t25 * pkin(8), -t28 * pkin(8), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t8;
