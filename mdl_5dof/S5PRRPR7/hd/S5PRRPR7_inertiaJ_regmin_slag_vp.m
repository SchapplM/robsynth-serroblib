% Calculate minimal parameter regressor of joint inertia matrix for
% S5PRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRPR7_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR7_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR7_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:37:34
% EndTime: 2019-12-05 16:37:36
% DurationCPUTime: 0.38s
% Computational Cost: add. (251->80), mult. (641->176), div. (0->0), fcn. (735->10), ass. (0->54)
t30 = cos(pkin(10));
t35 = cos(qJ(5));
t36 = cos(qJ(3));
t32 = sin(qJ(5));
t33 = sin(qJ(3));
t47 = t32 * t33;
t17 = t30 * t47 + t35 * t36;
t59 = -0.2e1 * t17;
t28 = sin(pkin(10));
t58 = -0.2e1 * t28;
t57 = 0.2e1 * t36;
t56 = pkin(7) * t28;
t24 = t33 * pkin(7);
t55 = t36 * pkin(7);
t25 = t28 ^ 2;
t54 = t25 * t35;
t23 = t28 * t33;
t29 = sin(pkin(5));
t53 = t29 * sin(qJ(2));
t52 = t29 * cos(qJ(2));
t21 = -t36 * pkin(3) - t33 * qJ(4) - pkin(2);
t51 = t30 * t21;
t50 = t30 * t33;
t49 = t32 * t28;
t48 = t32 * t30;
t46 = t35 * t28;
t45 = t35 * t30;
t13 = t28 * t21 + t30 * t55;
t26 = t30 ^ 2;
t44 = t25 + t26;
t43 = qJ(4) * t30;
t42 = t25 * qJ(4);
t31 = cos(pkin(5));
t16 = t31 * t33 + t36 * t53;
t5 = t16 * t28 + t30 * t52;
t7 = t16 * t30 - t28 * t52;
t41 = t5 * t28 + t7 * t30;
t40 = -pkin(3) * t33 + qJ(4) * t36;
t12 = -t28 * t55 + t51;
t39 = -t12 * t28 + t13 * t30;
t27 = t33 ^ 2;
t20 = -t30 * pkin(4) - t28 * pkin(8) - pkin(3);
t18 = -t32 * t36 + t33 * t45;
t15 = -t31 * t36 + t33 * t53;
t14 = t24 + (pkin(4) * t28 - pkin(8) * t30) * t33;
t11 = t32 * t20 + t35 * t43;
t10 = t35 * t20 - t32 * t43;
t9 = -t36 * pkin(8) + t13;
t8 = -t51 + (pkin(4) + t56) * t36;
t4 = t32 * t14 + t35 * t9;
t3 = t35 * t14 - t32 * t9;
t2 = t15 * t32 + t7 * t35;
t1 = t15 * t35 - t7 * t32;
t6 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15 ^ 2 + t5 ^ 2 + t7 ^ 2, 0, 0, 0, 0, 0, 0, 0; 0, 0, t52, -t53, 0, 0, 0, 0, 0, t36 * t52, -t33 * t52, t15 * t23 + t5 * t36, t15 * t50 + t7 * t36, (-t28 * t7 + t30 * t5) * t33, -t5 * t12 + t7 * t13 + t15 * t24, 0, 0, 0, 0, 0, t1 * t23 + t5 * t17, t5 * t18 - t2 * t23; 0, 1, 0, 0, t27, t33 * t57, 0, 0, 0, pkin(2) * t57, -0.2e1 * pkin(2) * t33, -0.2e1 * t12 * t36 + 0.2e1 * t27 * t56, 0.2e1 * t27 * pkin(7) * t30 + 0.2e1 * t13 * t36, 0.2e1 * (-t12 * t30 - t13 * t28) * t33, t27 * pkin(7) ^ 2 + t12 ^ 2 + t13 ^ 2, t18 ^ 2, t18 * t59, 0.2e1 * t18 * t23, t23 * t59, t25 * t27, 0.2e1 * t8 * t17 + 0.2e1 * t3 * t23, 0.2e1 * t8 * t18 - 0.2e1 * t4 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, -t16, -t15 * t30, t15 * t28, t41, -t15 * pkin(3) + t41 * qJ(4), 0, 0, 0, 0, 0, -t1 * t30 + t5 * t49, t2 * t30 + t5 * t46; 0, 0, 0, 0, 0, 0, t33, t36, 0, -t24, -t55, -pkin(7) * t50 + t40 * t28, pkin(7) * t23 + t40 * t30, t39, -pkin(3) * t24 + t39 * qJ(4), t18 * t46, (-t17 * t35 - t18 * t32) * t28, -t18 * t30 + t33 * t54, t17 * t30 - t25 * t47, -t28 * t50, -t3 * t30 + (qJ(4) * t17 + t10 * t33 + t32 * t8) * t28, t4 * t30 + (qJ(4) * t18 - t11 * t33 + t35 * t8) * t28; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3) * t30, pkin(3) * t58, 0.2e1 * t44 * qJ(4), t44 * qJ(4) ^ 2 + pkin(3) ^ 2, t35 ^ 2 * t25, -0.2e1 * t32 * t54, t45 * t58, 0.2e1 * t28 * t48, t26, -0.2e1 * t10 * t30 + 0.2e1 * t32 * t42, 0.2e1 * t11 * t30 + 0.2e1 * t35 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, t50, 0, t24, 0, 0, 0, 0, 0, t33 * t46, -t28 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, t28, 0, -pkin(3), 0, 0, 0, 0, 0, -t45, t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t17, t23, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, -t49, -t30, t10, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, -t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t6;
