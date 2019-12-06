% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRPPR3_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR3_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR3_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR3_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:26:55
% EndTime: 2019-12-05 15:26:56
% DurationCPUTime: 0.32s
% Computational Cost: add. (451->79), mult. (1004->113), div. (0->0), fcn. (681->6), ass. (0->67)
t28 = sin(pkin(8));
t29 = cos(pkin(8));
t31 = sin(qJ(2));
t33 = cos(qJ(2));
t19 = t28 * t33 + t29 * t31;
t15 = t19 * qJD(1);
t53 = t33 * qJD(1);
t21 = qJD(2) * pkin(2) + t53;
t58 = qJD(1) * t31;
t11 = t28 * t21 + t29 * t58;
t7 = qJD(2) * qJ(4) + t11;
t69 = -t7 + t15;
t14 = t19 * qJD(2);
t12 = qJD(1) * t14;
t30 = sin(qJ(5));
t32 = cos(qJ(5));
t22 = t28 * t58;
t10 = t29 * t21 - t22;
t39 = qJD(4) - t10;
t5 = (-pkin(3) - pkin(6)) * qJD(2) + t39;
t3 = -t30 * qJD(3) + t32 * t5;
t1 = t3 * qJD(5) + t30 * t12;
t4 = t32 * qJD(3) + t30 * t5;
t9 = t32 * t12;
t2 = -t4 * qJD(5) + t9;
t36 = -(t3 * t30 - t32 * t4) * qJD(5) + t1 * t30 + t2 * t32;
t66 = t29 * t33;
t67 = t28 * t31;
t18 = -t66 + t67;
t68 = t12 * t18;
t65 = t30 * t32;
t34 = qJD(5) ^ 2;
t64 = t34 * t30;
t63 = t34 * t32;
t26 = t30 ^ 2;
t27 = t32 ^ 2;
t62 = t26 - t27;
t61 = t26 + t27;
t35 = qJD(2) ^ 2;
t60 = -t34 - t35;
t59 = t7 * qJD(2);
t57 = qJD(5) * t30;
t56 = qJD(5) * t32;
t55 = t14 * qJD(2);
t48 = qJD(2) * t66;
t16 = -qJD(2) * t67 + t48;
t54 = t16 * qJD(2);
t49 = t29 * t53;
t17 = -t22 + t49;
t52 = qJD(4) - t17;
t51 = qJD(2) * qJD(5);
t50 = t35 * t65;
t47 = -t29 * pkin(2) - pkin(3);
t46 = t51 * t65;
t25 = t28 * pkin(2) + qJ(4);
t45 = qJD(2) * t25 - t69;
t44 = t17 - t49;
t20 = qJD(2) * t22;
t8 = -t20 + (qJD(4) + t49) * qJD(2);
t42 = t7 * t16 + t8 * t19;
t41 = t3 * t32 + t30 * t4;
t38 = t8 * t25 + t52 * t7;
t24 = -pkin(6) + t47;
t37 = t52 * qJD(2) - t24 * t34 + t8;
t13 = qJD(1) * t48 - t20;
t6 = -qJD(2) * pkin(3) + t39;
t23 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35 * t31, -t35 * t33, 0, 0, 0, 0, 0, 0, 0, 0, -t55, -t54, 0, -t10 * t14 + t11 * t16 + t13 * t19 + t68, 0, 0, 0, 0, 0, 0, 0, t55, t54, t6 * t14 + t42 + t68, 0, 0, 0, 0, 0, 0, t14 * t56 - t18 * t64 + (t16 * t30 + t19 * t56) * qJD(2), -t14 * t57 - t18 * t63 + (t16 * t32 - t19 * t57) * qJD(2), -t61 * t55, t41 * t14 + t36 * t18 + t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44 * qJD(2) + t20, 0, t10 * t15 - t11 * t17 + (-t12 * t29 + t13 * t28) * pkin(2), 0, 0, 0, 0, 0, 0, 0, 0, -t20 + (0.2e1 * qJD(4) - t44) * qJD(2), t12 * t47 - t6 * t15 + t38, -0.2e1 * t46, 0.2e1 * t62 * t51, -t64, 0.2e1 * t46, -t63, 0, t37 * t30 + t45 * t56, t37 * t32 - t45 * t57, t61 * t15 * qJD(2) - t36, -t41 * t15 + t36 * t24 + t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63, t64, 0, -t41 * qJD(5) + t1 * t32 - t2 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, t69 * qJD(2), 0, 0, 0, 0, 0, 0, t60 * t30, t60 * t32, 0, t36 - t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, -t62 * t35, 0, -t50, 0, 0, -t32 * t59 + t9, (-t12 + t59) * t30, 0, 0;];
tauc_reg = t23;
