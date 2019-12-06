% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRPRP4_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP4_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP4_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP4_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:36:06
% EndTime: 2019-12-05 15:36:09
% DurationCPUTime: 0.48s
% Computational Cost: add. (669->117), mult. (1650->157), div. (0->0), fcn. (1119->6), ass. (0->85)
t49 = cos(qJ(4));
t71 = t49 * qJD(3);
t50 = cos(qJ(2));
t70 = t50 * qJD(1);
t36 = qJD(2) * pkin(2) + t70;
t45 = sin(pkin(8));
t46 = cos(pkin(8));
t48 = sin(qJ(2));
t77 = qJD(1) * t48;
t19 = t45 * t36 + t46 * t77;
t14 = qJD(2) * pkin(6) + t19;
t47 = sin(qJ(4));
t86 = t47 * t14;
t11 = t71 - t86;
t92 = qJD(5) - t11;
t8 = -qJD(4) * pkin(4) + t92;
t12 = t47 * qJD(3) + t49 * t14;
t9 = qJD(4) * qJ(5) + t12;
t29 = t45 * t50 + t46 * t48;
t24 = t29 * qJD(1);
t91 = t46 * pkin(2);
t28 = t45 * t48 - t46 * t50;
t25 = t28 * qJD(2);
t21 = qJD(1) * t25;
t85 = t47 * t21;
t5 = t12 * qJD(4) - t85;
t90 = t5 * t47;
t89 = t5 * t49;
t23 = t29 * qJD(2);
t20 = qJD(1) * t23;
t88 = t20 * t28;
t39 = t45 * pkin(2) + pkin(6);
t51 = qJD(4) ^ 2;
t87 = t39 * t51;
t84 = t47 * t49;
t83 = t51 * t47;
t42 = t51 * t49;
t82 = qJD(4) * t71 - t49 * t21;
t37 = t45 * t77;
t26 = t46 * t70 - t37;
t73 = t24 * qJD(2);
t75 = qJD(4) * t47;
t81 = t26 * t75 + t49 * t73;
t58 = pkin(4) * t47 - qJ(5) * t49;
t22 = t58 * qJD(4) - t47 * qJD(5);
t80 = t22 - t24;
t43 = t47 ^ 2;
t44 = t49 ^ 2;
t79 = -t43 + t44;
t78 = t43 + t44;
t18 = t46 * t36 - t37;
t13 = -qJD(2) * pkin(3) - t18;
t76 = qJD(2) * t13;
t74 = qJD(4) * t49;
t72 = t25 * qJD(2);
t68 = qJD(2) * qJD(4);
t7 = (t24 + t22) * qJD(2);
t67 = -t7 - t87;
t66 = t11 + t86;
t65 = -t20 - t87;
t56 = -t49 * pkin(4) - t47 * qJ(5) - pkin(3);
t10 = t56 * qJD(2) - t18;
t27 = t56 - t91;
t64 = qJD(2) * t27 + t10;
t40 = -pkin(3) - t91;
t63 = qJD(2) * t40 + t13;
t61 = t68 * t84;
t60 = t78 * t26 * qJD(2);
t59 = t47 * t8 + t49 * t9;
t57 = t11 * t47 - t12 * t49;
t3 = (qJD(5) - t86) * qJD(4) + t82;
t54 = t3 * t49 + t90 + (-t47 * t9 + t49 * t8) * qJD(4);
t4 = -t14 * t75 + t82;
t53 = t4 * t49 + t90 + (-t11 * t49 - t12 * t47) * qJD(4);
t52 = qJD(2) ^ 2;
t38 = t52 * t84;
t35 = -0.2e1 * t61;
t34 = 0.2e1 * t61;
t33 = t79 * t52;
t31 = t58 * qJD(2);
t30 = t79 * t68;
t6 = t78 * t72;
t2 = t25 * t75 - t29 * t42 + (-t23 * t49 + t28 * t75) * qJD(2);
t1 = t25 * t74 + t29 * t83 + (t23 * t47 + t28 * t74) * qJD(2);
t15 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52 * t48, -t52 * t50, 0, 0, 0, 0, 0, 0, 0, 0, -t23 * qJD(2), t72, 0, -t18 * t23 - t19 * t25 - t21 * t29 + t88, 0, 0, 0, 0, 0, 0, t2, t1, -t6, t13 * t23 + t57 * t25 + t53 * t29 + t88, 0, 0, 0, 0, 0, 0, t2, -t6, -t1, t10 * t23 - t59 * t25 + t7 * t28 + t54 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (t28 * qJD(1) + t26) * qJD(2), 0, t18 * t24 - t19 * t26 + (-t20 * t46 - t21 * t45) * pkin(2), t34, 0.2e1 * t30, t42, t35, -t83, 0, t65 * t49 + t63 * t75 + t81, (-t65 - t73) * t47 + (t26 + t63) * t74, -t60 + t53, -t13 * t24 + t20 * t40 + t57 * t26 + t53 * t39, t34, t42, -0.2e1 * t30, 0, t83, t35, t64 * t75 + (-qJD(2) * t22 + t67) * t49 + t81, -t60 + t54, (-t26 - t64) * t74 + (-t80 * qJD(2) + t67) * t47, t80 * t10 - t59 * t26 + t7 * t27 + t54 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t83, -t42, 0, -t57 * qJD(4) + t4 * t47 - t89, 0, 0, 0, 0, 0, 0, -t83, 0, t42, t59 * qJD(4) + t3 * t47 - t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, -t33, 0, t38, 0, 0, (t21 - t76) * t47, t66 * qJD(4) - t49 * t76 - t82, 0, 0, -t38, 0, t33, 0, 0, t38, t85 + (-t10 * t47 + t31 * t49) * qJD(2), 0, (t10 * t49 + t31 * t47) * qJD(2) + (0.2e1 * qJD(5) - t66) * qJD(4) + t82, -t5 * pkin(4) + t3 * qJ(5) - t10 * t31 - t8 * t12 + t92 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, 0, -t43 * t52 - t51, (qJD(2) * t10 - t21) * t47 + (t12 - t9) * qJD(4);];
tauc_reg = t15;
