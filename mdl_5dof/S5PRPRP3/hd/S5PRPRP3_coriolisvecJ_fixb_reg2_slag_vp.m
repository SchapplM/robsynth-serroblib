% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PRPRP3
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
% Datum: 2019-12-05 15:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRPRP3_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP3_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP3_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP3_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:33:41
% EndTime: 2019-12-05 15:33:44
% DurationCPUTime: 0.60s
% Computational Cost: add. (678->132), mult. (1689->184), div. (0->0), fcn. (1140->6), ass. (0->94)
t53 = sin(pkin(8));
t54 = cos(pkin(8));
t56 = sin(qJ(2));
t58 = cos(qJ(2));
t35 = t53 * t56 - t54 * t58;
t30 = t35 * qJD(2);
t27 = qJD(1) * t30;
t101 = -qJD(3) * qJD(4) + t27;
t57 = cos(qJ(4));
t80 = t58 * qJD(1);
t42 = qJD(2) * pkin(2) + t80;
t88 = qJD(1) * t56;
t25 = t53 * t42 + t54 * t88;
t17 = qJD(2) * pkin(6) + t25;
t68 = qJ(5) * qJD(2) + t17;
t64 = t68 * t57;
t55 = sin(qJ(4));
t82 = t55 * qJD(3);
t10 = t64 + t82;
t36 = t53 * t58 + t54 * t56;
t29 = t36 * qJD(1);
t49 = t57 * qJD(3);
t9 = -t68 * t55 + t49;
t90 = qJD(4) * pkin(4);
t8 = t9 + t90;
t100 = t8 - t9;
t51 = t55 ^ 2;
t99 = pkin(4) * t51;
t98 = t54 * pkin(2);
t97 = t57 * pkin(4);
t28 = t36 * qJD(2);
t26 = qJD(1) * t28;
t96 = t26 * t35;
t60 = qJD(2) ^ 2;
t95 = t57 * t60;
t59 = qJD(4) ^ 2;
t94 = t59 * t55;
t50 = t59 * t57;
t43 = t53 * t88;
t31 = t54 * t80 - t43;
t84 = t29 * qJD(2);
t86 = qJD(4) * t55;
t93 = t31 * t86 + t57 * t84;
t52 = t57 ^ 2;
t92 = t51 - t52;
t91 = t51 + t52;
t46 = t53 * pkin(2) + pkin(6);
t89 = qJ(5) + t46;
t87 = qJD(2) * t57;
t85 = qJD(4) * t57;
t83 = t30 * qJD(2);
t81 = t57 * qJD(5);
t24 = t54 * t42 - t43;
t75 = -pkin(3) - t97;
t13 = t75 * qJD(2) + qJD(5) - t24;
t79 = -qJD(5) - t13;
t78 = qJD(2) * qJD(4);
t76 = t101 * t57 + t17 * t86;
t74 = qJ(5) * t86;
t73 = t91 * t31;
t72 = t55 * t78;
t71 = t57 * t78;
t70 = -t46 * t59 - t26;
t69 = qJD(4) * t89;
t67 = t55 * t71;
t66 = t10 * t57 - t55 * t8;
t11 = -t55 * t17 + t49;
t12 = t57 * t17 + t82;
t65 = t11 * t55 - t12 * t57;
t16 = -qJD(2) * pkin(3) - t24;
t47 = -pkin(3) - t98;
t63 = qJD(4) * (qJD(2) * t47 + t16);
t3 = (-t74 + t81) * qJD(2) - t76;
t4 = (-qJD(2) * qJD(5) + t27) * t55 - t10 * qJD(4);
t62 = t3 * t57 - t4 * t55 + (-t10 * t55 - t57 * t8) * qJD(4);
t6 = -t12 * qJD(4) + t55 * t27;
t61 = -t76 * t57 - t6 * t55 + (-t11 * t57 - t12 * t55) * qJD(4);
t45 = pkin(4) * t72;
t44 = t55 * t95;
t41 = -0.2e1 * t67;
t40 = 0.2e1 * t67;
t39 = t92 * t60;
t38 = t75 - t98;
t34 = t89 * t57;
t33 = t89 * t55;
t32 = -0.2e1 * t92 * t78;
t22 = t31 * t85;
t19 = -t55 * qJD(5) - t57 * t69;
t18 = -t55 * t69 + t81;
t15 = t45 + t26;
t7 = t91 * t83;
t2 = t30 * t86 - t36 * t50 + (-t28 * t57 + t35 * t86) * qJD(2);
t1 = t30 * t85 + t36 * t94 + (t28 * t55 + t35 * t85) * qJD(2);
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60 * t56, -t60 * t58, 0, 0, 0, 0, 0, 0, 0, 0, -t28 * qJD(2), t83, 0, -t24 * t28 - t25 * t30 - t27 * t36 + t96, 0, 0, 0, 0, 0, 0, t2, t1, -t7, t16 * t28 + t65 * t30 + t61 * t36 + t96, 0, 0, 0, 0, 0, 0, t2, t1, -t7, t13 * t28 + t15 * t35 - t66 * t30 + t62 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (t35 * qJD(1) + t31) * qJD(2), 0, t24 * t29 - t25 * t31 + (-t26 * t54 - t27 * t53) * pkin(2), t40, t32, t50, t41, -t94, 0, t55 * t63 + t70 * t57 + t93, t22 + t57 * t63 + (-t70 - t84) * t55, -qJD(2) * t73 + t61, -t16 * t29 + t26 * t47 + t65 * t31 + t61 * t46, t40, t32, t50, t41, -t94, 0, -t15 * t57 + (t19 + (t13 + (t38 - t97) * qJD(2)) * t55) * qJD(4) + t93, t22 + (t15 - t84) * t55 + (t13 * t57 - t18 + (t38 * t57 + t99) * qJD(2)) * qJD(4), (t18 * t57 - t19 * t55 - t73 + (t33 * t57 - t34 * t55) * qJD(4)) * qJD(2) + t62, t15 * t38 + t3 * t34 - t4 * t33 + (t31 * t55 + t19) * t8 + (pkin(4) * t86 - t29) * t13 + (-t31 * t57 + t18) * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t94, -t50, 0, -t65 * qJD(4) - t55 * t76 + t6 * t57, 0, 0, 0, 0, 0, 0, -t94, -t50, 0, t66 * qJD(4) + t3 * t55 + t4 * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, t39, 0, t44, 0, 0, (-qJD(2) * t16 + t27) * t55, t11 * qJD(4) - t16 * t87 + t76, 0, 0, -t44, t39, 0, t44, 0, 0, (t10 - t64) * qJD(4) + (pkin(4) * t95 + t79 * qJD(2) + t101) * t55, -t60 * t99 + t9 * qJD(4) + (t79 * t57 + t74) * qJD(2) + t76, (-t90 + t100) * t87, t100 * t10 + (-t13 * t55 * qJD(2) + t4) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t72, 0.2e1 * t71, -t91 * t60, t45 + (t29 - t66) * qJD(2);];
tauc_reg = t5;
