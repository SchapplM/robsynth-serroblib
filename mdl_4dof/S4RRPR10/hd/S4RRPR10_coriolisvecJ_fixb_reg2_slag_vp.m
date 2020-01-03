% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4RRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% 
% Output:
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RRPR10_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR10_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR10_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR10_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:11:56
% EndTime: 2019-12-31 17:11:59
% DurationCPUTime: 0.95s
% Computational Cost: add. (978->182), mult. (2350->275), div. (0->0), fcn. (1248->4), ass. (0->115)
t64 = cos(qJ(2));
t105 = qJD(1) * t64;
t61 = sin(qJ(4));
t63 = cos(qJ(4));
t96 = t63 * qJD(2);
t30 = -t61 * t105 + t96;
t93 = qJD(1) * qJD(2);
t128 = -0.2e1 * t93;
t62 = sin(qJ(2));
t86 = t62 * t93;
t52 = pkin(2) * t86;
t79 = pkin(6) * t62 - qJ(3) * t64;
t97 = t62 * qJD(3);
t69 = t79 * qJD(2) - t97;
t11 = t69 * qJD(1) + t52;
t85 = t64 * t93;
t51 = pkin(5) * t85;
t26 = pkin(3) * t85 + t51;
t65 = -pkin(2) - pkin(6);
t84 = -t62 * qJ(3) - pkin(1);
t27 = t65 * t64 + t84;
t16 = t27 * qJD(1);
t98 = t62 * qJD(1);
t54 = pkin(5) * t98;
t95 = pkin(3) * t98 + qJD(3) + t54;
t18 = t65 * qJD(2) + t95;
t77 = t61 * t16 - t63 * t18;
t1 = -t77 * qJD(4) + t63 * t11 + t61 * t26;
t53 = qJD(4) + t98;
t127 = t53 * t77 + t1;
t99 = t61 * qJD(2);
t28 = t63 * t105 + t99;
t12 = t28 * qJD(4) - t61 * t86;
t73 = t28 * t53;
t126 = t12 - t73;
t6 = t63 * t16 + t61 * t18;
t2 = -qJD(4) * t6 - t61 * t11 + t63 * t26;
t125 = t6 * t53 + t2;
t122 = pkin(3) + pkin(5);
t124 = t122 * t62;
t115 = t30 * t53;
t13 = t30 * qJD(4) - t63 * t86;
t123 = -t13 + t115;
t119 = t13 * t61;
t19 = (-qJD(1) * t124 + qJD(3)) * qJD(2);
t118 = t19 * t61;
t117 = t19 * t63;
t116 = t30 * t28;
t114 = t30 * t64;
t113 = t53 * t62;
t112 = t53 * t65;
t111 = t63 * t12;
t67 = qJD(1) ^ 2;
t110 = t64 * t67;
t66 = qJD(2) ^ 2;
t109 = t66 * t62;
t108 = t66 * t64;
t59 = t62 ^ 2;
t60 = t64 ^ 2;
t107 = t59 - t60;
t106 = qJD(2) * pkin(2);
t38 = -t64 * pkin(2) + t84;
t25 = qJD(1) * t38;
t104 = qJD(2) * t62;
t103 = qJD(2) * t64;
t102 = qJD(4) * t61;
t101 = qJD(4) * t63;
t100 = qJD(4) * t64;
t94 = qJD(2) * qJ(3);
t92 = t63 * t113;
t55 = pkin(5) * t105;
t44 = t122 * t64;
t90 = t61 * t100;
t89 = t53 * t101;
t88 = t63 * t100;
t83 = pkin(1) * t128;
t82 = qJD(3) - t106;
t40 = -t55 - t94;
t81 = t62 * t85;
t80 = t6 * t63 + t61 * t77;
t10 = t124 * t61 + t63 * t27;
t9 = t124 * t63 - t61 * t27;
t76 = -qJD(1) * t60 + t113;
t75 = -0.2e1 * qJD(2) * t25;
t74 = t53 * t61;
t70 = -t64 * t94 - t97;
t15 = t70 * qJD(1) + t52;
t57 = pkin(2) * t104;
t22 = t57 + t70;
t72 = pkin(5) * t66 + qJD(1) * t22 + t15;
t56 = pkin(3) * t105;
t24 = -t40 + t56;
t71 = t65 * t103 + t24 * t62;
t36 = (-qJD(3) + t54) * qJD(2);
t37 = t54 + t82;
t68 = -t36 * t64 + (t37 * t64 + (t40 + t55) * t62) * qJD(2);
t58 = pkin(2) * t98;
t50 = t62 * t110;
t46 = t63 * t85;
t42 = -0.2e1 * t81;
t41 = 0.2e1 * t81;
t39 = t107 * t67;
t35 = qJD(2) * t44;
t34 = t55 + t56;
t33 = qJD(2) * t124;
t31 = -qJ(3) * t105 + t58;
t23 = t107 * t128;
t21 = t79 * qJD(1) + t58;
t17 = t25 * t98;
t14 = t57 + t69;
t8 = t63 * t21 + t61 * t34;
t7 = -t61 * t21 + t63 * t34;
t4 = -t10 * qJD(4) - t61 * t14 + t63 * t35;
t3 = t9 * qJD(4) + t63 * t14 + t61 * t35;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, t23, t108, t42, -t109, 0, -pkin(5) * t108 + t62 * t83, pkin(5) * t109 + t64 * t83, 0, 0, 0, -t108, t109, t41, t23, t42, t68, t62 * t75 + t72 * t64, -t72 * t62 + t64 * t75, t68 * pkin(5) + t15 * t38 + t25 * t22, t12 * t61 * t64 + (t62 * t99 - t88) * t30, (-t28 * t61 + t30 * t63) * t104 + (t111 + t119 + (t28 * t63 + t30 * t61) * qJD(4)) * t64, -t53 * t88 - t12 * t62 + (t76 * t61 + t114) * qJD(2), t13 * t63 * t64 + (-t62 * t96 - t90) * t28, t53 * t90 - t13 * t62 + (-t28 * t64 + t76 * t63) * qJD(2), (t53 + t98) * t103, t44 * t13 - t33 * t28 + t4 * t53 + (-t24 * t96 + t2) * t62 + (-t24 * t102 + t117 + (qJD(1) * t9 - t77) * qJD(2)) * t64, -t44 * t12 - t3 * t53 - t33 * t30 + (t24 * t99 - t1) * t62 + (-t24 * t101 - t118 + (-qJD(1) * t10 - t6) * qJD(2)) * t64, -t10 * t13 + t9 * t12 - t3 * t28 - t4 * t30 + t80 * t104 + (-t1 * t63 + t2 * t61 + (t6 * t61 - t63 * t77) * qJD(4)) * t64, t1 * t10 + t19 * t44 + t2 * t9 - t24 * t33 + t6 * t3 - t4 * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, t39, 0, t50, 0, 0, t67 * pkin(1) * t62, pkin(1) * t110, 0, 0, 0, 0, 0, -t50, t39, t50, ((-t40 - t94) * t62 + (-t37 + t82) * t64) * qJD(1), -t31 * t105 + t17, 0.2e1 * qJD(2) * qJD(3) + (t25 * t64 + t31 * t62) * qJD(1), -t36 * qJ(3) - t40 * qJD(3) - t25 * t31 + (-t40 * t62 + (-t37 - t106) * t64) * qJD(1) * pkin(5), -t30 * t74 - t111, (-t13 - t115) * t63 + (t12 + t73) * t61, -t53 * t102 + t46 + (-t61 * t113 - t114) * qJD(1), t63 * t73 + t119, -t89 + (-t92 + (t28 - t99) * t64) * qJD(1), -t53 * t105, qJ(3) * t13 + t118 - t7 * t53 + t95 * t28 + (-t61 * t112 + t24 * t63) * qJD(4) + (t71 * t63 + t64 * t77) * qJD(1), -qJ(3) * t12 + t117 + t8 * t53 + t95 * t30 + (-t63 * t112 - t24 * t61) * qJD(4) + (t6 * t64 - t71 * t61) * qJD(1), t8 * t28 + t7 * t30 + (-t6 * t98 + t12 * t65 - t2 + (-t28 * t65 - t6) * qJD(4)) * t63 + (-t77 * t98 - t13 * t65 - t1 + (t30 * t65 - t77) * qJD(4)) * t61, t19 * qJ(3) + t77 * t7 - t6 * t8 + t95 * t24 + (t80 * qJD(4) + t1 * t61 + t2 * t63) * t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, -t59 * t67 - t66, t40 * qJD(2) + t17 + t51, 0, 0, 0, 0, 0, 0, -qJD(2) * t28 - t53 * t74 + t46, -t89 - qJD(2) * t30 + (-t64 * t99 - t92) * qJD(1), t123 * t61 + t126 * t63, -t24 * qJD(2) + t125 * t63 + t127 * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t116, -t28 ^ 2 + t30 ^ 2, -t126, -t116, t123, t85, -t24 * t30 + t125, t24 * t28 - t127, 0, 0;];
tauc_reg = t5;
