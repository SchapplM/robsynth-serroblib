% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPPRR9_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR9_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR9_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR9_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:02:48
% EndTime: 2019-12-31 18:02:52
% DurationCPUTime: 1.32s
% Computational Cost: add. (1819->215), mult. (3549->328), div. (0->0), fcn. (2001->6), ass. (0->129)
t102 = qJ(2) * qJD(1);
t60 = -pkin(1) - pkin(2);
t42 = t60 * qJD(1) + qJD(2);
t53 = sin(pkin(8));
t54 = cos(pkin(8));
t30 = t54 * t102 + t53 * t42;
t24 = -qJD(1) * pkin(6) + t30;
t57 = sin(qJ(4));
t59 = cos(qJ(4));
t16 = t59 * qJD(3) - t57 * t24;
t101 = qJD(1) * qJD(2);
t90 = t54 * t101;
t11 = t16 * qJD(4) + t59 * t90;
t81 = -pkin(4) * t57 + pkin(7) * t59;
t31 = t53 * qJD(2) + t81 * qJD(4);
t25 = t31 * qJD(1);
t56 = sin(qJ(5));
t58 = cos(qJ(5));
t17 = t57 * qJD(3) + t59 * t24;
t14 = qJD(4) * pkin(7) + t17;
t29 = -t53 * t102 + t54 * t42;
t23 = qJD(1) * pkin(3) - t29;
t82 = t59 * pkin(4) + t57 * pkin(7);
t15 = t82 * qJD(1) + t23;
t75 = t56 * t14 - t58 * t15;
t1 = -t75 * qJD(5) + t58 * t11 + t56 * t25;
t103 = t59 * qJD(1);
t43 = qJD(5) + t103;
t151 = t75 * t43 + t1;
t51 = t57 ^ 2;
t150 = (qJD(1) * t51 - t43 * t59) * t56;
t73 = t16 * t57 - t17 * t59;
t149 = t73 * t54;
t6 = t58 * t14 + t56 * t15;
t2 = -qJD(5) * t6 - t56 * t11 + t58 * t25;
t148 = -t6 * t43 - t2;
t12 = t17 * qJD(4) + t57 * t90;
t144 = t12 * t57;
t63 = -(t16 * t59 + t17 * t57) * qJD(4) + t11 * t59 + t144;
t145 = t12 * t56;
t143 = t12 * t58;
t13 = -qJD(4) * pkin(4) - t16;
t142 = t13 * t56;
t141 = t13 * t58;
t104 = t58 * qJD(4);
t112 = qJD(1) * t57;
t34 = t56 * t112 + t104;
t108 = qJD(5) * t34;
t100 = qJD(1) * qJD(4);
t91 = t58 * t100;
t21 = -t59 * t91 + t108;
t140 = t21 * t56;
t139 = t21 * t59;
t105 = t56 * qJD(4);
t106 = qJD(5) * t58;
t67 = t59 * t105 + t57 * t106;
t22 = t67 * qJD(1) - qJD(5) * t105;
t138 = t22 * t58;
t137 = t22 * t59;
t35 = t58 * t112 - t105;
t136 = t34 * t35;
t135 = t34 * t43;
t134 = t34 * t56;
t133 = t34 * t57;
t132 = t34 * t58;
t131 = t35 * t43;
t130 = t35 * t56;
t129 = t35 * t57;
t128 = t35 * t58;
t127 = t43 * t56;
t126 = t43 * t58;
t125 = t53 * t57;
t62 = qJD(1) ^ 2;
t124 = t54 * t62;
t123 = t56 * t59;
t122 = t58 * t59;
t61 = qJD(4) ^ 2;
t121 = t61 * t57;
t120 = t61 * t59;
t32 = -t53 * t123 - t58 * t54;
t110 = qJD(4) * t57;
t97 = t53 * t110;
t119 = t32 * qJD(5) - t58 * t97 - (t54 * t122 + t53 * t56) * qJD(1);
t33 = t53 * t122 - t56 * t54;
t118 = -t33 * qJD(5) + t56 * t97 - (-t54 * t123 + t53 * t58) * qJD(1);
t117 = t54 * qJ(2) + t53 * t60;
t52 = t59 ^ 2;
t116 = t51 - t52;
t115 = t51 + t52;
t114 = t61 + t62;
t88 = -t53 * qJ(2) + t54 * t60;
t36 = pkin(3) - t88;
t113 = qJD(1) * t36;
t111 = qJD(2) * t54;
t109 = qJD(4) * t59;
t107 = qJD(5) * t56;
t99 = t43 * t122;
t98 = t57 * t62 * t59;
t96 = t53 * t109;
t95 = t57 * t107;
t94 = t43 * t106;
t93 = 0.2e1 * t101;
t92 = 0.2e1 * t100;
t89 = t57 * t100;
t87 = -t23 - t111;
t86 = t57 * t94;
t85 = t54 * t92;
t84 = t53 * t93;
t83 = t59 * t89;
t37 = -pkin(6) + t117;
t80 = -qJD(5) * t37 * t59 + t31;
t79 = t56 * t6 - t58 * t75;
t78 = t56 * t75 + t58 * t6;
t72 = t29 * t53 - t30 * t54;
t70 = qJD(1) * (t23 - t111);
t69 = t43 * t95 + t51 * t91;
t68 = -t37 * t61 + t84;
t66 = qJD(4) * (t87 - t113);
t26 = t36 + t82;
t65 = qJD(5) * t26 - t37 * t110 + t59 * t111;
t64 = -t79 * qJD(5) + t1 * t58 - t2 * t56;
t38 = t81 * qJD(1);
t10 = t37 * t122 + t56 * t26;
t9 = -t37 * t123 + t58 * t26;
t8 = t58 * t16 + t56 * t38;
t7 = -t56 * t16 + t58 * t38;
t4 = -t65 * t56 + t80 * t58;
t3 = t80 * t56 + t65 * t58;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t93, qJ(2) * t93, 0, 0, 0, 0, 0, 0, t84, 0.2e1 * t90, 0, ((t54 * t117 - t53 * t88) * qJD(1) - t72) * qJD(2), 0.2e1 * t83, -t116 * t92, -t120, -0.2e1 * t83, t121, 0, t57 * t66 + t68 * t59, -t68 * t57 + t59 * t66, -t115 * t90 - t63, t63 * t37 + (-t149 + (t23 + t113) * t53) * qJD(2), -t21 * t58 * t57 + (t59 * t104 - t95) * t35, (-t130 - t132) * t109 + (t140 - t138 + (-t128 + t134) * qJD(5)) * t57, t139 + (-t99 + t129) * qJD(4) + t69, t22 * t56 * t57 + t67 * t34, t86 + t137 + (-t133 - t150) * qJD(4), (-t43 - t103) * t110, t4 * t43 + (t2 + (-t34 * t37 - t142) * qJD(4)) * t59 + (-t34 * t111 - t13 * t106 - t145 - t22 * t37 + (-qJD(1) * t9 + t75) * qJD(4)) * t57, -t3 * t43 + (-t1 + (-t35 * t37 - t141) * qJD(4)) * t59 + (-t35 * t111 + t13 * t107 - t143 + t21 * t37 + (qJD(1) * t10 + t6) * qJD(4)) * t57, t10 * t22 - t9 * t21 + t3 * t34 + t4 * t35 + t79 * t109 + (t78 * qJD(5) + t1 * t56 + t2 * t58) * t57, t37 * t144 + t1 * t10 + t2 * t9 + t6 * t3 - t75 * t4 + (t37 * t109 + t57 * t111) * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62, -t62 * qJ(2), 0, 0, 0, 0, 0, 0, -t53 * t62, -t124, 0, t72 * qJD(1), 0, 0, 0, 0, 0, 0, -t114 * t59 * t53 + t57 * t85, t114 * t125 + t59 * t85, t115 * t124, qJD(1) * t149 + (t87 * qJD(1) + t63) * t53, 0, 0, 0, 0, 0, 0, -t34 * t96 + t118 * t43 + (-t22 * t53 + (-qJD(4) * t32 + t34 * t54) * qJD(1)) * t57, -t35 * t96 - t119 * t43 + (t21 * t53 + (qJD(4) * t33 + t35 * t54) * qJD(1)) * t57, t118 * t35 + t119 * t34 - t32 * t21 + t33 * t22, t12 * t125 + t1 * t33 + t2 * t32 + t119 * t6 - t118 * t75 + (-t54 * t112 + t96) * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t121, -t120, 0, -t73 * qJD(4) + t11 * t57 - t12 * t59, 0, 0, 0, 0, 0, 0, -t86 + t137 + (-t133 + t150) * qJD(4), -t139 + (-t99 - t129) * qJD(4) + t69, (-t130 + t132) * t109 + (t140 + t138 + (-t128 - t134) * qJD(5)) * t57, (t78 * qJD(4) - t12) * t59 + (qJD(4) * t13 + t64) * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t98, t116 * t62, 0, t98, 0, 0, t57 * t70, t59 * t70, 0, 0, -t35 * t126 + t140, (t21 + t135) * t58 + (t22 + t131) * t56, t94 + (t99 + (-t35 - t105) * t57) * qJD(1), -t34 * t127 + t138, -t43 * t107 + (-t43 * t123 + (t34 - t104) * t57) * qJD(1), t43 * t112, pkin(4) * t22 - t143 + t17 * t34 - t7 * t43 + (-pkin(7) * t126 + t142) * qJD(5) + (-t75 * t57 + (pkin(7) * t110 + t13 * t59) * t56) * qJD(1), -pkin(4) * t21 + t145 + t17 * t35 + t8 * t43 + (pkin(7) * t127 + t141) * qJD(5) + (t13 * t122 + (pkin(7) * t104 - t6) * t57) * qJD(1), -t8 * t34 - t7 * t35 + ((-qJD(5) * t35 + t22) * pkin(7) + t151) * t58 + ((t21 - t108) * pkin(7) + t148) * t56, -t12 * pkin(4) + pkin(7) * t64 - t13 * t17 - t6 * t8 + t7 * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t136, -t34 ^ 2 + t35 ^ 2, t21 - t135, -t136, t22 - t131, -t89, t13 * t35 - t148, -t13 * t34 - t151, 0, 0;];
tauc_reg = t5;
