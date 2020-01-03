% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPRPR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRPR16_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR16_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR16_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPR16_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:39:32
% EndTime: 2019-12-31 18:39:37
% DurationCPUTime: 1.42s
% Computational Cost: add. (1435->222), mult. (2862->310), div. (0->0), fcn. (1444->4), ass. (0->135)
t75 = -pkin(1) - pkin(6);
t53 = t75 * qJD(1) + qJD(2);
t73 = cos(qJ(3));
t26 = (pkin(4) * qJD(1) - t53) * t73;
t113 = qJD(4) + t26;
t70 = sin(qJ(5));
t117 = t70 * qJD(3);
t71 = sin(qJ(3));
t123 = qJD(1) * t71;
t72 = cos(qJ(5));
t35 = -t72 * t123 + t117;
t132 = t73 * t53;
t110 = qJD(3) * qJ(4);
t40 = t71 * t53;
t31 = -t40 - t110;
t142 = t31 * t73;
t27 = (-qJD(4) - t132) * qJD(3);
t126 = qJD(3) * pkin(3);
t93 = -qJD(4) + t126;
t28 = -t93 - t132;
t160 = ((-t28 + t132) * t71 + t142) * qJD(3) + t27 * t71;
t109 = qJD(1) * qJD(3);
t100 = t73 * t109;
t66 = qJD(1) * qJD(2);
t99 = t71 * t109;
t107 = pkin(3) * t100 + qJ(4) * t99 + t66;
t83 = (qJD(3) * pkin(7) - qJD(4)) * t73;
t11 = qJD(1) * t83 + t107;
t122 = qJD(3) * t71;
t38 = t53 * t122;
t21 = -pkin(4) * t99 + t38;
t74 = -pkin(3) - pkin(7);
t12 = qJD(3) * t74 + t113;
t129 = pkin(3) * t123 + qJD(1) * qJ(2);
t125 = t73 * qJ(4);
t87 = t71 * pkin(7) - t125;
t19 = t87 * qJD(1) + t129;
t5 = t72 * t12 - t70 * t19;
t1 = qJD(5) * t5 + t72 * t11 + t70 * t21;
t115 = t73 * qJD(1);
t58 = qJD(5) + t115;
t159 = -t5 * t58 + t1;
t14 = t35 * qJD(5) - t70 * t100;
t80 = t35 * t58;
t158 = t14 - t80;
t6 = t70 * t12 + t72 * t19;
t2 = -qJD(5) * t6 - t70 * t11 + t72 * t21;
t157 = t6 * t58 + t2;
t156 = t71 * pkin(3) + qJ(2);
t116 = t72 * qJD(3);
t37 = t70 * t123 + t116;
t139 = t37 * t58;
t15 = t37 * qJD(5) - t72 * t100;
t155 = -t15 + t139;
t88 = t5 * t70 - t6 * t72;
t153 = -qJD(5) * t88 + t1 * t70 + t2 * t72;
t152 = 0.2e1 * t66;
t149 = pkin(4) - t75;
t148 = t14 * t73;
t147 = t15 * t70;
t146 = t15 * t73;
t17 = (qJD(4) - t26) * qJD(3);
t145 = t17 * t70;
t144 = t17 * t72;
t141 = t35 * t71;
t140 = t37 * t35;
t138 = t58 * t73;
t137 = t58 * t74;
t136 = t71 * t14;
t135 = t71 * t72;
t134 = t72 * t14;
t133 = t72 * t73;
t76 = qJD(3) ^ 2;
t131 = t76 * t71;
t130 = t76 * t73;
t111 = qJ(4) * qJD(1);
t39 = pkin(3) * t115 + t71 * t111;
t68 = t71 ^ 2;
t69 = t73 ^ 2;
t128 = t68 - t69;
t77 = qJD(1) ^ 2;
t127 = t76 + t77;
t124 = t77 * qJ(2);
t121 = qJD(3) * t73;
t120 = qJD(5) * t70;
t119 = qJD(5) * t72;
t25 = -pkin(4) * t123 + t40;
t20 = t25 + t110;
t118 = t20 * qJD(3);
t114 = t73 * qJD(4);
t112 = qJ(2) * qJD(3);
t108 = t70 * t138;
t106 = 0.2e1 * qJD(1);
t105 = t58 * t120;
t104 = t71 * t120;
t103 = t58 * t119;
t101 = pkin(3) * t121 + t71 * t110 + qJD(2);
t97 = qJD(3) * t149;
t29 = -t73 * t111 + t129;
t41 = -t125 + t156;
t95 = qJD(1) * t41 + t29;
t94 = t58 + t115;
t92 = qJD(5) * t73 + qJD(1);
t91 = t73 * t99;
t89 = t5 * t72 + t6 * t70;
t32 = t87 + t156;
t43 = t149 * t73;
t10 = t72 * t32 + t70 * t43;
t9 = -t70 * t32 + t72 * t43;
t85 = t94 * t71;
t84 = -qJD(1) * t68 + t138;
t82 = t58 * t92;
t13 = -qJD(1) * t114 + t107;
t23 = t101 - t114;
t79 = -qJD(1) * t23 + t75 * t76 - t13;
t64 = qJ(2) * t152;
t56 = t73 * t77 * t71;
t49 = t70 * t99;
t48 = -0.2e1 * t91;
t47 = 0.2e1 * t91;
t46 = t127 * t73;
t45 = t127 * t71;
t44 = t128 * t77;
t42 = t149 * t71;
t34 = t73 * t97;
t33 = t71 * t97;
t30 = 0.2e1 * t128 * t109;
t24 = pkin(7) * t115 + t39;
t22 = t29 * t115;
t16 = t83 + t101;
t8 = t72 * t24 + t70 * t25;
t7 = -t70 * t24 + t72 * t25;
t4 = -t10 * qJD(5) - t70 * t16 - t72 * t33;
t3 = t9 * qJD(5) + t72 * t16 - t70 * t33;
t18 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t152, t64, t48, t30, -t131, t47, -t130, 0, -t75 * t131 + (qJD(2) * t71 + t112 * t73) * t106, -t75 * t130 + (qJD(2) * t73 - t112 * t71) * t106, 0, t64, 0, t131, t130, t48, t30, t47, t160, -t95 * t121 + t71 * t79, t95 * t122 + t73 * t79, t13 * t41 - t160 * t75 + t29 * t23, -t70 * t136 + (t73 * t117 + t71 * t119) * t37, (-t35 * t70 + t37 * t72) * t121 + (-t134 - t147 + (-t35 * t72 - t37 * t70) * qJD(5)) * t71, t71 * t103 - t148 + (-t37 * t71 + t70 * t84) * qJD(3), -t15 * t135 + (-t116 * t73 + t104) * t35, -t58 * t104 - t146 + (t72 * t84 + t141) * qJD(3), -qJD(3) * t85, -t42 * t15 - t34 * t35 + t4 * t58 + (-t116 * t20 + t2) * t73 + (t20 * t120 - t144 + (-qJD(1) * t9 - t5) * qJD(3)) * t71, t42 * t14 - t3 * t58 - t34 * t37 + (t20 * t117 - t1) * t73 + (t20 * t119 + t145 + (qJD(1) * t10 + t6) * qJD(3)) * t71, -t10 * t15 + t9 * t14 - t3 * t35 - t4 * t37 - t88 * t121 + (-qJD(5) * t89 + t1 * t72 - t2 * t70) * t71, t1 * t10 - t17 * t42 + t2 * t9 - t20 * t34 + t6 * t3 + t5 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t77, -t124, 0, 0, 0, 0, 0, 0, -t45, -t46, 0, -t124, 0, 0, 0, 0, 0, 0, 0, t45, t46, -t29 * qJD(1) - t160, 0, 0, 0, 0, 0, 0, t71 * t15 + t70 * t82 + (t94 * t135 + t35 * t73) * qJD(3), -t136 + t72 * t82 + (t37 * t73 - t70 * t85) * qJD(3), (-t37 * t122 + t35 * t92 - t148) * t72 + (-t35 * t122 - t37 * t92 + t146) * t70, t88 * qJD(1) + (t89 * qJD(3) + t17) * t71 + (t118 - t153) * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, -t44, 0, -t56, 0, 0, -t73 * t124, t71 * t124, 0, 0, 0, 0, 0, t56, -t44, -t56, ((-t31 - t110) * t73 + (t28 + t93) * t71) * qJD(1), t39 * t123 + t22, 0.2e1 * qJD(3) * qJD(4) + (-t29 * t71 + t39 * t73) * qJD(1), -t27 * qJ(4) - t31 * qJD(4) - t29 * t39 + (t142 + (-t28 - t126) * t71) * t53, -t139 * t70 - t134, (-t15 - t139) * t72 + (t14 + t80) * t70, -t105 + (-t108 + (t37 - t116) * t71) * qJD(1), t72 * t80 + t147, -t103 + t49 + (-t58 * t133 - t141) * qJD(1), t58 * t123, qJ(4) * t15 + t145 - t7 * t58 + t113 * t35 + (-t70 * t137 + t20 * t72) * qJD(5) + (t20 * t133 + (-t116 * t74 + t5) * t71) * qJD(1), -qJ(4) * t14 + t144 + t8 * t58 + t113 * t37 + (-t72 * t137 - t20 * t70) * qJD(5) + (-t6 * t71 + (t74 * t122 - t20 * t73) * t70) * qJD(1), t8 * t35 + t7 * t37 + (-t6 * t115 + t14 * t74 - t2 + (-t35 * t74 - t6) * qJD(5)) * t72 + (t5 * t115 - t15 * t74 - t1 + (t37 * t74 + t5) * qJD(5)) * t70, t17 * qJ(4) + t113 * t20 + t153 * t74 - t5 * t7 - t6 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t56, -t69 * t77 - t76, t31 * qJD(3) + t22 + t38, 0, 0, 0, 0, 0, 0, -t105 - qJD(3) * t35 + (-t116 * t71 - t108) * qJD(1), -t58 ^ 2 * t72 - qJD(3) * t37 + t49, t155 * t70 + t158 * t72, t157 * t72 + t159 * t70 - t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t140, -t35 ^ 2 + t37 ^ 2, -t158, -t140, t155, -t99, -t20 * t37 + t157, t20 * t35 - t159, 0, 0;];
tauc_reg = t18;
