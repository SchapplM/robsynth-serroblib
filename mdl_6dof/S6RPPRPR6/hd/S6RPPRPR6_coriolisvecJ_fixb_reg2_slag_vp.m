% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RPPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPPRPR6_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR6_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR6_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRPR6_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:51:33
% EndTime: 2019-03-09 01:51:38
% DurationCPUTime: 1.67s
% Computational Cost: add. (1824->262), mult. (3394->348), div. (0->0), fcn. (1649->4), ass. (0->154)
t80 = cos(qJ(6));
t141 = qJD(4) * t80;
t79 = sin(qJ(4));
t144 = qJD(1) * t79;
t78 = sin(qJ(6));
t39 = t78 * t144 + t141;
t81 = cos(qJ(4));
t143 = qJD(1) * t81;
t61 = qJD(6) + t143;
t158 = t39 * t61;
t128 = qJD(1) * qJD(4);
t118 = t81 * t128;
t20 = qJD(6) * t39 - t80 * t118;
t178 = t20 - t158;
t68 = qJD(1) * qJ(2) + qJD(3);
t55 = -pkin(7) * qJD(1) + t68;
t30 = (pkin(5) * qJD(1) - t55) * t81;
t132 = qJD(5) + t30;
t131 = pkin(4) * t144 - qJD(2);
t77 = pkin(1) + qJ(3);
t108 = -qJ(5) * t81 + t77;
t87 = pkin(8) * t79 + t108;
t15 = qJD(1) * t87 + t131;
t82 = -pkin(4) - pkin(8);
t139 = qJD(4) * t82;
t16 = t139 + t132;
t100 = t15 * t78 - t16 * t80;
t119 = t79 * t128;
t72 = qJD(3) * qJD(1);
t125 = pkin(4) * t118 + qJ(5) * t119 + t72;
t93 = (qJD(4) * pkin(8) - qJD(5)) * t81;
t11 = qJD(1) * t93 + t125;
t133 = t81 * qJD(2);
t142 = qJD(4) * t79;
t36 = t55 * t142;
t14 = t36 + (-pkin(5) * t142 - t133) * qJD(1);
t1 = -qJD(6) * t100 + t80 * t11 + t78 * t14;
t177 = t100 * t61 + t1;
t121 = t80 * t144;
t137 = qJD(6) * t78;
t19 = qJD(4) * t137 - qJD(6) * t121 - t78 * t118;
t37 = qJD(4) * t78 - t121;
t90 = t37 * t61;
t86 = t19 - t90;
t6 = t15 * t80 + t16 * t78;
t2 = -qJD(6) * t6 - t78 * t11 + t80 * t14;
t176 = t6 * t61 + t2;
t174 = qJD(1) * t77;
t74 = t79 ^ 2;
t75 = t81 ^ 2;
t146 = t74 + t75;
t173 = t146 * qJD(2);
t149 = t81 * t55;
t129 = qJD(1) * qJD(2);
t62 = t79 * t129;
t22 = -t62 + (-qJD(5) - t149) * qJD(4);
t25 = -t81 * t129 + t36;
t112 = qJD(4) * pkin(4) - qJD(5);
t32 = -t112 - t149;
t130 = qJD(4) * qJ(5);
t41 = t79 * t55;
t34 = -t41 - t130;
t85 = -t22 * t79 - t25 * t81 + (t32 * t79 - t34 * t81) * qJD(4);
t103 = -t100 * t78 - t6 * t80;
t172 = -qJD(6) * t103 + t1 * t78 + t2 * t80;
t104 = -t100 * t80 + t6 * t78;
t171 = -t104 * qJD(6) + t1 * t80 - t2 * t78;
t170 = 0.2e1 * t72;
t76 = -pkin(7) + qJ(2);
t167 = pkin(5) - t76;
t12 = t62 + (qJD(5) - t30) * qJD(4);
t166 = t12 * t78;
t165 = t12 * t80;
t164 = t19 * t81;
t163 = t20 * t78;
t162 = t20 * t81;
t161 = t34 * t79;
t160 = t37 * t79;
t159 = t39 * t37;
t157 = t61 * t81;
t156 = t61 * t82;
t84 = qJD(1) ^ 2;
t155 = t74 * t84;
t83 = qJD(4) ^ 2;
t154 = t76 * t83;
t153 = t79 * t19;
t152 = t79 * t80;
t151 = t80 * t19;
t150 = t80 * t81;
t148 = t83 * t79;
t147 = t83 * t81;
t40 = pkin(4) * t143 + qJ(5) * t144;
t145 = t83 + t84;
t140 = qJD(4) * t81;
t138 = qJD(5) * t81;
t136 = qJD(6) * t80;
t29 = -pkin(5) * t144 + t41;
t23 = t29 + t130;
t135 = t23 * qJD(4);
t24 = t108 * qJD(1) + t131;
t134 = t24 * qJD(1);
t127 = t78 * t157;
t126 = t61 * t150;
t124 = t61 * t137;
t123 = t79 * t137;
t122 = t61 * t136;
t70 = 0.2e1 * t129;
t120 = pkin(4) * t140 + t79 * t130 + qJD(3);
t117 = qJD(4) * t167;
t115 = qJD(1) * t146;
t56 = -qJD(2) + t174;
t114 = t56 + t174;
t113 = t61 + t143;
t111 = (qJD(2) - t56) * qJD(1);
t110 = qJD(6) * t81 + qJD(1);
t109 = t79 * t118;
t51 = t78 * t119;
t107 = t51 - t122;
t102 = qJD(2) * t115;
t71 = t79 * pkin(4);
t33 = t71 + t87;
t44 = t167 * t81;
t10 = t33 * t80 + t44 * t78;
t9 = -t33 * t78 + t44 * t80;
t42 = t108 + t71;
t97 = qJD(1) * t42 + qJD(2) + t24;
t96 = qJD(2) + t114;
t95 = t113 * t79;
t94 = -qJD(1) * t74 + t157;
t92 = t61 * t110;
t89 = t170 - t154;
t17 = -qJD(1) * t138 + t125;
t28 = t120 - t138;
t88 = -qJD(1) * t28 + t154 - t17;
t69 = t75 * t84;
t59 = t81 * t84 * t79;
t50 = -0.2e1 * t109;
t49 = 0.2e1 * t109;
t48 = t145 * t81;
t47 = t145 * t79;
t46 = t69 - t155;
t45 = t69 + t155;
t43 = t167 * t79;
t35 = 0.2e1 * (t74 - t75) * t128;
t31 = pkin(8) * t143 + t40;
t27 = t79 * qJD(2) - t81 * t117;
t26 = -t79 * t117 - t133;
t21 = t93 + t120;
t18 = t81 * t134;
t8 = t29 * t78 + t31 * t80;
t7 = t29 * t80 - t31 * t78;
t4 = -qJD(6) * t10 - t78 * t21 + t80 * t26;
t3 = qJD(6) * t9 + t80 * t21 + t78 * t26;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, qJ(2) * t70, 0, 0, 0, 0, 0, 0, 0, t70, t170, t68 * qJD(2) + t56 * qJD(3) + (qJ(2) * qJD(2) + qJD(3) * t77) * qJD(1), t50, t35, -t148, t49, -t147, 0, t96 * t140 + t79 * t89, -t96 * t142 + t81 * t89, -0.2e1 * t102, t114 * qJD(3) + (qJD(1) * t76 + t55) * t173, 0, t148, t147, t50, t35, t49, -t102 - t85, -t97 * t140 + t79 * t88, t97 * t142 + t81 * t88, t17 * t42 + t24 * t28 + (-t32 * t81 - t161) * qJD(2) + t85 * t76, -t78 * t153 + (t136 * t79 + t140 * t78) * t39 (-t37 * t78 + t39 * t80) * t140 + (-t151 - t163 + (-t37 * t80 - t39 * t78) * qJD(6)) * t79, t79 * t122 - t164 + (-t39 * t79 + t78 * t94) * qJD(4), -t20 * t152 + (-t140 * t80 + t123) * t37, -t61 * t123 - t162 + (t80 * t94 + t160) * qJD(4), -qJD(4) * t95, -t43 * t20 + t27 * t37 + t4 * t61 + (-t135 * t80 + t2) * t81 + (t23 * t137 - t165 + (-qJD(1) * t9 + t100) * qJD(4)) * t79, t43 * t19 + t27 * t39 - t3 * t61 + (t135 * t78 - t1) * t81 + (t23 * t136 + t166 + (qJD(1) * t10 + t6) * qJD(4)) * t79, -t10 * t20 - t103 * t140 + t171 * t79 + t9 * t19 - t3 * t37 - t4 * t39, t1 * t10 - t100 * t4 - t12 * t43 + t2 * t9 + t23 * t27 + t3 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t84, -t84 * qJ(2), 0, 0, 0, 0, 0, 0, 0, -t84, 0, -qJD(1) * t68 - t72, 0, 0, 0, 0, 0, 0, -0.2e1 * t118, 0.2e1 * t119, t45, -t55 * t115 - t72, 0, 0, 0, 0, 0, 0, t45, 0.2e1 * t118, -0.2e1 * t119 (t161 + (qJD(5) + t32) * t81) * qJD(1) - t125, 0, 0, 0, 0, 0, 0 (t126 - t160) * qJD(1) - t107, -t124 + (-t127 + (-t39 - t141) * t79) * qJD(1), t178 * t80 + t86 * t78 (t104 * t81 - t23 * t79) * qJD(1) - t171; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t84, t111, 0, 0, 0, 0, 0, 0, -t47, -t48, 0 (-t56 + t173) * qJD(1), 0, 0, 0, 0, 0, 0, 0, t47, t48, t85 - t134, 0, 0, 0, 0, 0, 0, t79 * t20 + t78 * t92 + (t113 * t152 + t37 * t81) * qJD(4), -t153 + t80 * t92 + (t39 * t81 - t78 * t95) * qJD(4) (t110 * t37 - t39 * t142 - t164) * t80 + (-t110 * t39 - t37 * t142 + t162) * t78, t103 * qJD(1) + (t104 * qJD(4) + t12) * t79 + (t135 - t172) * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, t46, 0, -t59, 0, 0, t81 * t111, t56 * t144 - t62, 0, 0, 0, 0, 0, t59, t46, -t59 ((-t34 - t130) * t81 + (t112 + t32) * t79) * qJD(1), t18 + (t40 * t79 - t133) * qJD(1), 0.2e1 * qJD(4) * qJD(5) + t62 + (-t24 * t79 + t40 * t81) * qJD(1), -t32 * t41 - t25 * pkin(4) - t22 * qJ(5) - t24 * t40 + (-qJD(5) + t149) * t34, -t158 * t78 - t151 (-t20 - t158) * t80 + (t19 + t90) * t78, -t124 + (-t127 + (t39 - t141) * t79) * qJD(1), t80 * t90 + t163 (-t126 - t160) * qJD(1) + t107, t61 * t144, qJ(5) * t20 + t166 - t7 * t61 + t132 * t37 + (-t78 * t156 + t23 * t80) * qJD(6) + (t23 * t150 + (-t139 * t80 - t100) * t79) * qJD(1), -qJ(5) * t19 + t165 + t8 * t61 + t132 * t39 + (-t80 * t156 - t23 * t78) * qJD(6) + (-t6 * t79 + (t139 * t79 - t23 * t81) * t78) * qJD(1), t8 * t37 + t7 * t39 + (-t6 * t143 + t19 * t82 - t2 + (-t37 * t82 - t6) * qJD(6)) * t80 + (-t100 * t143 - t20 * t82 - t1 + (t39 * t82 - t100) * qJD(6)) * t78, t12 * qJ(5) + t100 * t7 + t132 * t23 + t172 * t82 - t6 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, -t69 - t83, qJD(4) * t34 + t18 + t25, 0, 0, 0, 0, 0, 0, -t124 - qJD(4) * t37 + (-t141 * t79 - t127) * qJD(1), -t61 ^ 2 * t80 - qJD(4) * t39 + t51, -t178 * t78 + t86 * t80, t176 * t80 + t177 * t78 - t135; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t159, -t37 ^ 2 + t39 ^ 2, -t86, -t159, -t178, -t119, -t23 * t39 + t176, t23 * t37 - t177, 0, 0;];
tauc_reg  = t5;
