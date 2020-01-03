% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRPR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5]';
% 
% Output:
% tau_reg [5x24]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRPR16_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR16_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR16_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR16_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR16_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPR16_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:39:35
% EndTime: 2019-12-31 18:39:39
% DurationCPUTime: 1.29s
% Computational Cost: add. (918->253), mult. (1703->317), div. (0->0), fcn. (956->6), ass. (0->148)
t83 = -pkin(1) - pkin(6);
t49 = t83 * qJD(1) + qJD(2);
t80 = cos(qJ(3));
t20 = (pkin(4) * qJD(1) - t49) * t80;
t137 = qJD(4) + t20;
t129 = qJD(1) * qJD(3);
t120 = t80 * t129;
t77 = sin(qJ(3));
t131 = t77 * qJDD(1);
t181 = t120 + t131;
t81 = cos(qJ(1));
t69 = g(2) * t81;
t78 = sin(qJ(1));
t70 = g(1) * t78;
t179 = t70 - t69;
t136 = qJ(4) * qJD(1);
t148 = qJD(1) * t77;
t157 = pkin(3) * t148 + qJD(1) * qJ(2);
t23 = -t80 * t136 + t157;
t180 = -qJD(1) * t23 - t179;
t121 = t77 * t129;
t64 = t80 * qJDD(1);
t178 = -t121 + t64;
t76 = sin(qJ(5));
t140 = qJD(5) * t76;
t34 = -qJDD(5) - t178;
t79 = cos(qJ(5));
t21 = t79 * t34;
t147 = qJD(1) * t80;
t55 = qJD(5) + t147;
t177 = -t55 * t140 - t21;
t125 = t79 * t148;
t144 = qJD(3) * t76;
t35 = -t125 + t144;
t146 = qJD(3) * t35;
t176 = t146 + t21;
t151 = qJ(4) * t80;
t67 = t77 * pkin(3);
t116 = -t67 + t151;
t72 = qJD(1) * qJD(2);
t124 = 0.2e1 * t72;
t71 = qJDD(1) * qJ(2);
t175 = t124 + 0.2e1 * t71;
t135 = qJD(3) * qJ(4);
t41 = t77 * t49;
t19 = -pkin(4) * t148 + t41;
t15 = t19 + t135;
t82 = -pkin(3) - pkin(7);
t174 = -t82 * t34 + (t15 - t19) * t55;
t128 = qJDD(3) * qJ(4);
t48 = t83 * qJDD(1) + qJDD(2);
t40 = t77 * t48;
t109 = -t40 - t128;
t165 = t49 * t80;
t10 = (-qJD(4) - t165) * qJD(3) + t109;
t143 = qJD(3) * t77;
t38 = t49 * t143;
t102 = -t48 * t80 + qJDD(4) + t38;
t138 = qJDD(3) * pkin(3);
t11 = t102 - t138;
t111 = qJD(3) * pkin(3) - qJD(4);
t22 = -t111 - t165;
t24 = -t41 - t135;
t89 = -t10 * t77 - t11 * t80 + (t22 * t77 - t24 * t80) * qJD(3);
t142 = qJD(3) * t79;
t37 = t76 * t148 + t142;
t9 = qJD(5) * t37 + qJDD(3) * t76 - t181 * t79;
t132 = qJDD(3) * t83;
t42 = qJ(2) - t116;
t173 = (qJD(1) * t42 + t23) * qJD(3) + t132;
t171 = g(3) * t77;
t170 = g(3) * t80;
t8 = -qJD(3) * t140 + qJD(5) * t125 + t79 * qJDD(3) + t181 * t76;
t169 = t8 * t79;
t168 = pkin(4) - t83;
t167 = t35 * t55;
t166 = t37 * t55;
t164 = t55 * t76;
t163 = t76 * t34;
t162 = t76 * t77;
t161 = t77 * t79;
t160 = t78 * t80;
t159 = t79 * t80;
t158 = t80 * t81;
t39 = pkin(3) * t147 + t77 * t136;
t156 = t81 * pkin(1) + t78 * qJ(2);
t74 = t77 ^ 2;
t75 = t80 ^ 2;
t154 = t74 - t75;
t84 = qJD(3) ^ 2;
t85 = qJD(1) ^ 2;
t153 = t84 + t85;
t152 = qJ(2) * t85;
t150 = pkin(1) * qJDD(1);
t145 = qJD(3) * t37;
t141 = qJD(3) * t80;
t139 = qJD(5) * t79;
t134 = qJDD(3) * t77;
t133 = qJDD(3) * t80;
t130 = qJ(4) * qJDD(1);
t127 = t80 * t85 * t77;
t123 = pkin(3) * t141 + t77 * t135 + qJD(2);
t106 = pkin(3) * t181 + qJ(4) * t121 + t71 + t72;
t110 = qJD(3) * pkin(7) - qJD(4);
t3 = pkin(7) * t131 + (qJD(1) * t110 - t130) * t80 + t106;
t5 = t178 * pkin(4) + t82 * qJDD(3) + t102;
t122 = -t76 * t3 + t79 * t5;
t119 = qJD(3) * t168;
t117 = qJD(1) * t39 - g(3);
t12 = t82 * qJD(3) + t137;
t115 = qJD(5) * t12 + t3;
t100 = pkin(7) * t77 - t151;
t14 = qJD(1) * t100 + t157;
t114 = qJD(5) * t14 - t5;
t112 = (-t74 - t75) * qJDD(1);
t108 = qJD(5) * t80 + qJD(1);
t107 = qJDD(2) - t150;
t105 = g(1) * t81 + g(2) * t78;
t101 = pkin(3) * t80 + qJ(4) * t77;
t99 = -t152 + t69;
t2 = t12 * t76 + t14 * t79;
t95 = t145 - t163;
t93 = -t55 * t139 + t163;
t92 = 0.2e1 * qJ(2) * t129 + t132;
t91 = -t83 * t84 - t105;
t56 = g(1) * t160;
t90 = qJDD(4) + t23 * t147 + t56 + (-t48 - t69) * t80;
t88 = t91 + t175;
t17 = -qJD(4) * t80 + t123;
t7 = (-qJD(1) * qJD(4) - t130) * t80 + t106;
t87 = -qJD(1) * t17 - qJDD(1) * t42 - t7 - t91;
t6 = -pkin(4) * t131 + (qJD(4) - t20) * qJD(3) - t109;
t86 = -t170 + t6 - t179 * t77 + (pkin(7) * t147 - qJD(5) * t82 + t39) * t55;
t66 = t81 * qJ(2);
t44 = t168 * t80;
t43 = t168 * t77;
t33 = t80 * t119;
t32 = t77 * t119;
t31 = qJ(2) + t100 + t67;
t30 = -t153 * t77 + t133;
t29 = t153 * t80 + t134;
t28 = -t76 * t160 + t79 * t81;
t27 = -t78 * t159 - t76 * t81;
t26 = -t76 * t158 - t78 * t79;
t25 = -t79 * t158 + t76 * t78;
t13 = t110 * t80 + t123;
t1 = t12 * t79 - t14 * t76;
t4 = [qJDD(1), t179, t105, qJDD(2) - 0.2e1 * t150 - t179, -t105 + t175, -t107 * pkin(1) - g(1) * (-pkin(1) * t78 + t66) - g(2) * t156 + (t124 + t71) * qJ(2), qJDD(1) * t75 - 0.2e1 * t120 * t77, 0.2e1 * t154 * t129 - 0.2e1 * t77 * t64, -t77 * t84 + t133, -t80 * t84 - t134, 0, t77 * t88 + t80 * t92, -t77 * t92 + t80 * t88, t112 * t83 + t179 - t89, -t173 * t80 + t87 * t77, t173 * t77 + t87 * t80, t7 * t42 + t23 * t17 - g(1) * (-t116 * t81 + t66) - g(2) * (pkin(6) * t81 + t156) + (-g(1) * t83 + g(2) * t116) * t78 + t89 * t83, t8 * t162 + (t77 * t139 + t76 * t141) * t37, (-t35 * t76 + t37 * t79) * t141 + (-t76 * t9 + t169 + (-t35 * t79 - t37 * t76) * qJD(5)) * t77, (t55 * t144 + t8) * t80 + (-t93 - t145) * t77, (t55 * t142 - t9) * t80 + (t146 + t177) * t77, -t55 * t143 - t34 * t80, (-t76 * t13 - t79 * t32) * t55 - (-t31 * t76 + t44 * t79) * t34 + t122 * t80 - t33 * t35 - t43 * t9 - t6 * t161 - g(1) * t26 - g(2) * t28 + (-t1 * t77 - t15 * t159) * qJD(3) + ((-t31 * t79 - t44 * t76) * t55 - t2 * t80 + t15 * t162) * qJD(5), t2 * t143 - g(1) * t25 - g(2) * t27 - t33 * t37 - t43 * t8 + (-(qJD(5) * t44 + t13) * t55 + t31 * t34 - t115 * t80 + t15 * qJD(5) * t77) * t79 + (-(-qJD(5) * t31 - t32) * t55 + t44 * t34 + t6 * t77 + (qJD(3) * t15 + t114) * t80) * t76; 0, 0, 0, qJDD(1), -t85, t107 - t152 - t179, 0, 0, 0, 0, 0, t30, -t29, t112, -t30, t29, t89 + t180, 0, 0, 0, 0, 0, t77 * t9 + t176 * t80 + (t108 * t76 + t77 * t142) * t55, t77 * t8 + t95 * t80 + (t108 * t79 - t76 * t143) * t55; 0, 0, 0, 0, 0, 0, t127, -t154 * t85, t64, -t131, qJDD(3), t171 - t56 + (t48 + t99) * t80, t170 - t40 + (-t99 + t70) * t77, -t101 * qJDD(1) + ((-t24 - t135) * t80 + (t111 + t22) * t77) * qJD(1), t117 * t77 - 0.2e1 * t138 + t90, 0.2e1 * qJD(3) * qJD(4) + t117 * t80 + t180 * t77 + 0.2e1 * t128 + t40, -t10 * qJ(4) - t11 * pkin(3) - t23 * t39 - t22 * t41 - g(3) * t116 + (-qJD(4) + t165) * t24 - t179 * t101, -t164 * t37 + t169, (-t9 - t166) * t79 + (-t8 + t167) * t76, (-t80 * t164 + t37 * t77) * qJD(1) + t177, (-t55 * t159 - t35 * t77) * qJD(1) + t93, t55 * t148, qJ(4) * t9 + t1 * t148 + t137 * t35 + t174 * t79 + t86 * t76, qJ(4) * t8 + t137 * t37 - t2 * t148 - t174 * t76 + t86 * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, qJDD(3) - t127, -t75 * t85 - t84, qJD(3) * t24 - t138 - t171 + t38 + t90, 0, 0, 0, 0, 0, -t164 * t55 - t176, -t55 ^ 2 * t79 - t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37 * t35, -t35 ^ 2 + t37 ^ 2, t8 + t167, t166 - t9, -t34, -g(1) * t27 + g(2) * t25 - g(3) * t161 - t15 * t37 + t122 + (-qJD(5) + t55) * t2, g(1) * t28 - g(2) * t26 + t1 * t55 + t15 * t35 - t115 * t79 + (t114 + t171) * t76;];
tau_reg = t4;
