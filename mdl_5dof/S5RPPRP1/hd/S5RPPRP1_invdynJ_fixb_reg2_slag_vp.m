% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPPRP1
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:13
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPRP1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP1_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:12:45
% EndTime: 2022-01-23 09:12:49
% DurationCPUTime: 1.52s
% Computational Cost: add. (1502->260), mult. (3023->344), div. (0->0), fcn. (1920->10), ass. (0->162)
t90 = sin(pkin(8));
t92 = cos(pkin(8));
t107 = -t92 * pkin(3) - t90 * pkin(6) - pkin(2);
t93 = cos(pkin(7));
t174 = t93 * pkin(1);
t53 = t107 - t174;
t97 = cos(qJ(4));
t165 = t92 * t97;
t91 = sin(pkin(7));
t75 = t91 * pkin(1) + qJ(3);
t57 = t75 * t165;
t95 = sin(qJ(4));
t26 = t95 * t53 + t57;
t185 = qJD(4) * t26;
t138 = qJD(1) * qJD(4);
t123 = t97 * t138;
t144 = qJDD(1) * t95;
t104 = t123 + t144;
t184 = t104 * t90;
t159 = pkin(1) * qJDD(1);
t56 = qJD(1) * qJD(3) + qJDD(1) * t75;
t166 = t92 * t95;
t87 = qJ(1) + pkin(7);
t82 = sin(t87);
t83 = cos(t87);
t41 = t82 * t166 + t83 * t97;
t43 = -t83 * t166 + t82 * t97;
t183 = -g(1) * t43 + g(2) * t41;
t176 = g(3) * t95;
t182 = t90 * t176 + t183;
t148 = t92 * qJD(1);
t70 = -qJD(4) + t148;
t158 = qJD(1) * t90;
t127 = qJ(5) * t158;
t38 = t53 * qJD(1) + qJD(3);
t62 = t75 * qJD(1);
t47 = t90 * qJD(2) + t92 * t62;
t14 = t97 * t38 - t95 * t47;
t8 = -t97 * t127 + t14;
t5 = -t70 * pkin(4) + t8;
t181 = -t8 + t5;
t180 = pkin(4) * t95;
t178 = g(1) * t82;
t77 = g(2) * t83;
t141 = t92 * qJDD(1);
t69 = -qJDD(4) + t141;
t175 = t69 * pkin(4);
t81 = t92 * qJD(2);
t46 = t90 * t62 - t81;
t173 = t46 * t90;
t172 = t69 * t92;
t171 = t83 * t92;
t170 = t83 * t95;
t85 = t90 ^ 2;
t99 = qJD(1) ^ 2;
t169 = t85 * t99;
t168 = t90 * t56;
t167 = t90 * (-qJ(5) - pkin(6));
t150 = qJD(4) * t97;
t154 = qJD(3) * t97;
t164 = t53 * t150 + t92 * t154;
t86 = t92 ^ 2;
t163 = t85 + t86;
t88 = t95 ^ 2;
t89 = t97 ^ 2;
t162 = -t88 - t89;
t161 = t88 - t89;
t160 = qJ(5) * t90;
t157 = qJD(1) * t95;
t156 = qJD(1) * t97;
t155 = qJD(3) * t95;
t15 = t95 * t38 + t97 * t47;
t153 = qJD(4) * t15;
t152 = qJD(4) * t47;
t151 = qJD(4) * t95;
t149 = qJD(5) * t90;
t28 = qJD(5) - t81 + (pkin(4) * t157 + t62) * t90;
t147 = qJD(5) + t28;
t78 = -pkin(2) - t174;
t145 = qJDD(1) * t78;
t143 = qJDD(1) * t97;
t142 = t90 * qJDD(1);
t140 = qJ(5) * qJDD(1);
t137 = qJD(1) * qJD(5);
t136 = t95 * t169;
t135 = t75 * t166;
t98 = cos(qJ(1));
t134 = t98 * pkin(1) + t83 * pkin(2) + t82 * qJ(3);
t133 = t97 * t160;
t132 = t90 * t157;
t131 = t75 * t151;
t130 = t70 * t151;
t129 = t46 * t158;
t128 = t92 * t155;
t96 = sin(qJ(1));
t126 = -t96 * pkin(1) + t83 * qJ(3);
t125 = t77 - t178;
t80 = t92 * qJDD(2);
t36 = -t80 + t168;
t124 = -t36 * t92 - g(3);
t35 = t53 * qJDD(1) + qJDD(3);
t37 = t90 * qJDD(2) + t92 * t56;
t3 = t38 * t150 - t47 * t151 + t95 * t35 + t97 * t37;
t50 = (t75 + t180) * t90;
t122 = qJD(1) * t50 + t28;
t121 = t69 - t141;
t120 = t69 + t141;
t119 = t90 * pkin(4) * t123 + t142 * t180 + qJDD(5) - t80;
t118 = qJD(4) * t132;
t117 = t95 * t123;
t116 = -g(1) * t41 - g(2) * t43;
t42 = -t82 * t165 + t170;
t44 = t83 * t165 + t82 * t95;
t115 = -g(1) * t42 - g(2) * t44;
t114 = -g(1) * t83 - g(2) * t82;
t113 = g(1) * t96 - g(2) * t98;
t9 = -t95 * t127 + t15;
t112 = t5 * t97 + t9 * t95;
t111 = t5 * t95 - t9 * t97;
t110 = t14 * t95 - t15 * t97;
t109 = t36 * t90 + t37 * t92;
t108 = t47 * t92 + t173;
t59 = qJDD(3) + t145;
t106 = t145 + t59 + t77;
t19 = t119 + t168;
t55 = (pkin(4) * t150 + qJD(3)) * t90;
t105 = qJD(1) * t55 + qJDD(1) * t50 + t19;
t31 = t97 * t35;
t103 = t31 + qJ(5) * t118 + (-qJD(4) * t38 - t37) * t95;
t102 = -t70 ^ 2 - t169;
t101 = g(3) * t90 * t97 + g(1) * t44 - g(2) * t42 - t3;
t4 = -t95 * t37 - t153 + t31;
t79 = t97 * pkin(4) + pkin(3);
t67 = t97 * t142;
t64 = t90 * t178;
t61 = t97 * t136;
t58 = t92 * t118;
t54 = t90 * t70 * t156;
t52 = t161 * t169;
t51 = t162 * t142;
t49 = t97 * t53;
t40 = (qJDD(1) * t89 - 0.2e1 * t117) * t85;
t39 = (qJDD(1) * t88 + 0.2e1 * t117) * t85;
t27 = 0.2e1 * (t161 * t138 - t95 * t143) * t85;
t25 = t49 - t135;
t24 = -t54 - t184;
t23 = t67 + (-qJD(4) - t70) * t132;
t22 = -t95 * t160 + t26;
t21 = -t128 - t185;
t20 = -t92 * t131 + t164;
t18 = -t133 + t49 + (-t75 * t95 - pkin(4)) * t92;
t17 = t102 * t97 + t95 * t69;
t16 = t102 * t95 - t97 * t69;
t13 = (t121 * t95 + (t70 - t148) * t150) * t90;
t12 = (t120 * t95 + (t70 + t148) * t150) * t90;
t11 = t58 + (-t120 * t97 + t130) * t90;
t10 = t58 + (t121 * t97 - t130) * t90;
t7 = -t128 - t97 * t149 + (-t57 + (-t53 + t160) * t95) * qJD(4);
t6 = -t95 * t149 + (-t133 - t135) * qJD(4) + t164;
t2 = (-t104 * qJ(5) - t95 * t137) * t90 + t3;
t1 = -t175 + (-t152 + (-t137 - t140) * t90) * t97 + t103;
t29 = [0, 0, 0, 0, 0, qJDD(1), t113, g(1) * t98 + g(2) * t96, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0.2e1 * t93 * t159 - t125, -0.2e1 * t91 * t159 - t114, 0, (t113 + (t91 ^ 2 + t93 ^ 2) * t159) * pkin(1), t85 * qJDD(1), 0.2e1 * t90 * t141, 0, t86 * qJDD(1), 0, 0, (-t106 + t178) * t92, t106 * t90 - t64, t56 * t163 + t109 + t114, t59 * t78 - g(1) * (-t82 * pkin(2) + t126) - g(2) * t134 + t109 * t75 + t108 * qJD(3), t40, t27, t11, t39, t12, t172, -t21 * t70 - t25 * t69 - t4 * t92 + (t46 * t150 + t36 * t95) * t90 + (t75 * t144 + (t75 * t150 + t155) * qJD(1)) * t85 + t115, t20 * t70 + t26 * t69 + t3 * t92 + (-t46 * t151 + t36 * t97) * t90 + (t75 * t143 + (-t131 + t154) * qJD(1)) * t85 + t116, t64 + (-t77 + (-t153 - qJDD(1) * t25 - t4 + (-t21 - t185) * qJD(1)) * t97 + (qJD(4) * t14 - qJDD(1) * t26 - t3 + (qJD(4) * t25 - t20) * qJD(1)) * t95) * t90, t3 * t26 + t15 * t20 + t4 * t25 + t14 * t21 - g(1) * t126 - g(2) * (pkin(3) * t171 + t134) + (-pkin(6) * t77 + t46 * qJD(3) + t36 * t75) * t90 - t107 * t178, t40, t27, t11, t39, t12, t172, -t1 * t92 - t18 * t69 - t7 * t70 + (t105 * t95 + t122 * t150) * t90 + t115, t2 * t92 + t22 * t69 + t6 * t70 + (t105 * t97 - t122 * t151) * t90 + t116, t64 + (-t77 + (-qJD(4) * t9 - qJDD(1) * t18 - t1 + (-qJD(4) * t22 - t7) * qJD(1)) * t97 + (qJD(4) * t5 - qJDD(1) * t22 - t2 + (qJD(4) * t18 - t6) * qJD(1)) * t95) * t90, t2 * t22 + t9 * t6 + t1 * t18 + t5 * t7 + t19 * t50 + t28 * t55 - g(1) * (pkin(4) * t170 + t126) - g(2) * (-t83 * t167 + t79 * t171 + t134) + (-g(1) * (-t79 * t92 - pkin(2) + t167) - g(2) * t180) * t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2) - g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t37 * t90 + t124, 0, 0, 0, 0, 0, 0, t13, t10, 0, (t3 * t97 - t4 * t95 + (-t14 * t97 - t15 * t95) * qJD(4)) * t90 + t124, 0, 0, 0, 0, 0, 0, t13, t10, 0, -t19 * t92 - g(3) + (-qJD(4) * t112 - t1 * t95 + t2 * t97) * t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t141, t142, -t163 * t99, -qJD(1) * t108 + t125 + t59, 0, 0, 0, 0, 0, 0, t16, t17, t51, t3 * t95 + t4 * t97 - t110 * qJD(4) + (t110 * t92 - t173) * qJD(1) + t125, 0, 0, 0, 0, 0, 0, t16, t17, t51, t1 * t97 + t2 * t95 - t111 * qJD(4) + (t111 * t92 - t28 * t90) * qJD(1) + t125; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, -t52, t23, -t61, t24, -t69, -t129 * t97 - t15 * t70 + t182 + t4, t129 * t95 - t14 * t70 + t101, 0, 0, t61, -t52, t23, -t61, t24, -t69, -0.2e1 * t175 - t9 * t70 + (-pkin(4) * t136 - t152 + (-t147 * qJD(1) - t140) * t90) * t97 + t103 + t182, -t89 * pkin(4) * t169 - t8 * t70 + (t95 * t140 + (qJ(5) * t150 + t147 * t95) * qJD(1)) * t90 + t101, (-pkin(4) * t143 + (pkin(4) * qJD(4) - t181) * t157) * t90, t181 * t9 + (t1 + (-t28 * t156 + t176) * t90 + t183) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54 + t184, t67 + (-qJD(4) + t70) * t132, t162 * t169, g(3) * t92 + (qJD(1) * t112 + t114 + t56) * t90 + t119;];
tau_reg = t29;
