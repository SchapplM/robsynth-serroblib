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
% Datum: 2020-01-03 11:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:25:58
% EndTime: 2020-01-03 11:26:01
% DurationCPUTime: 1.60s
% Computational Cost: add. (1502->255), mult. (3023->337), div. (0->0), fcn. (1920->10), ass. (0->159)
t90 = qJ(1) + pkin(7);
t84 = sin(t90);
t85 = cos(t90);
t166 = g(2) * t85 + g(3) * t84;
t96 = cos(pkin(7));
t80 = -pkin(1) * t96 - pkin(2);
t145 = qJDD(1) * t80;
t59 = qJDD(3) + t145;
t188 = t166 + t59;
t95 = cos(pkin(8));
t182 = pkin(3) * t95;
t93 = sin(pkin(8));
t53 = -pkin(6) * t93 - t182 + t80;
t100 = cos(qJ(4));
t160 = t100 * t95;
t94 = sin(pkin(7));
t75 = pkin(1) * t94 + qJ(3);
t57 = t75 * t160;
t98 = sin(qJ(4));
t26 = t98 * t53 + t57;
t187 = qJD(4) * t26;
t139 = qJD(1) * qJD(4);
t127 = t100 * t139;
t144 = qJDD(1) * t98;
t108 = t127 + t144;
t186 = t108 * t93;
t157 = pkin(1) * qJDD(1);
t62 = t75 * qJD(1);
t56 = qJD(1) * qJD(3) + qJDD(1) * t75;
t169 = t95 * t98;
t41 = -t100 * t85 - t84 * t169;
t159 = t84 * t100;
t43 = t85 * t169 - t159;
t185 = -g(2) * t41 - g(3) * t43;
t180 = g(1) * t98;
t184 = t93 * t180 + t185;
t156 = qJD(1) * t95;
t71 = -qJD(4) + t156;
t148 = qJD(1) * t100;
t129 = t93 * t148;
t38 = t53 * qJD(1) + qJD(3);
t47 = qJD(2) * t93 + t62 * t95;
t14 = t100 * t38 - t47 * t98;
t8 = -qJ(5) * t129 + t14;
t5 = -pkin(4) * t71 + t8;
t183 = -t8 + t5;
t181 = pkin(4) * t98;
t142 = t95 * qJDD(1);
t70 = -qJDD(4) + t142;
t177 = t70 * pkin(4);
t171 = t93 * t56;
t82 = t95 * qJDD(2);
t36 = -t82 + t171;
t176 = t36 * t93;
t83 = t95 * qJD(2);
t46 = t62 * t93 - t83;
t175 = t46 * t93;
t174 = t70 * t95;
t173 = t84 * t95;
t172 = t84 * t98;
t170 = t93 * (-qJ(5) - pkin(6));
t147 = qJD(4) * t100;
t168 = qJD(3) * t160 + t53 * t147;
t99 = sin(qJ(1));
t167 = t99 * pkin(1) + t84 * pkin(2);
t88 = t93 ^ 2;
t89 = t95 ^ 2;
t165 = t88 + t89;
t91 = t98 ^ 2;
t92 = t100 ^ 2;
t164 = -t91 - t92;
t163 = t91 - t92;
t162 = qJ(5) * t93;
t161 = t100 * t93;
t102 = qJD(1) ^ 2;
t158 = t88 * t102;
t155 = qJD(1) * t98;
t154 = qJD(3) * t98;
t15 = t100 * t47 + t38 * t98;
t153 = qJD(4) * t15;
t152 = qJD(4) * t47;
t151 = qJD(4) * t98;
t150 = qJD(5) * t93;
t28 = qJD(5) - t83 + (pkin(4) * t155 + t62) * t93;
t149 = qJD(5) + t28;
t143 = t93 * qJDD(1);
t141 = qJ(5) * qJDD(1);
t138 = qJD(1) * qJD(5);
t137 = qJDD(1) * t100;
t136 = t75 * t169;
t101 = cos(qJ(1));
t135 = t101 * pkin(1) + t85 * pkin(2) + t84 * qJ(3);
t134 = t98 * t158;
t133 = t93 * t155;
t132 = t71 * t151;
t131 = qJ(5) * t161;
t130 = t95 * t154;
t128 = -t36 * t95 - g(1);
t35 = t53 * qJDD(1) + qJDD(3);
t37 = qJDD(2) * t93 + t56 * t95;
t3 = t100 * t37 + t38 * t147 - t47 * t151 + t98 * t35;
t50 = (t75 + t181) * t93;
t126 = qJD(1) * t50 + t28;
t125 = t70 - t142;
t124 = t70 + t142;
t123 = t93 * pkin(4) * t127 + t143 * t181 + qJDD(5) - t82;
t122 = qJD(4) * t133;
t121 = t98 * t127;
t120 = -t85 * qJ(3) + t167;
t119 = g(2) * t43 - g(3) * t41;
t42 = t95 * t159 - t85 * t98;
t44 = t85 * t160 + t172;
t118 = -g(2) * t44 - g(3) * t42;
t116 = -g(2) * t84 + g(3) * t85;
t115 = -g(2) * t101 - g(3) * t99;
t9 = -qJ(5) * t133 + t15;
t114 = t100 * t9 - t5 * t98;
t113 = t100 * t5 + t9 * t98;
t112 = t37 * t95 + t176;
t111 = t47 * t95 + t175;
t110 = t100 * t15 - t14 * t98;
t19 = t123 + t171;
t55 = (pkin(4) * t147 + qJD(3)) * t93;
t109 = qJD(1) * t55 + qJDD(1) * t50 + t19;
t31 = t100 * t35;
t107 = t31 + qJ(5) * t122 + (-qJD(4) * t38 - t37) * t98;
t106 = t145 + t188;
t105 = -t71 ^ 2 - t158;
t104 = g(1) * t161 + g(2) * t42 - g(3) * t44 - t3;
t4 = -t98 * t37 - t153 + t31;
t81 = pkin(4) * t100 + pkin(3);
t68 = t93 * t137;
t61 = t100 * t134;
t58 = t95 * t122;
t54 = t71 * t129;
t52 = t163 * t158;
t51 = t164 * t143;
t49 = t100 * t53;
t40 = (qJDD(1) * t92 - 0.2e1 * t121) * t88;
t39 = (qJDD(1) * t91 + 0.2e1 * t121) * t88;
t27 = 0.2e1 * (-t98 * t137 + t163 * t139) * t88;
t25 = t49 - t136;
t24 = -t54 - t186;
t23 = t68 + (-qJD(4) - t71) * t133;
t22 = -t98 * t162 + t26;
t21 = -t130 - t187;
t20 = -qJD(4) * t136 + t168;
t18 = -t131 + t49 + (-t75 * t98 - pkin(4)) * t95;
t17 = t105 * t100 + t98 * t70;
t16 = -t100 * t70 + t105 * t98;
t13 = (t125 * t98 + (t71 - t156) * t147) * t93;
t12 = (t124 * t98 + (t71 + t156) * t147) * t93;
t11 = t58 + (-t124 * t100 + t132) * t93;
t10 = t58 + (t125 * t100 - t132) * t93;
t7 = -t130 - t100 * t150 + (-t57 + (-t53 + t162) * t98) * qJD(4);
t6 = -t98 * t150 + (-t131 - t136) * qJD(4) + t168;
t2 = (-t108 * qJ(5) - t98 * t138) * t93 + t3;
t1 = -t177 + (-t152 + (-t138 - t141) * t93) * t100 + t107;
t29 = [0, 0, 0, 0, 0, qJDD(1), t115, g(2) * t99 - g(3) * t101, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0.2e1 * t96 * t157 - t166, -0.2e1 * t94 * t157 - t116, 0, (t115 + (t94 ^ 2 + t96 ^ 2) * t157) * pkin(1), t88 * qJDD(1), 0.2e1 * t93 * t142, 0, t89 * qJDD(1), 0, 0, -t106 * t95, t106 * t93, t56 * t165 + t112 + t116, -g(2) * t135 - g(3) * t120 + qJD(3) * t111 + t112 * t75 + t59 * t80, t40, t27, t11, t39, t12, t174, -t21 * t71 - t25 * t70 - t4 * t95 + (t46 * t147 + t36 * t98) * t93 + (t75 * t144 + (t75 * t147 + t154) * qJD(1)) * t88 + t118, t20 * t71 + t26 * t70 + t3 * t95 + (-t88 * t62 - t175) * t151 + (t56 * t88 + t176) * t100 + t119, ((qJD(4) * t14 - qJDD(1) * t26 - t3 + (qJD(4) * t25 - t20) * qJD(1)) * t98 + (-t153 - qJDD(1) * t25 - t4 + (-t21 - t187) * qJD(1)) * t100 - t166) * t93, t3 * t26 + t15 * t20 + t4 * t25 + t14 * t21 - g(2) * (t85 * t182 + t135) - g(3) * (pkin(3) * t173 + t120) + (-pkin(6) * t166 + t46 * qJD(3) + t36 * t75) * t93, t40, t27, t11, t39, t12, t174, -t1 * t95 - t18 * t70 - t7 * t71 + (t109 * t98 + t126 * t147) * t93 + t118, t2 * t95 + t22 * t70 + t6 * t71 + (t100 * t109 - t126 * t151) * t93 + t119, ((qJD(4) * t5 - qJDD(1) * t22 - t2 + (qJD(4) * t18 - t6) * qJD(1)) * t98 + (-qJD(4) * t9 - qJDD(1) * t18 - t1 + (-qJD(4) * t22 - t7) * qJD(1)) * t100 - t166) * t93, t2 * t22 + t9 * t6 + t1 * t18 + t5 * t7 + t19 * t50 + t28 * t55 - g(2) * (pkin(4) * t172 + t135) - g(3) * (-t84 * t170 + t81 * t173 + t167) + (-g(2) * (t81 * t95 - t170) - g(3) * (-qJ(3) - t181)) * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2) - g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t37 * t93 + t128, 0, 0, 0, 0, 0, 0, t13, t10, 0, (t100 * t3 - t4 * t98 + (-t100 * t14 - t15 * t98) * qJD(4)) * t93 + t128, 0, 0, 0, 0, 0, 0, t13, t10, 0, -t19 * t95 - g(1) + (-qJD(4) * t113 - t1 * t98 + t100 * t2) * t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t142, t143, -t165 * t102, -qJD(1) * t111 + t188, 0, 0, 0, 0, 0, 0, t16, t17, t51, t4 * t100 + t3 * t98 + t110 * qJD(4) + (-t110 * t95 - t175) * qJD(1) + t166, 0, 0, 0, 0, 0, 0, t16, t17, t51, t1 * t100 + t2 * t98 + t114 * qJD(4) + (-t114 * t95 - t28 * t93) * qJD(1) + t166; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, -t52, t23, -t61, t24, -t70, -t129 * t46 - t15 * t71 + t184 + t4, t133 * t46 - t14 * t71 + t104, 0, 0, t61, -t52, t23, -t61, t24, -t70, -0.2e1 * t177 - t9 * t71 + (-pkin(4) * t134 - t152 + (-t149 * qJD(1) - t141) * t93) * t100 + t107 + t184, -t92 * pkin(4) * t158 - t8 * t71 + (t98 * t141 + (qJ(5) * t147 + t149 * t98) * qJD(1)) * t93 + t104, (-pkin(4) * t137 + (pkin(4) * qJD(4) - t183) * t155) * t93, t183 * t9 + (t1 + (-t28 * t148 + t180) * t93 + t185) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54 + t186, t68 + (-qJD(4) + t71) * t133, t164 * t158, g(1) * t95 + (qJD(1) * t113 + t116 + t56) * t93 + t123;];
tau_reg = t29;
