% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4RRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RRRP6_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP6_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP6_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP6_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP6_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP6_invdynJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:19:13
% EndTime: 2019-12-31 17:19:17
% DurationCPUTime: 1.99s
% Computational Cost: add. (1803->328), mult. (4270->427), div. (0->0), fcn. (2662->6), ass. (0->172)
t108 = sin(qJ(2));
t164 = qJD(1) * t108;
t216 = qJD(3) * t164 - qJDD(2);
t107 = sin(qJ(3));
t110 = cos(qJ(3));
t154 = t108 * qJDD(1);
t111 = cos(qJ(2));
t163 = qJD(1) * t111;
t28 = ((qJD(3) + t163) * qJD(2) + t154) * t107 + t216 * t110;
t155 = qJD(1) * qJD(2);
t99 = t111 * qJDD(1);
t215 = -t108 * t155 + t99;
t150 = t107 * t164;
t157 = t110 * qJD(2);
t67 = t150 - t157;
t86 = -qJD(3) + t163;
t194 = t67 * t86;
t145 = t111 * t155;
t27 = -qJD(3) * t157 + (-t145 - t154) * t110 + t216 * t107;
t214 = -t27 + t194;
t149 = t110 * t164;
t162 = qJD(2) * t107;
t69 = t149 + t162;
t192 = t69 * t86;
t213 = t28 - t192;
t109 = sin(qJ(1));
t112 = cos(qJ(1));
t132 = g(1) * t112 + g(2) * t109;
t122 = t132 * t108;
t175 = t107 * t108;
t169 = t110 * t112;
t174 = t107 * t111;
t51 = t109 * t174 + t169;
t168 = t111 * t112;
t171 = t109 * t110;
t53 = -t107 * t168 + t171;
t212 = -g(1) * t53 + g(2) * t51 + g(3) * t175;
t210 = -2 * pkin(1);
t209 = t69 ^ 2;
t208 = pkin(3) * t67;
t64 = qJDD(3) - t215;
t205 = t64 * pkin(3);
t204 = pkin(3) * t107;
t203 = pkin(5) * t107;
t202 = g(1) * t109;
t199 = g(2) * t112;
t198 = g(3) * t108;
t197 = g(3) * t111;
t136 = pkin(2) * t111 + pkin(6) * t108;
t74 = -pkin(1) - t136;
t59 = t74 * qJD(1);
t96 = pkin(5) * t163;
t80 = qJD(2) * pkin(6) + t96;
t33 = -t107 * t80 + t110 * t59;
t196 = t33 * t86;
t34 = t107 * t59 + t110 * t80;
t195 = t34 * t86;
t193 = t69 * t67;
t191 = qJ(4) + pkin(6);
t23 = -qJ(4) * t69 + t33;
t20 = -pkin(3) * t86 + t23;
t190 = -t23 + t20;
t142 = qJD(3) * t191;
t156 = t110 * qJD(4);
t172 = t108 * t110;
t135 = pkin(2) * t108 - pkin(6) * t111;
t71 = t135 * qJD(1);
t55 = t107 * t71;
t189 = -t107 * t142 + t156 - t55 - (-pkin(5) * t172 - qJ(4) * t174) * qJD(1);
t170 = t110 * t111;
t125 = pkin(3) * t108 - qJ(4) * t170;
t39 = pkin(5) * t150 + t110 * t71;
t188 = -t125 * qJD(1) - t107 * qJD(4) - t110 * t142 - t39;
t158 = qJD(3) * t110;
t72 = t135 * qJD(2);
t187 = t107 * t72 + t74 * t158;
t161 = qJD(2) * t108;
t186 = t110 * t72 + t161 * t203;
t185 = (g(1) * t169 + g(2) * t171) * t108;
t87 = pkin(5) * t170;
t43 = t107 * t74 + t87;
t184 = t107 * t67;
t183 = t110 * t69;
t182 = t27 * qJ(4);
t181 = t27 * t107;
t180 = t28 * qJ(4);
t179 = t28 * t110;
t178 = pkin(5) * qJDD(1);
t177 = qJD(3) * t69;
t176 = t191 * t108;
t173 = t107 * t112;
t167 = t112 * pkin(1) + t109 * pkin(5);
t104 = t108 ^ 2;
t105 = t111 ^ 2;
t166 = t104 - t105;
t165 = t104 + t105;
t160 = qJD(2) * t111;
t159 = qJD(3) * t107;
t114 = qJD(1) ^ 2;
t153 = t108 * t114 * t111;
t143 = -qJD(4) - t208;
t79 = -qJD(2) * pkin(2) + pkin(5) * t164;
t38 = -t143 + t79;
t152 = t38 * t158;
t151 = t86 * t164;
t148 = t108 * t158;
t94 = pkin(5) * t154;
t50 = -qJDD(2) * pkin(2) + pkin(5) * t145 + t94;
t19 = t28 * pkin(3) + qJDD(4) + t50;
t147 = -t19 - t197;
t35 = qJD(1) * t72 + t74 * qJDD(1);
t49 = t215 * pkin(5) + qJDD(2) * pkin(6);
t7 = t107 * t35 + t110 * t49 + t59 * t158 - t80 * t159;
t140 = -pkin(6) * qJD(3) * t86 + t50;
t139 = t108 * t145;
t138 = -g(1) * t51 - g(2) * t53;
t52 = -t109 * t170 + t173;
t54 = t107 * t109 + t110 * t168;
t137 = -g(1) * t52 - g(2) * t54;
t134 = pkin(5) * t67 + t107 * t79;
t133 = pkin(5) * t69 + t110 * t79;
t131 = -pkin(6) * t64 + qJD(3) * t79;
t24 = -qJ(4) * t67 + t34;
t130 = t107 * t24 + t110 * t20;
t129 = -t107 * t34 - t110 * t33;
t93 = pkin(3) * t110 + pkin(2);
t127 = t111 * t93 + t176;
t124 = t107 * t64 - t86 * t158;
t123 = t110 * t64 + t86 * t159;
t121 = -pkin(5) * qJDD(2) + t155 * t210;
t120 = pkin(1) * t114 + t132;
t113 = qJD(2) ^ 2;
t119 = pkin(5) * t113 + qJDD(1) * t210 + t199;
t118 = g(1) * t54 - g(2) * t52 + g(3) * t172 - t7;
t117 = -t132 * t111 - t198;
t8 = -qJD(3) * t34 - t107 * t49 + t110 * t35;
t116 = t8 + t212;
t102 = t112 * pkin(5);
t91 = t108 * t202;
t90 = g(3) * t174;
t76 = t191 * t110;
t75 = t191 * t107;
t73 = (pkin(5) + t204) * t108;
t66 = t110 * t74;
t63 = t67 ^ 2;
t61 = t163 * t204 + t96;
t42 = -pkin(5) * t174 + t66;
t41 = pkin(5) * t160 + (t107 * t160 + t148) * pkin(3);
t40 = -pkin(5) * t149 + t55;
t37 = -qJ(4) * t175 + t43;
t36 = -t111 * t64 - t86 * t161;
t32 = -qJ(4) * t172 + t66 + (-pkin(3) - t203) * t111;
t25 = -t63 + t209;
t22 = -t43 * qJD(3) + t186;
t21 = (-t108 * t157 - t111 * t159) * pkin(5) + t187;
t18 = -t192 - t28;
t17 = -t27 - t194;
t16 = (-t108 * t69 + t86 * t170) * qJD(1) + t124;
t15 = (t108 * t67 - t86 * t174) * qJD(1) + t123;
t14 = (-pkin(5) * qJD(2) - qJ(4) * qJD(3)) * t172 + (-qJD(4) * t108 + (-pkin(5) * qJD(3) - qJ(4) * qJD(2)) * t111) * t107 + t187;
t13 = -t108 * t156 + t125 * qJD(2) + (-t87 + (qJ(4) * t108 - t74) * t107) * qJD(3) + t186;
t12 = -t86 * t184 - t179;
t11 = -t86 * t183 - t181;
t10 = t67 * t148 + (t108 * t28 + t67 * t160) * t107;
t9 = t69 * t111 * t157 + (-t27 * t110 - t69 * t159) * t108;
t6 = (t86 * t162 + t28) * t111 + (-qJD(2) * t67 - t124) * t108;
t5 = (-t86 * t157 + t27) * t111 + (qJD(2) * t69 + t123) * t108;
t4 = -qJD(4) * t67 - t180 + t7;
t3 = -t69 * qJD(4) + t182 + t205 + t8;
t2 = -t213 * t107 + t214 * t110;
t1 = (-t107 * t69 - t110 * t67) * t160 + (t181 - t179 + (-t183 + t184) * qJD(3)) * t108;
t26 = [0, 0, 0, 0, 0, qJDD(1), -t199 + t202, t132, 0, 0, qJDD(1) * t104 + 0.2e1 * t139, 0.2e1 * t108 * t99 - 0.2e1 * t166 * t155, qJDD(2) * t108 + t111 * t113, qJDD(1) * t105 - 0.2e1 * t139, qJDD(2) * t111 - t108 * t113, 0, t121 * t108 + (-t119 + t202) * t111, t108 * t119 + t111 * t121 - t91, 0.2e1 * t165 * t178 - t132, -g(1) * (-pkin(1) * t109 + t102) - g(2) * t167 + (t165 * pkin(5) ^ 2 + (pkin(1) ^ 2)) * qJDD(1), t9, t1, t5, t10, t6, t36, -t22 * t86 + t42 * t64 + (qJD(2) * t134 - t8) * t111 + (pkin(5) * t28 + qJD(2) * t33 + t50 * t107 + t79 * t158) * t108 + t137, t21 * t86 - t43 * t64 + (qJD(2) * t133 + t7) * t111 + (-pkin(5) * t27 - qJD(2) * t34 + t50 * t110 - t79 * t159) * t108 + t138, -t21 * t67 - t22 * t69 + t42 * t27 - t43 * t28 + t91 + t129 * t160 + (-t199 - t107 * t7 - t110 * t8 + (t107 * t33 - t110 * t34) * qJD(3)) * t108, t7 * t43 + t34 * t21 + t8 * t42 + t33 * t22 - g(1) * t102 - g(2) * (t112 * t136 + t167) - t74 * t202 + (t108 * t50 + t79 * t160) * pkin(5), t9, t1, t5, t10, t6, t36, -t13 * t86 + t73 * t28 + t32 * t64 + t41 * t67 + (t38 * t162 - t3) * t111 + (qJD(2) * t20 + t19 * t107 + t152) * t108 + t137, t14 * t86 - t73 * t27 - t37 * t64 + t41 * t69 + (t38 * t157 + t4) * t111 + (-qJD(2) * t24 + t19 * t110 - t38 * t159) * t108 + t138, -t13 * t69 - t14 * t67 + t32 * t27 - t37 * t28 + t91 - t130 * t160 + (-t199 - t107 * t4 - t110 * t3 + (t107 * t20 - t110 * t24) * qJD(3)) * t108, t4 * t37 + t24 * t14 + t3 * t32 + t20 * t13 + t19 * t73 + t38 * t41 - g(1) * (pkin(3) * t173 + t102) - g(2) * (t112 * t176 + t93 * t168 + t167) + (-g(1) * (-pkin(1) - t127) - g(2) * t204) * t109; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t153, t166 * t114, t154, t153, t99, qJDD(2), t108 * t120 - t197 - t94, t198 + (t120 - t178) * t111, 0, 0, t11, t2, t16, t12, t15, t151, -pkin(2) * t28 + t39 * t86 + t131 * t107 + (-t140 - t197) * t110 + (-t108 * t33 - t111 * t134) * qJD(1) + t185, pkin(2) * t27 - t40 * t86 + t90 + t131 * t110 + (t108 * t34 - t111 * t133) * qJD(1) + (-t122 + t140) * t107, t39 * t69 + t40 * t67 + (t7 + t196 + (-t28 + t177) * pkin(6)) * t110 + (-t8 + t195 + (qJD(3) * t67 - t27) * pkin(6)) * t107 + t117, -t79 * t96 - t33 * t39 - t34 * t40 + (-t197 - t50 + t122) * pkin(2) + (qJD(3) * t129 - t8 * t107 + t7 * t110 + t117) * pkin(6), t11, t2, t16, t12, t15, t151, -t20 * t164 - t93 * t28 - t61 * t67 - t75 * t64 - t188 * t86 + t147 * t110 + (-t38 * t163 + (t38 + t208) * qJD(3)) * t107 + t185, t152 + t93 * t27 - t61 * t69 - t76 * t64 + t90 + t189 * t86 + (t108 * t24 - t38 * t170) * qJD(1) + (pkin(3) * t177 - t122 + t19) * t107, -t198 - t107 * t3 + t110 * t4 - t27 * t75 - t28 * t76 - t188 * t69 - t189 * t67 - t130 * qJD(3) + (qJD(1) * t130 - t132) * t111, t4 * t76 - t3 * t75 - t19 * t93 - g(3) * t127 + (pkin(3) * t159 - t61) * t38 + t189 * t24 + t188 * t20 + t132 * (t108 * t93 - t111 * t191); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t193, t25, t17, -t193, t18, t64, -t79 * t69 + t116 - t195, t67 * t79 + t118 - t196, 0, 0, t193, t25, t17, -t193, t18, t64, 0.2e1 * t205 + t182 - t24 * t86 + (t143 - t38) * t69 + t116, -t209 * pkin(3) + t180 - t23 * t86 + (qJD(4) + t38) * t67 + t118, t27 * pkin(3) - t190 * t67, t190 * t24 + (-t38 * t69 + t212 + t3) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t213, t214, -t63 - t209, t20 * t69 + t24 * t67 - t122 - t147;];
tau_reg = t26;
