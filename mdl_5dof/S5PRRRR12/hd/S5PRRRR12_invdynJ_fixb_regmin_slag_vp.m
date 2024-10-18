% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PRRRR12
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d2,d3,d4,d5,theta1]';
% 
% Output:
% tau_reg [5x17]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-28 18:09
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRRR12_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR12_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR12_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR12_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR12_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR12_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-28 18:08:20
% EndTime: 2024-09-28 18:08:23
% DurationCPUTime: 1.47s
% Computational Cost: add. (2422->283), mult. (4755->439), div. (0->0), fcn. (3914->26), ass. (0->194)
t116 = qJD(2) + qJD(3);
t107 = qJD(4) + t116;
t124 = cos(pkin(6));
t202 = t107 * t124;
t84 = qJD(5) + t202;
t188 = qJD(5) - t84;
t128 = sin(qJ(3));
t132 = cos(qJ(3));
t122 = sin(pkin(5));
t129 = sin(qJ(2));
t176 = qJDD(1) * t129;
t133 = cos(qJ(2));
t183 = qJD(2) * t133;
t141 = (qJD(1) * t183 + t176) * t122;
t185 = qJD(1) * t122;
t81 = qJD(2) * pkin(2) + t133 * t185;
t210 = qJD(3) * t81;
t167 = t129 * t185;
t199 = t122 * t133;
t91 = qJDD(1) * t199;
t58 = qJDD(2) * pkin(2) - qJD(2) * t167 + t91;
t245 = t128 * t58 + (t141 + t210) * t132;
t224 = g(3) * t122;
t121 = sin(pkin(6));
t203 = t107 * t121;
t127 = sin(qJ(4));
t131 = cos(qJ(4));
t53 = t128 * t81 + t132 * t167;
t211 = t131 * t53;
t52 = -t128 * t167 + t132 * t81;
t39 = pkin(3) * t116 + t52;
t21 = -t127 * t39 - t211;
t16 = pkin(10) * t203 - t21;
t243 = t188 * t16;
t67 = (-t128 * t129 + t132 * t133) * t122;
t189 = t129 * t132;
t151 = t128 * t133 + t189;
t68 = t151 * t122;
t119 = qJ(2) + qJ(3);
t110 = sin(t119);
t120 = sin(pkin(11));
t123 = cos(pkin(11));
t108 = pkin(5) + t119;
t93 = cos(t108) / 0.2e1;
t109 = pkin(5) - t119;
t98 = cos(t109);
t79 = t98 / 0.2e1 + t93;
t92 = sin(t109) / 0.2e1;
t97 = sin(t108);
t242 = -g(3) * (t97 / 0.2e1 + t92) - g(2) * (-t110 * t120 + t123 * t79) - g(1) * (-t110 * t123 - t120 * t79);
t113 = qJ(4) + t119;
t101 = sin(t113);
t99 = qJ(4) + t108;
t90 = cos(t99) / 0.2e1;
t100 = -qJ(4) + t109;
t95 = cos(t100);
t73 = t95 / 0.2e1 + t90;
t89 = sin(t100) / 0.2e1;
t94 = sin(t99);
t241 = -g(3) * (t94 / 0.2e1 + t89) - g(2) * (-t101 * t120 + t123 * t73) - g(1) * (-t101 * t123 - t120 * t73);
t213 = t127 * t53;
t20 = t131 * t39 - t213;
t17 = pkin(4) * t107 + t20;
t125 = cos(pkin(5));
t184 = qJD(1) * t125;
t13 = -t121 * t17 + t124 * t184;
t149 = t121 * t184 + t124 * t17;
t114 = qJDD(2) + qJDD(3);
t106 = qJDD(4) + t114;
t207 = t106 * t121;
t54 = t132 * t58;
t135 = t54 + (-t128 * t176 + (-qJD(3) * t189 - t128 * t183) * qJD(1)) * t122 - t128 * t210;
t14 = pkin(3) * t114 + t135;
t156 = qJD(3) * t167;
t80 = t128 * t156;
t15 = t245 - t80;
t236 = -(qJD(4) * t39 + t15) * t131 - t127 * t14;
t182 = qJD(4) * t127;
t44 = t53 * t182;
t6 = pkin(10) * t207 - t236 - t44;
t240 = t101 * t224 - t13 * t203 - t188 * t149 - t6;
t31 = -t127 * t68 + t131 * t67;
t152 = t121 * t125 + t124 * t31;
t24 = -t121 * t31 + t124 * t125;
t239 = -t152 * t84 + t24 * t203;
t226 = pkin(3) * t131;
t103 = pkin(4) + t226;
t22 = -t127 * t52 - t211;
t155 = pkin(3) * t182 + t22;
t206 = t106 * t124;
t82 = qJDD(5) + t206;
t238 = -t103 * t82 + t155 * t84;
t104 = pkin(2) * t132 + pkin(3);
t181 = qJD(4) * t131;
t190 = t128 * t131;
t55 = t151 * t185;
t56 = qJD(1) * t67;
t218 = t127 * t56 + t131 * t55 - t104 * t182 + (-t128 * t181 + (-t127 * t132 - t190) * qJD(3)) * pkin(2);
t191 = t127 * t128;
t157 = -pkin(2) * t191 + t131 * t104;
t69 = pkin(4) + t157;
t237 = t218 * t84 + t69 * t82;
t34 = t116 * t67;
t35 = t116 * t68;
t8 = t31 * qJD(4) - t127 * t35 + t131 * t34;
t229 = t8 * t84;
t130 = cos(qJ(5));
t175 = t125 * qJDD(1);
t162 = t121 * t175;
t197 = t124 * t130;
t161 = -t127 * t15 + t131 * t14;
t140 = t21 * qJD(4) + t161;
t225 = pkin(4) * t106;
t7 = t140 + t225;
t228 = t130 * t162 + t7 * t197;
t227 = pkin(2) * t114;
t112 = t121 * pkin(10);
t4 = -t121 * t7 + t124 * t175;
t223 = t4 * t121;
t217 = pkin(2) * t190 + t127 * t104;
t57 = t112 + t217;
t222 = t57 * t82;
t85 = pkin(3) * t127 + t112;
t220 = t85 * t82;
t219 = -t127 * t55 + t131 * t56 - t104 * t181 - (-t128 * t182 + (t131 * t132 - t191) * qJD(3)) * pkin(2);
t215 = t107 * t21;
t74 = t82 * t124;
t209 = qJD(5) * t84;
t115 = t121 ^ 2;
t208 = t107 ^ 2 * t115;
t126 = sin(qJ(5));
t205 = t106 * t126;
t204 = t107 * t115;
t201 = t107 * t126;
t200 = t121 * t126;
t198 = t124 * t126;
t196 = t125 * t126;
t195 = t125 * t129;
t194 = t125 * t130;
t193 = t125 * t133;
t192 = t126 * t130;
t187 = qJDD(1) - g(3);
t117 = t126 ^ 2;
t186 = -t130 ^ 2 + t117;
t180 = qJD(5) * t107;
t179 = qJD(5) * t124;
t178 = qJD(5) * t126;
t177 = qJD(5) * t130;
t174 = pkin(4) * t204;
t172 = qJD(5) * t121 * t13;
t171 = t84 * t178;
t170 = t69 * t177;
t169 = t124 * t196;
t168 = t124 * t194;
t166 = t103 * t177;
t165 = t115 * t180;
t164 = -pkin(3) * t107 - t39;
t163 = t218 * t107;
t159 = t82 + t206;
t158 = t84 + t202;
t23 = t131 * t52 - t213;
t154 = -pkin(3) * t181 + t23;
t153 = g(1) * t120 - g(2) * t123;
t32 = t127 * t67 + t131 * t68;
t148 = pkin(4) * t198 + t130 * t112;
t9 = -t32 * qJD(4) - t127 * t34 - t131 * t35;
t146 = t9 * t204 - t24 * t207;
t102 = cos(t113);
t71 = t89 - t94 / 0.2e1;
t145 = -g(1) * (-t102 * t123 - t120 * t71) - g(2) * (-t102 * t120 + t123 * t71) - g(3) * (t90 - t95 / 0.2e1) + t44;
t111 = cos(t119);
t77 = t92 - t97 / 0.2e1;
t144 = -g(1) * (-t111 * t123 - t120 * t77) - g(2) * (-t111 * t120 + t123 * t77) - g(3) * (t93 - t98 / 0.2e1) + t80;
t59 = t120 * t169 - t123 * t130;
t61 = t120 * t130 + t123 * t169;
t63 = t120 * t198 - t123 * t194;
t65 = t120 * t194 + t123 * t198;
t143 = -g(1) * (t101 * t59 - t102 * t65) - g(2) * (-t101 * t61 - t102 * t63) - (-t101 * t198 + t102 * t130) * t224 + (-t126 * t6 + (-t149 * t126 - t130 * t16) * qJD(5) + t228) * t124 + t126 * t172;
t60 = -t120 * t168 - t123 * t126;
t62 = -t120 * t126 + t123 * t168;
t64 = -t120 * t197 - t123 * t196;
t66 = t120 * t196 - t123 * t197;
t142 = -g(1) * (-t101 * t60 + t102 * t66) - g(2) * (-t101 * t62 + t102 * t64) - (-t101 * t197 - t102 * t126) * t224 + t130 * t172 + t4 * t200;
t139 = -t141 + (-pkin(2) * t116 - t81) * qJD(3);
t138 = t145 + t236;
t136 = t140 + t241;
t134 = qJD(2) ^ 2;
t86 = t121 * t194;
t47 = (t106 * t117 + 0.2e1 * t177 * t201) * t115;
t33 = 0.2e1 * (t106 * t192 - t186 * t180) * t115;
t19 = (t159 * t126 + t158 * t177) * t121;
t18 = (t159 * t130 - t158 * t178) * t121;
t2 = t130 * t6 + (t124 * t7 + t162) * t126 + (-t126 * t16 + t149 * t130) * qJD(5);
t1 = [t187, 0, (qJDD(2) * t133 - t129 * t134) * t122, (-qJDD(2) * t129 - t133 * t134) * t122, 0, t114 * t67 - t116 * t35, -t114 * t68 - t116 * t34, 0, t106 * t31 + t107 * t9, -t106 * t32 - t107 * t8, 0, 0, 0, 0, 0, -t126 * t229 + (-t126 * t32 + t86) * t82 + ((t31 * t82 + t84 * t9) * t124 + t146) * t130 + (-t130 * t32 * t84 + t239 * t126) * qJD(5), (t239 * qJD(5) - t32 * t82 - t229) * t130 + (-(-qJD(5) * t32 + t124 * t9) * t84 - t152 * t82 - t146) * t126; 0, qJDD(2), t91 - g(1) * (-t120 * t193 - t123 * t129) - g(2) * (-t120 * t129 + t123 * t193) - g(3) * t199, -g(1) * (t120 * t195 - t123 * t133) - g(2) * (-t120 * t133 - t123 * t195) - t187 * t129 * t122, t114, t116 * t55 + t54 + (-t156 + t227) * t132 + t139 * t128 + t242, t116 * t56 + (-t58 - t227) * t128 + t139 * t132 + t144, t106, t157 * t106 + t136 + t163, -t217 * t106 + t219 * t107 + t138, t47, t33, t19, t18, t74, (-t69 * t165 - t222 + (-t69 * t179 + t219) * t84) * t126 + (-t57 * t209 - t223 + t237 * t124 + (t106 * t69 + t163) * t115) * t130 + t143, t57 * t171 + (t219 * t84 - t222) * t130 + (-t237 * t126 - t84 * t170 - t2) * t124 + (-t69 * t205 + (-t218 * t126 - t170) * t107) * t115 + t142; 0, 0, 0, 0, t114, t116 * t53 + t135 + t242, t116 * t52 + t144 - t245, t106, t106 * t226 - t107 * t22 + (t164 * t127 - t211) * qJD(4) + t161 + t241, t107 * t23 + (-pkin(3) * t106 - t14) * t127 + (t164 * qJD(4) - t15) * t131 + t145, t47, t33, t19, t18, t74, (-t103 * t165 - t220 + (-t103 * t179 + t154) * t84) * t126 + (-t85 * t209 - t223 - t238 * t124 + (t103 * t106 - t155 * t107) * t115) * t130 + t143, t85 * t171 + (t154 * t84 - t220) * t130 + (t238 * t126 - t84 * t166 - t2) * t124 + (-t103 * t205 + (t155 * t126 - t166) * t107) * t115 + t142; 0, 0, 0, 0, 0, 0, 0, t106, t136 - t215, t107 * t20 + t138, t47, t33, t19, t18, t74, (-t82 * t112 + t20 * t84) * t126 + (-t223 + (-t124 * t84 - t204) * t21 + (t106 * t115 + t74) * pkin(4)) * t130 + (-t126 * t174 - t148 * t84) * qJD(5) + t143, -t148 * t82 - t2 * t124 + (t130 * t20 + t21 * t198) * t84 + (t215 - t225) * t126 * t115 + (-(pkin(4) * t197 - pkin(10) * t200) * t84 - t130 * t174) * qJD(5) + t142; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t192 * t208, t186 * t208, (t188 * t130 * t107 + t205) * t121, (t106 * t130 - t188 * t201) * t121, t82, -g(1) * (t101 * t66 + t102 * t60) - g(2) * (t101 * t64 + t102 * t62) - g(3) * t86 + (-t243 + (-g(3) * t102 * t124 - t153 * t121) * t122) * t130 + t240 * t126 + t228, -g(1) * (t101 * t65 + t102 * t59) - g(2) * (t101 * t63 - t102 * t61) + (t243 + (t102 * t224 - t7) * t124 + (t153 * t122 - t187 * t125) * t121) * t126 + t240 * t130;];
tau_reg = t1;
