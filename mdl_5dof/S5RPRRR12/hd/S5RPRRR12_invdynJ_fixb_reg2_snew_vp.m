% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPRRR12
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPRRR12_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR12_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR12_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR12_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR12_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR12_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:13:09
% EndTime: 2019-12-31 19:13:16
% DurationCPUTime: 2.64s
% Computational Cost: add. (7569->286), mult. (15096->383), div. (0->0), fcn. (10017->8), ass. (0->193)
t159 = sin(qJ(4));
t153 = qJDD(3) + qJDD(4);
t160 = sin(qJ(3));
t163 = cos(qJ(4));
t164 = cos(qJ(3));
t205 = t164 * t159;
t132 = (t160 * t163 + t205) * qJD(1);
t202 = qJD(1) * t164;
t134 = -t159 * t160 * qJD(1) + t163 * t202;
t213 = t134 * t132;
t235 = t153 - t213;
t238 = t159 * t235;
t237 = t163 * t235;
t154 = qJD(3) + qJD(4);
t212 = t154 * t132;
t199 = qJD(1) * qJD(3);
t189 = t164 * t199;
t197 = t160 * qJDD(1);
t140 = -t189 - t197;
t151 = t164 * qJDD(1);
t188 = t160 * t199;
t141 = t151 - t188;
t99 = -t132 * qJD(4) + t159 * t140 + t163 * t141;
t236 = t99 - t212;
t158 = sin(qJ(5));
t162 = cos(qJ(5));
t115 = t134 * t158 - t162 * t154;
t117 = t134 * t162 + t154 * t158;
t96 = t117 * t115;
t185 = -t163 * t140 + t141 * t159;
t98 = -qJD(4) * t134 - t185;
t97 = qJDD(5) - t98;
t232 = -t96 + t97;
t234 = t158 * t232;
t233 = t162 * t232;
t198 = qJD(2) * qJD(1);
t152 = 0.2e1 * t198;
t156 = t160 ^ 2;
t167 = qJD(1) ^ 2;
t177 = qJD(3) * pkin(3) - pkin(7) * t202;
t155 = qJDD(1) * qJ(2);
t161 = sin(qJ(1));
t165 = cos(qJ(1));
t183 = g(1) * t165 + g(2) * t161;
t178 = -t155 + t183;
t228 = pkin(6) + pkin(1);
t230 = -pkin(3) * t140 - (pkin(7) * t156 + t228) * t167 + t177 * t202 - t178;
t35 = t152 - t236 * pkin(8) + (t134 * t154 - t98) * pkin(4) + t230;
t108 = pkin(4) * t132 - pkin(8) * t134;
t229 = t154 ^ 2;
t182 = g(1) * t161 - g(2) * t165;
t176 = qJDD(2) - t182;
t219 = qJ(2) * t167;
t174 = t176 - t219;
t171 = -t228 * qJDD(1) + t174;
t112 = t164 * g(3) - t160 * t171;
t209 = t156 * t167;
t93 = -pkin(3) * t209 + t140 * pkin(7) - qJD(3) * t177 - t112;
t222 = t163 * t93;
t170 = t164 * t171;
t169 = -t141 * pkin(7) + t170;
t204 = t164 * t167;
t231 = qJDD(3) * pkin(3) + t169 + (-pkin(3) * t204 - pkin(7) * t199 + g(3)) * t160;
t69 = t231 * t159 + t222;
t40 = -t229 * pkin(4) + t153 * pkin(8) - t132 * t108 + t69;
t16 = t158 * t40 - t162 * t35;
t17 = t158 * t35 + t162 * t40;
t7 = t158 * t16 + t162 * t17;
t129 = qJD(5) + t132;
t186 = -t162 * t153 + t158 * t99;
t61 = (qJD(5) - t129) * t117 + t186;
t113 = t115 ^ 2;
t114 = t117 ^ 2;
t127 = t129 ^ 2;
t130 = t132 ^ 2;
t131 = t134 ^ 2;
t68 = t159 * t93 - t163 * t231;
t39 = -t153 * pkin(4) - t229 * pkin(8) + t108 * t134 + t68;
t227 = -pkin(4) * t39 + pkin(8) * t7;
t226 = t160 * g(3);
t36 = t158 * t39;
t71 = t96 + t97;
t225 = t158 * t71;
t194 = -0.2e1 * t198;
t94 = t194 - t230;
t224 = t159 * t94;
t37 = t162 * t39;
t223 = t162 * t71;
t221 = t163 * t94;
t32 = t159 * t69 - t163 * t68;
t220 = t164 * t32;
t218 = qJDD(1) * pkin(1);
t106 = t153 + t213;
t217 = t106 * t159;
t216 = t106 * t163;
t215 = t129 * t158;
t214 = t129 * t162;
t211 = t154 * t159;
t210 = t154 * t163;
t157 = t164 ^ 2;
t208 = t157 * t167;
t191 = t160 * t204;
t207 = t160 * (qJDD(3) + t191);
t145 = qJDD(3) - t191;
t206 = t164 * t145;
t203 = t156 + t157;
t200 = qJD(5) + t129;
t92 = -t114 - t127;
t46 = -t158 * t92 - t223;
t180 = -t153 * t158 - t162 * t99;
t66 = t200 * t115 + t180;
t196 = pkin(4) * t66 + pkin(8) * t46 + t36;
t81 = -t127 - t113;
t43 = t162 * t81 - t234;
t62 = -t200 * t117 - t186;
t195 = pkin(4) * t62 + pkin(8) * t43 - t37;
t193 = t159 * t96;
t192 = t163 * t96;
t190 = -pkin(4) * t163 - pkin(3);
t33 = t159 * t68 + t163 * t69;
t103 = t129 * t115;
t76 = -qJD(5) * t115 - t180;
t65 = t103 + t76;
t31 = t158 * t65 - t162 * t61;
t78 = t113 + t114;
t187 = pkin(4) * t78 + pkin(8) * t31 + t7;
t3 = t159 * t7 - t163 * t39;
t1 = t160 * (t159 * t39 + t163 * t7) + t164 * t3;
t181 = t158 * t17 - t16 * t162;
t111 = t170 + t226;
t90 = t164 * t111 - t160 * t112;
t175 = (-qJD(4) + t154) * t134 - t185;
t166 = qJD(3) ^ 2;
t143 = t203 * qJDD(1);
t142 = t151 - 0.2e1 * t188;
t139 = 0.2e1 * t189 + t197;
t128 = -t174 + t218;
t124 = -t131 + t229;
t123 = t130 - t229;
t122 = t228 * t167 + t178 + t194;
t120 = -t131 - t229;
t119 = -t207 + t164 * (-t166 - t208);
t118 = t160 * (-t166 - t209) + t206;
t109 = t131 - t130;
t104 = -t229 - t130;
t102 = -t114 + t127;
t101 = t113 - t127;
t100 = -t130 - t131;
t95 = t114 - t113;
t89 = -t120 * t159 - t216;
t88 = t120 * t163 - t217;
t87 = t99 + t212;
t82 = (qJD(4) + t154) * t134 + t185;
t80 = t104 * t163 - t238;
t79 = t104 * t159 + t237;
t75 = -qJD(5) * t117 - t186;
t74 = (-t115 * t162 + t117 * t158) * t129;
t73 = (-t115 * t158 - t117 * t162) * t129;
t64 = -t103 + t76;
t58 = -t117 * t215 + t162 * t76;
t57 = t117 * t214 + t158 * t76;
t56 = t115 * t214 - t158 * t75;
t55 = t115 * t215 + t162 * t75;
t54 = t160 * t89 + t164 * t88;
t53 = t159 * t87 + t163 * t175;
t52 = t159 * t175 - t163 * t87;
t51 = t101 * t162 - t225;
t50 = -t102 * t158 + t233;
t49 = t101 * t158 + t223;
t48 = t102 * t162 + t234;
t47 = t160 * t80 + t164 * t79;
t45 = t162 * t92 - t225;
t42 = t158 * t81 + t233;
t30 = -t158 * t64 + t162 * t62;
t29 = -t158 * t61 - t162 * t65;
t28 = t158 * t62 + t162 * t64;
t26 = t160 * t53 + t164 * t52;
t25 = -t159 * t66 + t163 * t46;
t24 = t159 * t46 + t163 * t66;
t23 = -t159 * t62 + t163 * t43;
t22 = t159 * t43 + t163 * t62;
t21 = -t159 * t78 + t163 * t31;
t20 = t159 * t31 + t163 * t78;
t19 = -pkin(8) * t45 + t37;
t18 = -pkin(8) * t42 + t36;
t13 = t160 * t33 + t220;
t12 = -pkin(4) * t45 + t17;
t11 = -pkin(4) * t42 + t16;
t10 = t160 * t25 + t164 * t24;
t9 = t160 * t23 + t164 * t22;
t8 = t160 * t21 + t164 * t20;
t2 = -pkin(8) * t29 - t181;
t4 = [0, 0, 0, 0, 0, qJDD(1), t182, t183, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t176 - 0.2e1 * t218, t152 + 0.2e1 * t155 - t183, pkin(1) * t128 + qJ(2) * (-pkin(1) * t167 + t152 - t178), (t141 - t188) * t164, -t139 * t164 - t142 * t160, t206 - t160 * (t166 - t208), (-t140 + t189) * t160, t164 * (-t166 + t209) - t207, 0, qJ(2) * t139 - t228 * t118 - t122 * t160, qJ(2) * t142 - t228 * t119 - t122 * t164, t228 * t143 - t203 * t219 - t90, -qJ(2) * t122 - t228 * t90, t164 * (-t134 * t211 + t163 * t99) - t160 * (t134 * t210 + t159 * t99), t164 * (-t159 * t236 - t163 * t82) - t160 * (-t159 * t82 + t163 * t236), t164 * (-t124 * t159 + t237) - t160 * (t124 * t163 + t238), t164 * (t132 * t210 - t159 * t98) - t160 * (t132 * t211 + t163 * t98), t164 * (t123 * t163 - t217) - t160 * (t123 * t159 + t216), (t164 * (-t132 * t163 + t134 * t159) - t160 * (-t132 * t159 - t134 * t163)) * t154, t164 * (-pkin(7) * t79 - t224) - t160 * (-pkin(3) * t82 + pkin(7) * t80 + t221) + qJ(2) * t82 - t228 * t47, t164 * (-pkin(7) * t88 - t221) - t160 * (-pkin(3) * t236 + pkin(7) * t89 - t224) + qJ(2) * t236 - t228 * t54, t164 * (-pkin(7) * t52 - t32) - t160 * (-pkin(3) * t100 + pkin(7) * t53 + t33) + qJ(2) * t100 - t228 * t26, -pkin(7) * t220 - t160 * (pkin(3) * t94 + pkin(7) * t33) - qJ(2) * t94 - t228 * t13, t164 * (t163 * t58 + t193) - t160 * (t159 * t58 - t192), t164 * (t159 * t95 + t163 * t30) - t160 * (t159 * t30 - t163 * t95), t164 * (t159 * t65 + t163 * t50) - t160 * (t159 * t50 - t163 * t65), t164 * (t163 * t56 - t193) - t160 * (t159 * t56 + t192), t164 * (-t159 * t61 + t163 * t51) - t160 * (t159 * t51 + t163 * t61), t164 * (t159 * t97 + t163 * t74) - t160 * (t159 * t74 - t163 * t97), t164 * (-pkin(7) * t22 - t11 * t159 + t163 * t18) - t160 * (-pkin(3) * t42 + pkin(7) * t23 + t11 * t163 + t159 * t18) + qJ(2) * t42 - t228 * t9, t164 * (-pkin(7) * t24 - t12 * t159 + t163 * t19) - t160 * (-pkin(3) * t45 + pkin(7) * t25 + t12 * t163 + t159 * t19) + qJ(2) * t45 - t228 * t10, t164 * (-pkin(7) * t20 + t163 * t2) - t160 * (pkin(7) * t21 + t159 * t2) - t228 * t8 + (pkin(4) * t205 - t160 * t190 + qJ(2)) * t29, (t164 * (pkin(4) * t159 - pkin(8) * t163) - t160 * (-pkin(8) * t159 + t190) + qJ(2)) * t181 + (-t228 - pkin(7)) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t167, -t128, 0, 0, 0, 0, 0, 0, t118, t119, -t143, t90, 0, 0, 0, 0, 0, 0, t47, t54, t26, t13, 0, 0, 0, 0, 0, 0, t9, t10, t8, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t191, (-t156 + t157) * t167, t151, -t191, -t197, qJDD(3), t111, t112, 0, 0, t213, t109, t87, -t213, t175, t153, pkin(3) * t79 - t68, -t222 - t159 * (-pkin(7) * t188 + t169 + t226) + (-t145 * t159 + t88) * pkin(3), pkin(3) * t52, pkin(3) * t32, t57, t28, t48, t55, t49, t73, pkin(3) * t22 + t195, pkin(3) * t24 + t196, pkin(3) * t20 + t187, pkin(3) * t3 + t227; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t213, t109, t87, -t213, t175, t153, -t68, -t69, 0, 0, t57, t28, t48, t55, t49, t73, t195, t196, t187, t227; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96, t95, t65, -t96, -t61, t97, -t16, -t17, 0, 0;];
tauJ_reg = t4;
