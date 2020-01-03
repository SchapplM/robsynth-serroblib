% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S4RRRR6
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,d2,d3,d4]';
% 
% Output:
% tauJ_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S4RRRR6_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR6_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR6_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR6_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR6_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4RRRR6_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:30:52
% EndTime: 2019-12-31 17:30:59
% DurationCPUTime: 2.84s
% Computational Cost: add. (10559->329), mult. (23195->489), div. (0->0), fcn. (17838->10), ass. (0->216)
t168 = sin(qJ(4));
t167 = cos(pkin(4));
t161 = t167 * qJD(1) + qJD(2);
t169 = sin(qJ(3));
t172 = cos(qJ(3));
t166 = sin(pkin(4));
t170 = sin(qJ(2));
t213 = qJD(1) * t170;
t198 = t166 * t213;
t136 = t169 * t161 + t172 * t198;
t173 = cos(qJ(2));
t214 = qJD(1) * t166;
t196 = qJD(2) * t214;
t207 = qJDD(1) * t166;
t143 = t170 * t207 + t173 * t196;
t160 = t167 * qJDD(1) + qJDD(2);
t193 = t169 * t143 - t172 * t160;
t107 = -t136 * qJD(3) - t193;
t106 = qJDD(4) - t107;
t212 = qJD(1) * t173;
t197 = t166 * t212;
t153 = -qJD(3) + t197;
t171 = cos(qJ(4));
t117 = t168 * t136 + t171 * t153;
t119 = t171 * t136 - t168 * t153;
t97 = t119 * t117;
t246 = t106 - t97;
t248 = t168 * t246;
t247 = t171 * t246;
t163 = t166 ^ 2;
t174 = qJD(1) ^ 2;
t245 = t163 * (qJD(1) * t161 - t167 * t174);
t154 = t170 * t196;
t206 = qJDD(1) * t173;
t190 = t166 * t206 - t154;
t138 = -qJDD(3) + t190;
t134 = -t172 * t161 + t169 * t198;
t227 = t136 * t134;
t179 = -t138 - t227;
t244 = t169 * t179;
t243 = t172 * t179;
t187 = -t172 * t143 - t169 * t160;
t108 = -t134 * qJD(3) - t187;
t125 = t134 * t153;
t92 = t108 + t125;
t239 = sin(qJ(1));
t240 = cos(qJ(1));
t182 = t240 * g(1) + t239 * g(2);
t139 = -t174 * pkin(1) + pkin(6) * t207 - t182;
t181 = t239 * g(1) - t240 * g(2);
t221 = t166 * t174;
t178 = qJDD(1) * pkin(1) + pkin(6) * t221 + t181;
t177 = t167 * t178;
t194 = t170 * t139 - t173 * t177;
t222 = t166 * t173;
t110 = g(3) * t222 + t194;
t223 = t166 * t170;
t176 = -g(3) * t223 + t170 * t177;
t111 = t173 * t139 + t176;
t242 = t170 * t110 + t173 * t111;
t130 = qJD(4) + t134;
t195 = t168 * t108 + t171 * t138;
t58 = (qJD(4) - t130) * t119 + t195;
t89 = (qJD(3) + t153) * t136 + t193;
t115 = t117 ^ 2;
t116 = t119 ^ 2;
t128 = t130 ^ 2;
t132 = t134 ^ 2;
t133 = t136 ^ 2;
t241 = t153 ^ 2;
t159 = t161 ^ 2;
t238 = pkin(3) * t169;
t237 = pkin(6) * t166;
t236 = t167 * g(3);
t175 = t154 * pkin(2) - t143 * pkin(7) - t236 + (-t161 * pkin(7) * t212 + (t161 * t213 - t206) * pkin(2) - t178) * t166;
t191 = -t173 * pkin(2) - t170 * pkin(7);
t142 = t191 * t214;
t86 = t160 * pkin(7) - t159 * pkin(2) + (t142 * t214 + t139) * t173 + t176;
t52 = t169 * t175 + t172 * t86;
t112 = t134 * pkin(3) - t136 * pkin(8);
t51 = t169 * t86 - t172 * t175;
t36 = t138 * pkin(3) - t241 * pkin(8) + t136 * t112 + t51;
t235 = t168 * t36;
t68 = t106 + t97;
t234 = t168 * t68;
t85 = -t160 * pkin(2) - t159 * pkin(7) + (g(3) * t173 + t142 * t213) * t166 + t194;
t233 = t169 * t85;
t232 = t171 * t36;
t231 = t171 * t68;
t230 = t172 * t85;
t229 = t130 * t168;
t228 = t130 * t171;
t226 = t153 * t169;
t225 = t153 * t172;
t224 = t163 * t174;
t102 = t138 - t227;
t220 = t169 * t102;
t152 = t173 * t170 * t224;
t140 = t152 + t160;
t218 = t170 * t140;
t217 = t172 * t102;
t141 = -t152 + t160;
t215 = t173 * t141;
t211 = qJD(3) - t153;
t208 = qJD(4) + t130;
t164 = t170 ^ 2;
t205 = t164 * t224;
t165 = t173 ^ 2;
t204 = t165 * t224;
t203 = t169 * t97;
t202 = t172 * t97;
t201 = t173 * t227;
t147 = t161 * t197;
t200 = t147 + t143;
t199 = -pkin(3) * t172 - pkin(2);
t37 = -t241 * pkin(3) - t138 * pkin(8) - t134 * t112 + t52;
t39 = -t92 * pkin(8) + (-t153 * t136 - t107) * pkin(3) + t85;
t14 = t168 * t37 - t171 * t39;
t15 = t168 * t39 + t171 * t37;
t9 = t168 * t14 + t171 * t15;
t27 = t169 * t51 + t172 * t52;
t192 = -pkin(3) * t36 + pkin(8) * t9;
t8 = -t171 * t14 + t168 * t15;
t189 = t169 * t52 - t172 * t51;
t188 = -t171 * t108 + t168 * t138;
t185 = -pkin(1) + t191;
t95 = -t116 - t128;
t45 = -t168 * t95 - t231;
t63 = t208 * t117 + t188;
t184 = pkin(3) * t63 + pkin(8) * t45 + t235;
t81 = -t128 - t115;
t43 = t171 * t81 - t248;
t59 = -t208 * t119 - t195;
t183 = pkin(3) * t59 + pkin(8) * t43 - t232;
t105 = t130 * t117;
t75 = -t117 * qJD(4) - t188;
t62 = t105 + t75;
t33 = t168 * t62 - t171 * t58;
t76 = t115 + t116;
t180 = pkin(3) * t76 + pkin(8) * t33 + t9;
t146 = t161 * t198;
t145 = (t164 - t165) * t224;
t144 = -t159 - t204;
t129 = -t205 - t159;
t126 = t166 * t178 + t236;
t124 = -t146 + t190;
t123 = t146 + t190;
t122 = -t147 + t143;
t121 = -t133 + t241;
t120 = t132 - t241;
t114 = -t133 - t241;
t113 = t133 - t132;
t109 = -t241 - t132;
t101 = -t116 + t128;
t100 = t115 - t128;
t99 = t132 + t133;
t98 = (t134 * t169 + t136 * t172) * t153;
t96 = t116 - t115;
t94 = t211 * t134 + t187;
t93 = t108 - t125;
t90 = -t211 * t136 - t193;
t88 = t169 * t108 - t136 * t225;
t87 = t172 * t107 - t134 * t226;
t80 = t169 * t120 - t217;
t79 = t172 * t121 + t244;
t78 = -t169 * t114 + t217;
t77 = t172 * t114 + t220;
t74 = -t119 * qJD(4) - t195;
t73 = t172 * t109 - t244;
t72 = t169 * t109 + t243;
t71 = (-t117 * t171 + t119 * t168) * t130;
t70 = (-t117 * t168 - t119 * t171) * t130;
t66 = t169 * t93 - t172 * t89;
t64 = t169 * t90 + t172 * t92;
t61 = -t105 + t75;
t57 = -t119 * t229 + t171 * t75;
t56 = t119 * t228 + t168 * t75;
t55 = t117 * t228 - t168 * t74;
t54 = -t117 * t229 - t171 * t74;
t53 = -t172 * t106 + t169 * t71;
t49 = t171 * t100 - t234;
t48 = -t168 * t101 + t247;
t47 = t168 * t100 + t231;
t46 = t171 * t101 + t248;
t44 = t171 * t95 - t234;
t42 = t168 * t81 + t247;
t41 = t169 * t57 - t202;
t40 = t169 * t55 + t202;
t35 = pkin(2) * t94 + pkin(7) * t78 + t233;
t34 = pkin(2) * t90 + pkin(7) * t73 - t230;
t32 = -t168 * t61 + t171 * t59;
t31 = -t168 * t58 - t171 * t62;
t30 = t168 * t59 + t171 * t61;
t29 = t169 * t49 + t172 * t58;
t28 = t169 * t48 - t172 * t62;
t25 = -t169 * t63 + t172 * t45;
t24 = t169 * t45 + t172 * t63;
t23 = -t169 * t59 + t172 * t43;
t22 = t169 * t43 + t172 * t59;
t21 = t169 * t32 - t172 * t96;
t20 = -t169 * t76 + t172 * t33;
t19 = t169 * t33 + t172 * t76;
t18 = -pkin(2) * t85 + pkin(7) * t27;
t17 = -pkin(8) * t44 + t232;
t16 = -pkin(8) * t42 + t235;
t12 = pkin(2) * t99 + pkin(7) * t66 + t27;
t11 = -pkin(3) * t44 + t15;
t10 = -pkin(3) * t42 + t14;
t7 = t169 * t36 + t172 * t9;
t6 = t169 * t9 - t172 * t36;
t5 = -pkin(8) * t31 - t8;
t4 = -pkin(2) * t44 + pkin(7) * t25 + t172 * t11 + t169 * t17;
t3 = -pkin(2) * t42 + pkin(7) * t23 + t172 * t10 + t169 * t16;
t2 = pkin(7) * t20 + t169 * t5 + t199 * t31;
t1 = pkin(7) * t7 + (-pkin(8) * t169 + t199) * t8;
t13 = [0, 0, 0, 0, 0, qJDD(1), t181, t182, 0, 0, (t143 * t166 + t173 * t245) * t170, t167 * t145 + (t170 * t124 + t173 * t200) * t166, t167 * t122 + (t218 + t173 * (t159 - t205)) * t166, (t166 * t190 - t170 * t245) * t173, t167 * t123 + (t170 * (-t159 + t204) + t215) * t166, t167 * t160, (-t110 + pkin(1) * (t140 * t173 + t144 * t170)) * t167 + (t173 * t126 + pkin(1) * t124 + pkin(6) * (t173 * t144 - t218)) * t166, -t126 * t223 - t167 * t111 + pkin(1) * (-t166 * t200 + (t173 * t129 - t170 * t141) * t167) + (-t170 * t129 - t215) * t237, pkin(1) * ((-t173 * t122 + t170 * t123) * t167 - (-t164 - t165) * t163 * t221) + (t170 * t122 + t173 * t123) * t237 + t242 * t166, pkin(1) * (t166 * t126 + (-t110 * t173 + t111 * t170) * t167) + t242 * t237, t167 * t88 + (t170 * (t172 * t108 + t136 * t226) - t201) * t166, t167 * t64 + (t170 * (-t169 * t92 + t172 * t90) - t173 * t113) * t166, t167 * t79 + (t170 * (-t169 * t121 + t243) - t173 * t93) * t166, t167 * t87 + (t170 * (-t169 * t107 - t134 * t225) + t201) * t166, t167 * t80 + (t170 * (t172 * t120 + t220) + t173 * t89) * t166, t138 * t222 + t167 * t98 + (t134 * t172 - t136 * t169) * t153 * t223, (t34 + pkin(1) * (t170 * t73 + t173 * t90)) * t167 + (t170 * (-pkin(7) * t72 + t233) + t173 * (-pkin(2) * t72 + t51) - pkin(1) * t72 + pkin(6) * (-t170 * t90 + t173 * t73)) * t166, (t35 + pkin(1) * (t170 * t78 + t173 * t94)) * t167 + (t170 * (-pkin(7) * t77 + t230) + t173 * (-pkin(2) * t77 + t52) - pkin(1) * t77 + pkin(6) * (-t170 * t94 + t173 * t78)) * t166, (t12 + pkin(1) * (t170 * t66 + t173 * t99)) * t167 + (-t170 * t189 + pkin(6) * (-t170 * t99 + t173 * t66) + t185 * (-t169 * t89 - t172 * t93)) * t166, (t18 + pkin(1) * (t170 * t27 - t173 * t85)) * t167 + (pkin(6) * (t170 * t85 + t173 * t27) + t185 * t189) * t166, t167 * t41 + (t170 * (t172 * t57 + t203) - t173 * t56) * t166, t167 * t21 + (t170 * (t169 * t96 + t172 * t32) - t173 * t30) * t166, t167 * t28 + (t170 * (t169 * t62 + t172 * t48) - t173 * t46) * t166, t167 * t40 + (t170 * (t172 * t55 - t203) + t173 * t54) * t166, t167 * t29 + (t170 * (-t169 * t58 + t172 * t49) - t173 * t47) * t166, t167 * t53 + (t170 * (t169 * t106 + t172 * t71) - t173 * t70) * t166, (t3 + pkin(1) * (t170 * t23 - t173 * t42)) * t167 + (t170 * (-pkin(7) * t22 - t169 * t10 + t172 * t16) + t173 * (-pkin(2) * t22 - t183) - pkin(1) * t22 + pkin(6) * (t170 * t42 + t173 * t23)) * t166, (t4 + pkin(1) * (t170 * t25 - t173 * t44)) * t167 + (t170 * (-pkin(7) * t24 - t169 * t11 + t172 * t17) + t173 * (-pkin(2) * t24 - t184) - pkin(1) * t24 + pkin(6) * (t170 * t44 + t173 * t25)) * t166, (t2 + pkin(1) * (t170 * t20 - t173 * t31)) * t167 + (t170 * (-pkin(7) * t19 + t172 * t5 + t31 * t238) + t173 * (-pkin(2) * t19 - t180) - pkin(1) * t19 + pkin(6) * (t170 * t31 + t173 * t20)) * t166, (t1 + pkin(1) * (t170 * t7 - t173 * t8)) * t167 + (t170 * (-pkin(7) * t6 + (-pkin(8) * t172 + t238) * t8) + t173 * (-pkin(2) * t6 - t192) - pkin(1) * t6 + pkin(6) * (t170 * t8 + t173 * t7)) * t166; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t152, t145, t122, t152, t123, t160, -t110, -t111, 0, 0, t88, t64, t79, t87, t80, t98, t34, t35, t12, t18, t41, t21, t28, t40, t29, t53, t3, t4, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t227, t113, t93, -t227, -t89, -t138, -t51, -t52, 0, 0, t56, t30, t46, -t54, t47, t70, t183, t184, t180, t192; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97, t96, t62, -t97, -t58, t106, -t14, -t15, 0, 0;];
tauJ_reg = t13;
