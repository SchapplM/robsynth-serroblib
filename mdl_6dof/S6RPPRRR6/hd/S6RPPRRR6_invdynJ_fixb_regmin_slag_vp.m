% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
% 
% Output:
% tau_reg [6x30]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPPRRR6_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR6_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR6_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR6_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR6_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR6_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:31:32
% EndTime: 2019-03-09 02:31:41
% DurationCPUTime: 3.20s
% Computational Cost: add. (2369->385), mult. (4432->499), div. (0->0), fcn. (2956->10), ass. (0->209)
t140 = sin(qJ(4));
t217 = t140 * qJD(1);
t110 = qJD(5) + t217;
t144 = cos(qJ(4));
t141 = sin(qJ(1));
t145 = cos(qJ(1));
t171 = g(1) * t145 + g(2) * t141;
t265 = g(3) * t140;
t152 = t171 * t144 - t265;
t114 = qJ(2) * qJD(1) + qJD(3);
t101 = -pkin(7) * qJD(1) + t114;
t227 = qJD(4) * t140;
t129 = qJD(1) * qJD(2);
t130 = qJ(2) * qJDD(1);
t186 = qJDD(3) + t129 + t130;
t82 = -pkin(7) * qJDD(1) + t186;
t36 = -qJDD(4) * pkin(4) + t101 * t227 - t144 * t82;
t280 = qJD(5) * pkin(8) * t110 + t152 + t36;
t138 = sin(qJ(6));
t142 = cos(qJ(6));
t139 = sin(qJ(5));
t143 = cos(qJ(5));
t216 = t143 * qJD(4);
t229 = qJD(1) * t144;
t77 = t139 * t229 - t216;
t218 = t139 * qJD(4);
t79 = t143 * t229 + t218;
t170 = t138 * t77 - t142 * t79;
t32 = t138 * t79 + t142 * t77;
t279 = t170 * t32;
t278 = t170 ^ 2 - t32 ^ 2;
t103 = qJD(6) + t110;
t220 = qJD(6) * t142;
t221 = qJD(6) * t138;
t190 = t140 * t216;
t222 = qJD(5) * t144;
t194 = t139 * t222;
t209 = t144 * qJDD(1);
t27 = t143 * t209 + qJD(5) * t216 + t139 * qJDD(4) + (-t190 - t194) * qJD(1);
t215 = qJD(1) * qJD(4);
t189 = t140 * t215;
t28 = t79 * qJD(5) - t143 * qJDD(4) + (-t189 + t209) * t139;
t6 = -t138 * t28 + t142 * t27 - t77 * t220 - t79 * t221;
t277 = t32 * t103 + t6;
t137 = pkin(1) + qJ(3);
t172 = t140 * pkin(4) - t144 * pkin(8);
t88 = t172 + t137;
t48 = t88 * qJD(1) - qJD(2);
t87 = t140 * t101;
t61 = qJD(4) * pkin(8) + t87;
t19 = t139 * t48 + t143 * t61;
t16 = -t77 * pkin(9) + t19;
t12 = t16 * t221;
t134 = qJ(5) + qJ(6);
t120 = cos(t134);
t264 = g(3) * t144;
t238 = t144 * t101;
t62 = -qJD(4) * pkin(4) - t238;
t41 = t77 * pkin(5) + t62;
t119 = sin(t134);
t237 = t145 * t119;
t243 = t141 * t120;
t50 = -t140 * t243 - t237;
t236 = t145 * t120;
t244 = t141 * t119;
t52 = t140 * t236 - t244;
t276 = g(1) * t52 - g(2) * t50 + t120 * t264 + t41 * t32 + t12;
t128 = qJDD(1) * qJ(3);
t135 = qJDD(1) * pkin(1);
t207 = t135 - qJDD(2);
t185 = t128 + t207;
t173 = pkin(4) * t144 + pkin(8) * t140;
t75 = t173 * qJD(4) + qJD(3);
t26 = t75 * qJD(1) + t172 * qJDD(1) + t185;
t21 = t143 * t26;
t226 = qJD(4) * t144;
t37 = qJDD(4) * pkin(8) + t101 * t226 + t140 * t82;
t188 = t144 * t215;
t210 = t140 * qJDD(1);
t74 = qJDD(5) + t188 + t210;
t2 = t74 * pkin(5) - t27 * pkin(9) - t19 * qJD(5) - t139 * t37 + t21;
t223 = qJD(5) * t143;
t206 = -t139 * t26 - t143 * t37 - t48 * t223;
t225 = qJD(5) * t139;
t161 = -t61 * t225 - t206;
t3 = -t28 * pkin(9) + t161;
t200 = -t138 * t3 + t142 * t2;
t49 = t140 * t244 - t236;
t18 = -t139 * t61 + t143 * t48;
t15 = -t79 * pkin(9) + t18;
t11 = t110 * pkin(5) + t15;
t256 = t142 * t16;
t5 = t138 * t11 + t256;
t51 = -t140 * t237 - t243;
t275 = -g(1) * t51 + g(2) * t49 - t5 * qJD(6) + t119 * t264 + t41 * t170 + t200;
t150 = t170 * qJD(6) - t138 * t27 - t142 * t28;
t274 = -t103 * t170 + t150;
t81 = t138 * t143 + t142 * t139;
t53 = t81 * t144;
t232 = g(1) * t141 - g(2) * t145;
t272 = -t138 * t225 - t139 * t221;
t271 = -qJD(6) * t143 - t223;
t270 = qJD(1) * t137;
t208 = qJD(5) + qJD(6);
t158 = t81 * qJD(1);
t260 = t140 * t158 + t208 * t81;
t69 = qJDD(6) + t74;
t240 = t142 * t143;
t246 = t138 * t139;
t80 = -t240 + t246;
t263 = t80 * t69;
t269 = -t103 * t260 - t263;
t102 = -qJD(2) + t270;
t136 = -pkin(7) + qJ(2);
t268 = (qJD(2) + t102 + t270) * qJD(4) + qJDD(4) * t136;
t123 = 0.2e1 * t129;
t267 = pkin(8) + pkin(9);
t266 = pkin(9) * t144;
t262 = t81 * t69;
t156 = t271 * t142;
t198 = t139 * t217;
t261 = t138 * t198 + t208 * t246 - t217 * t240 + t156;
t84 = t173 * qJD(1);
t259 = t139 * t84 + t143 * t238;
t245 = t140 * t143;
t258 = t136 * t245 + t139 * t88;
t257 = t139 * t74;
t255 = t143 * t74;
t254 = t144 * t27;
t253 = t144 * t77;
t252 = t144 * t79;
t251 = t27 * t139;
t250 = t77 * t110;
t249 = t79 * t110;
t248 = qJD(4) * t77;
t247 = qJD(4) * t79;
t242 = t141 * t139;
t241 = t141 * t143;
t239 = t143 * t144;
t235 = t145 * t139;
t234 = t145 * t143;
t233 = t145 * pkin(1) + t141 * qJ(2);
t133 = t144 ^ 2;
t231 = t140 ^ 2 - t133;
t146 = qJD(4) ^ 2;
t147 = qJD(1) ^ 2;
t230 = -t146 - t147;
t228 = qJD(2) * t140;
t224 = qJD(5) * t140;
t214 = qJD(3) * qJD(1);
t212 = qJDD(4) * t140;
t211 = t137 * qJDD(1);
t204 = t110 * t245;
t203 = t139 * t140 * t110;
t202 = t32 * t229;
t201 = t170 * t229;
t199 = qJD(5) * t267;
t197 = t136 * t226;
t196 = t144 * t216;
t195 = t136 * t224;
t191 = t140 * t218;
t187 = qJDD(2) - t232;
t184 = qJD(6) * t11 + t3;
t182 = -qJD(5) * t48 - t37;
t181 = t110 * t136 + t61;
t180 = t136 * t196 + t139 * t75 + t143 * t228 + t88 * t223;
t179 = -0.2e1 * t188;
t178 = qJD(1) + t224;
t177 = -t135 + t187;
t176 = -t87 + (t198 + t225) * pkin(5);
t167 = pkin(5) * t144 + pkin(9) * t245;
t68 = t143 * t84;
t96 = t267 * t143;
t175 = t167 * qJD(1) + qJD(6) * t96 - t139 * t238 + t143 * t199 + t68;
t95 = t267 * t139;
t174 = pkin(9) * t198 + qJD(6) * t95 + t139 * t199 + t259;
t168 = -t128 + t177;
t165 = -t261 * t103 + t262;
t164 = t110 * t223 + t257;
t163 = t110 * t225 - t255;
t160 = t123 + 0.2e1 * t130 - t171;
t159 = t103 * t80;
t157 = t102 * qJD(1) + t171;
t155 = -t143 * t222 + t191;
t154 = -pkin(8) * t74 + t110 * t62;
t153 = t157 - t82;
t83 = t185 + t214;
t149 = -t136 * t146 + t211 + t214 + t232 + t83;
t122 = t145 * qJ(2);
t117 = qJDD(4) * t144;
t113 = -t143 * pkin(5) - pkin(4);
t73 = (pkin(5) * t139 - t136) * t144;
t71 = t143 * t88;
t66 = t140 * t234 - t242;
t65 = -t140 * t235 - t241;
t64 = -t140 * t241 - t235;
t63 = t140 * t242 - t234;
t57 = t143 * t75;
t54 = t80 * t144;
t42 = -pkin(5) * t155 - t144 * qJD(2) + t136 * t227;
t38 = -t139 * t266 + t258;
t30 = -pkin(9) * t239 + t71 + (-t136 * t139 + pkin(5)) * t140;
t14 = -t138 * t190 + t272 * t144 + (t208 * t239 - t191) * t142;
t13 = -t208 * t53 + t80 * t227;
t10 = t28 * pkin(5) + t36;
t9 = pkin(9) * t155 - t139 * t195 + t180;
t8 = -t143 * t195 + t57 + t167 * qJD(4) + (-t197 - t228 + (-t88 + t266) * qJD(5)) * t139;
t4 = t142 * t11 - t138 * t16;
t1 = [qJDD(1), t232, t171, -0.2e1 * t135 + t187, t160, t207 * pkin(1) - g(1) * (-t141 * pkin(1) + t122) - g(2) * t233 + (t123 + t130) * qJ(2), qJDD(3) + t160, -t168 + t211 + 0.2e1 * t214, t83 * t137 + t102 * qJD(3) + t186 * qJ(2) + t114 * qJD(2) - g(1) * (-t137 * t141 + t122) - g(2) * (t145 * qJ(3) + t233) t133 * qJDD(1) + t140 * t179, -0.2e1 * t140 * t209 + 0.2e1 * t231 * t215, -t146 * t140 + t117, -t146 * t144 - t212, 0, t149 * t140 + t144 * t268, -t140 * t268 + t149 * t144, -t79 * t194 + (-t227 * t79 + t254) * t143 (t139 * t79 + t143 * t77) * t227 + (-t251 - t143 * t28 + (t139 * t77 - t143 * t79) * qJD(5)) * t144 (-t110 * t216 + t27) * t140 + (-t163 + t247) * t144 (t110 * t218 - t28) * t140 + (-t164 - t248) * t144, t110 * t226 + t74 * t140, -g(1) * t64 - g(2) * t66 + t57 * t110 + t71 * t74 + (t136 * t248 - t181 * t223 + t21) * t140 + (-qJD(2) * t77 + t18 * qJD(4) - t136 * t28 + t223 * t62) * t144 + ((-qJD(5) * t88 - t197) * t110 + t36 * t144 + (-qJD(2) * t110 - t62 * qJD(4) - t136 * t74 + t182) * t140) * t139, -t180 * t110 - t258 * t74 - g(1) * t63 - g(2) * t65 + (t181 * t225 + (t136 * t79 - t143 * t62) * qJD(4) + t206) * t140 + (-qJD(2) * t79 - t19 * qJD(4) - t136 * t27 + t36 * t143 - t225 * t62) * t144, -t13 * t170 - t6 * t54, -t13 * t32 + t14 * t170 - t150 * t54 - t6 * t53, t13 * t103 + t6 * t140 - t170 * t226 - t54 * t69, -t14 * t103 + t140 * t150 - t226 * t32 - t53 * t69, t103 * t226 + t69 * t140 (-t138 * t9 + t142 * t8) * t103 + (-t138 * t38 + t142 * t30) * t69 + t200 * t140 + t4 * t226 + t42 * t32 - t73 * t150 + t10 * t53 + t41 * t14 - g(1) * t50 - g(2) * t52 + ((-t138 * t30 - t142 * t38) * t103 - t5 * t140) * qJD(6), -t5 * t226 - g(1) * t49 - g(2) * t51 - t10 * t54 + t12 * t140 + t41 * t13 - t42 * t170 + t73 * t6 + (-(-qJD(6) * t38 + t8) * t103 - t30 * t69 - t2 * t140) * t138 + (-(qJD(6) * t30 + t9) * t103 - t38 * t69 - t184 * t140) * t142; 0, 0, 0, qJDD(1), -t147, -t147 * qJ(2) + t177, -t147, -qJDD(1) (-qJD(3) - t114) * qJD(1) + t168, 0, 0, 0, 0, 0, t179 - t210, 0.2e1 * t189 - t209, 0, 0, 0, 0, 0 (t203 + t253) * qJD(1) + t163 (t204 + t252) * qJD(1) + t164, 0, 0, 0, 0, 0, t202 - t269, t165 - t201; 0, 0, 0, 0, 0, 0, qJDD(1), -t147, -t157 + t186, 0, 0, 0, 0, 0, t230 * t140 + t117, t144 * t230 - t212, 0, 0, 0, 0, 0, -t144 * t28 + (t248 - t257) * t140 + (-t143 * t178 - t144 * t218) * t110, -t254 + (t247 - t255) * t140 + (t139 * t178 - t196) * t110, 0, 0, 0, 0, 0, qJD(1) * t159 + (-qJD(4) * t103 * t81 + t150) * t144 + ((t156 - t272) * t103 - t262 + qJD(4) * t32) * t140, t103 * t158 + (qJD(4) * t159 - t6) * t144 + (-(t138 * t271 - t139 * t220 - t142 * t225) * t103 + t263 - qJD(4) * t170) * t140; 0, 0, 0, 0, 0, 0, 0, 0, 0, t144 * t147 * t140, -t231 * t147, t209, -t210, qJDD(4), -t144 * t153 + t265, t140 * t153 + t264, t143 * t249 + t251 (t27 - t250) * t143 + (-t28 - t249) * t139 (t204 - t252) * qJD(1) + t164 (-t203 + t253) * qJD(1) - t163, -t110 * t229, -t18 * t229 - t77 * t87 - pkin(4) * t28 - t68 * t110 + (t110 * t238 + t154) * t139 - t280 * t143, -pkin(4) * t27 + t259 * t110 + t280 * t139 + t154 * t143 + t19 * t229 - t79 * t87, t170 * t261 + t6 * t81, t150 * t81 + t170 * t260 + t261 * t32 - t6 * t80, t165 + t201, t202 + t269, -t103 * t229 (-t138 * t96 - t142 * t95) * t69 - t113 * t150 + t10 * t80 - t4 * t229 + t260 * t41 + t176 * t32 + (t138 * t174 - t142 * t175) * t103 - t152 * t120 -(-t138 * t95 + t142 * t96) * t69 + t113 * t6 + t10 * t81 + t5 * t229 - t261 * t41 - t176 * t170 + (t138 * t175 + t142 * t174) * t103 + t152 * t119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79 * t77, -t77 ^ 2 + t79 ^ 2, t27 + t250, t249 - t28, t74, -t61 * t223 - g(1) * t65 + g(2) * t63 + t19 * t110 - t62 * t79 + t21 + (t182 + t264) * t139, g(1) * t66 - g(2) * t64 + g(3) * t239 + t18 * t110 + t62 * t77 - t161, -t279, t278, t277, t274, t69 -(-t138 * t15 - t256) * t103 + (-t103 * t221 + t142 * t69 - t79 * t32) * pkin(5) + t275 (-t16 * t103 - t2) * t138 + (t15 * t103 - t184) * t142 + (-t103 * t220 - t138 * t69 + t170 * t79) * pkin(5) + t276; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t279, t278, t277, t274, t69, t5 * t103 + t275, t4 * t103 - t138 * t2 - t142 * t184 + t276;];
tau_reg  = t1;
