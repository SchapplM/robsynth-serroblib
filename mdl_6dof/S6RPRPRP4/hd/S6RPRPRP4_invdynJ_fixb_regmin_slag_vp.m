% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRPRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
% 
% Output:
% tau_reg [6x26]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRPRP4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP4_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP4_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP4_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:13:04
% EndTime: 2019-03-09 03:13:12
% DurationCPUTime: 3.16s
% Computational Cost: add. (3153->414), mult. (5964->493), div. (0->0), fcn. (3487->10), ass. (0->215)
t139 = sin(pkin(9));
t114 = pkin(1) * t139 + pkin(7);
t92 = t114 * qJDD(1);
t290 = qJD(2) * qJD(3) + t92;
t142 = sin(qJ(3));
t218 = qJD(1) * qJD(3);
t203 = t142 * t218;
t145 = cos(qJ(3));
t213 = t145 * qJDD(1);
t286 = -t203 + t213;
t128 = t142 * qJ(4);
t289 = pkin(2) + t128;
t262 = pkin(4) + t114;
t94 = t114 * qJD(1);
t60 = -t145 * qJD(2) + t142 * t94;
t288 = qJD(4) + t60;
t224 = qJD(3) * t145;
t144 = cos(qJ(5));
t233 = t144 * t145;
t229 = qJD(1) * t142;
t109 = qJD(5) + t229;
t141 = sin(qJ(5));
t225 = qJD(3) * t144;
t208 = t142 * t225;
t221 = qJD(5) * t145;
t257 = (t141 * t221 + t208) * t109;
t228 = qJD(1) * t145;
t204 = t141 * t228;
t31 = -qJD(5) * t204 + t141 * qJDD(3) + (qJD(3) * qJD(5) + t286) * t144;
t28 = t142 * t31;
t202 = t145 * t218;
t214 = t142 * qJDD(1);
t169 = t202 + t214;
t78 = qJDD(5) + t169;
t227 = qJD(3) * t141;
t83 = t144 * t228 + t227;
t287 = t83 * t224 - t78 * t233 + t257 + t28;
t268 = g(3) * t142;
t133 = qJ(1) + pkin(9);
t122 = sin(t133);
t123 = cos(t133);
t280 = -g(1) * t123 - g(2) * t122;
t159 = t145 * t280 - t268;
t147 = -pkin(3) - pkin(8);
t220 = qJD(5) * t147;
t151 = -t109 * t220 + t159;
t134 = qJDD(3) * qJ(4);
t135 = qJD(3) * qJD(4);
t226 = qJD(3) * t142;
t197 = -t142 * qJDD(2) - t290 * t145 + t94 * t226;
t22 = -t134 - t135 + t197;
t13 = pkin(4) * t286 - t22;
t30 = qJD(5) * t83 - t144 * qJDD(3) + t141 * t286;
t85 = -t204 + t225;
t3 = pkin(5) * t31 + qJ(6) * t30 - qJD(6) * t85 + t13;
t285 = t3 + t151;
t136 = qJD(3) * qJ(4);
t61 = t142 * qJD(2) + t145 * t94;
t48 = -t136 - t61;
t45 = -qJD(3) * pkin(3) + t288;
t255 = t109 * t83;
t284 = t30 - t255;
t254 = t109 * t85;
t283 = -t31 + t254;
t178 = t109 * t141;
t243 = qJD(3) * t83;
t69 = t144 * t78;
t282 = -t109 * t178 - t243 + t69;
t232 = pkin(4) * t229 + t288;
t279 = t290 * t142 + t94 * t224;
t132 = g(3) * t145;
t216 = qJDD(2) * t145;
t168 = qJDD(4) - t216 + t279;
t12 = t169 * pkin(4) + t147 * qJDD(3) + t168;
t108 = pkin(3) * t203;
t246 = qJ(4) * t145;
t185 = pkin(8) * t142 - t246;
t219 = t142 * qJD(4);
t156 = qJD(3) * t185 - t219;
t140 = cos(pkin(9));
t265 = t140 * pkin(1);
t158 = t147 * t145 - t265 - t289;
t19 = qJD(1) * t156 + qJDD(1) * t158 + t108;
t222 = qJD(5) * t144;
t35 = t147 * qJD(3) + t232;
t212 = t141 * t12 + t144 * t19 + t35 * t222;
t40 = t158 * qJD(1);
t241 = qJD(5) * t40;
t235 = t141 * t142;
t54 = t122 * t144 + t123 * t235;
t56 = -t122 * t235 + t123 * t144;
t278 = -g(1) * t54 + g(2) * t56 + (-t241 + t132) * t141 + t212;
t242 = qJD(3) * t85;
t252 = t141 * t78;
t274 = t109 ^ 2;
t277 = t274 * t144 + t242 + t252;
t223 = qJD(5) * t141;
t256 = qJ(6) * t78;
t1 = qJD(6) * t109 - t40 * t223 + t212 + t256;
t8 = -t141 * t40 + t144 * t35;
t247 = qJD(6) - t8;
t6 = -pkin(5) * t109 + t247;
t9 = t141 * t35 + t144 * t40;
t7 = qJ(6) * t109 + t9;
t189 = t141 * t6 + t144 * t7;
t199 = -t144 * t12 + t141 * t19 + t40 * t222 + t35 * t223;
t273 = pkin(5) * t78;
t2 = qJDD(6) + t199 - t273;
t276 = qJD(5) * t189 + t1 * t141 - t2 * t144;
t275 = t85 ^ 2;
t272 = g(1) * t122;
t269 = g(2) * t123;
t267 = t109 * t7;
t266 = t109 * t9;
t130 = t145 * pkin(3);
t264 = t145 * pkin(8);
t263 = t85 * t83;
t115 = -pkin(2) - t265;
t231 = t130 + t128;
t70 = t115 - t231;
t57 = t70 - t264;
t71 = t262 * t142;
t261 = t141 * t71 + t144 * t57;
t44 = pkin(4) * t228 + t61;
t121 = pkin(3) * t229;
t62 = qJD(1) * t185 + t121;
t260 = t141 * t44 + t144 * t62;
t187 = pkin(5) * t144 + qJ(6) * t141;
t176 = -pkin(4) - t187;
t259 = qJD(5) * t187 - t144 * qJD(6) - t176 * t229 + t288;
t29 = t142 * t30;
t258 = t85 * t224 - t29;
t253 = t141 * t31;
t251 = t141 * t85;
t250 = t144 * t83;
t249 = t147 * t78;
t248 = t30 * t144;
t72 = t262 * t145;
t87 = -qJ(4) * t228 + t121;
t245 = qJD(1) * t87;
t38 = t136 + t44;
t18 = pkin(5) * t83 - qJ(6) * t85 + t38;
t244 = qJD(3) * t18;
t240 = qJDD(3) * pkin(3);
t239 = t122 * t142;
t238 = t122 * t145;
t237 = t123 * t142;
t236 = t123 * t145;
t234 = t142 * t144;
t137 = t142 ^ 2;
t138 = t145 ^ 2;
t230 = t137 - t138;
t95 = qJD(1) * t115;
t215 = qJDD(3) * t114;
t149 = qJD(1) ^ 2;
t211 = t142 * t149 * t145;
t210 = -g(1) * t237 - g(2) * t239 + t132;
t209 = t141 * t226;
t206 = t144 * t221;
t143 = sin(qJ(1));
t201 = -t143 * pkin(1) + t123 * pkin(7);
t196 = t85 * t208;
t53 = t122 * t141 - t123 * t234;
t55 = t122 * t234 + t123 * t141;
t195 = -g(1) * t55 - g(2) * t53;
t194 = -g(1) * t56 - g(2) * t54;
t146 = cos(qJ(1));
t192 = g(1) * t143 - g(2) * t146;
t190 = t141 * t7 - t144 * t6;
t188 = t146 * pkin(1) + pkin(3) * t236 + t122 * pkin(7) + t289 * t123;
t186 = pkin(5) * t141 - qJ(6) * t144;
t148 = qJD(3) ^ 2;
t184 = t114 * t148 + t269;
t177 = -t289 - t130;
t173 = -t109 * t222 - t252;
t120 = pkin(3) * t226;
t46 = t120 + t156;
t67 = t262 * t224;
t172 = t141 * t67 + t144 * t46 + t71 * t222 - t57 * t223;
t171 = -qJD(3) * t60 + t197;
t170 = -qJ(4) * t224 - t219;
t166 = t177 - t265;
t165 = t109 * t18 + t249;
t164 = -qJD(3) * t61 + t210 + t279;
t49 = t166 * qJD(1);
t163 = t215 + (-qJD(1) * t70 - t49) * qJD(3);
t162 = -0.2e1 * qJDD(1) * t115 - t184;
t161 = 0.2e1 * t95 * qJD(3) - t215;
t160 = g(1) * t53 - g(2) * t55 + g(3) * t233 - t199;
t25 = t168 - t240;
t26 = qJD(1) * t170 + qJDD(1) * t166 + t108;
t68 = t120 + t170;
t154 = qJD(1) * t68 + qJDD(1) * t70 + t184 + t26;
t153 = t25 * t142 - t22 * t145 + (t142 * t48 + t145 * t45) * qJD(3);
t152 = t18 * t85 + qJDD(6) - t160;
t102 = g(1) * t238;
t98 = qJ(4) * t236;
t96 = qJ(4) * t238;
t91 = qJDD(3) * t145 - t142 * t148;
t90 = qJDD(3) * t142 + t145 * t148;
t89 = qJ(4) + t186;
t75 = t109 * t209;
t66 = t262 * t226;
t52 = t83 * t209;
t41 = t49 * t229;
t37 = t145 * t187 + t72;
t36 = pkin(5) * t85 + qJ(6) * t83;
t24 = -pkin(5) * t142 + t141 * t57 - t144 * t71;
t23 = qJ(6) * t142 + t261;
t21 = -pkin(5) * t228 + t141 * t62 - t144 * t44;
t20 = qJ(6) * t228 + t260;
t16 = (-qJD(5) * t186 + qJD(6) * t141) * t145 + (-t114 + t176) * t226;
t5 = -pkin(5) * t224 + t261 * qJD(5) + t141 * t46 - t144 * t67;
t4 = qJ(6) * t224 + qJD(6) * t142 + t172;
t10 = [qJDD(1), t192, g(1) * t146 + g(2) * t143 (t192 + (t139 ^ 2 + t140 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), qJDD(1) * t137 + 0.2e1 * t142 * t202, 0.2e1 * t142 * t213 - 0.2e1 * t230 * t218, t90, t91, 0, t161 * t142 + t145 * t162 + t102, t161 * t145 + (-t162 - t272) * t142 (t137 + t138) * t92 + t153 + t280, t163 * t142 + t145 * t154 - t102, t163 * t145 + (-t154 + t272) * t142, -g(1) * t201 - g(2) * t188 + t114 * t153 - t177 * t272 + t26 * t70 + t49 * t68, -t85 * t206 + (t145 * t30 + t85 * t226) * t141, t196 - t52 + (t253 + t248 + (t250 + t251) * qJD(5)) * t145, t145 * t173 + t258 + t75, -t28 + (-t243 - t69) * t145 + t257, t109 * t224 + t142 * t78, -t199 * t142 + t8 * t224 - t66 * t83 + t72 * t31 + ((-qJD(5) * t71 - t46) * t109 - t57 * t78 - t38 * t221) * t141 + ((-qJD(5) * t57 + t67) * t109 + t71 * t78 + t13 * t145 - t38 * t226) * t144 + t194, -t172 * t109 - t261 * t78 - t66 * t85 - t72 * t30 + ((qJD(3) * t38 + t241) * t141 - t212) * t142 + (-qJD(3) * t9 - t13 * t141 - t222 * t38) * t145 - t195, -t109 * t5 + t16 * t83 - t24 * t78 + t31 * t37 + (-t18 * t225 - t2) * t142 + (-qJD(3) * t6 + t144 * t3 - t18 * t223) * t145 + t194, -t23 * t31 - t24 * t30 - t4 * t83 + t5 * t85 + t102 + t189 * t226 + (qJD(5) * t190 - t1 * t144 - t141 * t2 - t269) * t145, t109 * t4 - t16 * t85 + t23 * t78 + t30 * t37 + (-t18 * t227 + t1) * t142 + (qJD(3) * t7 + t141 * t3 + t18 * t222) * t145 + t195, t1 * t23 + t7 * t4 + t3 * t37 + t18 * t16 + t2 * t24 + t6 * t5 - g(1) * (pkin(4) * t123 + pkin(5) * t56 + qJ(6) * t55 + t201) - g(2) * (pkin(5) * t54 + pkin(8) * t236 + qJ(6) * t53 + t188) + (-g(1) * (t177 - t264) - g(2) * pkin(4)) * t122; 0, 0, 0, qJDD(2) - g(3), 0, 0, 0, 0, 0, t91, -t90, 0, -t91, t90, -t22 * t142 - t25 * t145 - g(3) + (t142 * t45 - t145 * t48) * qJD(3), 0, 0, 0, 0, 0, t287, t145 * t252 + (t206 - t209) * t109 + t258, t287, -t196 - t52 + (t253 - t248 + (t250 - t251) * qJD(5)) * t145, t29 + t75 + (t173 - t242) * t145, -g(3) + (qJD(3) * t190 + t3) * t142 + (t244 - t276) * t145; 0, 0, 0, 0, -t211, t230 * t149, t214, t213, qJDD(3), -t95 * t229 - t164 + t216, t268 + (-qJD(1) * t95 - t280) * t145 + t171 (-pkin(3) * t142 + t246) * qJDD(1), -0.2e1 * t240 + qJDD(4) + t41 + (-qJDD(2) - t245) * t145 + t164, 0.2e1 * t134 + 0.2e1 * t135 + (-g(3) + t245) * t142 + (qJD(1) * t49 + t280) * t145 - t171, -t22 * qJ(4) - t25 * pkin(3) - t49 * t87 - t45 * t61 - g(1) * (-pkin(3) * t237 + t98) - g(2) * (-pkin(3) * t239 + t96) - g(3) * t231 - t288 * t48, -t178 * t85 - t248 (-t31 - t254) * t144 + (t30 + t255) * t141, -t109 * t223 + t69 + (-t109 * t235 - t145 * t85) * qJD(1) (-t109 * t234 + t145 * t83) * qJD(1) + t173, -t109 * t228, -t8 * t228 + qJ(4) * t31 + t232 * t83 + (t249 + (t38 - t44) * t109) * t144 + (t13 + (t62 - t220) * t109 + t159) * t141, -qJ(4) * t30 + t260 * t109 + t9 * t228 + t232 * t85 + (-t109 * t38 - t249) * t141 + (t13 + t151) * t144, t109 * t21 + t141 * t285 + t165 * t144 + t6 * t228 + t259 * t83 + t31 * t89, t20 * t83 - t21 * t85 + (-t7 * t229 + t147 * t30 + t2 + (-t147 * t83 - t7) * qJD(5)) * t144 + (-t6 * t229 - t147 * t31 - t1 + (t147 * t85 - t6) * qJD(5)) * t141 - t210, -t109 * t20 + t165 * t141 - t144 * t285 - t7 * t228 - t259 * t85 + t30 * t89, t3 * t89 - t7 * t20 - t6 * t21 - g(1) * t98 - g(2) * t96 - g(3) * (t142 * t186 + t231 + t264) + t259 * t18 + t276 * t147 + t280 * (t147 * t142 + t186 * t145); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t214, qJDD(3) + t211, -t137 * t149 - t148, qJD(3) * t48 + t210 + t25 + t41, 0, 0, 0, 0, 0, t282, -t277, t282, t141 * t283 + t144 * t284, t277, -t244 + (-t2 + t267) * t144 + (t109 * t6 + t1) * t141 + t210; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t263, -t83 ^ 2 + t275, -t284, t283, t78, -t38 * t85 + t160 + t266, t109 * t8 + t38 * t83 - t278, -t36 * t83 - t152 + t266 + 0.2e1 * t273, pkin(5) * t30 - qJ(6) * t31 + (t7 - t9) * t85 + (t6 - t247) * t83, 0.2e1 * t256 - t18 * t83 + t36 * t85 + (0.2e1 * qJD(6) - t8) * t109 + t278, t1 * qJ(6) - t2 * pkin(5) - t18 * t36 - t6 * t9 - g(1) * (-pkin(5) * t53 + qJ(6) * t54) - g(2) * (pkin(5) * t55 - qJ(6) * t56) + t247 * t7 + t187 * t132; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t78 + t263, -t284, -t274 - t275, t152 - t267 - t273;];
tau_reg  = t10;
