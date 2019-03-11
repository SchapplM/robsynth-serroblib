% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRRPRP2
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% tau_reg [6x24]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRRPRP2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP2_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP2_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRP2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP2_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:32:19
% EndTime: 2019-03-08 21:32:29
% DurationCPUTime: 4.43s
% Computational Cost: add. (5524->481), mult. (12877->642), div. (0->0), fcn. (10207->14), ass. (0->243)
t180 = sin(pkin(11));
t187 = sin(qJ(3));
t190 = cos(qJ(3));
t293 = cos(pkin(11));
t237 = t293 * t190;
t212 = -t180 * t187 + t237;
t191 = cos(qJ(2));
t182 = sin(pkin(6));
t266 = qJD(1) * t182;
t245 = t191 * t266;
t104 = t212 * t245;
t308 = qJ(4) + pkin(8);
t239 = qJD(3) * t308;
t130 = t190 * qJD(4) - t187 * t239;
t131 = -t187 * qJD(4) - t190 * t239;
t68 = t130 * t293 + t180 * t131;
t324 = t68 - t104;
t162 = qJD(2) * t237;
t264 = qJD(2) * t187;
t137 = -t180 * t264 + t162;
t129 = qJD(5) - t137;
t152 = t308 * t187;
t153 = t308 * t190;
t106 = -t180 * t152 + t153 * t293;
t186 = sin(qJ(5));
t189 = cos(qJ(5));
t188 = sin(qJ(2));
t248 = t188 * t266;
t260 = qJD(5) * t189;
t261 = qJD(5) * t186;
t146 = t180 * t190 + t187 * t293;
t138 = t146 * qJD(3);
t141 = t212 * qJD(3);
t303 = qJD(3) * pkin(3);
t253 = t187 * t303;
t70 = pkin(4) * t138 - pkin(9) * t141 + t253;
t171 = t190 * pkin(3) + pkin(2);
t83 = -pkin(4) * t212 - pkin(9) * t146 - t171;
t325 = t106 * t261 - t83 * t260 - t324 * t189 + (t248 - t70) * t186;
t294 = t180 * t130 - t293 * t131 - t146 * t245;
t277 = t182 * t191;
t159 = t186 * t277;
t279 = t182 * t188;
t251 = t187 * t279;
t184 = cos(pkin(6));
t275 = t184 * t190;
t142 = -t251 + t275;
t278 = t182 * t190;
t143 = t184 * t187 + t188 * t278;
t75 = t180 * t142 + t143 * t293;
t59 = t189 * t75 - t159;
t181 = sin(pkin(10));
t183 = cos(pkin(10));
t274 = t184 * t191;
t133 = t181 * t188 - t183 * t274;
t135 = t181 * t274 + t183 * t188;
t232 = g(1) * t135 + g(2) * t133;
t323 = -g(3) * t277 + t232;
t258 = qJD(2) * qJD(3);
t243 = t187 * t258;
t322 = pkin(3) * t243 + qJDD(4);
t276 = t184 * t188;
t134 = t181 * t191 + t183 * t276;
t136 = -t181 * t276 + t183 * t191;
t177 = qJ(3) + pkin(11);
t172 = sin(t177);
t173 = cos(t177);
t280 = t182 * t183;
t281 = t181 * t182;
t209 = -g(3) * (-t172 * t279 + t173 * t184) - g(2) * (-t134 * t172 - t173 * t280) - g(1) * (-t136 * t172 + t173 * t281);
t139 = t146 * qJD(2);
t107 = -t189 * qJD(3) + t139 * t186;
t109 = qJD(3) * t186 + t139 * t189;
t234 = qJD(2) * t308 + t248;
t265 = qJD(1) * t184;
t102 = t187 * t265 + t190 * t234;
t92 = t180 * t102;
t101 = -t187 * t234 + t190 * t265;
t96 = t101 + t303;
t43 = t293 * t96 - t92;
t38 = -qJD(3) * pkin(4) - t43;
t21 = t107 * pkin(5) - t109 * qJ(6) + t38;
t311 = pkin(3) * t180;
t168 = pkin(9) + t311;
t255 = t187 * qJDD(2);
t226 = qJDD(2) * t237 - t180 * t255;
t81 = qJD(2) * t138 + qJDD(5) - t226;
t300 = t168 * t81;
t321 = t129 * t21 - t300;
t241 = qJDD(1) * t277;
t259 = qJD(1) * qJD(2);
t244 = t188 * t259;
t227 = t182 * t244 - t241;
t292 = qJDD(2) * pkin(2);
t116 = t227 - t292;
t192 = qJD(3) ^ 2;
t320 = -pkin(8) * t192 + t182 * (-g(3) * t191 + t244) - t116 + t232 + t292;
t319 = qJD(3) * t162 + qJDD(2) * t146 - t180 * t243;
t318 = t109 ^ 2;
t317 = t129 ^ 2;
t316 = pkin(5) * t81;
t313 = qJ(6) * t138 - qJD(6) * t212 - t325;
t305 = t189 * t106 + t186 * t83;
t71 = t104 * t186 - t189 * t248;
t312 = -t138 * pkin(5) + qJD(5) * t305 + t186 * t68 - t189 * t70 - t71;
t310 = pkin(3) * t187;
t256 = t184 * qJDD(1);
t161 = t190 * t256;
t117 = qJDD(2) * pkin(8) + (qJDD(1) * t188 + t191 * t259) * t182;
t201 = qJ(4) * qJDD(2) + qJD(2) * qJD(4) + qJD(3) * t265 + t117;
t224 = t234 * qJD(3);
t34 = qJDD(3) * pkin(3) - t187 * t201 - t190 * t224 + t161;
t35 = (-t224 + t256) * t187 + t201 * t190;
t11 = t180 * t34 + t293 * t35;
t194 = -t189 * qJDD(3) + t186 * t319;
t41 = qJD(5) * t109 + t194;
t307 = -t107 * t260 - t186 * t41;
t48 = t101 * t293 - t92;
t69 = pkin(3) * t264 + pkin(4) * t139 - pkin(9) * t137;
t306 = t186 * t69 + t189 * t48;
t238 = t293 * t102;
t44 = t180 * t96 + t238;
t304 = qJD(2) * pkin(2);
t39 = qJD(3) * pkin(9) + t44;
t128 = -qJD(2) * t171 + qJD(4) - t245;
t57 = -t137 * pkin(4) - t139 * pkin(9) + t128;
t19 = t186 * t57 + t189 * t39;
t13 = qJ(6) * t129 + t19;
t302 = t129 * t13;
t301 = t129 * t19;
t76 = t186 * t81;
t77 = t189 * t81;
t257 = qJD(3) * qJD(5);
t40 = -t186 * qJDD(3) + t139 * t261 + (-t319 - t257) * t189;
t298 = t40 * t186;
t297 = t81 * qJ(6);
t228 = pkin(5) * t186 - qJ(6) * t189;
t46 = t180 * t101 + t238;
t296 = -t186 * qJD(6) + t129 * t228 - t46;
t229 = t189 * pkin(5) + t186 * qJ(6);
t295 = t228 * t141 + (qJD(5) * t229 - qJD(6) * t189) * t146 + t294;
t291 = t107 * t137;
t290 = t109 * t107;
t236 = t109 * t129;
t289 = t109 * t139;
t288 = t129 * t186;
t287 = t139 * t107;
t286 = t141 * t186;
t285 = t141 * t189;
t284 = t146 * t189;
t283 = t173 * t186;
t282 = t173 * t189;
t273 = t308 * t188;
t272 = t189 * t191;
t18 = -t186 * t39 + t189 * t57;
t271 = qJD(6) - t18;
t270 = qJDD(1) - g(3);
t269 = -t133 * t171 + t134 * t308;
t268 = -t135 * t171 + t136 * t308;
t178 = t187 ^ 2;
t267 = -t190 ^ 2 + t178;
t263 = qJD(2) * t188;
t262 = qJD(5) * t168;
t254 = t190 * qJDD(2);
t252 = t183 * t278;
t250 = t182 * t272;
t249 = t293 * pkin(3);
t247 = t182 * t263;
t246 = qJD(2) * t277;
t242 = t191 * t258;
t85 = -qJD(3) * t139 + t226;
t26 = -pkin(3) * t254 - t85 * pkin(4) - pkin(9) * t319 + t116 + t322;
t9 = qJDD(3) * pkin(9) + t11;
t240 = t186 * t9 - t189 * t26 + t39 * t260 + t57 * t261;
t105 = t293 * t152 + t180 * t153;
t233 = pkin(4) * t173 + pkin(9) * t172;
t231 = g(1) * t136 + g(2) * t134;
t230 = g(1) * t181 - g(2) * t183;
t10 = -t180 * t35 + t293 * t34;
t12 = -pkin(5) * t129 + t271;
t225 = t12 * t189 - t13 * t186;
t193 = qJD(2) ^ 2;
t223 = qJDD(2) * t191 - t188 * t193;
t222 = -t129 * t261 + t137 * t288 + t77;
t221 = t76 + (-t137 * t189 + t260) * t129;
t220 = pkin(4) + t229;
t58 = t186 * t75 + t250;
t218 = t186 * t26 + t189 * t9 + t57 * t260 - t261 * t39;
t216 = -t146 * t261 + t285;
t100 = -qJD(3) * t143 - t187 * t246;
t99 = qJD(3) * t142 + t190 * t246;
t47 = t180 * t100 + t293 * t99;
t16 = qJD(5) * t59 + t186 * t47 - t189 * t247;
t45 = -t100 * t293 + t180 * t99;
t74 = -t142 * t293 + t180 * t143;
t214 = t45 * t107 - t129 * t16 + t74 * t41 - t58 * t81;
t213 = t129 * t38 - t300;
t111 = t159 * t173 - t189 * t279;
t63 = -t133 * t283 - t134 * t189;
t65 = -t135 * t283 - t136 * t189;
t211 = g(1) * t65 + g(2) * t63 + g(3) * t111;
t112 = (t173 * t272 + t186 * t188) * t182;
t64 = -t133 * t282 + t134 * t186;
t66 = -t135 * t282 + t136 * t186;
t210 = -g(1) * t66 - g(2) * t64 - g(3) * t112;
t8 = -qJDD(3) * pkin(4) - t10;
t207 = -g(3) * t279 - t231;
t150 = -t245 - t304;
t206 = -qJD(2) * t150 - t117 + t231;
t15 = -qJD(5) * t58 + t186 * t247 + t189 * t47;
t205 = t109 * t45 - t129 * t15 - t40 * t74 - t59 * t81;
t87 = t134 * t173 - t172 * t280;
t52 = -t133 * t189 + t186 * t87;
t89 = t136 * t173 + t172 * t281;
t54 = -t135 * t189 + t186 * t89;
t120 = t172 * t184 + t173 * t279;
t90 = t120 * t186 + t250;
t204 = g(1) * t54 + g(2) * t52 + g(3) * t90 - t240;
t203 = -t129 * t262 + t209;
t3 = t41 * pkin(5) + t40 * qJ(6) - t109 * qJD(6) + t8;
t200 = t203 - t3;
t199 = -pkin(8) * qJDD(3) + (t150 + t245 - t304) * qJD(3);
t53 = t133 * t186 + t189 * t87;
t55 = t135 * t186 + t189 * t89;
t91 = t120 * t189 - t159;
t198 = -g(1) * t55 - g(2) * t53 - g(3) * t91 + t218;
t84 = -qJDD(2) * t171 + t227 + t322;
t197 = t109 * t21 + qJDD(6) - t204;
t169 = -t249 - pkin(4);
t164 = pkin(3) * t275;
t151 = t181 * pkin(3) * t278;
t148 = t171 * t277;
t144 = -t249 - t220;
t56 = pkin(5) * t109 + qJ(6) * t107;
t49 = t146 * t228 + t105;
t30 = pkin(5) * t212 + t106 * t186 - t189 * t83;
t29 = -qJ(6) * t212 + t305;
t22 = t107 * t129 - t40;
t20 = -pkin(5) * t139 + t186 * t48 - t189 * t69;
t17 = qJ(6) * t139 + t306;
t2 = qJDD(6) + t240 - t316;
t1 = qJD(6) * t129 + t218 + t297;
t4 = [t270, 0, t223 * t182 (-qJDD(2) * t188 - t191 * t193) * t182, 0, 0, 0, 0, 0, t100 * qJD(3) + t142 * qJDD(3) + (-t187 * t242 + t190 * t223) * t182, -t99 * qJD(3) - t143 * qJDD(3) + (-t187 * t223 - t190 * t242) * t182, t47 * t137 + t45 * t139 + t319 * t74 + t75 * t85, -t10 * t74 + t11 * t75 - t43 * t45 + t44 * t47 - g(3) + (t128 * t263 - t191 * t84) * t182, 0, 0, 0, 0, 0, t214, t205, t214, -t107 * t15 + t109 * t16 - t40 * t58 - t41 * t59, -t205, t1 * t59 + t12 * t16 + t13 * t15 + t2 * t58 + t21 * t45 + t3 * t74 - g(3); 0, qJDD(2), t323 + t241, -t270 * t279 + t231, qJDD(2) * t178 + 0.2e1 * t190 * t243, 0.2e1 * t187 * t254 - 0.2e1 * t258 * t267, qJDD(3) * t187 + t190 * t192, qJDD(3) * t190 - t187 * t192, 0, t199 * t187 + t190 * t320, -t187 * t320 + t199 * t190, -t10 * t146 + t105 * t319 + t106 * t85 + t11 * t212 + t137 * t324 - t44 * t138 + t139 * t294 - t43 * t141 + t207, t11 * t106 - t10 * t105 - t84 * t171 - g(1) * t268 - g(2) * t269 - g(3) * (t182 * t273 + t148) + t324 * t44 - t294 * t43 + (-t248 + t253) * t128, t109 * t216 - t284 * t40 (-t107 * t189 - t109 * t186) * t141 + (t298 - t189 * t41 + (t107 * t186 - t109 * t189) * qJD(5)) * t146, t109 * t138 + t129 * t216 + t212 * t40 + t284 * t81, -t146 * t76 - t107 * t138 + t41 * t212 + (-t146 * t260 - t286) * t129, t129 * t138 - t212 * t81, t240 * t212 + t18 * t138 + t105 * t41 + t71 * t129 + t294 * t107 + ((-qJD(5) * t106 + t70) * t129 + t83 * t81 + t38 * qJD(5) * t146) * t189 + ((-qJD(5) * t83 - t68) * t129 - t106 * t81 + t8 * t146 + t38 * t141) * t186 + t210, -t305 * t81 + t218 * t212 - t19 * t138 - t105 * t40 + t38 * t285 + (t8 * t189 - t261 * t38) * t146 + t325 * t129 + t294 * t109 + t211, t21 * t286 - t12 * t138 + t2 * t212 - t30 * t81 + t49 * t41 + (t3 * t186 + t21 * t260) * t146 - t312 * t129 + t295 * t107 + t210, -t29 * t41 - t30 * t40 + t225 * t141 + t312 * t109 - t313 * t107 + t323 * t172 + (-t1 * t186 + t2 * t189 + (-t12 * t186 - t13 * t189) * qJD(5)) * t146, -t21 * t285 - t1 * t212 + t13 * t138 + t29 * t81 + t49 * t40 + (-t3 * t189 + t21 * t261) * t146 + t313 * t129 - t295 * t109 - t211, t1 * t29 + t3 * t49 + t2 * t30 - g(1) * (t66 * pkin(5) + t65 * qJ(6) - t135 * t233 + t268) - g(2) * (t64 * pkin(5) + t63 * qJ(6) - t133 * t233 + t269) + t295 * t21 + t313 * t13 + t312 * t12 + (-t112 * pkin(5) - t111 * qJ(6) - t148 - (t191 * t233 + t273) * t182) * g(3); 0, 0, 0, 0, -t187 * t193 * t190, t267 * t193, t255, t254, qJDD(3), -g(3) * t142 + t187 * t206 - t230 * t278 + t161, g(3) * t143 + (t182 * t230 - t256) * t187 + t206 * t190, -t319 * t249 + t85 * t311 - (-t44 + t46) * t139 + (-t48 + t43) * t137, -g(1) * t151 - g(3) * t164 + t43 * t46 - t44 * t48 + (g(2) * t252 + t10 * t293 + t11 * t180 + (-qJD(2) * t128 - t207) * t187) * pkin(3), t189 * t236 - t298 (-t40 + t291) * t189 - t109 * t288 + t307, t221 - t289, t222 + t287, -t129 * t139, -t46 * t107 - t18 * t139 + t169 * t41 + (t48 * t129 + t213) * t186 + (-t8 + (-t69 - t262) * t129 + t209) * t189, -t169 * t40 + t306 * t129 + t19 * t139 - t46 * t109 + t213 * t189 + (-t203 + t8) * t186, t296 * t107 + t12 * t139 + t20 * t129 + t144 * t41 + t186 * t321 + t200 * t189, -g(1) * t89 - g(2) * t87 - g(3) * t120 + t17 * t107 - t20 * t109 + (-t12 * t137 - t168 * t41 + t1 + (t109 * t168 + t12) * qJD(5)) * t189 + (t13 * t137 - t168 * t40 + t2 + (t107 * t168 - t13) * qJD(5)) * t186, -t296 * t109 - t17 * t129 - t13 * t139 + t144 * t40 + t200 * t186 - t189 * t321, t3 * t144 - t13 * t17 - t12 * t20 - g(1) * (t89 * pkin(9) - t136 * t310 + t151) - g(2) * (-pkin(3) * t252 + t87 * pkin(9) - t134 * t310) - g(3) * (-pkin(3) * t251 + t120 * pkin(9) + t164) + t296 * t21 + (qJD(5) * t225 + t1 * t189 + t2 * t186) * t168 + t209 * t220; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t137 ^ 2 - t139 ^ 2, -t44 * t137 + t43 * t139 - t323 + t84, 0, 0, 0, 0, 0, t222 - t287, -t317 * t189 - t289 - t76, -t129 * t288 - t287 + t77 (t40 + t291) * t189 + t186 * t236 + t307, t221 + t289, -t21 * t139 + (-t2 + t302) * t189 + (t12 * t129 + t1) * t186 - t323; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t290, -t107 ^ 2 + t318, t22, -t139 * t260 - t186 * t257 - t194 + t236, t81, -t109 * t38 + t204 + t301, t107 * t38 + t129 * t18 - t198, -t107 * t56 - t197 + t301 + 0.2e1 * t316, pkin(5) * t40 - t41 * qJ(6) + (t13 - t19) * t109 + (t12 - t271) * t107, 0.2e1 * t297 - t21 * t107 + t56 * t109 + (0.2e1 * qJD(6) - t18) * t129 + t198, t1 * qJ(6) - t2 * pkin(5) - t21 * t56 - t12 * t19 - g(1) * (-pkin(5) * t54 + qJ(6) * t55) - g(2) * (-pkin(5) * t52 + qJ(6) * t53) - g(3) * (-pkin(5) * t90 + qJ(6) * t91) + t271 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(5) + t85 + t290, t22, -t317 - t318, t197 - t302 - t316;];
tau_reg  = t4;
