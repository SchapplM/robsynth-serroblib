% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RRPPPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta4]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 08:45
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RRPPPR4_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR4_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR4_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPPR4_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR4_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR4_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 08:44:17
% EndTime: 2019-05-06 08:44:31
% DurationCPUTime: 5.76s
% Computational Cost: add. (13302->401), mult. (29923->483), div. (0->0), fcn. (17834->8), ass. (0->245)
t203 = sin(pkin(9));
t207 = sin(qJ(2));
t210 = cos(qJ(2));
t278 = t207 * qJ(3);
t305 = -pkin(2) - qJ(4);
t225 = t210 * t305 - pkin(1) - t278;
t204 = cos(pkin(9));
t268 = qJD(1) * t210;
t161 = qJD(2) * t203 + t204 * t268;
t163 = qJD(2) * t204 - t203 * t268;
t127 = t163 * t161;
t262 = qJD(1) * qJD(2);
t251 = t210 * t262;
t259 = t207 * qJDD(1);
t172 = t251 + t259;
t323 = t127 + t172;
t290 = t323 * t204;
t194 = t210 * qJDD(1);
t253 = t207 * t262;
t238 = -t194 + t253;
t137 = t204 * qJDD(2) + t203 * t238;
t269 = qJD(1) * t207;
t254 = t161 * t269;
t320 = -t254 + t137;
t159 = t163 ^ 2;
t200 = t207 ^ 2;
t213 = qJD(1) ^ 2;
t195 = t200 * t213;
t322 = -t159 - t195;
t291 = t323 * t203;
t78 = -t204 * t322 + t291;
t362 = -pkin(7) * (t207 * t78 - t210 * t320) - t225 * (t203 * t322 + t290);
t361 = pkin(3) * t78;
t360 = t305 * t78;
t136 = qJDD(2) * t203 - t204 * t238;
t255 = t163 * t269;
t107 = -t136 + t255;
t313 = t161 ^ 2;
t142 = t313 - t195;
t357 = -t207 * t107 + t210 * (t142 * t203 + t290);
t292 = t320 * t204;
t105 = t136 + t255;
t295 = t105 * t203;
t356 = t210 * (t292 - t295) - t207 * (t159 - t313);
t100 = -t313 - t159;
t109 = t254 + t137;
t334 = t107 * t203 - t109 * t204;
t299 = t207 * t334;
t355 = t225 * (t107 * t204 + t203 * t109) + pkin(7) * (t100 * t210 + t299);
t354 = pkin(3) * t100;
t353 = qJ(5) * t320;
t143 = -t159 + t195;
t324 = -t127 + t172;
t289 = t324 * t204;
t351 = -t143 * t203 + t289;
t266 = qJD(4) * t163;
t317 = -t313 - t195;
t333 = t203 * t317 + t289;
t348 = pkin(3) * t333 - 0.2e1 * t266;
t347 = qJ(3) * t105 + t305 * t333;
t345 = t142 * t204 - t291;
t98 = t203 * t324;
t344 = t207 * t109 + t210 * (-t143 * t204 - t98);
t343 = qJ(3) * t100 + t305 * t334;
t342 = t225 * (t204 * t317 - t98) + pkin(7) * (t210 * t105 + t207 * t333);
t341 = pkin(3) * t320;
t212 = qJD(2) ^ 2;
t182 = -t195 - t212;
t274 = t207 * t213;
t256 = t210 * t274;
t178 = -qJDD(2) + t256;
t272 = t210 * t178;
t340 = pkin(7) * (-t182 * t207 + t272);
t294 = t105 * t204;
t332 = t203 * t320 + t294;
t312 = pkin(3) + pkin(7);
t311 = pkin(4) + pkin(5);
t206 = sin(qJ(6));
t209 = cos(qJ(6));
t119 = -t209 * t161 + t163 * t206;
t187 = -qJD(6) + t269;
t103 = t119 * t187;
t70 = -qJD(6) * t119 + t136 * t206 + t137 * t209;
t329 = t103 + t70;
t164 = -qJDD(6) + t172;
t121 = t161 * t206 + t163 * t209;
t288 = t121 * t119;
t221 = -t164 - t288;
t328 = t206 * t221;
t326 = t209 * t221;
t239 = t172 + t251;
t325 = t239 * qJ(3);
t208 = sin(qJ(1));
t211 = cos(qJ(1));
t241 = g(1) * t211 + g(2) * t208;
t296 = qJDD(1) * pkin(7);
t155 = -pkin(1) * t213 - t241 + t296;
t309 = pkin(2) * t210;
t240 = -t278 - t309;
t169 = t240 * qJD(1);
t233 = (qJD(1) * t169 + t155) * t210;
t179 = pkin(3) * t269 - qJD(2) * qJ(4);
t252 = qJD(3) * t269;
t191 = -0.2e1 * t252;
t201 = t210 ^ 2;
t249 = t208 * g(1) - t211 * g(2);
t228 = -qJDD(1) * pkin(1) - t249;
t217 = -t228 + (-t238 - t253) * pkin(2);
t62 = t191 - t179 * t269 + t238 * qJ(4) + (-pkin(3) * t201 - pkin(7)) * t213 - t325 - t217;
t276 = t207 * t155;
t226 = -qJDD(2) * pkin(2) - t212 * qJ(3) + t169 * t269 + qJDD(3) + t276;
t77 = t172 * pkin(3) - qJDD(2) * qJ(4) + (-pkin(3) * t262 - qJ(4) * t274 + g(3)) * t210 + t226;
t40 = -0.2e1 * qJD(4) * t161 + t203 * t77 + t204 * t62;
t279 = t201 * t213;
t319 = t272 - (-t212 + t279) * t207;
t173 = t194 - 0.2e1 * t253;
t183 = t212 + t279;
t177 = qJDD(2) + t256;
t284 = t177 * t207;
t316 = pkin(7) * (t183 * t210 + t284) - pkin(1) * t173;
t246 = -t209 * t136 + t206 * t137;
t55 = (qJD(6) + t187) * t121 + t246;
t273 = t210 * t155;
t285 = t163 * t207;
t315 = -(pkin(4) * t285 + t169 * t210) * qJD(1) - t273;
t314 = t136 * pkin(4) - 0.2e1 * qJD(5) * t163 - t353;
t117 = t119 ^ 2;
t118 = t121 ^ 2;
t185 = t187 ^ 2;
t308 = pkin(3) * t105;
t197 = t207 * g(3);
t307 = t210 * g(3);
t306 = t213 * pkin(7);
t260 = qJDD(2) * qJ(3);
t270 = -t212 * pkin(2) - t197;
t223 = -t238 * pkin(3) - qJ(4) * t279 + qJDD(4) + t260 + t270;
t219 = (0.2e1 * qJD(3) + t179) * qJD(2) + t223;
t73 = t233 + t219;
t303 = t203 * t73;
t302 = t204 * t73;
t231 = -pkin(5) * t269 - pkin(8) * t163;
t41 = t219 + t314 - t315;
t36 = t136 * pkin(5) + pkin(8) * t313 - t163 * t231 + t41;
t301 = t206 * t36;
t74 = t164 - t288;
t300 = t206 * t74;
t298 = t209 * t36;
t297 = t209 * t74;
t281 = t187 * t206;
t280 = t187 * t209;
t275 = t207 * t173;
t175 = t195 + t279;
t271 = (t200 + t201) * t296 + pkin(1) * t175;
t124 = pkin(4) * t161 - qJ(5) * t163;
t264 = 0.2e1 * qJD(4) + t124;
t261 = qJD(3) * qJD(2);
t258 = t207 * t288;
t257 = t161 * t285;
t250 = -pkin(4) * t203 - qJ(3);
t248 = t203 * t62 - t204 * t77;
t227 = -t172 * pkin(4) - qJ(5) * t195 + qJDD(5) + t248;
t22 = -t172 * pkin(5) - t109 * pkin(8) + (pkin(5) * t161 + t264) * t163 + t227;
t232 = t172 * qJ(5) + 0.2e1 * qJD(5) * t269 - t161 * t124 + t40;
t30 = -pkin(4) * t195 + t232;
t25 = -pkin(5) * t313 + t136 * pkin(8) + t231 * t269 + t30;
t9 = t206 * t25 - t209 * t22;
t132 = t276 + t307;
t133 = -t197 + t273;
t247 = t132 * t207 + t210 * t133;
t244 = t203 * t254;
t243 = t204 * t254;
t139 = t204 * t255;
t242 = t210 * (-t137 * t203 - t139) + t257;
t10 = t206 * t22 + t209 * t25;
t8 = t209 * t10 + t206 * t9;
t7 = t10 * t206 - t209 * t9;
t39 = t248 + 0.2e1 * t266;
t19 = t203 * t40 - t204 * t39;
t236 = t203 * t39 + t204 * t40;
t235 = (-t195 + t212) * t210 + t284;
t138 = t203 * t255;
t230 = t138 - t243;
t224 = t136 * t203 + t243;
t156 = t207 * t172;
t222 = t156 + t210 * (t139 + t244);
t92 = t226 + t307;
t220 = t233 + 0.2e1 * t261 + t270;
t31 = t264 * t163 + t227;
t218 = -t257 + t210 * (t136 * t204 - t244);
t89 = t220 + t260;
t216 = t217 + t306;
t215 = -qJD(2) * t179 - t223 - 0.2e1 * t261 - t314;
t214 = t216 + 0.2e1 * t252;
t176 = t195 - t279;
t171 = 0.2e1 * t251 + t259;
t154 = -t228 + t306;
t131 = t207 * t251 + t156;
t130 = t173 * t210;
t125 = t210 * t171 + t275;
t96 = -t118 + t185;
t95 = t117 - t185;
t93 = t137 * t204 - t138;
t91 = -t118 - t185;
t83 = t118 - t117;
t82 = -t185 - t117;
t69 = -qJD(6) * t121 - t246;
t64 = (t119 * t209 - t121 * t206) * t187;
t63 = (-t119 * t206 - t121 * t209) * t187;
t60 = -t117 - t118;
t59 = -t103 + t70;
t54 = (qJD(6) - t187) * t121 + t246;
t53 = t209 * t95 + t300;
t52 = -t206 * t96 + t326;
t51 = -t206 * t95 + t297;
t50 = -t209 * t96 - t328;
t49 = t121 * t281 + t209 * t70;
t48 = t121 * t280 - t206 * t70;
t47 = -t119 * t280 - t206 * t69;
t46 = t119 * t281 - t209 * t69;
t45 = -t206 * t91 + t297;
t44 = t209 * t91 + t300;
t43 = t209 * t82 - t328;
t42 = t206 * t82 + t326;
t38 = t215 - t233 + (-t105 - t255) * pkin(4);
t37 = t215 + t315 + t353;
t35 = t206 * t59 - t209 * t55;
t34 = -t206 * t329 - t209 * t54;
t33 = -t206 * t55 - t209 * t59;
t32 = t206 * t54 - t209 * t329;
t28 = t203 * t45 - t204 * t44;
t27 = -qJ(5) * t100 + t31;
t26 = (-t100 - t195) * pkin(4) + t232;
t23 = t203 * t43 - t204 * t42;
t18 = -pkin(8) * t44 + qJ(5) * t329 - t298;
t16 = t203 * t35 - t204 * t33;
t14 = t203 * t30 - t204 * t31;
t13 = -pkin(8) * t42 + qJ(5) * t54 - t301;
t12 = -pkin(8) * t45 + t311 * t329 + t301;
t11 = -pkin(8) * t43 + t311 * t54 - t298;
t6 = -pkin(8) * t7 - qJ(5) * t36;
t5 = -pkin(8) * t33 + qJ(5) * t60 - t7;
t4 = -pkin(8) * t35 + t311 * t60 - t8;
t3 = -pkin(8) * t8 - t311 * t36;
t1 = t203 * t8 - t204 * t7;
t2 = [0, 0, 0, 0, 0, qJDD(1), t249, t241, 0, 0, t131, t125, t235, t130, -t319, 0, t210 * t154 - t316, -pkin(1) * t171 - t207 * t154 + t340, t247 + t271, pkin(1) * t154 + pkin(7) * t247, 0, -t235, t319, t131, t125, t130, t207 * (qJ(3) * t175 + t226) + (pkin(2) * t175 + t197 + t89) * t210 + t271, t210 * (-pkin(2) * t173 + t191 - t216) + (-t210 * t239 - t275) * qJ(3) + t316, t207 * t214 - t340 + (pkin(1) + t309) * t171 + (t171 + t239) * t278, pkin(7) * (t207 * t92 + t210 * t89) + (pkin(1) - t240) * (t214 + t325), t242, -t356, t344, t218, -t357, t222, t207 * (-t248 + t348) + t210 * (t302 + t308) + t342, t207 * (-t40 - t361) + t210 * (-t303 + t341) + t362, pkin(3) * t299 + t210 * (-t236 + t354) + t355, t225 * t236 + t312 * (t19 * t207 + t210 * t73), t242, t344, t356, t222, t357, t218, t207 * (pkin(4) * t324 + qJ(5) * t317 - t163 * t124 - t227 + t348) + t210 * (qJ(5) * t295 - t204 * t38 + t308) + t342, t207 * (pkin(3) * t334 - pkin(4) * t109 + qJ(5) * t107) + t210 * (-t203 * t27 - t204 * t26 + t354) + t355, t207 * (-pkin(4) * t322 + qJ(5) * t323 + t30 + t361) + t210 * (-pkin(4) * t292 - t203 * t37 - t341) - t362, (-pkin(4) * t31 + qJ(5) * t30 + t14 * t312) * t207 + (pkin(4) * t204 + qJ(5) * t203 + t312) * t210 * t41 + t225 * (t203 * t31 + t204 * t30), -t258 + t210 * (-t203 * t49 - t204 * t48), -t207 * t83 + t210 * (-t203 * t34 - t204 * t32), -t207 * t59 + t210 * (-t203 * t52 - t204 * t50), t258 + t210 * (-t203 * t47 - t204 * t46), t207 * t55 + t210 * (-t203 * t53 - t204 * t51), t207 * t164 + t210 * (-t203 * t64 - t204 * t63), t207 * (pkin(3) * t23 + qJ(5) * t43 - t311 * t42 + t9) + t210 * (-pkin(3) * t54 - t204 * t11 - t203 * t13) + pkin(7) * (t207 * t23 - t210 * t54) + t225 * (t203 * t42 + t204 * t43), t207 * (pkin(3) * t28 + qJ(5) * t45 - t311 * t44 + t10) + t210 * (-pkin(3) * t329 - t204 * t12 - t203 * t18) + pkin(7) * (t207 * t28 - t210 * t329) + t225 * (t203 * t44 + t204 * t45), t207 * (pkin(3) * t16 + qJ(5) * t35 - t311 * t33) + t210 * (-pkin(3) * t60 - t203 * t5 - t204 * t4) + pkin(7) * (t16 * t207 - t210 * t60) + t225 * (t203 * t33 + t204 * t35), t207 * (pkin(3) * t1 + qJ(5) * t8 - t311 * t7) + t210 * (pkin(3) * t36 - t203 * t6 - t204 * t3) + pkin(7) * (t1 * t207 + t210 * t36) + t225 * (t203 * t7 + t204 * t8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t256, t176, t259, t256, t194, qJDD(2), -t132, -t133, 0, 0, qJDD(2), -t259, -t194, -t256, t176, t256, (-pkin(2) * t207 + qJ(3) * t210) * qJDD(1), -pkin(2) * t177 + qJ(3) * t183 + t92, -pkin(2) * t182 + (qJDD(2) - t178) * qJ(3) + t220, -pkin(2) * t92 + qJ(3) * t89, t93, -t332, t351, t224, t345, t230, t303 + t347, qJ(3) * t320 + t302 - t360, -t19 + t343, qJ(3) * t73 + t19 * t305, t93, t351, t332, t230, -t345, t224, -qJ(5) * t294 - t203 * t38 + t347, -t203 * t26 + t204 * t27 + t343, t204 * t37 + t250 * t320 + t360, t305 * t14 + (-qJ(5) * t204 - t250) * t41, -t203 * t48 + t204 * t49, -t203 * t32 + t204 * t34, -t203 * t50 + t204 * t52, -t203 * t46 + t204 * t47, -t203 * t51 + t204 * t53, -t203 * t63 + t204 * t64, -qJ(3) * t54 - t203 * t11 + t204 * t13 + t23 * t305, -qJ(3) * t329 - t203 * t12 + t204 * t18 + t28 * t305, -qJ(3) * t60 + t16 * t305 - t203 * t4 + t204 * t5, qJ(3) * t36 + t1 * t305 - t203 * t3 + t204 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t259, t177, t182, t92, 0, 0, 0, 0, 0, 0, t333, -t78, t334, t19, 0, 0, 0, 0, 0, 0, t333, t334, t78, t14, 0, 0, 0, 0, 0, 0, t23, t28, t16, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t105, t320, t100, t73, 0, 0, 0, 0, 0, 0, t105, t100, -t320, t41, 0, 0, 0, 0, 0, 0, -t54, -t329, -t60, t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t324, t109, t322, t31, 0, 0, 0, 0, 0, 0, t42, t44, t33, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t288, t83, t59, -t288, -t55, -t164, -t9, -t10, 0, 0;];
tauJ_reg  = t2;
