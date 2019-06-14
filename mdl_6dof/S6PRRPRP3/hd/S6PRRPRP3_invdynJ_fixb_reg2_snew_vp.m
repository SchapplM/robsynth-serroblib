% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6PRRPRP3
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
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 03:59
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6PRRPRP3_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP3_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP3_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRP3_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP3_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP3_invdynJ_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 03:56:06
% EndTime: 2019-05-05 03:56:25
% DurationCPUTime: 8.74s
% Computational Cost: add. (18315->412), mult. (38289->564), div. (0->0), fcn. (27979->12), ass. (0->262)
t239 = sin(pkin(6));
t241 = cos(pkin(6));
t244 = sin(qJ(3));
t245 = sin(qJ(2));
t247 = cos(qJ(3));
t248 = cos(qJ(2));
t238 = sin(pkin(11));
t240 = cos(pkin(11));
t289 = qJD(2) * t244;
t206 = -t240 * qJD(3) + t238 * t289;
t207 = qJD(3) * t238 + t240 * t289;
t243 = sin(qJ(5));
t246 = cos(qJ(5));
t179 = t246 * t206 + t207 * t243;
t286 = qJD(2) * qJD(3);
t280 = t247 * t286;
t285 = t244 * qJDD(2);
t211 = t280 + t285;
t188 = qJDD(3) * t238 + t211 * t240;
t264 = qJDD(3) * t240 - t211 * t238;
t254 = -t179 * qJD(5) + t246 * t188 + t243 * t264;
t288 = qJD(2) * t247;
t226 = -qJD(5) + t288;
t300 = t179 * t226;
t334 = t254 + t300;
t229 = t244 * t286;
t284 = t247 * qJDD(2);
t212 = -t229 + t284;
t208 = -qJDD(5) + t212;
t181 = -t206 * t243 + t207 * t246;
t299 = t181 * t179;
t256 = t208 - t299;
t306 = t256 * t243;
t178 = t181 ^ 2;
t322 = t226 ^ 2;
t336 = -t178 - t322;
t81 = t246 * t336 + t306;
t305 = t256 * t246;
t83 = -t243 * t336 + t305;
t50 = t238 * t81 - t240 * t83;
t38 = -t244 * t334 + t247 * t50;
t56 = t238 * t83 + t240 * t81;
t397 = (t245 * t38 + t248 * t56) * t239 + t241 * (t244 * t50 + t247 * t334);
t396 = pkin(2) * t56 + pkin(8) * t38;
t392 = qJ(4) * t56;
t390 = -pkin(3) * t56 - pkin(4) * t81;
t389 = pkin(3) * t334 + qJ(4) * t50;
t164 = t181 * t226;
t277 = -t243 * t188 + t246 * t264;
t261 = qJD(5) * t181 - t277;
t103 = t164 + t261;
t324 = t179 ^ 2;
t156 = t324 - t322;
t94 = t156 * t243 - t305;
t98 = t156 * t246 + t306;
t387 = -t247 * t103 + t244 * (t238 * t94 - t240 * t98);
t335 = t178 - t324;
t337 = -t164 + t261;
t66 = -t337 * t243 + t246 * t334;
t310 = t334 * t243;
t70 = -t337 * t246 - t310;
t386 = t244 * (-t238 * t66 + t240 * t70) - t247 * t335;
t384 = pkin(9) * t81;
t383 = pkin(9) * t83;
t133 = t208 + t299;
t304 = t133 * t243;
t332 = -t322 - t324;
t340 = t246 * t332 + t304;
t131 = t246 * t133;
t341 = t243 * t332 - t131;
t353 = t238 * t340 + t240 * t341;
t354 = -t238 * t341 + t240 * t340;
t369 = t244 * t337 + t247 * t354;
t382 = -pkin(2) * t353 + pkin(8) * t369;
t381 = t238 * t70 + t240 * t66;
t379 = t238 * t98 + t240 * t94;
t377 = t241 * (t244 * t354 - t247 * t337) + (t245 * t369 - t248 * t353) * t239;
t373 = qJ(4) * t353;
t371 = -pkin(3) * t353 - pkin(4) * t341;
t370 = -pkin(3) * t337 + qJ(4) * t354;
t157 = -t178 + t322;
t355 = t246 * t157 - t304;
t356 = -t157 * t243 - t131;
t367 = t238 * t356 + t240 * t355;
t333 = -t300 + t254;
t366 = t244 * (-t238 * t355 + t240 * t356) - t247 * t333;
t113 = -t324 - t178;
t365 = pkin(3) * t113;
t364 = pkin(4) * t113;
t362 = pkin(9) * t340;
t361 = pkin(9) * t341;
t359 = qJ(6) * t334;
t314 = sin(pkin(10));
t315 = cos(pkin(10));
t257 = t314 * g(1) - t315 * g(2);
t291 = -g(3) + qJDD(1);
t253 = -t239 * t257 + t241 * t291;
t186 = t247 * t253;
t216 = -t315 * g(1) - t314 * g(2);
t255 = t241 * t257;
t342 = t239 * t291 + t255;
t174 = t248 * t216 + t245 * t342;
t249 = qJD(2) ^ 2;
t154 = -t249 * pkin(2) + qJDD(2) * pkin(8) + t174;
t267 = -pkin(3) * t247 - qJ(4) * t244;
t275 = t249 * t267 + t154;
t321 = qJD(3) ^ 2;
t116 = -qJDD(3) * pkin(3) - t321 * qJ(4) + t244 * t275 + qJDD(4) - t186;
t189 = -pkin(4) * t288 - pkin(9) * t207;
t323 = t206 ^ 2;
t75 = -t264 * pkin(4) - t323 * pkin(9) + t207 * t189 + t116;
t360 = pkin(5) * t261 - t359 + t75;
t358 = t113 * t244;
t357 = t113 * t247;
t298 = t206 * t207;
t259 = -t212 - t298;
t339 = t238 * t259;
t338 = t240 * t259;
t193 = t206 * t288;
t169 = -t188 + t193;
t194 = t207 * t288;
t167 = t264 - t194;
t144 = pkin(5) * t179 - qJ(6) * t181;
t251 = t244 * t253;
t117 = -t321 * pkin(3) + qJDD(3) * qJ(4) + t247 * t275 + t251;
t269 = t245 * t216 - t248 * t342;
t153 = -qJDD(2) * pkin(2) - t249 * pkin(8) + t269;
t265 = t211 + t280;
t129 = -t265 * qJ(4) + (-t212 + t229) * pkin(3) + t153;
t73 = 0.2e1 * qJD(4) * t207 + t238 * t117 - t240 * t129;
t55 = t259 * pkin(4) + pkin(9) * t169 - t73;
t74 = -0.2e1 * qJD(4) * t206 + t240 * t117 + t238 * t129;
t63 = -t323 * pkin(4) + pkin(9) * t264 + t189 * t288 + t74;
t32 = t243 * t55 + t246 * t63;
t276 = -t208 * qJ(6) - t179 * t144 + t32;
t329 = -pkin(5) * (t336 + t322) - qJ(6) * t256 + t276;
t258 = (t179 * t243 + t181 * t246) * t226;
t296 = t226 * t243;
t155 = t181 * t296;
t295 = t226 * t246;
t283 = t179 * t295;
t270 = -t155 + t283;
t328 = t238 * t270 + t240 * t258;
t262 = t243 * t261 - t283;
t271 = -t179 * t296 - t246 * t261;
t327 = t238 * t262 + t240 * t271;
t326 = t244 * (-t238 * t258 + t240 * t270) + t247 * t208;
t282 = t247 * t299;
t325 = t244 * (-t238 * t271 + t240 * t262) + t282;
t205 = t207 ^ 2;
t320 = pkin(5) * t246;
t31 = t243 * t63 - t246 * t55;
t16 = t243 * t32 - t246 * t31;
t319 = t16 * t238;
t318 = t16 * t240;
t317 = t243 * t75;
t316 = t246 * t75;
t313 = qJ(6) * t246;
t309 = t333 * t243;
t308 = t116 * t238;
t307 = t116 * t240;
t171 = t212 - t298;
t302 = t171 * t238;
t301 = t171 * t240;
t225 = t244 * t249 * t247;
t218 = qJDD(3) + t225;
t294 = t244 * t218;
t217 = -t225 + qJDD(3);
t292 = t247 * t217;
t287 = qJD(6) * t226;
t281 = t247 * t298;
t279 = -qJ(6) * t243 - pkin(4);
t46 = t238 * t73 + t240 * t74;
t17 = t243 * t31 + t246 * t32;
t142 = t244 * t154 - t186;
t143 = t247 * t154 + t251;
t80 = t142 * t244 + t247 * t143;
t90 = -t181 * t295 + t243 * t254;
t91 = t246 * t254 + t155;
t273 = t244 * (-t238 * t90 + t240 * t91) - t282;
t219 = -0.2e1 * t287;
t268 = t219 + t276;
t24 = -pkin(5) * t322 + t268;
t25 = t208 * pkin(5) - qJ(6) * t322 + t144 * t181 + qJDD(6) + t31;
t272 = -pkin(5) * t25 + qJ(6) * t24;
t266 = -pkin(5) * t333 - qJ(6) * t103;
t45 = t238 * t74 - t240 * t73;
t263 = -pkin(2) + t267;
t252 = -pkin(5) * t133 + qJ(6) * t332 - t25;
t250 = 0.2e1 * qJD(6) * t181 - t360;
t235 = t247 ^ 2;
t234 = t244 ^ 2;
t232 = t235 * t249;
t231 = t234 * t249;
t223 = -t232 - t321;
t222 = -t231 - t321;
t215 = t231 + t232;
t214 = (t234 + t235) * qJDD(2);
t213 = -0.2e1 * t229 + t284;
t210 = 0.2e1 * t280 + t285;
t201 = t247 * t212;
t192 = -t205 - t232;
t191 = -t205 + t232;
t190 = -t232 + t323;
t185 = -t222 * t244 - t292;
t184 = t223 * t247 - t294;
t182 = -t232 - t323;
t168 = t188 + t193;
t166 = -t194 - t264;
t160 = -t205 - t323;
t149 = -t192 * t238 + t301;
t148 = t192 * t240 + t302;
t138 = t182 * t240 - t339;
t137 = t182 * t238 + t338;
t128 = t167 * t240 - t169 * t238;
t127 = t167 * t238 + t169 * t240;
t111 = t149 * t247 + t168 * t244;
t110 = t138 * t247 + t166 * t244;
t104 = (-qJD(5) - t226) * t181 + t277;
t100 = t246 * t333;
t85 = t128 * t247 + t160 * t244;
t71 = t104 * t246 + t309;
t69 = -t103 * t246 + t309;
t67 = t104 * t243 - t100;
t65 = -t103 * t243 - t100;
t61 = t238 * t91 + t240 * t90;
t52 = t316 - t384;
t47 = t317 - t361;
t42 = -pkin(4) * t334 + t317 + t383;
t41 = (-pkin(5) * t226 - 0.2e1 * qJD(6)) * t181 + t360;
t40 = -pkin(4) * t337 - t316 + t362;
t37 = t116 * t244 + t247 * t46;
t36 = -t238 * t67 + t240 * t71;
t35 = -t238 * t65 + t240 * t69;
t34 = t238 * t71 + t240 * t67;
t33 = t238 * t69 + t240 * t65;
t29 = (-t337 + t164) * pkin(5) + t250;
t28 = pkin(5) * t164 + t250 + t359;
t27 = t247 * t36 + t358;
t26 = t247 * t35 + t358;
t23 = -qJ(6) * t113 + t25;
t22 = (-t113 - t322) * pkin(5) + t268;
t21 = -t243 * t29 - t313 * t337 - t361;
t20 = -pkin(5) * t310 + t246 * t28 + t384;
t19 = t246 * t29 + t279 * t337 + t362;
t18 = -t383 + t243 * t28 + (pkin(4) + t320) * t334;
t15 = -pkin(4) * t75 + pkin(9) * t17;
t14 = -pkin(9) * t67 - t16;
t13 = t24 * t246 + t243 * t25;
t12 = t24 * t243 - t246 * t25;
t11 = pkin(9) * t71 + t17 - t364;
t10 = -pkin(9) * t65 - t22 * t243 + t23 * t246;
t9 = pkin(9) * t69 + t22 * t246 + t23 * t243 - t364;
t8 = t17 * t240 - t319;
t7 = t17 * t238 + t318;
t6 = -pkin(9) * t12 + (pkin(5) * t243 - t313) * t41;
t5 = t244 * t75 + t247 * t8;
t4 = pkin(9) * t13 + (t279 - t320) * t41;
t3 = -t12 * t238 + t13 * t240;
t2 = t12 * t240 + t13 * t238;
t1 = t244 * t41 + t247 * t3;
t30 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t291, 0, 0, 0, 0, 0, 0, (qJDD(2) * t248 - t245 * t249) * t239, (-qJDD(2) * t245 - t248 * t249) * t239, 0, t241 ^ 2 * t291 + (t245 * t174 - t248 * t269 - t255) * t239, 0, 0, 0, 0, 0, 0, t241 * (t218 * t247 + t223 * t244) + (t184 * t245 + t213 * t248) * t239, t241 * (-t217 * t244 + t222 * t247) + (t185 * t245 - t210 * t248) * t239, (t214 * t245 + t215 * t248) * t239, t241 * (-t142 * t247 + t143 * t244) + (-t153 * t248 + t245 * t80) * t239, 0, 0, 0, 0, 0, 0, t241 * (t138 * t244 - t166 * t247) + (t110 * t245 - t137 * t248) * t239, t241 * (t149 * t244 - t168 * t247) + (t111 * t245 - t148 * t248) * t239, t241 * (t128 * t244 - t160 * t247) + (-t127 * t248 + t245 * t85) * t239, t241 * (-t116 * t247 + t244 * t46) + (t245 * t37 - t248 * t45) * t239, 0, 0, 0, 0, 0, 0, t377, -t397, t241 * (t244 * t36 - t357) + (t245 * t27 - t248 * t34) * t239, t241 * (t244 * t8 - t247 * t75) + (t245 * t5 - t248 * t7) * t239, 0, 0, 0, 0, 0, 0, t377, t241 * (t244 * t35 - t357) + (t245 * t26 - t248 * t33) * t239, t397, t241 * (t244 * t3 - t247 * t41) + (t1 * t245 - t2 * t248) * t239; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t269, -t174, 0, 0, t265 * t244, t210 * t247 + t213 * t244, t294 + t247 * (-t231 + t321), -t244 * t280 + t201, t244 * (t232 - t321) + t292, 0, pkin(2) * t213 + pkin(8) * t184 - t153 * t247, -pkin(2) * t210 + pkin(8) * t185 + t153 * t244, pkin(2) * t215 + pkin(8) * t214 + t80, -pkin(2) * t153 + pkin(8) * t80, t244 * (t188 * t240 + t238 * t194) - t281, t244 * (-t166 * t240 - t168 * t238) + t247 * (-t205 + t323), t244 * (-t191 * t238 + t338) + t247 * t169, t244 * (-t193 * t240 - t238 * t264) + t281, t244 * (t190 * t240 + t302) - t247 * t167, t201 + t244 * (t206 * t240 - t207 * t238) * t288, t244 * (-qJ(4) * t137 + t308) + t247 * (-pkin(3) * t137 + t73) - pkin(2) * t137 + pkin(8) * t110, t244 * (-qJ(4) * t148 + t307) + t247 * (-pkin(3) * t148 + t74) - pkin(2) * t148 + pkin(8) * t111, pkin(8) * t85 + t127 * t263 - t244 * t45, pkin(8) * t37 + t263 * t45, t273, t386, t366, t325, -t387, t326, t244 * (-t238 * t40 + t240 * t47 - t373) + t247 * (t31 + t371) + t382, t244 * (-t238 * t42 + t240 * t52 - t392) + t247 * (t32 + t390) - t396, t244 * (-qJ(4) * t34 - t11 * t238 + t14 * t240) + t247 * (-pkin(3) * t34 - pkin(4) * t67) - pkin(2) * t34 + pkin(8) * t27, t244 * (-pkin(9) * t318 - qJ(4) * t7 - t15 * t238) + t247 * (-pkin(3) * t7 - pkin(4) * t16) - pkin(2) * t7 + pkin(8) * t5, t273, t366, -t386, t326, t387, t325, t244 * (-t19 * t238 + t21 * t240 - t373) + t247 * (-t252 + t371) + t382, t244 * (-qJ(4) * t33 + t10 * t240 - t238 * t9) + t247 * (-pkin(3) * t33 - pkin(4) * t65 - t266) - pkin(2) * t33 + pkin(8) * t26, t244 * (-t18 * t238 + t20 * t240 + t392) + t247 * (0.2e1 * t287 - t329 - t390) + t396, t244 * (-qJ(4) * t2 - t238 * t4 + t240 * t6) + t247 * (-pkin(3) * t2 - pkin(4) * t12 - t272) - pkin(2) * t2 + pkin(8) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t225, t231 - t232, t285, t225, t284, qJDD(3), -t142, -t143, 0, 0, t188 * t238 - t194 * t240, -t166 * t238 + t168 * t240, t191 * t240 + t339, -t193 * t238 + t240 * t264, t190 * t238 - t301, (t206 * t238 + t207 * t240) * t288, -pkin(3) * t166 + qJ(4) * t138 - t307, -pkin(3) * t168 + qJ(4) * t149 + t308, -pkin(3) * t160 + qJ(4) * t128 + t46, -pkin(3) * t116 + qJ(4) * t46, t61, t381, t367, t327, t379, t328, t238 * t47 + t240 * t40 + t370, t238 * t52 + t240 * t42 - t389, qJ(4) * t36 + t11 * t240 + t14 * t238 - t365, -pkin(3) * t75 - pkin(9) * t319 + qJ(4) * t8 + t15 * t240, t61, t367, -t381, t328, -t379, t327, t19 * t240 + t21 * t238 + t370, qJ(4) * t35 + t10 * t238 + t240 * t9 - t365, t18 * t240 + t20 * t238 + t389, -pkin(3) * t41 + qJ(4) * t3 + t238 * t6 + t240 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t166, t168, t160, t116, 0, 0, 0, 0, 0, 0, t337, t334, t113, t75, 0, 0, 0, 0, 0, 0, t337, t113, -t334, t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t299, t335, t333, -t299, -t103, -t208, -t31, -t32, 0, 0, t299, t333, -t335, -t208, t103, -t299, t252, t266, t219 + t329, t272; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t133, t333, t336, t25;];
tauJ_reg  = t30;
