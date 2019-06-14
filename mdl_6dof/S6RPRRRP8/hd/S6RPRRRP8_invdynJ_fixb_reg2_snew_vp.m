% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPRRRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 01:49
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPRRRP8_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP8_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP8_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP8_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP8_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP8_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:48:07
% EndTime: 2019-05-06 01:48:24
% DurationCPUTime: 6.87s
% Computational Cost: add. (17461->396), mult. (34575->495), div. (0->0), fcn. (23670->8), ass. (0->269)
t252 = sin(qJ(3));
t255 = cos(qJ(4));
t251 = sin(qJ(4));
t256 = cos(qJ(3));
t322 = t251 * t256;
t225 = (t252 * t255 + t322) * qJD(1);
t222 = qJD(5) + t225;
t357 = t222 ^ 2;
t314 = qJD(1) * t256;
t227 = -t251 * t252 * qJD(1) + t255 * t314;
t246 = qJD(3) + qJD(4);
t250 = sin(qJ(5));
t254 = cos(qJ(5));
t206 = t227 * t250 - t254 * t246;
t358 = t206 ^ 2;
t185 = t358 - t357;
t312 = qJD(1) * qJD(3);
t299 = t256 * t312;
t310 = t252 * qJDD(1);
t233 = -t299 - t310;
t244 = t256 * qJDD(1);
t300 = t252 * t312;
t234 = t244 - t300;
t295 = -t255 * t233 + t251 * t234;
t180 = -qJD(4) * t227 - t295;
t178 = qJDD(5) - t180;
t208 = t227 * t254 + t246 * t250;
t331 = t208 * t206;
t137 = -t331 - t178;
t341 = t137 * t250;
t102 = -t185 * t254 - t341;
t189 = t222 * t208;
t181 = -t225 * qJD(4) + t251 * t233 + t255 * t234;
t309 = qJDD(3) + qJDD(4);
t296 = -t250 * t181 + t254 * t309;
t278 = qJD(5) * t208 - t296;
t121 = -t189 + t278;
t421 = t252 * (t102 * t251 - t121 * t255) - t256 * (t102 * t255 + t121 * t251);
t205 = t208 ^ 2;
t370 = -t205 - t357;
t91 = t254 * t370 + t341;
t420 = pkin(3) * t91;
t419 = pkin(4) * t91;
t418 = pkin(9) * t91;
t340 = t137 * t254;
t93 = -t250 * t370 + t340;
t417 = pkin(9) * t93;
t416 = qJ(2) * t91;
t415 = t251 * t93;
t414 = t255 * t93;
t369 = t205 - t358;
t274 = -t254 * t181 - t250 * t309;
t268 = -t206 * qJD(5) - t274;
t332 = t206 * t222;
t374 = t332 - t268;
t343 = t374 * t250;
t371 = t189 + t278;
t68 = t371 * t254 - t343;
t411 = t252 * (t251 * t68 + t255 * t369) - t256 * (-t251 * t369 + t255 * t68);
t355 = pkin(7) + pkin(1);
t366 = -t331 + t178;
t132 = t254 * t366;
t363 = -t357 - t358;
t373 = t250 * t363 + t132;
t339 = t366 * t250;
t372 = t254 * t363 - t339;
t387 = t251 * t371 + t255 * t372;
t388 = t251 * t372 - t255 * t371;
t401 = t252 * t387 + t256 * t388;
t408 = qJ(2) * t373 - t355 * t401;
t407 = pkin(3) * t388;
t406 = pkin(8) * t388;
t328 = t227 * t225;
t389 = t309 - t328;
t405 = t251 * t389;
t404 = t255 * t389;
t403 = t374 * qJ(6);
t98 = -t185 * t250 + t340;
t402 = -pkin(3) * t373 + pkin(8) * t387;
t364 = t332 + t268;
t186 = -t205 + t357;
t392 = -t186 * t250 + t132;
t400 = t256 * (t251 * t364 + t255 * t392) - t252 * (t251 * t392 - t255 * t364);
t398 = pkin(4) * t373;
t397 = pkin(9) * t372;
t396 = pkin(9) * t373;
t391 = t254 * t186 + t339;
t327 = t246 * t225;
t390 = t181 - t327;
t342 = t374 * t254;
t66 = -t371 * t250 - t342;
t368 = t205 + t358;
t386 = pkin(4) * t368;
t385 = -qJ(6) * t250 - pkin(4);
t382 = t251 * t368;
t377 = t255 * t368;
t311 = qJD(2) * qJD(1);
t245 = 0.2e1 * t311;
t248 = t252 ^ 2;
t259 = qJD(1) ^ 2;
t282 = qJD(3) * pkin(3) - pkin(8) * t314;
t247 = qJDD(1) * qJ(2);
t253 = sin(qJ(1));
t257 = cos(qJ(1));
t291 = t257 * g(1) + t253 * g(2);
t283 = -t247 + t291;
t361 = -t233 * pkin(3) - (pkin(8) * t248 + t355) * t259 + t282 * t314 - t283;
t76 = t245 - t390 * pkin(9) + (t227 * t246 - t180) * pkin(4) + t361;
t290 = t253 * g(1) - t257 * g(2);
t281 = qJDD(2) - t290;
t317 = t259 * qJ(2);
t271 = t281 - t317;
t265 = -qJDD(1) * t355 + t271;
t202 = t256 * g(3) - t252 * t265;
t324 = t248 * t259;
t167 = -pkin(3) * t324 + t233 * pkin(8) - qJD(3) * t282 - t202;
t320 = t255 * t167;
t263 = t256 * t265;
t261 = -t234 * pkin(8) + t263;
t318 = t256 * t259;
t367 = qJDD(3) * pkin(3) + t261 + (-pkin(3) * t318 - pkin(8) * t312 + g(3)) * t252;
t131 = t251 * t367 + t320;
t196 = pkin(4) * t225 - pkin(9) * t227;
t356 = t246 ^ 2;
t82 = -pkin(4) * t356 + pkin(9) * t309 - t225 * t196 + t131;
t45 = t250 * t82 - t254 * t76;
t46 = t250 * t76 + t254 * t82;
t20 = t250 * t45 + t254 * t46;
t168 = pkin(5) * t206 - qJ(6) * t208;
t294 = t178 * qJ(6) - t206 * t168 + t46;
t362 = -pkin(5) * (t370 + t357) - qJ(6) * t137 + t294;
t329 = t222 * t254;
t303 = t206 * t329;
t279 = t250 * t278 + t303;
t304 = t255 * t331;
t305 = t251 * t331;
t360 = t256 * (t255 * t279 - t305) - t252 * (t251 * t279 + t304);
t330 = t222 * t250;
t183 = t208 * t330;
t288 = t183 - t303;
t359 = t256 * (t178 * t251 + t255 * t288) - t252 * (-t255 * t178 + t251 * t288);
t223 = t225 ^ 2;
t224 = t227 ^ 2;
t354 = t252 * g(3);
t130 = t251 * t167 - t255 * t367;
t81 = -t309 * pkin(4) - t356 * pkin(9) + t227 * t196 + t130;
t353 = -pkin(4) * t81 + pkin(9) * t20;
t77 = t250 * t81;
t78 = t254 * t81;
t72 = -t130 * t255 + t131 * t251;
t352 = t256 * t72;
t350 = qJ(6) * t254;
t349 = qJDD(1) * pkin(1);
t344 = t364 * t250;
t306 = -0.2e1 * t311;
t169 = t306 - t361;
t338 = t169 * t251;
t337 = t169 * t255;
t194 = t328 + t309;
t334 = t194 * t251;
t333 = t194 * t255;
t326 = t246 * t251;
t325 = t246 * t255;
t249 = t256 ^ 2;
t323 = t249 * t259;
t302 = t252 * t318;
t321 = t252 * (qJDD(3) + t302);
t238 = qJDD(3) - t302;
t319 = t256 * t238;
t315 = t248 + t249;
t313 = qJD(6) * t222;
t128 = (qJD(5) + t222) * t206 + t274;
t308 = pkin(4) * t128 + t417 + t77;
t307 = -pkin(4) * t371 + t397 - t78;
t301 = -pkin(4) * t255 - pkin(3);
t216 = 0.2e1 * t313;
t287 = t216 + t294;
t26 = (t368 - t357) * pkin(5) + t287;
t38 = -t178 * pkin(5) - qJ(6) * t357 + t168 * t208 + qJDD(6) + t45;
t30 = qJ(6) * t368 + t38;
t69 = -t121 * t254 + t344;
t298 = pkin(9) * t69 + t250 * t30 + t254 * t26 + t386;
t123 = (-qJD(5) + t222) * t208 + t296;
t71 = t123 * t254 + t344;
t297 = pkin(9) * t71 + t20 + t386;
t73 = t130 * t251 + t255 * t131;
t37 = -pkin(5) * t357 + t287;
t292 = -pkin(5) * t38 + qJ(6) * t37;
t289 = t206 * t330 - t254 * t278;
t286 = -pkin(5) * t364 - qJ(6) * t121;
t12 = t20 * t251 - t255 * t81;
t4 = t12 * t256 + (t20 * t255 + t251 * t81) * t252;
t285 = t250 * t46 - t254 * t45;
t201 = t263 + t354;
t163 = t256 * t201 - t252 * t202;
t267 = t278 * pkin(5) + t403 + t81;
t266 = 0.2e1 * qJD(6) * t208 - t267;
t33 = -pkin(5) * t189 + t266 - t403;
t280 = -pkin(4) * t374 - pkin(5) * t342 + t250 * t33 - t417;
t34 = (-t371 - t189) * pkin(5) + t266;
t277 = t254 * t34 + t385 * t371 + t397;
t275 = (-t206 * t250 - t208 * t254) * t222;
t273 = (-qJD(4) + t246) * t227 - t295;
t10 = t250 * t38 + t254 * t37;
t42 = (pkin(5) * t222 - 0.2e1 * qJD(6)) * t208 + t267;
t272 = pkin(9) * t10 + (-pkin(5) * t254 + t385) * t42;
t114 = t254 * t268 - t183;
t264 = t256 * (t255 * t114 + t305) - t252 * (t251 * t114 - t304);
t262 = pkin(5) * t366 + qJ(6) * t363 - t38;
t258 = qJD(3) ^ 2;
t236 = t315 * qJDD(1);
t235 = t244 - 0.2e1 * t300;
t232 = 0.2e1 * t299 + t310;
t221 = -t271 + t349;
t215 = -t224 + t356;
t214 = t223 - t356;
t213 = t259 * t355 + t283 + t306;
t211 = -t224 - t356;
t210 = -t321 + t256 * (-t258 - t323);
t209 = t252 * (-t258 - t324) + t319;
t198 = t224 - t223;
t192 = -t356 - t223;
t182 = -t223 - t224;
t162 = -t211 * t251 - t333;
t161 = t211 * t255 - t334;
t160 = t327 + t181;
t155 = (qJD(4) + t246) * t227 + t295;
t152 = t192 * t255 - t405;
t151 = t192 * t251 + t404;
t115 = t254 * t364;
t113 = t208 * t329 + t250 * t268;
t106 = t161 * t256 + t162 * t252;
t105 = t160 * t251 + t255 * t273;
t104 = -t160 * t255 + t251 * t273;
t95 = t151 * t256 + t152 * t252;
t67 = t123 * t250 - t115;
t65 = -t121 * t250 - t115;
t61 = t104 * t256 + t105 * t252;
t59 = -t128 * t251 + t414;
t57 = t128 * t255 + t415;
t55 = t251 * t374 - t414;
t53 = -t255 * t374 - t415;
t52 = t255 * t71 - t382;
t51 = t255 * t69 - t382;
t50 = t251 * t71 + t377;
t49 = t251 * t69 + t377;
t48 = t78 - t418;
t47 = t77 - t396;
t40 = t252 * t73 + t352;
t39 = -pkin(4) * t65 - t286;
t36 = t46 - t419;
t35 = t45 - t398;
t28 = t252 * t59 + t256 * t57;
t24 = t252 * t55 + t256 * t53;
t22 = t252 * t52 + t256 * t50;
t21 = t252 * t51 + t256 * t49;
t17 = -t262 - t398;
t16 = -0.2e1 * t313 - t362 + t419;
t15 = -t250 * t34 - t350 * t371 - t396;
t14 = pkin(5) * t343 + t254 * t33 + t418;
t11 = -pkin(9) * t67 - t285;
t9 = t250 * t37 - t254 * t38;
t7 = -pkin(9) * t65 - t250 * t26 + t254 * t30;
t6 = t10 * t255 + t251 * t42;
t5 = t10 * t251 - t255 * t42;
t3 = -pkin(9) * t9 + (pkin(5) * t250 - t350) * t42;
t2 = -pkin(4) * t9 - t292;
t1 = t252 * t6 + t256 * t5;
t8 = [0, 0, 0, 0, 0, qJDD(1), t290, t291, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t281 - 0.2e1 * t349, t245 + 0.2e1 * t247 - t291, pkin(1) * t221 + qJ(2) * (-t259 * pkin(1) + t245 - t283), (t234 - t300) * t256, -t232 * t256 - t235 * t252, t319 - t252 * (t258 - t323), (-t233 + t299) * t252, t256 * (-t258 + t324) - t321, 0, qJ(2) * t232 - t209 * t355 - t252 * t213, qJ(2) * t235 - t210 * t355 - t256 * t213, t236 * t355 - t315 * t317 - t163, -qJ(2) * t213 - t163 * t355, t256 * (t181 * t255 - t227 * t326) - t252 * (t181 * t251 + t227 * t325), t256 * (-t155 * t255 - t251 * t390) - t252 * (-t155 * t251 + t255 * t390), t256 * (-t215 * t251 + t404) - t252 * (t215 * t255 + t405), t256 * (-t180 * t251 + t225 * t325) - t252 * (t180 * t255 + t225 * t326), t256 * (t214 * t255 - t334) - t252 * (t214 * t251 + t333), (t256 * (-t225 * t255 + t227 * t251) - t252 * (-t225 * t251 - t227 * t255)) * t246, t256 * (-pkin(8) * t151 - t338) - t252 * (-pkin(3) * t155 + pkin(8) * t152 + t337) + qJ(2) * t155 - t355 * t95, t256 * (-pkin(8) * t161 - t337) - t252 * (-pkin(3) * t390 + pkin(8) * t162 - t338) + qJ(2) * t390 - t355 * t106, t256 * (-pkin(8) * t104 - t72) - t252 * (-pkin(3) * t182 + pkin(8) * t105 + t73) + qJ(2) * t182 - t355 * t61, -pkin(8) * t352 - t252 * (pkin(3) * t169 + pkin(8) * t73) - qJ(2) * t169 - t355 * t40, t264, t411, t400, t360, t421, t359, t256 * (-t251 * t35 + t255 * t47 - t406) - t252 * (t251 * t47 + t255 * t35 + t402) + t408, t256 * (-pkin(8) * t57 - t251 * t36 + t255 * t48) - t252 * (pkin(8) * t59 + t251 * t48 + t255 * t36 - t420) + t416 - t355 * t28, t256 * (-pkin(8) * t50 + t11 * t255) - t252 * (pkin(8) * t52 + t251 * t11) + (pkin(4) * t322 - t252 * t301 + qJ(2)) * t67 - t355 * t22, (t256 * (pkin(4) * t251 - pkin(9) * t255) - t252 * (-pkin(9) * t251 + t301) + qJ(2)) * t285 + (-t355 - pkin(8)) * t4, t264, t400, -t411, t359, -t421, t360, t256 * (t15 * t255 - t17 * t251 - t406) - t252 * (t15 * t251 + t17 * t255 + t402) + t408, t256 * (-pkin(8) * t49 - t251 * t39 + t255 * t7) - t252 * (-pkin(3) * t65 + pkin(8) * t51 + t251 * t7 + t255 * t39) + qJ(2) * t65 - t355 * t21, t256 * (-pkin(8) * t53 + t14 * t255 - t16 * t251) - t252 * (pkin(8) * t55 + t14 * t251 + t16 * t255 + t420) - t416 - t355 * t24, t256 * (-pkin(8) * t5 - t2 * t251 + t255 * t3) - t252 * (-pkin(3) * t9 + pkin(8) * t6 + t2 * t255 + t251 * t3) + qJ(2) * t9 - t355 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t259, -t221, 0, 0, 0, 0, 0, 0, t209, t210, -t236, t163, 0, 0, 0, 0, 0, 0, t95, t106, t61, t40, 0, 0, 0, 0, 0, 0, t401, t28, t22, t4, 0, 0, 0, 0, 0, 0, t401, t21, t24, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t302, (-t248 + t249) * t259, t244, -t302, -t310, qJDD(3), t201, t202, 0, 0, t328, t198, t160, -t328, t273, t309, pkin(3) * t151 - t130, -t320 - t251 * (-pkin(8) * t300 + t261 + t354) + (-t238 * t251 + t161) * pkin(3), pkin(3) * t104, pkin(3) * t72, t113, t66, t391, t289, -t98, t275, t307 + t407, pkin(3) * t57 + t308, pkin(3) * t50 + t297, pkin(3) * t12 + t353, t113, t391, -t66, t275, t98, t289, t277 + t407, pkin(3) * t49 + t298, pkin(3) * t53 + t280, pkin(3) * t5 + t272; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t328, t198, t160, -t328, t273, t309, -t130, -t131, 0, 0, t113, t66, t391, t289, -t98, t275, t307, t308, t297, t353, t113, t391, -t66, t275, t98, t289, t277, t298, t280, t272; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t331, t369, t364, -t331, -t121, t178, -t45, -t46, 0, 0, t331, t364, -t369, t178, t121, -t331, t262, t286, t216 + t362, t292; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t366, t364, t370, t38;];
tauJ_reg  = t8;
