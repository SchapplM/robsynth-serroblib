% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRRPRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% tau_reg [6x33]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRRPRR2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR2_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR2_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR2_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:08:44
% EndTime: 2019-03-09 18:08:58
% DurationCPUTime: 6.65s
% Computational Cost: add. (12014->522), mult. (28740->687), div. (0->0), fcn. (22647->18), ass. (0->309)
t273 = qJ(2) + qJ(3);
t260 = pkin(11) + t273;
t248 = sin(t260);
t279 = sin(qJ(1));
t284 = cos(qJ(1));
t325 = g(1) * t284 + g(2) * t279;
t455 = t325 * t248;
t277 = sin(qJ(3));
t278 = sin(qJ(2));
t282 = cos(qJ(3));
t283 = cos(qJ(2));
t210 = t277 * t278 - t282 * t283;
t192 = t210 * qJD(1);
t371 = qJD(1) * t278;
t381 = t277 * t283;
t194 = -qJD(1) * t381 - t282 * t371;
t274 = sin(pkin(11));
t401 = cos(pkin(11));
t336 = -t401 * t192 + t194 * t274;
t431 = qJD(5) + qJD(6);
t454 = t336 - t431;
t429 = pkin(7) + pkin(8);
t230 = t429 * t283;
t222 = qJD(1) * t230;
t199 = t282 * t222;
t229 = t429 * t278;
t220 = qJD(1) * t229;
t334 = t220 * t277 - t199;
t400 = qJ(4) * t192;
t138 = t334 + t400;
t186 = t194 * qJ(4);
t195 = t277 * t222;
t375 = -t282 * t220 - t195;
t139 = t186 + t375;
t386 = t274 * t277;
t414 = pkin(2) * qJD(3);
t403 = -t274 * t138 - t139 * t401 + (t282 * t401 - t386) * t414;
t276 = sin(qJ(5));
t368 = qJD(5) * t276;
t394 = t336 * t276;
t453 = t368 - t394;
t266 = t283 * pkin(2);
t418 = pkin(1) + t266;
t330 = pkin(3) * t210 - t418;
t262 = sin(t273);
t264 = cos(t273);
t452 = -g(3) * t264 + t262 * t325;
t228 = qJD(1) * t418;
t451 = qJDD(1) * t418;
t281 = cos(qJ(5));
t304 = -t274 * t192 - t194 * t401;
t428 = pkin(3) * t194;
t102 = pkin(4) * t304 - pkin(9) * t336 - t428;
t258 = pkin(2) * t371;
t99 = t102 + t258;
t450 = -t276 * t403 - t281 * t99;
t269 = qJD(2) + qJD(3);
t135 = -t281 * t269 + t276 * t304;
t137 = t269 * t276 + t281 * t304;
t275 = sin(qJ(6));
t280 = cos(qJ(6));
t316 = t135 * t275 - t280 * t137;
t91 = t280 * t135 + t137 * t275;
t449 = t316 * t91;
t360 = t283 * qJDD(1);
t363 = qJD(1) * qJD(2);
t347 = t283 * t363;
t361 = t278 * qJDD(1);
t303 = -t347 - t361;
t362 = qJD(1) * qJD(3);
t443 = -t283 * t362 + t303;
t132 = (-t269 * t371 + t360) * t277 - t443 * t282;
t267 = qJDD(2) + qJDD(3);
t162 = qJDD(2) * pkin(2) + t303 * t429;
t348 = t278 * t363;
t302 = -t348 + t360;
t166 = t429 * t302;
t413 = qJD(2) * pkin(2);
t201 = -t220 + t413;
t315 = -t201 * t277 - t199;
t296 = qJD(3) * t315 + t282 * t162 - t277 * t166;
t57 = pkin(3) * t267 - qJ(4) * t132 + qJD(4) * t194 + t296;
t370 = qJD(3) * t277;
t190 = t222 * t370;
t331 = -qJD(3) * t201 - t166;
t63 = -t192 * qJD(4) - t190 + (qJ(4) * t443 + t162) * t277 + ((-t278 * t362 + t302) * qJ(4) - t331) * t282;
t22 = -t274 * t63 + t401 * t57;
t20 = -t267 * pkin(4) - t22;
t249 = cos(t260);
t448 = g(3) * t249 + t20;
t432 = qJD(5) - t336;
t447 = t135 * t432;
t212 = t278 * t282 + t381;
t161 = t269 * t212;
t385 = t275 * t281;
t211 = t276 * t280 + t385;
t445 = t454 * t211;
t209 = t275 * t276 - t280 * t281;
t444 = t454 * t209;
t332 = t432 * t281;
t290 = qJD(1) * t161;
t288 = -t210 * qJDD(1) - t290;
t83 = -t132 * t274 + t401 * t288;
t80 = qJDD(5) - t83;
t441 = -t276 * t80 - t432 * t332;
t440 = t316 ^ 2 - t91 ^ 2;
t143 = qJD(6) + t432;
t365 = qJD(6) * t280;
t366 = qJD(6) * t275;
t367 = qJD(5) * t281;
t84 = t401 * t132 + t274 * t288;
t60 = t276 * t267 + t269 * t367 + t281 * t84 - t304 * t368;
t337 = -t281 * t267 + t276 * t84;
t61 = qJD(5) * t137 + t337;
t16 = -t135 * t365 - t137 * t366 - t275 * t61 + t280 * t60;
t439 = t143 * t91 + t16;
t272 = qJ(5) + qJ(6);
t263 = cos(t272);
t388 = t263 * t279;
t261 = sin(t272);
t389 = t261 * t284;
t169 = -t249 * t388 + t389;
t387 = t263 * t284;
t390 = t261 * t279;
t171 = t249 * t387 + t390;
t335 = t282 * t201 - t195;
t130 = t186 + t335;
t121 = pkin(3) * t269 + t130;
t131 = -t315 - t400;
t340 = t401 * t131;
t73 = t274 * t121 + t340;
t71 = pkin(9) * t269 + t73;
t163 = t192 * pkin(3) + qJD(4) - t228;
t86 = -pkin(4) * t336 - pkin(9) * t304 + t163;
t43 = t276 * t86 + t281 * t71;
t38 = -pkin(10) * t135 + t43;
t35 = t38 * t366;
t421 = g(3) * t263;
t122 = t274 * t131;
t72 = t121 * t401 - t122;
t70 = -t269 * pkin(4) - t72;
t62 = t135 * pkin(5) + t70;
t438 = g(1) * t171 - g(2) * t169 + t248 * t421 + t62 * t91 + t35;
t42 = -t276 * t71 + t281 * t86;
t37 = -pkin(10) * t137 + t42;
t27 = pkin(5) * t432 + t37;
t408 = t280 * t38;
t10 = t27 * t275 + t408;
t168 = t249 * t390 + t387;
t170 = -t249 * t389 + t388;
t23 = t274 * t57 + t401 * t63;
t21 = pkin(9) * t267 + t23;
t246 = pkin(2) * t348;
t287 = pkin(3) * t290 + qJDD(1) * t330 + qJDD(4) + t246;
t36 = -t83 * pkin(4) - t84 * pkin(9) + t287;
t34 = t281 * t36;
t3 = pkin(5) * t80 - pkin(10) * t60 - qJD(5) * t43 - t276 * t21 + t34;
t309 = t281 * t21 + t276 * t36 + t86 * t367 - t368 * t71;
t4 = -pkin(10) * t61 + t309;
t354 = -t275 * t4 + t280 * t3;
t422 = g(3) * t261;
t437 = -g(1) * t170 + g(2) * t168 - qJD(6) * t10 + t248 * t422 + t62 * t316 + t354;
t297 = qJD(6) * t316 - t275 * t60 - t280 * t61;
t436 = -t143 * t316 + t297;
t339 = t401 * t277;
t404 = t401 * t138 - t139 * t274 + (t274 * t282 + t339) * t414;
t156 = -t274 * t210 + t212 * t401;
t113 = t211 * t156;
t374 = -t277 * t229 + t282 * t230;
t434 = t453 * pkin(5);
t433 = t276 * t99 - t281 * t403;
t78 = qJDD(6) + t80;
t430 = t143 * t444 + t211 * t78;
t424 = g(3) * t248;
t419 = t281 * pkin(5);
t265 = t281 * pkin(10);
t256 = pkin(2) * t282 + pkin(3);
t185 = pkin(2) * t339 + t274 * t256;
t180 = pkin(9) + t185;
t417 = -pkin(10) - t180;
t250 = pkin(3) * t274 + pkin(9);
t416 = -pkin(10) - t250;
t412 = t304 * t91;
t411 = t336 * t70;
t407 = t60 * t276;
t406 = t316 * t304;
t405 = t434 + t404;
t82 = t130 * t401 - t122;
t402 = t276 * t102 + t281 * t82;
t160 = t269 * t210;
t119 = -t160 * t401 - t274 * t161;
t399 = t119 * t281;
t398 = t135 * t304;
t397 = t137 * t304;
t396 = t143 * t304;
t395 = t432 * t304;
t393 = t156 * t276;
t392 = t156 * t281;
t391 = t194 * t192;
t384 = t276 * t119;
t383 = t276 * t279;
t382 = t276 * t284;
t380 = t279 * t281;
t333 = -t282 * t229 - t230 * t277;
t148 = -qJ(4) * t212 + t333;
t149 = -qJ(4) * t210 + t374;
t111 = t274 * t148 + t149 * t401;
t106 = t281 * t111;
t378 = t281 * t284;
t155 = t210 * t401 + t212 * t274;
t109 = pkin(4) * t155 - pkin(9) * t156 + t330;
t377 = t276 * t109 + t106;
t373 = pkin(3) * t264 + t266;
t270 = t278 ^ 2;
t372 = -t283 ^ 2 + t270;
t369 = qJD(3) * t282;
t359 = pkin(10) * t394;
t259 = t278 * t413;
t67 = t70 * t368;
t353 = qJD(2) * t429;
t352 = t156 * t368;
t351 = t156 * t367;
t350 = qJD(5) * t250 * t432;
t345 = pkin(3) * t161 + t259;
t344 = qJD(6) * t27 + t4;
t221 = t278 * t353;
t223 = t283 * t353;
t301 = -t282 * t221 - t277 * t223 - t229 * t369 - t230 * t370;
t95 = -qJ(4) * t161 - qJD(4) * t210 + t301;
t295 = -qJD(3) * t374 + t277 * t221 - t282 * t223;
t96 = qJ(4) * t160 - qJD(4) * t212 + t295;
t46 = t274 * t95 - t401 * t96;
t342 = qJD(5) * t417;
t341 = qJD(5) * t416;
t338 = -qJD(5) * t86 - t21;
t81 = t130 * t274 + t340;
t110 = -t401 * t148 + t149 * t274;
t329 = -t81 + t434;
t328 = pkin(5) * t304 - t265 * t336;
t251 = -pkin(3) * t401 - pkin(4);
t327 = -t367 * t71 + t34;
t326 = t143 * t445 - t209 * t78;
t324 = g(1) * t279 - g(2) * t284;
t164 = t417 * t276;
t323 = -qJD(6) * t164 - t276 * t342 - t359 + t433;
t165 = t180 * t281 + t265;
t322 = qJD(6) * t165 - t281 * t342 + t328 - t450;
t204 = t416 * t276;
t321 = -qJD(6) * t204 - t276 * t341 - t359 + t402;
t101 = t281 * t102;
t205 = t250 * t281 + t265;
t320 = qJD(6) * t205 - t276 * t82 - t281 * t341 + t101 + t328;
t319 = t304 * t73 + t336 * t72;
t318 = -t180 * t80 - t411;
t317 = -t250 * t80 - t411;
t184 = -pkin(2) * t386 + t256 * t401;
t314 = t276 * t448 + t304 * t43 + t70 * t367;
t313 = t281 * t455 - t304 * t42 + t67;
t312 = t281 * t80 - t432 * t453;
t179 = -pkin(4) - t184;
t310 = -0.2e1 * pkin(1) * t363 - pkin(7) * qJDD(2);
t308 = t351 + t384;
t307 = -t352 + t399;
t47 = t274 * t96 + t401 * t95;
t118 = -t160 * t274 + t161 * t401;
t53 = pkin(4) * t118 - pkin(9) * t119 + t345;
t305 = t109 * t367 - t111 * t368 + t276 * t53 + t281 * t47;
t285 = qJD(2) ^ 2;
t299 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t285 + t324;
t286 = qJD(1) ^ 2;
t298 = pkin(1) * t286 - pkin(7) * qJDD(1) + t325;
t13 = t61 * pkin(5) + t20;
t9 = t27 * t280 - t275 * t38;
t294 = t13 * t209 - t249 * t421 + t263 * t455 - t304 * t9 - t445 * t62;
t293 = g(3) * t262 - t277 * t162 - t228 * t192 + t264 * t325 + t282 * t331 + t190;
t292 = t10 * t304 + t13 * t211 + t249 * t422 - t261 * t455 + t444 * t62;
t289 = -t228 * t194 + t296 + t452;
t268 = -qJ(4) - t429;
t225 = t251 - t419;
t219 = pkin(1) + t373;
t189 = t246 - t451;
t178 = t249 * t378 + t383;
t177 = -t249 * t382 + t380;
t176 = -t249 * t380 + t382;
t175 = t249 * t383 + t378;
t172 = t179 - t419;
t140 = -t192 ^ 2 + t194 ^ 2;
t117 = -t194 * t269 + t288;
t116 = t192 * t269 + t132;
t114 = t209 * t156;
t105 = t281 * t109;
t69 = pkin(5) * t393 + t110;
t51 = t281 * t53;
t44 = -pkin(10) * t393 + t377;
t41 = pkin(5) * t155 - pkin(10) * t392 - t111 * t276 + t105;
t30 = t119 * t385 - t275 * t352 - t366 * t393 + (t392 * t431 + t384) * t280;
t29 = -t113 * t431 - t209 * t119;
t28 = t137 * t332 + t407;
t26 = pkin(5) * t308 + t46;
t25 = -t397 - t441;
t24 = t312 + t398;
t15 = t326 + t412;
t14 = t406 + t430;
t8 = (t60 - t447) * t281 + (-t137 * t432 - t61) * t276;
t7 = -pkin(10) * t308 + t305;
t6 = -pkin(10) * t399 + pkin(5) * t118 - t276 * t47 + t51 + (-t106 + (pkin(10) * t156 - t109) * t276) * qJD(5);
t5 = t16 * t211 - t316 * t444;
t1 = -t16 * t209 + t211 * t297 - t316 * t445 - t444 * t91;
t2 = [qJDD(1), t324, t325, qJDD(1) * t270 + 0.2e1 * t278 * t347, 0.2e1 * t278 * t360 - 0.2e1 * t363 * t372, qJDD(2) * t278 + t283 * t285, qJDD(2) * t283 - t278 * t285, 0, t278 * t310 + t283 * t299, -t278 * t299 + t283 * t310, t132 * t212 + t160 * t194, -t132 * t210 + t160 * t192 + t194 * t161 + t212 * t288, -t160 * t269 + t212 * t267, -t161 * t269 - t210 * t267, 0, t192 * t259 + t264 * t324 + t267 * t333 + t269 * t295 + (t189 - t451) * t210 - 0.2e1 * t228 * t161, -t132 * t418 + t228 * t160 + t189 * t212 - t194 * t259 - t262 * t324 - t267 * t374 - t269 * t301, t110 * t84 + t111 * t83 - t118 * t73 - t119 * t72 - t155 * t23 - t156 * t22 + t304 * t46 + t336 * t47 - t325, t23 * t111 + t73 * t47 - t22 * t110 - t72 * t46 + t287 * t330 + t163 * t345 - g(1) * (-t219 * t279 - t268 * t284) - g(2) * (t219 * t284 - t268 * t279) t137 * t307 + t392 * t60 (-t135 * t281 - t137 * t276) * t119 + (-t407 - t281 * t61 + (t135 * t276 - t137 * t281) * qJD(5)) * t156, t137 * t118 + t60 * t155 + t307 * t432 + t392 * t80, -t135 * t118 - t61 * t155 - t308 * t432 - t393 * t80, t118 * t432 + t155 * t80 (-t111 * t367 + t51) * t432 + t105 * t80 + t327 * t155 + t42 * t118 + t46 * t135 + t110 * t61 + t70 * t351 - g(1) * t176 - g(2) * t178 + ((-qJD(5) * t109 - t47) * t432 - t111 * t80 + t338 * t155 + t20 * t156 + t70 * t119) * t276, -t305 * t432 - t377 * t80 - t309 * t155 - t43 * t118 + t46 * t137 + t110 * t60 + t70 * t399 - g(1) * t175 - g(2) * t177 + (t20 * t281 - t67) * t156, -t114 * t16 - t29 * t316, -t113 * t16 - t114 * t297 - t29 * t91 + t30 * t316, -t114 * t78 - t118 * t316 + t143 * t29 + t155 * t16, -t113 * t78 - t118 * t91 - t143 * t30 + t155 * t297, t118 * t143 + t155 * t78 (-t275 * t7 + t280 * t6) * t143 + (-t275 * t44 + t280 * t41) * t78 + t354 * t155 + t9 * t118 + t26 * t91 - t69 * t297 + t13 * t113 + t62 * t30 - g(1) * t169 - g(2) * t171 + ((-t275 * t41 - t280 * t44) * t143 - t10 * t155) * qJD(6), -g(1) * t168 - g(2) * t170 - t10 * t118 - t13 * t114 + t35 * t155 + t69 * t16 - t26 * t316 + t62 * t29 + (-(-qJD(6) * t44 + t6) * t143 - t41 * t78 - t3 * t155) * t275 + (-(qJD(6) * t41 + t7) * t143 - t44 * t78 - t344 * t155) * t280; 0, 0, 0, -t278 * t286 * t283, t372 * t286, t361, t360, qJDD(2), -g(3) * t283 + t278 * t298, g(3) * t278 + t283 * t298, -t391, t140, t116, t117, t267, -t334 * t269 + (-t192 * t371 + t267 * t282 - t269 * t370) * pkin(2) + t289, t375 * t269 + (t194 * t371 - t267 * t277 - t269 * t369) * pkin(2) + t293, -t184 * t84 + t185 * t83 + t304 * t404 + t336 * t403 + t319, t23 * t185 + t22 * t184 - t163 * (t258 - t428) - g(3) * t373 + t403 * t73 - t404 * t72 - t325 * (-pkin(2) * t278 - pkin(3) * t262) t28, t8, t25, t24, -t395, t179 * t61 - t448 * t281 + t318 * t276 + t404 * t135 + (-t180 * t367 + t450) * t432 + t313, t179 * t60 + t318 * t281 - t276 * t455 + t404 * t137 + (t180 * t368 + t433) * t432 + t314, t5, t1, t14, t15, -t396 (t164 * t280 - t165 * t275) * t78 - t172 * t297 + t405 * t91 + (t275 * t323 - t280 * t322) * t143 + t294 -(t164 * t275 + t165 * t280) * t78 + t172 * t16 - t405 * t316 + (t275 * t322 + t280 * t323) * t143 + t292; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t391, t140, t116, t117, t267, -t269 * t315 + t289, t269 * t335 + t293, -t82 * t336 - t81 * t304 + (t274 * t83 - t401 * t84) * pkin(3) + t319, t72 * t81 - t73 * t82 + (t163 * t194 + t22 * t401 + t23 * t274 + t452) * pkin(3), t28, t8, t25, t24, -t395, -t101 * t432 - t81 * t135 + t251 * t61 + (t432 * t82 + t317) * t276 + (-t448 - t350) * t281 + t313, t251 * t60 + t402 * t432 - t81 * t137 + t317 * t281 + (-t455 + t350) * t276 + t314, t5, t1, t14, t15, -t396 (t204 * t280 - t205 * t275) * t78 - t225 * t297 + t329 * t91 + (t275 * t321 - t280 * t320) * t143 + t294 -(t204 * t275 + t205 * t280) * t78 + t225 * t16 - t329 * t316 + (t275 * t320 + t280 * t321) * t143 + t292; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t304 ^ 2 - t336 ^ 2, t304 * t72 - t336 * t73 + t287 - t324, 0, 0, 0, 0, 0, t312 - t398, -t397 + t441, 0, 0, 0, 0, 0, t326 - t412, t406 - t430; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t137 * t135, -t135 ^ 2 + t137 ^ 2, t60 + t447, -t337 + (-qJD(5) + t432) * t137, t80, -g(1) * t177 + g(2) * t175 - t137 * t70 + t432 * t43 + (t338 + t424) * t276 + t327, g(1) * t178 - g(2) * t176 + t135 * t70 + t281 * t424 + t42 * t432 - t309, -t449, t440, t439, t436, t78 -(-t275 * t37 - t408) * t143 + (-t137 * t91 - t143 * t366 + t280 * t78) * pkin(5) + t437 (-t143 * t38 - t3) * t275 + (t143 * t37 - t344) * t280 + (t137 * t316 - t143 * t365 - t275 * t78) * pkin(5) + t438; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t449, t440, t439, t436, t78, t10 * t143 + t437, t143 * t9 - t275 * t3 - t280 * t344 + t438;];
tau_reg  = t2;
