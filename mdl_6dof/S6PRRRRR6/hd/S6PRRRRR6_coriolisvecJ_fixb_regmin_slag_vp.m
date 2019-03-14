% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% tauc_reg [6x32]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRRRRR6_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR6_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR6_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRRRRR6_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:25:04
% EndTime: 2019-03-09 01:25:38
% DurationCPUTime: 13.27s
% Computational Cost: add. (14590->609), mult. (43605->910), div. (0->0), fcn. (37529->16), ass. (0->302)
t256 = cos(pkin(7));
t266 = cos(qJ(3));
t267 = cos(qJ(2));
t396 = t266 * t267;
t261 = sin(qJ(3));
t262 = sin(qJ(2));
t402 = t261 * t262;
t287 = -t256 * t402 + t396;
t380 = qJD(3) * t266;
t353 = t256 * t380;
t254 = sin(pkin(6));
t386 = qJD(1) * t254;
t457 = -pkin(2) * t353 + t287 * t386;
t253 = sin(pkin(7));
t255 = cos(pkin(8));
t434 = pkin(11) * t255;
t318 = t253 * (-pkin(10) - t434);
t298 = t261 * t318;
t456 = qJD(3) * t298 - t457;
t406 = t256 * t261;
t369 = pkin(2) * t406;
t180 = (t266 * t318 - t369) * qJD(3);
t399 = t262 * t266;
t400 = t261 * t267;
t289 = -t256 * t399 - t400;
t195 = t289 * t386;
t252 = sin(pkin(8));
t435 = pkin(11) * t252;
t297 = pkin(3) * t261 - t266 * t435;
t281 = t297 * qJD(3);
t206 = t253 * t281;
t358 = t262 * t386;
t332 = t253 * t358;
t393 = (t206 - t332) * t255 + (-t180 + t195) * t252;
t265 = cos(qJ(4));
t398 = t265 * t266;
t260 = sin(qJ(4));
t404 = t260 * t261;
t291 = -t255 * t404 + t398;
t384 = qJD(2) * t253;
t196 = t291 * t384;
t378 = qJD(4) * t265;
t351 = t252 * t378;
t455 = -t196 + t351;
t224 = pkin(10) * t384 + t358;
t367 = t253 * t434;
t347 = t267 * t386;
t429 = qJD(2) * pkin(2);
t234 = t347 + t429;
t257 = cos(pkin(6));
t385 = qJD(1) * t257;
t359 = t253 * t385;
t236 = t266 * t359;
t405 = t256 * t266;
t389 = t234 * t405 + t236;
t133 = (-qJD(2) * t367 - t224) * t261 + t389;
t407 = t255 * t266;
t356 = qJD(2) * t407;
t397 = t266 * t224;
t134 = -t234 * t406 - t397 + (-pkin(11) * t356 - t261 * t385) * t253;
t205 = t297 * t384;
t409 = t255 * t260;
t414 = t252 * t260;
t243 = pkin(11) * t414;
t408 = t255 * t265;
t443 = pkin(3) * t408 - t243;
t454 = t443 * qJD(4) - t265 * t133 - t134 * t409 - t205 * t414;
t412 = t253 * t266;
t169 = t369 + pkin(10) * t412 + (t252 * t256 + t253 * t407) * pkin(11);
t178 = (pkin(2) * t266 + pkin(3)) * t256 + t298;
t296 = -pkin(3) * t266 - t261 * t435;
t202 = (-pkin(2) + t296) * t253;
t277 = t195 * t255 + t252 * t332;
t349 = t255 * t378;
t379 = qJD(4) * t260;
t453 = -t169 * t379 + t178 * t349 + t180 * t409 + t202 * t351 + t206 * t414 - t260 * t277 + t265 * t456;
t381 = qJD(3) * t253;
t355 = t261 * t381;
t327 = t252 * t355;
t452 = -pkin(12) * t327 - t453;
t441 = t255 * t398 - t404;
t131 = t256 * t351 + (t291 * qJD(3) + qJD(4) * t441) * t253;
t401 = t261 * t265;
t403 = t260 * t266;
t292 = t255 * t403 + t401;
t293 = t255 * t401 + t403;
t436 = qJD(3) * t293 + qJD(4) * t292;
t269 = t436 * t253;
t352 = t252 * t379;
t132 = t256 * t352 + t269;
t451 = -pkin(4) * t132 + pkin(12) * t131 - t393;
t194 = t293 * t384;
t94 = -t134 * t252 + t255 * t205;
t450 = pkin(4) * t194 - pkin(12) * t196 - (pkin(4) * t260 - pkin(12) * t265) * t252 * qJD(4) + t94;
t357 = t261 * t384;
t330 = t252 * t357;
t449 = -pkin(12) * t330 + t454;
t350 = t255 * t379;
t448 = t169 * t378 + t178 * t350 + t202 * t352 + t260 * t456 + t265 * t277;
t259 = sin(qJ(5));
t264 = cos(qJ(5));
t218 = -t264 * t255 + t259 * t414;
t392 = -qJD(5) * t218 - t259 * t330 + t264 * t455;
t219 = t255 * t259 + t264 * t414;
t391 = qJD(5) * t219 + t259 * t455 + t264 * t330;
t447 = t194 - t352;
t413 = t252 * t265;
t388 = pkin(3) * t409 + pkin(11) * t413;
t446 = t388 * qJD(4) - t260 * t133;
t326 = t253 * t356;
t228 = t265 * t326;
t383 = qJD(2) * t256;
t337 = qJD(3) + t383;
t299 = t337 * t252;
t328 = t260 * t357;
t157 = -t265 * t299 - t228 + t328;
t156 = qJD(5) + t157;
t279 = t292 * t253;
t159 = qJD(2) * t279 + t260 * t299;
t366 = t252 * t412;
t235 = qJD(2) * t366;
t371 = t235 - qJD(4);
t278 = -t255 * t337 + t371;
t190 = t264 * t278;
t118 = t159 * t259 + t190;
t116 = qJD(6) + t118;
t423 = t134 * t408 - (-pkin(4) * t357 - t205 * t265) * t252 + t446;
t374 = qJD(5) * t264;
t376 = qJD(5) * t259;
t128 = -t178 * t252 + t255 * t202;
t172 = -t253 * t441 - t256 * t413;
t175 = t256 * t414 + t279;
t72 = pkin(4) * t172 - pkin(12) * t175 + t128;
t410 = t255 * t256;
t216 = t366 - t410;
t361 = t265 * t169 + t178 * t409 + t202 * t414;
t77 = -pkin(12) * t216 + t361;
t445 = t451 * t259 + t264 * t452 - t72 * t374 + t376 * t77;
t421 = -t180 * t408 + (-pkin(4) * t355 - t206 * t265) * t252 + t448;
t444 = t261 * t266;
t209 = pkin(12) * t255 + t388;
t210 = (-pkin(4) * t265 - pkin(12) * t260 - pkin(3)) * t252;
t390 = t264 * t209 + t259 * t210;
t442 = t209 * t376 - t210 * t374 + t259 * t450 - t264 * t449;
t120 = t264 * t159 - t259 * t278;
t283 = qJD(4) * t299;
t280 = t260 * t283;
t122 = qJD(2) * t269 + t280;
t285 = t234 * t256 + t359;
t117 = t397 + t285 * t261 + (t299 + t326) * pkin(11);
t123 = pkin(3) * t337 + t133;
t242 = t256 * t385;
t146 = t242 + (qJD(2) * t296 - t234) * t253;
t166 = (t281 + t358) * t384;
t331 = t256 * t358;
t325 = qJD(2) * t347;
t360 = qJD(3) * t236 + t234 * t353 + t266 * t325;
t382 = qJD(3) * t224;
t92 = (-t382 + (-qJD(3) * t367 - t331) * qJD(2)) * t261 + t360;
t93 = qJD(2) * t195 + qJD(3) * t134;
t275 = t117 * t379 - t123 * t349 - t146 * t351 - t166 * t414 - t265 * t92 - t93 * t409;
t370 = qJD(2) * qJD(3);
t346 = t253 * t370;
t324 = t261 * t346;
t300 = t252 * t324;
t19 = pkin(12) * t300 - t275;
t53 = t265 * t117 + t123 * t409 + t146 * t414;
t45 = -pkin(12) * t278 + t53;
t78 = -t123 * t252 + t255 * t146;
t51 = pkin(4) * t157 - pkin(12) * t159 + t78;
t24 = t259 * t51 + t264 * t45;
t323 = t266 * t346;
t335 = -t255 * qJD(3) - qJD(4);
t121 = qJD(4) * t228 + t328 * t335 + (t283 + t323) * t265;
t67 = t255 * t166 - t252 * t93;
t39 = pkin(4) * t122 - pkin(12) * t121 + t67;
t6 = -qJD(5) * t24 - t19 * t259 + t264 * t39;
t4 = -pkin(5) * t122 - t6;
t440 = t116 * (pkin(5) * t120 + pkin(13) * t116) + t4;
t258 = sin(qJ(6));
t263 = cos(qJ(6));
t60 = -qJD(5) * t190 + t264 * t121 - t159 * t376 + t259 * t300;
t81 = t120 * t263 + t156 * t258;
t27 = qJD(6) * t81 - t263 * t122 + t258 * t60;
t52 = -t260 * t117 + t265 * (t123 * t255 + t146 * t252);
t439 = -t260 * t169 + t265 * (t178 * t255 + t202 * t252);
t288 = t256 * t400 + t399;
t174 = t253 * t257 * t261 + t254 * t288;
t290 = t256 * t396 - t402;
t173 = t254 * t290 + t257 * t412;
t217 = -t253 * t254 * t267 + t256 * t257;
t304 = t173 * t255 + t217 * t252;
t438 = -t174 * t260 + t265 * t304;
t61 = qJD(5) * t120 + t259 * t121 - t264 * t300;
t22 = -t117 * t378 - t123 * t350 - t146 * t352 + t166 * t413 - t260 * t92 + t93 * t408;
t20 = -pkin(4) * t300 - t22;
t10 = pkin(5) * t61 - pkin(13) * t60 + t20;
t5 = t264 * t19 + t259 * t39 + t51 * t374 - t376 * t45;
t3 = pkin(13) * t122 + t5;
t17 = pkin(13) * t156 + t24;
t44 = pkin(4) * t278 - t52;
t25 = t118 * pkin(5) - t120 * pkin(13) + t44;
t312 = t17 * t258 - t25 * t263;
t1 = -qJD(6) * t312 + t10 * t258 + t263 * t3;
t307 = t259 * t72 + t264 * t77;
t437 = -qJD(5) * t307 + t259 * t452 - t451 * t264;
t268 = qJD(2) ^ 2;
t433 = -pkin(5) * t132 - t437;
t98 = pkin(4) * t159 + pkin(12) * t157;
t432 = t259 * t98 + t264 * t52;
t430 = t447 * pkin(5) + t390 * qJD(5) + t449 * t259 + t450 * t264;
t79 = t120 * t258 - t263 * t156;
t428 = t116 * t79;
t427 = t116 * t81;
t372 = qJD(6) * t263;
t373 = qJD(6) * t258;
t26 = -t120 * t373 + t258 * t122 + t156 * t372 + t263 * t60;
t426 = t258 * t26;
t425 = t258 * t61;
t424 = t263 * t61;
t422 = -t53 + t156 * (pkin(5) * t259 - pkin(13) * t264);
t420 = t118 * t156;
t419 = t120 * t156;
t418 = t157 * t264;
t199 = -t234 * t253 + t242;
t416 = t199 * t253;
t249 = t253 ^ 2;
t415 = t249 * t268;
t411 = t254 * t268;
t183 = t219 * t258 + t263 * t413;
t395 = -qJD(6) * t183 - t258 * t447 + t263 * t392;
t365 = t258 * t413;
t394 = -qJD(6) * t365 + t219 * t372 + t258 * t392 + t263 * t447;
t387 = t261 ^ 2 - t266 ^ 2;
t377 = qJD(5) * t258;
t375 = qJD(5) * t263;
t364 = t262 * t411;
t354 = t253 * t380;
t348 = t253 * t256 * t268;
t239 = -pkin(5) * t264 - pkin(13) * t259 - pkin(4);
t342 = pkin(13) * t159 - qJD(6) * t239 + t432;
t340 = t116 * t263;
t339 = t264 * t156;
t338 = 0.2e1 * t249 * t370;
t336 = qJD(3) + 0.2e1 * t383;
t334 = t249 * t364;
t329 = t254 * t262 * t384;
t137 = t175 * t259 + t264 * t216;
t138 = t175 * t264 - t216 * t259;
t76 = pkin(4) * t216 - t439;
t43 = pkin(5) * t137 - pkin(13) * t138 + t76;
t322 = -pkin(13) * t132 - qJD(6) * t43 + t445;
t35 = pkin(13) * t172 + t307;
t68 = -qJD(5) * t137 + t131 * t264 + t259 * t327;
t69 = qJD(5) * t138 + t131 * t259 - t264 * t327;
t321 = -pkin(5) * t69 + pkin(13) * t68 + qJD(6) * t35 - t421;
t208 = t243 + (-pkin(3) * t265 - pkin(4)) * t255;
t139 = pkin(5) * t218 - pkin(13) * t219 + t208;
t320 = pkin(13) * t447 - qJD(6) * t139 + t442;
t141 = -pkin(13) * t413 + t390;
t319 = -pkin(5) * t391 + pkin(13) * t392 + qJD(6) * t141 - t423;
t96 = t159 * t258 - t263 * t418;
t316 = t263 * t374 - t96;
t315 = -t224 + t358;
t8 = t17 * t263 + t25 * t258;
t136 = -t173 * t252 + t217 * t255;
t87 = t174 * t265 + t260 * t304;
t63 = t136 * t259 + t264 * t87;
t311 = -t258 * t438 + t263 * t63;
t310 = -t258 * t63 - t263 * t438;
t23 = -t259 * t45 + t264 * t51;
t308 = -t259 * t77 + t264 * t72;
t306 = t136 * t264 - t259 * t87;
t91 = t138 * t263 + t172 * t258;
t90 = t138 * t258 - t263 * t172;
t302 = -t209 * t259 + t210 * t264;
t301 = t252 ^ 2 * t324;
t295 = -t116 * t372 - t425;
t294 = -t116 * t373 + t424;
t284 = -pkin(12) * t122 + t156 * t44;
t129 = -t257 * t355 + (qJD(2) * t289 - qJD(3) * t288) * t254;
t276 = t129 * t255 + t252 * t329;
t274 = -qJD(2) * t331 - t382;
t16 = -pkin(5) * t156 - t23;
t273 = -pkin(13) * t61 + (t16 + t23) * t116;
t271 = qJD(4) * t278;
t2 = -qJD(6) * t8 + t263 * t10 - t258 * t3;
t270 = -t256 * t285 + t416;
t184 = t219 * t263 - t365;
t140 = pkin(5) * t413 - t302;
t130 = t257 * t354 + (qJD(2) * t287 + qJD(3) * t290) * t254;
t100 = -t129 * t252 + t255 * t329;
t95 = -t263 * t159 - t258 * t418;
t42 = qJD(4) * t438 + t130 * t265 + t276 * t260;
t41 = qJD(4) * t87 + t130 * t260 - t265 * t276;
t34 = -pkin(5) * t172 - t308;
t33 = qJD(6) * t91 - t263 * t132 + t258 * t68;
t32 = -qJD(6) * t90 + t132 * t258 + t263 * t68;
t30 = -pkin(5) * t159 + t259 * t52 - t264 * t98;
t14 = qJD(5) * t306 + t100 * t259 + t264 * t42;
t13 = qJD(5) * t63 - t100 * t264 + t259 * t42;
t7 = [0, 0, -t364, -t267 * t411, 0, 0, 0, 0, 0, t129 * t337 + t217 * t324 - t266 * t334, -t130 * t337 + t217 * t323 + t261 * t334, 0, 0, 0, 0, 0, t100 * t157 + t136 * t122 + t278 * t41 + t300 * t438, t100 * t159 + t136 * t121 + t278 * t42 - t300 * t87, 0, 0, 0, 0, 0, t118 * t41 + t122 * t306 - t13 * t156 - t438 * t61, t120 * t41 - t122 * t63 - t14 * t156 - t438 * t60, 0, 0, 0, 0, 0 (-qJD(6) * t311 - t14 * t258 + t263 * t41) * t116 + t310 * t61 + t13 * t79 - t306 * t27 -(qJD(6) * t310 + t14 * t263 + t258 * t41) * t116 - t311 * t61 + t13 * t81 - t306 * t26; 0, 0, 0, 0, t338 * t444, -t387 * t338, t336 * t354, -t336 * t355, 0, -t195 * t337 + (-pkin(10) * t337 * t381 + t256 * t274) * t266 + (-t256 * t325 + t270 * qJD(3) + (-qJD(3) * t256 + (-t256 ^ 2 - t249) * qJD(2)) * qJD(3) * pkin(2)) * t261 -(t261 * t274 + t360) * t256 + (-t249 * t429 + t416) * t380 + (pkin(10) * t355 + t457) * t337, t121 * t175 + t131 * t159, -t121 * t172 - t122 * t175 - t131 * t157 - t132 * t159, -t216 * t121 - t278 * t131 + (qJD(2) * t175 + t159) * t327, t216 * t122 + t278 * t132 + (-qJD(2) * t172 - t157) * t327 (-t235 + (-t216 + t410) * qJD(2) - t335) * t327, t128 * t122 + t67 * t172 + t78 * t132 - t22 * t216 + (qJD(2) * t439 + t52) * t327 + t393 * t157 + (-(t180 * t255 + t206 * t252) * t265 + t448) * t278, t128 * t121 + t67 * t175 + t78 * t131 - t275 * t216 + (-qJD(2) * t361 - t53) * t327 + t393 * t159 + t453 * t278, t120 * t68 + t138 * t60, -t118 * t68 - t120 * t69 - t137 * t60 - t138 * t61, t120 * t132 + t122 * t138 + t156 * t68 + t172 * t60, -t118 * t132 - t122 * t137 - t156 * t69 - t172 * t61, t122 * t172 + t132 * t156, t421 * t118 + t308 * t122 + t23 * t132 + t20 * t137 + t156 * t437 + t6 * t172 + t44 * t69 + t76 * t61, t421 * t120 - t307 * t122 - t24 * t132 + t20 * t138 + t156 * t445 - t5 * t172 + t44 * t68 + t76 * t60, t26 * t91 + t32 * t81, -t26 * t90 - t27 * t91 - t32 * t79 - t33 * t81, t116 * t32 + t137 * t26 + t61 * t91 + t69 * t81, -t116 * t33 - t137 * t27 - t61 * t90 - t69 * t79, t116 * t69 + t137 * t61 (-t258 * t35 + t263 * t43) * t61 + t2 * t137 - t312 * t69 + t34 * t27 + t4 * t90 + t16 * t33 + t433 * t79 + (t258 * t322 - t263 * t321) * t116 -(t258 * t43 + t263 * t35) * t61 - t1 * t137 - t8 * t69 + t34 * t26 + t4 * t91 + t16 * t32 + t433 * t81 + (t258 * t321 + t263 * t322) * t116; 0, 0, 0, 0, -t415 * t444, t387 * t415, -t266 * t348, t261 * t348, 0 (-t315 * t405 + (-t270 - t347) * t261) * qJD(2) (-t199 * t412 + (t261 * t315 + t389) * t256) * qJD(2) + t389 * qJD(3) - t360, t121 * t414 + t159 * t455, t157 * t196 + t159 * t194 + (t121 * t265 - t122 * t260 + (-t157 * t265 - t159 * t260) * qJD(4)) * t252, t255 * t121 + t260 * t301 + t278 * t196 + (-t159 * t357 - t265 * t271) * t252, -t255 * t122 + t265 * t301 - t278 * t194 + (t157 * t357 + t260 * t271) * t252 -(t255 * t383 - t371) * t330, -t252 * pkin(3) * t122 - t67 * t413 + t22 * t255 - t94 * t157 - t447 * t78 + (qJD(3) * t443 - t52) * t330 + ((t134 * t255 + t205 * t252) * t265 + t446) * t278, t275 * t255 - t94 * t159 - t78 * t196 + (t78 * t378 - pkin(3) * t121 + t67 * t260 + (-qJD(3) * t388 + t53) * t357) * t252 + t454 * t278, t120 * t392 + t219 * t60, -t118 * t392 - t120 * t391 - t218 * t60 - t219 * t61, -t120 * t447 + t122 * t219 + t156 * t392 - t413 * t60, t118 * t447 - t122 * t218 - t156 * t391 + t413 * t61, -t122 * t413 - t156 * t447, t302 * t122 - t6 * t413 + t208 * t61 + t20 * t218 + t391 * t44 - t447 * t23 + ((-qJD(5) * t209 - t450) * t264 + (-qJD(5) * t210 - t449) * t259) * t156 + t423 * t118, t423 * t120 - t390 * t122 + t156 * t442 + t20 * t219 + t208 * t60 + t24 * t447 + t392 * t44 + t5 * t413, t184 * t26 + t395 * t81, -t183 * t26 - t184 * t27 - t394 * t81 - t395 * t79, t116 * t395 + t184 * t61 + t218 * t26 + t391 * t81, -t116 * t394 - t183 * t61 - t218 * t27 - t391 * t79, t116 * t391 + t218 * t61 (t139 * t263 - t141 * t258) * t61 + t2 * t218 + t140 * t27 + t4 * t183 + t430 * t79 - t391 * t312 + t394 * t16 + (t258 * t320 - t263 * t319) * t116 -(t139 * t258 + t141 * t263) * t61 - t1 * t218 + t140 * t26 + t4 * t184 + t430 * t81 - t391 * t8 + t395 * t16 + (t258 * t319 + t263 * t320) * t116; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t157 * t159, -t157 ^ 2 + t159 ^ 2, -t157 * t278 + t121, -t278 * t159 - t384 * t436 - t280, t300, -t78 * t159 - t278 * t53 + t22, t78 * t157 - t278 * t52 + t275, t120 * t339 + t259 * t60 (t60 - t420) * t264 + (-t61 - t419) * t259, -t120 * t159 + t259 * t122 + t156 * t339, -t156 ^ 2 * t259 + t118 * t159 + t264 * t122, -t156 * t159, -pkin(4) * t61 - t53 * t118 - t23 * t159 + (-t20 + (-pkin(12) * qJD(5) - t98) * t156) * t264 + (t52 * t156 + t284) * t259, -pkin(4) * t60 - t53 * t120 + t24 * t159 + t20 * t259 + (pkin(12) * t376 + t432) * t156 + t284 * t264, t259 * t26 * t263 + (-t259 * t373 + t316) * t81, t79 * t96 + t81 * t95 + (-t258 * t81 - t263 * t79) * t374 + (-t426 - t263 * t27 + (t258 * t79 - t263 * t81) * qJD(6)) * t259, -t26 * t264 + t316 * t116 + (t156 * t81 + t294) * t259, t264 * t27 + (-t258 * t374 + t95) * t116 + (-t156 * t79 + t295) * t259, t116 * t156 * t259 - t264 * t61, t239 * t424 - t16 * t95 - t30 * t79 + (t258 * t342 + t263 * t422) * t116 + (t16 * t377 - t2 + (qJD(5) * t79 + t295) * pkin(12)) * t264 + (t16 * t372 + t4 * t258 - t156 * t312 + (t116 * t377 + t27) * pkin(12)) * t259, -t239 * t425 - t16 * t96 - t30 * t81 + (-t258 * t422 + t263 * t342) * t116 + (t16 * t375 + t1 + (qJD(5) * t81 - t294) * pkin(12)) * t264 + (-t16 * t373 + t4 * t263 - t156 * t8 + (t116 * t375 + t26) * pkin(12)) * t259; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t120 * t118, -t118 ^ 2 + t120 ^ 2, t60 + t420, t419 - t61, t122, -t120 * t44 + t156 * t24 + t6, t118 * t44 + t156 * t23 - t5, t340 * t81 + t426 (t26 - t428) * t263 + (-t27 - t427) * t258, t116 * t340 - t120 * t81 + t425, -t116 ^ 2 * t258 + t120 * t79 + t424, -t116 * t120, -pkin(5) * t27 + t120 * t312 - t24 * t79 + t273 * t258 - t263 * t440, -pkin(5) * t26 + t8 * t120 - t24 * t81 + t258 * t440 + t273 * t263; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81 * t79, -t79 ^ 2 + t81 ^ 2, t26 + t428, -t27 + t427, t61, t116 * t8 - t16 * t81 + t2, -t116 * t312 + t16 * t79 - t1;];
tauc_reg  = t7;