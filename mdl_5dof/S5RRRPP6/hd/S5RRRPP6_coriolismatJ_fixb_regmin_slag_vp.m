% Calculate minimal parameter regressor of coriolis matrix for
% S5RRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x23]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RRRPP6_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP6_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP6_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP6_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:02:02
% EndTime: 2019-12-31 21:02:14
% DurationCPUTime: 5.16s
% Computational Cost: add. (6232->462), mult. (13034->637), div. (0->0), fcn. (12867->6), ass. (0->327)
t277 = sin(pkin(8));
t278 = cos(pkin(8));
t280 = sin(qJ(2));
t281 = cos(qJ(3));
t404 = t281 * t280;
t279 = sin(qJ(3));
t407 = t279 * t280;
t198 = -t277 * t407 + t278 * t404;
t461 = t198 / 0.2e1;
t408 = t278 * t281;
t216 = t277 * t279 - t408;
t460 = -t216 / 0.2e1;
t459 = t216 / 0.2e1;
t409 = t278 * t279;
t412 = t277 * t281;
t218 = t409 + t412;
t458 = -t218 / 0.2e1;
t482 = -qJ(4) - pkin(7);
t237 = t482 * t281;
t222 = t278 * t237;
t345 = t482 * t279;
t340 = t277 * t345;
t479 = -t222 + t340;
t489 = -t479 / 0.2e1;
t212 = t218 ^ 2;
t480 = t216 ^ 2 + t212;
t488 = qJD(4) * t480;
t194 = t198 ^ 2;
t195 = t218 * t280;
t481 = t195 ^ 2 + t194;
t487 = qJD(4) * t481;
t326 = t216 * t195 + t218 * t198;
t486 = qJD(1) * t481 + qJD(2) * t326;
t485 = qJD(1) * t326 + qJD(2) * t480;
t282 = cos(qJ(2));
t339 = -t282 * pkin(2) - t280 * pkin(7);
t236 = -pkin(1) + t339;
t406 = t279 * t282;
t253 = pkin(6) * t406;
t395 = -t281 * t236 + t253;
t165 = -qJ(4) * t404 - t395;
t148 = -t282 * pkin(3) + t165;
t137 = t277 * t148;
t403 = t281 * t282;
t366 = pkin(6) * t403;
t166 = t366 + (-qJ(4) * t280 + t236) * t279;
t410 = t278 * t166;
t84 = t410 + t137;
t484 = t84 / 0.2e1;
t420 = t216 * qJ(5);
t443 = t279 * pkin(3);
t444 = t218 * pkin(4);
t140 = t420 + t443 + t444;
t483 = t140 / 0.2e1;
t240 = t277 * t406;
t200 = t278 * t403 - t240;
t422 = t200 * qJ(5);
t197 = t218 * t282;
t446 = t197 * pkin(4);
t478 = -t422 / 0.2e1 + t446 / 0.2e1;
t171 = -t277 * t237 - t278 * t345;
t327 = t171 * t218 - t216 * t479;
t474 = qJD(2) * t327;
t472 = qJD(4) * t327;
t471 = t326 * qJD(4);
t470 = t171 * t461 + t195 * t489;
t414 = t277 * t166;
t83 = t278 * t148 - t414;
t469 = t83 * t458 + t84 * t460;
t402 = t282 * qJ(5);
t71 = t84 - t402;
t72 = t282 * pkin(4) - t83;
t468 = t72 * t458 + t71 * t459;
t273 = t279 ^ 2;
t275 = t281 ^ 2;
t247 = t275 - t273;
t374 = t280 * qJD(1);
t360 = t281 * t374;
t467 = t247 * qJD(2) - 0.2e1 * t279 * t360;
t415 = t277 * t165;
t95 = t410 + t415;
t466 = -t95 / 0.2e1;
t96 = t278 * t165 - t414;
t465 = t96 / 0.2e1;
t464 = t479 / 0.2e1;
t463 = -t171 / 0.2e1;
t457 = t218 / 0.2e1;
t350 = -t222 / 0.2e1;
t257 = t277 * pkin(3) + qJ(5);
t456 = -t257 / 0.2e1;
t455 = t257 / 0.2e1;
t259 = -t278 * pkin(3) - pkin(4);
t454 = t259 / 0.2e1;
t453 = -t277 / 0.2e1;
t452 = t277 / 0.2e1;
t451 = -t278 / 0.2e1;
t450 = -t279 / 0.2e1;
t449 = -t280 / 0.2e1;
t448 = -t282 / 0.2e1;
t447 = t282 / 0.2e1;
t445 = t198 * pkin(4);
t442 = t280 * pkin(2);
t441 = t280 * pkin(4);
t440 = t282 * pkin(7);
t439 = pkin(3) * qJD(3);
t228 = pkin(3) * t407 + t280 * pkin(6);
t108 = t195 * pkin(4) - t198 * qJ(5) + t228;
t254 = pkin(3) * t406;
t271 = t282 * pkin(6);
t229 = t271 + t254;
t109 = t229 - t422 + t446;
t238 = -t440 + t442;
t227 = t281 * t238;
t256 = pkin(6) * t407;
t158 = t280 * pkin(3) - qJ(4) * t403 + t227 + t256;
t141 = t277 * t158;
t226 = t279 * t238;
t368 = pkin(6) * t404;
t167 = -qJ(4) * t406 + t226 - t368;
t159 = t278 * t167;
t92 = t159 + t141;
t79 = t280 * qJ(5) + t92;
t411 = t278 * t158;
t413 = t277 * t167;
t91 = t411 - t413;
t80 = -t91 - t441;
t9 = t108 * t109 + t71 * t79 + t72 * t80;
t434 = t9 * qJD(1);
t433 = t95 * t171;
t432 = t96 * t479;
t367 = pkin(3) * t404;
t423 = t195 * qJ(5);
t111 = t367 + t423 + t445;
t10 = t108 * t111 + t71 * t96 + t72 * t95;
t431 = t10 * qJD(1);
t430 = t108 * t198;
t11 = -t79 * t195 - t71 * t197 + t80 * t198 + t72 * t200;
t429 = t11 * qJD(1);
t16 = (-t71 + t95) * t198 + (-t72 - t96) * t195;
t428 = t16 * qJD(1);
t17 = -t92 * t195 - t84 * t197 - t91 * t198 - t83 * t200;
t427 = t17 * qJD(1);
t18 = (-t84 + t95) * t198 + (t83 - t96) * t195;
t424 = t18 * qJD(1);
t21 = t228 * t229 + t83 * t91 + t84 * t92;
t421 = t21 * qJD(1);
t22 = t108 * t200 + t109 * t198 - t71 * t280 + t79 * t282;
t419 = t22 * qJD(1);
t23 = t108 * t197 + t109 * t195 - t72 * t280 + t80 * t282;
t418 = t23 * qJD(1);
t24 = t228 * t367 - t83 * t95 + t84 * t96;
t417 = t24 * qJD(1);
t264 = -t281 * pkin(3) - pkin(2);
t136 = t216 * pkin(4) - t218 * qJ(5) + t264;
t313 = t108 * t457 + t136 * t461;
t289 = t448 * t479 - t313;
t310 = t413 / 0.2e1 - t411 / 0.2e1;
t291 = -t441 / 0.2e1 + t310;
t27 = -t289 + t291;
t416 = t27 * qJD(1);
t274 = t280 ^ 2;
t405 = t281 * t274;
t29 = t111 * t195 + t95 * t282 + t430;
t401 = t29 * qJD(1);
t30 = -t108 * t195 + t111 * t198 + t96 * t282;
t400 = t30 * qJD(1);
t31 = -t71 * t195 + t72 * t198;
t399 = t31 * qJD(1);
t33 = -t84 * t195 - t83 * t198;
t398 = t33 * qJD(1);
t38 = t71 * t282 + t430;
t397 = t38 * qJD(1);
t190 = t279 * t236 + t366;
t99 = t226 * t282 + (-t190 + t366) * t280;
t396 = t99 * qJD(1);
t276 = t282 ^ 2;
t248 = t276 - t274;
t394 = qJD(2) * t279;
t393 = qJD(2) * t281;
t392 = qJD(3) * t279;
t391 = qJD(3) * t281;
t390 = qJD(5) * t195;
t100 = t395 * t280 + (-t256 + t227) * t282;
t389 = t100 * qJD(1);
t307 = t412 / 0.2e1 + t409 / 0.2e1;
t127 = (t458 + t307) * t282;
t388 = t127 * qJD(1);
t128 = -t240 / 0.2e1 + (t408 / 0.2e1 + t460) * t282;
t387 = t128 * qJD(1);
t129 = t240 / 0.2e1 + (-t408 / 0.2e1 + t460) * t282;
t386 = t129 * qJD(1);
t151 = t274 * pkin(6) * t279 + t395 * t282;
t385 = t151 * qJD(1);
t152 = -pkin(6) * t405 - t190 * t282;
t384 = t152 * qJD(1);
t383 = t195 * qJD(4);
t382 = t198 * qJD(1);
t381 = t198 * qJD(4);
t380 = t216 * qJD(3);
t379 = t218 * qJD(2);
t378 = t218 * qJD(5);
t224 = t248 * t279;
t377 = t224 * qJD(1);
t225 = t281 * t276 - t405;
t376 = t225 * qJD(1);
t375 = t248 * qJD(1);
t373 = t280 * qJD(2);
t372 = t282 * qJD(1);
t371 = t282 * qJD(2);
t370 = t282 * qJD(3);
t369 = t282 * qJD(5);
t365 = pkin(1) * t374;
t364 = pkin(1) * t372;
t265 = t443 / 0.2e1;
t362 = t454 - pkin(4) / 0.2e1;
t361 = -t432 / 0.2e1;
t359 = t279 * t370;
t358 = t281 * t370;
t357 = t195 * t382;
t356 = t279 * t391;
t355 = t279 * t393;
t354 = t280 * t371;
t353 = t280 * t372;
t352 = t281 * t373;
t349 = -t404 / 0.2e1;
t348 = qJ(5) / 0.2e1 + t455;
t347 = -t137 / 0.2e1 - t410 / 0.2e1;
t346 = t254 / 0.2e1 + t271 / 0.2e1;
t344 = t171 * t200 - t197 * t479;
t183 = t195 * t372;
t343 = t128 * qJD(2) - t183;
t243 = t367 / 0.2e1;
t251 = -qJD(3) + t372;
t341 = t279 * t352;
t338 = t410 / 0.2e1 + t415 / 0.2e1;
t285 = t108 * t483 + t111 * t136 / 0.2e1 + t71 * t463 + t72 * t464 + t433 / 0.2e1;
t317 = t80 * t454 + t79 * t455;
t1 = t361 - t285 + t317;
t32 = t136 * t140;
t337 = -t1 * qJD(1) + t32 * qJD(2);
t287 = (t463 + t171 / 0.2e1) * t195;
t311 = t197 * t456 + t200 * t454;
t3 = (t71 / 0.2e1 + t466) * t218 + (t465 + t72 / 0.2e1) * t216 + t287 + t311;
t336 = t3 * qJD(1);
t296 = (t197 * t453 + t200 * t451) * pkin(3);
t5 = (t484 + t466) * t218 + (t465 - t83 / 0.2e1) * t216 + t296 + t287;
t335 = t5 * qJD(1);
t39 = t264 * t443;
t292 = t83 * t464 + t171 * t484 - t433 / 0.2e1;
t316 = t92 * t452 + t91 * t278 / 0.2e1;
t7 = t361 + (t228 * t450 + t264 * t349 + t316) * pkin(3) + t292;
t334 = -t7 * qJD(1) + t39 * qJD(2);
t284 = t108 * t460 + t111 * t457 - t136 * t195 / 0.2e1 + t140 * t461 - t171 * t447;
t293 = t348 * t280 + t141 / 0.2e1 + t159 / 0.2e1;
t12 = t284 + t293;
t45 = t136 * t216 - t140 * t218;
t333 = -t12 * qJD(1) + t45 * qJD(2);
t283 = t111 * t459 + t195 * t483 + t447 * t479 + t313;
t14 = t362 * t280 + t283 + t310;
t44 = t136 * t218 + t140 * t216;
t332 = t14 * qJD(1) + t44 * qJD(2);
t288 = t346 - t470;
t19 = t288 + t468 + t478;
t331 = -t19 * qJD(1) + t474;
t25 = t288 - t469;
t330 = -t25 * qJD(1) + t474;
t40 = t348 * t195 - t362 * t198 + t243;
t67 = t348 * t216 - t362 * t218 + t265;
t329 = t40 * qJD(1) + t67 * qJD(2);
t321 = t251 * t280;
t309 = -t195 * t452 + t198 * t451;
t105 = (t349 + t309) * pkin(3);
t308 = t216 * t453 + t218 * t451;
t123 = (t450 + t308) * pkin(3);
t320 = t105 * qJD(1) + t123 * qJD(2);
t319 = -t195 * qJD(1) - t216 * qJD(2);
t132 = t379 + t382;
t318 = t442 / 0.2e1 - t440 / 0.2e1;
t301 = t318 * t279;
t160 = t226 / 0.2e1 + t301;
t315 = pkin(2) * t393 - t160 * qJD(1);
t302 = t318 * t281;
t161 = -t227 / 0.2e1 - t302;
t314 = pkin(2) * t394 - t161 * qJD(1);
t312 = -t195 * t458 + t198 * t459;
t81 = t449 + t312;
t306 = t81 * qJD(1) + t216 * t379;
t305 = t281 * t321;
t209 = (t273 / 0.2e1 - t275 / 0.2e1) * t280;
t304 = -t209 * qJD(1) + t355;
t303 = t129 * qJD(2) + t195 * qJD(3) - t183;
t300 = t346 + t470;
t299 = t279 * qJD(1) * t405 + t209 * qJD(2);
t223 = t247 * t274;
t298 = t223 * qJD(1) + 0.2e1 * t341;
t164 = t350 + t222 / 0.2e1;
t36 = t348 * t282 + t338 + t347;
t295 = t36 * qJD(1) - t164 * qJD(2) - t257 * qJD(3);
t182 = t276 + t194;
t290 = t182 * qJD(1) + t198 * t379 - t370;
t286 = t198 * t489 + t95 * t457 + t96 * t460 + t479 * t461;
t260 = t373 / 0.2e1;
t215 = (t372 - qJD(3) / 0.2e1) * t280;
t206 = t209 * qJD(3);
t144 = t198 * t378;
t139 = t256 + t227 / 0.2e1 - t302;
t138 = t368 - t226 / 0.2e1 + t301;
t130 = t218 * t447 + t307 * t282;
t124 = 0.2e1 * t350 + t340;
t122 = pkin(3) * t308 + t265;
t110 = t212 * qJD(2) + t218 * t382;
t104 = pkin(3) * t309 + t243;
t82 = t449 - t312;
t68 = t216 * t456 + t218 * t454 + t265 + t420 / 0.2e1 + t444 / 0.2e1;
t50 = t432 / 0.2e1;
t41 = -t195 * t455 + t198 * t454 + t243 + t423 / 0.2e1 + t445 / 0.2e1;
t37 = t257 * t448 - t402 / 0.2e1 + t338 - t347;
t28 = t289 + t291;
t26 = t300 + t469;
t20 = t300 - t468 + t478;
t15 = t259 * t449 + t283 - t291;
t13 = -t284 + t293;
t8 = t316 * pkin(3) + t228 * t265 + t264 * t243 - t292 + t50;
t6 = t84 * t458 + t83 * t459 + t286 + t296;
t4 = t71 * t458 + t72 * t460 + t286 + t311;
t2 = t50 + t285 + t317;
t34 = [0, 0, 0, t354, t248 * qJD(2), 0, 0, 0, -pkin(1) * t373, -pkin(1) * t371, -t274 * t356 + t275 * t354, -t223 * qJD(3) - 0.2e1 * t282 * t341, -t225 * qJD(2) + t280 * t359, t224 * qJD(2) + t280 * t358, -t354, -t100 * qJD(2) - t152 * qJD(3), t99 * qJD(2) - t151 * qJD(3), t17 * qJD(2) + t18 * qJD(3) + t487, t21 * qJD(2) + t24 * qJD(3) + t33 * qJD(4), t23 * qJD(2) + t29 * qJD(3) + (qJD(4) * t282 - t390) * t198, t11 * qJD(2) + t16 * qJD(3) + t195 * t369 + t487, -t22 * qJD(2) - t30 * qJD(3) + t182 * qJD(5) + t282 * t383, t9 * qJD(2) + t10 * qJD(3) + t31 * qJD(4) - t38 * qJD(5); 0, 0, 0, t353, t375, t371, -t373, 0, -pkin(6) * t371 - t365, pkin(6) * t373 - t364, -t206 + (t275 * t374 + t355) * t282, -0.2e1 * t280 * t356 + t467 * t282, t279 * t373 - t376, t352 + t377, -t215, -t389 + (t279 * t339 - t366) * qJD(2) + t139 * qJD(3), t396 + (t281 * t339 + t253) * qJD(2) + t138 * qJD(3), t427 + (-t92 * t216 - t91 * t218 + t344) * qJD(2) + t6 * qJD(3) + t471, t421 + (-t91 * t171 + t229 * t264 + t479 * t92) * qJD(2) + t8 * qJD(3) + t26 * qJD(4), t418 + (t109 * t216 + t136 * t197 - t171 * t280) * qJD(2) + t15 * qJD(3) + t130 * qJD(4) + t82 * qJD(5), t429 + (-t79 * t216 + t80 * t218 + t344) * qJD(2) + t4 * qJD(3) + t471 - t129 * qJD(5), -t419 + (-t109 * t218 - t136 * t200 + t280 * t479) * qJD(2) + t13 * qJD(3) - t128 * qJD(4) + t144, t434 + (t109 * t136 + t80 * t171 + t479 * t79) * qJD(2) + t2 * qJD(3) + t20 * qJD(4) + t28 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t299, -t298, t279 * t321, t305, t260, t139 * qJD(2) - t190 * qJD(3) - t384, t138 * qJD(2) + qJD(3) * t395 - t385, t424 + t6 * qJD(2) + (t195 * t278 - t198 * t277) * t439, t417 + t8 * qJD(2) + t104 * qJD(4) + (t277 * t96 - t278 * t95) * t439, t15 * qJD(2) - t95 * qJD(3) + t401, t428 + t4 * qJD(2) + (-t259 * t195 - t257 * t198) * qJD(3) - t390, t13 * qJD(2) + t96 * qJD(3) - t369 - t400, t431 + t2 * qJD(2) + (t96 * t257 + t95 * t259) * qJD(3) + t41 * qJD(4) + t37 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t486, t26 * qJD(2) + t104 * qJD(3) + t398, t130 * qJD(2) + t198 * t372, t486, -t343, t20 * qJD(2) + t41 * qJD(3) + t399; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82 * qJD(2) - t357, -t303, t290, t28 * qJD(2) + t37 * qJD(3) - t397; 0, 0, 0, -t353, -t375, 0, 0, 0, t365, t364, -t275 * t353 - t206, 0.2e1 * t279 * t305, -t358 + t376, t359 - t377, t215, t161 * qJD(3) + t389, t160 * qJD(3) - t396, -t5 * qJD(3) - t427 + t471, -t7 * qJD(3) - t25 * qJD(4) - t421, t14 * qJD(3) - t127 * qJD(4) - t81 * qJD(5) - t418, -t3 * qJD(3) - t128 * qJD(5) - t429 + t471, -t12 * qJD(3) - t129 * qJD(4) + t144 + t419, -t1 * qJD(3) - t19 * qJD(4) - t27 * qJD(5) - t434; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t356, t247 * qJD(3), 0, 0, 0, -pkin(2) * t392, -pkin(2) * t391, t488, t39 * qJD(3) + t472, t44 * qJD(3) - t216 * t378, t488, t45 * qJD(3) + t212 * qJD(5), t32 * qJD(3) - t136 * t378 + t472; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t304, t467, -t251 * t281, t251 * t279, -t374 / 0.2e1, -pkin(7) * t391 - t314, pkin(7) * t392 - t315, (t216 * t278 - t218 * t277) * t439 - t335, t122 * qJD(4) + (-t171 * t277 - t278 * t479) * t439 + t334, -qJD(3) * t479 + t332, (-t259 * t216 - t257 * t218) * qJD(3) - qJD(5) * t216 - t336, -qJD(3) * t171 + t333, (-t171 * t257 + t259 * t479) * qJD(3) + t68 * qJD(4) + t124 * qJD(5) + t337; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t485, t122 * qJD(3) + t330, -t388, t485, -t386, t68 * qJD(3) + t331; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t306, -t380 - t387, t110, t124 * qJD(3) - t136 * t379 - t416; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t299, t298, (-t279 * t374 + t393) * t282, (-t360 - t394) * t282, t260, -t161 * qJD(2) + t384, -t160 * qJD(2) + t385, t5 * qJD(2) - t424, t7 * qJD(2) + t105 * qJD(4) - t417, -t14 * qJD(2) - t381 - t401, t3 * qJD(2) - t428, t12 * qJD(2) - t369 - t383 + t400, t1 * qJD(2) - t40 * qJD(4) - t36 * qJD(5) - t431; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t304, -t467, t281 * t372, -t279 * t372, t374 / 0.2e1, t314, t315, t335, t123 * qJD(4) - t334, -t218 * qJD(4) - t332, t336, -t216 * qJD(4) - t333, -t67 * qJD(4) + t164 * qJD(5) - t337; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5), t257 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t320, -t132, 0, t319, -t329; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t251, -t295; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t486, t25 * qJD(2) - t105 * qJD(3) - t398, t127 * qJD(2) - t198 * t251, -t486, t303, t19 * qJD(2) + t40 * qJD(3) - t198 * qJD(5) - t399; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t485, -t123 * qJD(3) - t330, t218 * qJD(3) + t388, -t485, t380 + t386, t67 * qJD(3) - t331 - t378; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t320, t132, 0, -t319, t329; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t132; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81 * qJD(2) + t357, t343, -t290, t27 * qJD(2) + t36 * qJD(3) + t381 + t397; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t306, t387, -t110, t416 - t164 * qJD(3) + (qJD(2) * t136 + qJD(4)) * t218; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t251, t295; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t132; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t34;
