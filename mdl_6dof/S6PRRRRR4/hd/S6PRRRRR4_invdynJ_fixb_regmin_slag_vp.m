% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRRRRR4
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% tau_reg [6x32]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRRRRR4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR4_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR4_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRR4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR4_invdynJ_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:00:04
% EndTime: 2019-03-09 01:00:32
% DurationCPUTime: 12.84s
% Computational Cost: add. (11137->629), mult. (27574->929), div. (0->0), fcn. (24011->18), ass. (0->320)
t268 = sin(pkin(7));
t276 = sin(qJ(3));
t281 = cos(qJ(3));
t315 = (pkin(3) * t276 - pkin(10) * t281) * t268;
t269 = sin(pkin(6));
t277 = sin(qJ(2));
t420 = t269 * t277;
t371 = qJD(1) * t420;
t493 = qJD(3) * t315 - t268 * t371;
t271 = cos(pkin(7));
t282 = cos(qJ(2));
t409 = t281 * t282;
t414 = t276 * t277;
t312 = -t271 * t414 + t409;
t187 = t312 * t269;
t423 = t268 * t276;
t254 = pkin(9) * t423;
t415 = t271 * t281;
t463 = pkin(2) * t415 - t254;
t492 = qJD(1) * t187 - t463 * qJD(3);
t278 = cos(qJ(6));
t387 = qJD(6) * t278;
t396 = qJD(2) * t271;
t253 = qJD(3) + t396;
t280 = cos(qJ(4));
t275 = sin(qJ(4));
t397 = qJD(2) * t268;
t370 = t276 * t397;
t343 = t275 * t370;
t184 = t253 * t280 - t343;
t185 = t253 * t275 + t280 * t370;
t274 = sin(qJ(5));
t279 = cos(qJ(5));
t126 = -t279 * t184 + t185 * t274;
t480 = t126 * t278;
t491 = t387 + t480;
t416 = t271 * t276;
t422 = t268 * t281;
t401 = pkin(2) * t416 + pkin(9) * t422;
t196 = pkin(10) * t271 + t401;
t331 = -pkin(3) * t281 - pkin(10) * t276;
t197 = (-pkin(2) + t331) * t268;
t464 = t280 * t196 + t275 * t197;
t490 = qJD(4) * t464 - t492 * t275 - t280 * t493;
t391 = qJD(4) * t280;
t392 = qJD(4) * t275;
t489 = -t196 * t392 + t197 * t391 + t275 * t493 - t492 * t280;
t213 = -t280 * t271 + t275 * t423;
t393 = qJD(3) * t281;
t366 = t268 * t393;
t342 = t280 * t366;
t159 = -qJD(4) * t213 + t342;
t394 = qJD(3) * t276;
t367 = t268 * t394;
t488 = -pkin(4) * t367 + pkin(11) * t159 + t490;
t376 = t280 * t423;
t214 = t271 * t275 + t376;
t364 = t275 * t393;
t160 = qJD(4) * t214 + t268 * t364;
t487 = -pkin(11) * t160 + t489;
t395 = qJD(2) * t281;
t369 = t268 * t395;
t344 = t275 * t369;
t452 = pkin(10) + pkin(11);
t373 = qJD(4) * t452;
t217 = pkin(9) * t397 + t371;
t272 = cos(pkin(6));
t399 = qJD(1) * t272;
t372 = t268 * t399;
t398 = qJD(1) * t282;
t228 = qJD(2) * pkin(2) + t269 * t398;
t428 = t228 * t271;
t132 = -t276 * t217 + t281 * (t372 + t428);
t203 = qJD(2) * t315;
t407 = t280 * t132 + t275 * t203;
t486 = pkin(11) * t344 - t275 * t373 - t407;
t190 = t280 * t203;
t410 = t280 * t281;
t485 = t280 * t373 - t275 * t132 + t190 + (pkin(4) * t276 - pkin(11) * t410) * t397;
t273 = sin(qJ(6));
t319 = t184 * t274 + t279 * t185;
t382 = qJDD(2) * t271;
t252 = qJDD(3) + t382;
t111 = qJD(2) * t342 - qJD(4) * t343 + qJDD(2) * t376 + t275 * t252 + t253 * t391;
t381 = qJDD(2) * t276;
t112 = -t280 * t252 + t253 * t392 + t268 * (qJD(2) * (t276 * t391 + t364) + t275 * t381);
t44 = qJD(5) * t319 + t274 * t111 + t279 * t112;
t40 = qJDD(6) + t44;
t38 = t278 * t40;
t385 = -qJD(6) - t126;
t241 = -qJD(4) + t369;
t229 = -qJD(5) + t241;
t95 = t278 * t229 + t273 * t319;
t484 = -t385 ^ 2 * t273 + t95 * t319 + t38;
t380 = qJDD(2) * t281;
t250 = t268 * t380;
t384 = qJD(2) * qJD(3);
t363 = t276 * t384;
t198 = t268 * t363 + qJDD(4) - t250;
t193 = qJDD(5) + t198;
t388 = qJD(6) * t273;
t389 = qJD(5) * t279;
t390 = qJD(5) * t274;
t43 = t279 * t111 - t274 * t112 + t184 * t389 - t185 * t390;
t26 = t273 * t193 - t229 * t387 + t278 * t43 - t319 * t388;
t24 = t26 * t273;
t97 = -t229 * t273 + t278 * t319;
t483 = t491 * t97 + t24;
t482 = t273 * t40 - t319 * t97 - t385 * t491;
t439 = sin(pkin(13));
t356 = t439 * t282;
t270 = cos(pkin(13));
t418 = t270 * t277;
t211 = t272 * t418 + t356;
t357 = t439 * t277;
t417 = t270 * t282;
t296 = -t272 * t417 + t357;
t421 = t269 * t270;
t375 = t268 * t421;
t117 = t211 * t281 + (-t271 * t296 - t375) * t276;
t212 = -t272 * t357 + t417;
t297 = t272 * t356 + t418;
t358 = t269 * t439;
t335 = t268 * t358;
t119 = t212 * t281 + (-t271 * t297 + t335) * t276;
t412 = t277 * t281;
t413 = t276 * t282;
t313 = t271 * t413 + t412;
t154 = t269 * t313 + t272 * t423;
t155 = t268 * t296 - t271 * t421;
t156 = t268 * t297 + t271 * t358;
t419 = t269 * t282;
t210 = -t268 * t419 + t271 * t272;
t267 = qJ(4) + qJ(5);
t262 = sin(t267);
t263 = cos(t267);
t307 = -g(3) * (-t154 * t262 + t210 * t263) - g(2) * (-t117 * t262 + t155 * t263) - g(1) * (-t119 * t262 + t156 * t263);
t251 = qJDD(1) * t419;
t368 = qJD(2) * t420;
t341 = qJD(1) * t368;
t194 = qJDD(2) * pkin(2) + t251 - t341;
t383 = qJDD(1) * t272;
t361 = t268 * t383;
t477 = qJDD(2) * t268 * pkin(9) + (qJD(2) * t398 + qJDD(1) * t277) * t269 + qJD(3) * t372;
t294 = t194 * t416 - t217 * t394 + t276 * t361 + t281 * t477 + t393 * t428;
t59 = pkin(10) * t252 + t294;
t133 = t281 * t217 + t228 * t416 + t276 * t372;
t115 = pkin(10) * t253 + t133;
t249 = t271 * t399;
t143 = t249 + (qJD(2) * t331 - t228) * t268;
t70 = t115 * t280 + t143 * t275;
t248 = t271 * t383;
t305 = t363 - t380;
t362 = t281 * t384;
t306 = t362 + t381;
t93 = t248 + (pkin(3) * t305 - pkin(10) * t306 - t194) * t268;
t87 = t280 * t93;
t288 = -qJD(4) * t70 - t275 * t59 + t87;
t16 = pkin(4) * t198 - pkin(11) * t111 + t288;
t308 = t115 * t392 - t143 * t391 - t275 * t93 - t280 * t59;
t18 = -pkin(11) * t112 - t308;
t57 = pkin(11) * t184 + t70;
t442 = t279 * t57;
t69 = -t115 * t275 + t280 * t143;
t56 = -pkin(11) * t185 + t69;
t48 = -pkin(4) * t241 + t56;
t29 = t274 * t48 + t442;
t289 = -qJD(5) * t29 + t279 * t16 - t274 * t18;
t5 = -t193 * pkin(5) - t289;
t299 = t307 - t5;
t443 = t274 * t57;
t28 = t279 * t48 - t443;
t21 = pkin(5) * t229 - t28;
t481 = t126 * t21;
t479 = t319 * t126;
t218 = t274 * t275 - t279 * t280;
t461 = qJD(4) + qJD(5);
t157 = t461 * t218;
t166 = t218 * t369;
t405 = -t157 + t166;
t219 = t274 * t280 + t275 * t279;
t404 = (-t369 + t461) * t219;
t314 = t271 * t412 + t413;
t186 = t314 * t269;
t365 = t271 * t394;
t402 = pkin(2) * t365 + pkin(9) * t366 - qJD(1) * t186;
t476 = -t126 ^ 2 + t319 ^ 2;
t68 = pkin(5) * t319 + pkin(12) * t126;
t475 = -t126 * t229 + t43;
t109 = t154 * t263 + t210 * t262;
t355 = -t274 * t16 - t279 * t18 - t48 * t389 + t57 * t390;
t74 = t117 * t263 + t155 * t262;
t76 = t119 * t263 + t156 * t262;
t114 = -t253 * pkin(3) - t132;
t83 = -pkin(4) * t184 + t114;
t474 = g(1) * t76 + g(2) * t74 + g(3) * t109 + t126 * t83 + t355;
t27 = qJD(6) * t97 - t278 * t193 + t273 * t43;
t324 = t273 * t97 + t278 * t95;
t473 = -t126 * t324 + t26 * t278 - t273 * t27 - t95 * t387 - t388 * t97;
t351 = -t196 * t275 + t280 * t197;
t92 = -pkin(4) * t422 - pkin(11) * t214 + t351;
t99 = -pkin(11) * t213 + t464;
t470 = -t274 * t488 + t487 * t279 + t92 * t389 - t390 * t99;
t445 = t274 * t92 + t279 * t99;
t469 = qJD(5) * t445 + t487 * t274 + t279 * t488;
t468 = t385 * t319;
t244 = t452 * t275;
t245 = t452 * t280;
t317 = -t244 * t279 - t245 * t274;
t467 = qJD(5) * t317 - t274 * t485 + t279 * t486;
t169 = -t244 * t274 + t245 * t279;
t466 = qJD(5) * t169 + t274 * t486 + t279 * t485;
t334 = -t133 + (-t344 + t392) * pkin(4);
t406 = pkin(4) * t160 + t402;
t462 = t271 * t409 - t414;
t153 = -t269 * t462 - t272 * t422;
t150 = t153 * t278;
t120 = -t154 * t275 + t210 * t280;
t121 = t154 * t280 + t210 * t275;
t64 = t120 * t274 + t121 * t279;
t465 = -t273 * t64 + t150;
t22 = -pkin(12) * t229 + t29;
t36 = pkin(5) * t126 - pkin(12) * t319 + t83;
t325 = t22 * t273 - t278 * t36;
t460 = t21 * t388 + t319 * t325;
t9 = t22 * t278 + t273 * t36;
t459 = t21 * t387 - t273 * t299 + t9 * t319;
t457 = -t319 * t83 + t289 + t307;
t455 = -t229 * t319 - t44;
t4 = pkin(12) * t193 - t355;
t316 = -t194 * t415 + t217 * t393 + t228 * t365 + t276 * t477 - t281 * t361;
t60 = -pkin(3) * t252 + t316;
t35 = pkin(4) * t112 + t60;
t7 = pkin(5) * t44 - pkin(12) * t43 + t35;
t1 = -t325 * qJD(6) + t273 * t7 + t278 * t4;
t448 = -pkin(5) * t367 + t469;
t441 = pkin(5) * t370 + t466;
t434 = t153 * t273;
t433 = t184 * t241;
t432 = t185 * t241;
t431 = t210 * t268;
t430 = t219 * t273;
t429 = t219 * t278;
t427 = t262 * t268;
t426 = t263 * t273;
t425 = t263 * t278;
t264 = t268 ^ 2;
t283 = qJD(2) ^ 2;
t424 = t264 * t283;
t411 = t277 * t283;
t408 = qJDD(1) - g(3);
t265 = t276 ^ 2;
t400 = -t281 ^ 2 + t265;
t386 = qJD(3) - t253;
t377 = t273 * t422;
t261 = -pkin(4) * t280 - pkin(3);
t259 = pkin(4) * t274 + pkin(12);
t353 = pkin(4) * t185 + qJD(6) * t259 + t68;
t349 = t253 + t396;
t348 = t252 + t382;
t347 = t264 * t269 * t411;
t345 = t268 * t368;
t31 = t274 * t56 + t442;
t338 = pkin(4) * t390 - t31;
t32 = t279 * t56 - t443;
t337 = -pkin(4) * t389 + t32;
t146 = t279 * t213 + t214 * t274;
t147 = -t213 * t274 + t214 * t279;
t195 = t254 + (-pkin(2) * t281 - pkin(3)) * t271;
t148 = t213 * pkin(4) + t195;
t62 = pkin(5) * t146 - pkin(12) * t147 + t148;
t336 = -pkin(12) * t367 - qJD(6) * t62 - t470;
t333 = -pkin(5) * t404 + pkin(12) * t405 + qJD(6) * t169 - t334;
t46 = -pkin(12) * t422 + t445;
t66 = -qJD(5) * t146 + t279 * t159 - t274 * t160;
t67 = qJD(5) * t147 + t274 * t159 + t279 * t160;
t332 = -pkin(5) * t67 + pkin(12) * t66 + qJD(6) * t46 - t406;
t330 = g(1) * t212 + g(2) * t211;
t149 = pkin(5) * t218 - pkin(12) * t219 + t261;
t328 = pkin(12) * t370 - qJD(6) * t149 - t467;
t326 = -t259 * t40 + t481;
t322 = -t274 * t99 + t279 * t92;
t321 = t278 * t64 + t434;
t320 = t120 * t279 - t121 * t274;
t122 = t147 * t273 + t278 * t422;
t292 = t296 * t281;
t116 = t211 * t276 + t271 * t292 + t281 * t375;
t293 = t297 * t281;
t118 = t212 * t276 + t271 * t293 - t281 * t335;
t303 = g(1) * t118 + g(2) * t116 + g(3) * t153;
t137 = -t211 * t416 - t292;
t139 = -t212 * t416 - t293;
t302 = g(1) * t139 + g(2) * t137 + g(3) * t187;
t144 = -t166 * t273 - t278 * t370;
t301 = -t157 * t273 + t219 * t387 - t144;
t145 = -t166 * t278 + t273 * t370;
t300 = -t157 * t278 - t219 * t388 - t145;
t295 = -g(3) * t420 - t330;
t291 = -pkin(10) * t198 - t114 * t241;
t2 = -qJD(6) * t9 - t273 * t4 + t278 * t7;
t286 = pkin(10) * qJD(4) * t241 + t303 - t60;
t260 = -pkin(4) * t279 - pkin(5);
t178 = -t268 * t228 + t249;
t151 = -t194 * t268 + t248;
t142 = t187 * t263 + t420 * t427;
t138 = t212 * t415 - t276 * t297;
t136 = t211 * t415 - t276 * t296;
t123 = t147 * t278 - t377;
t107 = t272 * t366 + (t312 * qJD(2) + qJD(3) * t462) * t269;
t106 = t272 * t367 + (qJD(2) * t314 + qJD(3) * t313) * t269;
t91 = t139 * t263 + t212 * t427;
t90 = t137 * t263 + t211 * t427;
t51 = qJD(4) * t120 + t107 * t280 + t275 * t345;
t50 = -qJD(4) * t121 - t107 * t275 + t280 * t345;
t45 = pkin(5) * t422 - t322;
t42 = -qJD(6) * t377 + t147 * t387 + t273 * t66 - t278 * t367;
t41 = -qJD(6) * t122 + t273 * t367 + t278 * t66;
t15 = qJD(5) * t64 + t274 * t51 - t279 * t50;
t14 = qJD(5) * t320 + t274 * t50 + t279 * t51;
t3 = [t408, 0 (qJDD(2) * t282 - t411) * t269 (-qJDD(2) * t277 - t282 * t283) * t269, 0, 0, 0, 0, 0, -t106 * t253 - t153 * t252 - t281 * t347 + t305 * t431, -t107 * t253 - t154 * t252 + t276 * t347 + t306 * t431, 0, 0, 0, 0, 0, -t106 * t184 + t112 * t153 + t120 * t198 - t241 * t50, t106 * t185 + t111 * t153 - t121 * t198 + t241 * t51, 0, 0, 0, 0, 0, t106 * t126 + t15 * t229 + t153 * t44 + t193 * t320, t106 * t319 + t14 * t229 + t153 * t43 - t193 * t64, 0, 0, 0, 0, 0 -(-qJD(6) * t321 + t106 * t278 - t273 * t14) * t385 + t465 * t40 + t15 * t95 - t320 * t27 (qJD(6) * t465 + t106 * t273 + t278 * t14) * t385 - t321 * t40 + t15 * t97 - t320 * t26; 0, qJDD(2), g(1) * t297 + g(2) * t296 - g(3) * t419 + t251, -t408 * t420 + t330 (qJDD(2) * t265 + 0.2e1 * t276 * t362) * t264, 0.2e1 * (t276 * t380 - t384 * t400) * t264 (t276 * t348 + t349 * t393) * t268 (t281 * t348 - t349 * t394) * t268, t252 * t271, t463 * t252 - t316 * t271 + (-t151 * t281 + t178 * t394) * t268 - t402 * t253 + (-pkin(2) * t305 + t281 * t341) * t264 - t302, -t401 * t252 - t294 * t271 + g(1) * t138 + g(2) * t136 + g(3) * t186 + (t151 * t276 + t178 * t393) * t268 + t492 * t253 + (-pkin(2) * t306 - t276 * t341) * t264, t111 * t214 + t159 * t185, -t111 * t213 - t112 * t214 + t159 * t184 - t160 * t185, -t159 * t241 + t214 * t198 + (-t111 * t281 + t185 * t394) * t268, t160 * t241 - t213 * t198 + (t112 * t281 + t184 * t394) * t268 (-t198 * t281 - t241 * t394) * t268, t351 * t198 + t195 * t112 + t60 * t213 + t114 * t160 - t302 * t280 + (-(-t115 * t391 + t87) * t281 + t69 * t394 + (-(-qJD(4) * t143 - t59) * t281 + t295) * t275) * t268 + t490 * t241 - t402 * t184, t195 * t111 + t60 * t214 + t114 * t159 - t464 * t198 + t302 * t275 + (t280 * t295 - t281 * t308 - t394 * t70) * t268 + t489 * t241 + t402 * t185, t147 * t43 + t319 * t66, -t126 * t66 - t146 * t43 - t147 * t44 - t319 * t67, t147 * t193 - t66 * t229 + (-t281 * t43 + t319 * t394) * t268, -t146 * t193 + t67 * t229 + (-t126 * t394 + t281 * t44) * t268 (-t193 * t281 - t229 * t394) * t268, t322 * t193 + t148 * t44 + t35 * t146 + t83 * t67 - g(1) * t91 - g(2) * t90 - g(3) * t142 + (t28 * t394 - t281 * t289) * t268 + t469 * t229 + t406 * t126, -t445 * t193 + t148 * t43 + t35 * t147 + t83 * t66 + t302 * t262 + (t263 * t295 - t281 * t355 - t29 * t394) * t268 + t470 * t229 + t406 * t319, t123 * t26 + t41 * t97, -t122 * t26 - t123 * t27 - t41 * t95 - t42 * t97, t123 * t40 + t146 * t26 - t385 * t41 + t67 * t97, -t122 * t40 - t146 * t27 + t385 * t42 - t67 * t95, t146 * t40 - t385 * t67 (-t273 * t46 + t278 * t62) * t40 + t2 * t146 - t325 * t67 + t45 * t27 + t5 * t122 + t21 * t42 - g(1) * (t138 * t273 + t278 * t91) - g(2) * (t136 * t273 + t278 * t90) - g(3) * (t142 * t278 + t186 * t273) + t448 * t95 - (t273 * t336 - t278 * t332) * t385 -(t273 * t62 + t278 * t46) * t40 - t1 * t146 - t9 * t67 + t45 * t26 + t5 * t123 + t21 * t41 - g(1) * (t138 * t278 - t273 * t91) - g(2) * (t136 * t278 - t273 * t90) - g(3) * (-t142 * t273 + t186 * t278) + t448 * t97 - (t273 * t332 + t278 * t336) * t385; 0, 0, 0, 0, -t276 * t281 * t424, t400 * t424 (t386 * t395 + t381) * t268, -t370 * t386 + t250, t252, t133 * t253 - t178 * t370 + t303 - t316, g(1) * t119 + g(2) * t117 + g(3) * t154 + t132 * t253 - t178 * t369 - t294, t111 * t275 - t280 * t432 (t111 - t433) * t280 + (-t112 + t432) * t275, -t241 * t391 + t275 * t198 + (-t185 * t276 + t241 * t410) * t397, t241 * t392 + t280 * t198 + (-t241 * t275 * t281 - t184 * t276) * t397, t241 * t370, -t69 * t370 - pkin(3) * t112 + t133 * t184 + t190 * t241 + (-t132 * t241 + t291) * t275 + t286 * t280, -pkin(3) * t111 - t133 * t185 - t241 * t407 - t275 * t286 + t280 * t291 + t370 * t70, t43 * t219 + t319 * t405, -t126 * t405 - t218 * t43 - t219 * t44 - t319 * t404, t219 * t193 - t229 * t405 - t319 * t370, t126 * t370 - t218 * t193 + t229 * t404, t229 * t370, t334 * t126 + t193 * t317 + t35 * t218 + t229 * t466 + t261 * t44 + t303 * t263 - t28 * t370 + t404 * t83, -t169 * t193 + t35 * t219 + t229 * t467 + t261 * t43 - t303 * t262 + t29 * t370 + t334 * t319 + t405 * t83, t26 * t429 + t300 * t97, t144 * t97 + t145 * t95 + t324 * t157 + (-t24 - t27 * t278 + (t273 * t95 - t278 * t97) * qJD(6)) * t219, t218 * t26 - t300 * t385 + t40 * t429 + t404 * t97, -t218 * t27 + t301 * t385 - t40 * t430 - t404 * t95, t40 * t218 - t385 * t404 (t149 * t278 - t169 * t273) * t40 + t2 * t218 - t317 * t27 + t5 * t430 - g(1) * (-t118 * t425 + t119 * t273) - g(2) * (-t116 * t425 + t117 * t273) - g(3) * (-t153 * t425 + t154 * t273) + t441 * t95 - t404 * t325 - (t273 * t328 - t278 * t333) * t385 + t301 * t21 -(t149 * t273 + t169 * t278) * t40 - t1 * t218 - t317 * t26 + t5 * t429 - g(1) * (t118 * t426 + t119 * t278) - g(2) * (t116 * t426 + t117 * t278) - g(3) * (t153 * t426 + t154 * t278) + t441 * t97 - t404 * t9 - (t273 * t333 + t278 * t328) * t385 + t300 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t185 * t184, -t184 ^ 2 + t185 ^ 2, t111 + t433, -t112 - t432, t198, -t70 * t241 - t114 * t185 - g(1) * (-t119 * t275 + t156 * t280) - g(2) * (-t117 * t275 + t155 * t280) - g(3) * t120 + t288, -t114 * t184 - t69 * t241 - g(1) * (-t119 * t280 - t156 * t275) - g(2) * (-t117 * t280 - t155 * t275) + g(3) * t121 + t308, t479, t476, t475, t455, t193, -t229 * t31 + (-t126 * t185 + t193 * t279 + t229 * t390) * pkin(4) + t457, -t229 * t32 + (-t185 * t319 - t193 * t274 + t229 * t389) * pkin(4) + t474, t483, t473, t482, t484, t468, t260 * t27 + t338 * t95 + (-t337 * t385 + t326) * t273 + (t353 * t385 + t299) * t278 + t460, t260 * t26 + t338 * t97 + t326 * t278 - (t273 * t353 + t278 * t337) * t385 + t459; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t479, t476, t475, t455, t193, -t229 * t29 + t457, -t229 * t28 + t474, t483, t473, t482, t484, t468, -pkin(5) * t27 - t29 * t95 + (-pkin(12) * t40 - t28 * t385 + t481) * t273 + (-(-pkin(12) * qJD(6) - t68) * t385 + t299) * t278 + t460, -pkin(5) * t26 - (t273 * t68 + t278 * t28) * t385 - t29 * t97 + t21 * t480 + (-t385 * t388 - t38) * pkin(12) + t459; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97 * t95, -t95 ^ 2 + t97 ^ 2, -t385 * t95 + t26, -t385 * t97 - t27, t40, -t9 * t385 - t21 * t97 - g(1) * (t118 * t278 - t273 * t76) - g(2) * (t116 * t278 - t273 * t74) - g(3) * (-t109 * t273 + t150) + t2, t325 * t385 + t21 * t95 - g(1) * (-t118 * t273 - t278 * t76) - g(2) * (-t116 * t273 - t278 * t74) - g(3) * (-t109 * t278 - t434) - t1;];
tau_reg  = t3;
