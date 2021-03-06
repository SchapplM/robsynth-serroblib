% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RRPRRP5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 17:56
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RRPRRP5_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP5_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP5_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP5_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP5_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP5_invdynJ_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 17:53:58
% EndTime: 2019-05-06 17:54:25
% DurationCPUTime: 11.34s
% Computational Cost: add. (80309->640), mult. (211367->908), div. (0->0), fcn. (169067->12), ass. (0->401)
t333 = sin(pkin(11));
t335 = cos(pkin(11));
t342 = cos(qJ(2));
t334 = sin(pkin(6));
t413 = qJD(1) * t334;
t390 = t342 * t413;
t339 = sin(qJ(2));
t391 = t339 * t413;
t305 = t333 * t390 + t335 * t391;
t336 = cos(pkin(6));
t328 = qJD(1) * t336 + qJD(2);
t338 = sin(qJ(4));
t341 = cos(qJ(4));
t287 = t305 * t341 + t328 * t338;
t337 = sin(qJ(5));
t340 = cos(qJ(5));
t303 = t333 * t391 - t335 * t390;
t367 = qJD(4) + t303;
t262 = t287 * t337 - t340 * t367;
t264 = t340 * t287 + t337 * t367;
t218 = t264 * t262;
t405 = t342 * qJDD(1);
t408 = qJD(1) * qJD(2);
t370 = -t339 * t408 + t405;
t359 = t370 * t334;
t406 = qJDD(1) * t339;
t371 = t342 * t408 + t406;
t360 = t371 * t334;
t279 = t333 * t359 + t335 * t360;
t327 = t336 * qJDD(1) + qJDD(2);
t385 = t338 * t279 - t341 * t327;
t242 = -t287 * qJD(4) - t385;
t241 = qJDD(5) - t242;
t486 = -t218 + t241;
t494 = pkin(5) * t486;
t285 = t305 * t338 - t341 * t328;
t376 = -t341 * t279 - t338 * t327;
t243 = -t285 * qJD(4) - t376;
t350 = t333 * t360 - t335 * t359;
t349 = -qJDD(4) - t350;
t191 = -t262 * qJD(5) + t340 * t243 - t337 * t349;
t282 = qJD(5) + t285;
t235 = t282 * t262;
t168 = t235 + t191;
t493 = qJ(6) * t168;
t278 = t305 * t303;
t483 = -t278 + t327;
t492 = t333 * t483;
t491 = t335 * t483;
t248 = t287 * t285;
t484 = -t248 - t349;
t490 = t338 * t484;
t489 = t341 * t484;
t439 = t486 * t337;
t438 = t486 * t340;
t261 = t264 ^ 2;
t280 = t282 ^ 2;
t215 = -t261 - t280;
t260 = t262 ^ 2;
t386 = t243 * t337 + t340 * t349;
t190 = -qJD(5) * t264 - t386;
t230 = pkin(5) * t282 - qJ(6) * t264;
t343 = qJD(1) ^ 2;
t472 = sin(qJ(1));
t473 = cos(qJ(1));
t369 = g(1) * t473 + g(2) * t472;
t407 = qJDD(1) * t334;
t309 = -t343 * pkin(1) + pkin(8) * t407 - t369;
t425 = t334 * t339;
t329 = g(3) * t425;
t368 = g(1) * t472 - g(2) * t473;
t353 = qJDD(1) * pkin(1) + t368;
t423 = t334 * t343;
t348 = pkin(8) * t423 + t353;
t346 = t336 * t348;
t382 = qJD(1) * (-qJD(2) + t328);
t331 = t334 ^ 2;
t426 = t331 * t343;
t396 = t342 * t426;
t479 = t328 ^ 2;
t231 = -t479 * pkin(2) - t329 + (qJ(3) * t334 * t382 + t346) * t339 + (-pkin(2) * t396 + qJ(3) * t407 + t309) * t342;
t323 = t339 * t396;
t310 = t327 + t323;
t384 = -t339 * t309 + t342 * t346;
t345 = t310 * pkin(2) + (-t342 * g(3) + (t342 * t382 - t406) * qJ(3)) * t334 + t384;
t178 = -0.2e1 * qJD(3) * t303 + t335 * t231 + t333 * t345;
t275 = pkin(3) * t303 - pkin(9) * t305;
t157 = -pkin(3) * t479 + pkin(9) * t327 - t275 * t303 + t178;
t429 = t328 * t305;
t249 = t350 + t429;
t455 = t336 * g(3);
t294 = t334 * t348 + t455;
t478 = t342 ^ 2;
t395 = t478 * t426;
t480 = qJ(3) * t395 - qJDD(3);
t293 = t328 * t303;
t482 = -t293 + t279;
t344 = (pkin(2) * t328 - qJ(3) * t391) * t391 - pkin(2) * t359 - t294 - t482 * pkin(9) + t249 * pkin(3) - t480;
t125 = t341 * t157 + t338 * t344;
t246 = pkin(4) * t285 - pkin(10) * t287;
t363 = t367 ^ 2;
t103 = -pkin(4) * t363 - pkin(10) * t349 - t285 * t246 + t125;
t387 = t333 * t231 - t335 * t345;
t477 = 0.2e1 * qJD(3);
t156 = -t327 * pkin(3) - t479 * pkin(9) + (t477 + t275) * t305 + t387;
t267 = t367 * t285;
t210 = t243 - t267;
t358 = t367 * t287;
t116 = -t210 * pkin(10) + (-t242 + t358) * pkin(4) + t156;
t76 = t340 * t103 + t337 * t116;
t372 = t190 * qJ(6) - 0.2e1 * qJD(6) * t262 - t230 * t282 + t76;
t488 = -t372 + (t215 + t260) * pkin(5);
t485 = -t235 + t191;
t424 = t334 * t342;
t268 = g(3) * t424 - t384;
t269 = t342 * t309 + t339 * t346 - t329;
t481 = t339 * t268 + t342 * t269;
t165 = t264 * (qJD(5) - t282) + t386;
t283 = t285 ^ 2;
t284 = t287 ^ 2;
t301 = t303 ^ 2;
t302 = t305 ^ 2;
t332 = t339 ^ 2;
t122 = -t165 * t340 + t168 * t337;
t197 = -t260 - t261;
t95 = t122 * t338 - t197 * t341;
t476 = pkin(3) * t95;
t410 = qJD(6) * t264;
t256 = -0.2e1 * t410;
t75 = t103 * t337 - t340 * t116;
t362 = -t493 - t75 + t494;
t60 = t256 + t362;
t475 = pkin(5) * t60;
t474 = pkin(9) * t95;
t177 = t305 * t477 + t387;
t129 = -t335 * t177 + t333 * t178;
t471 = pkin(2) * t129;
t254 = t293 + t279;
t347 = -t350 + t429;
t213 = -t254 * t335 + t333 * t347;
t470 = pkin(2) * t213;
t202 = -t280 - t260;
t140 = t202 * t340 - t439;
t164 = (qJD(5) + t282) * t264 + t386;
t105 = t140 * t338 - t164 * t341;
t469 = pkin(3) * t105;
t182 = t218 + t241;
t440 = t182 * t340;
t146 = -t215 * t337 - t440;
t108 = t146 * t338 - t341 * t485;
t468 = pkin(3) * t108;
t467 = pkin(3) * t333;
t139 = t202 * t337 + t438;
t466 = pkin(4) * t139;
t441 = t182 * t337;
t145 = t215 * t340 - t441;
t465 = pkin(4) * t145;
t464 = pkin(4) * t338;
t463 = pkin(4) * t341;
t462 = pkin(5) * t168;
t461 = pkin(8) * t334;
t460 = pkin(9) * t105;
t459 = pkin(9) * t108;
t120 = -t165 * t337 - t168 * t340;
t458 = pkin(10) * t120;
t457 = pkin(10) * t139;
t456 = pkin(10) * t145;
t96 = t122 * t341 + t197 * t338;
t68 = -t120 * t335 + t333 * t96;
t69 = t120 * t333 + t335 * t96;
t454 = pkin(1) * (-t334 * t95 + (t339 * t69 + t342 * t68) * t336) + (-t339 * t68 + t342 * t69) * t461;
t106 = t140 * t341 + t164 * t338;
t82 = t106 * t333 - t139 * t335;
t83 = t106 * t335 + t139 * t333;
t453 = pkin(1) * (-t105 * t334 + (t339 * t83 + t342 * t82) * t336) + (-t339 * t82 + t342 * t83) * t461;
t109 = t146 * t341 + t338 * t485;
t90 = t109 * t333 - t145 * t335;
t91 = t109 * t335 + t145 * t333;
t452 = pkin(1) * (-t108 * t334 + (t339 * t91 + t342 * t90) * t336) + (-t339 * t90 + t342 * t91) * t461;
t451 = qJ(3) * t68;
t450 = qJ(3) * t82;
t449 = qJ(3) * t90;
t448 = t337 * t60;
t447 = t340 * t60;
t124 = t157 * t338 - t341 * t344;
t102 = t349 * pkin(4) - t363 * pkin(10) + t246 * t287 + t124;
t446 = t102 * t337;
t445 = t102 * t340;
t444 = t129 * t339;
t443 = t156 * t338;
t442 = t156 * t341;
t222 = t248 - t349;
t437 = t222 * t338;
t436 = t222 * t341;
t383 = qJD(1) * (-qJD(2) - t328);
t245 = t455 + ((t339 * t383 + t405) * pkin(2) + t353 + (qJ(3) * t332 + pkin(8)) * t423) * t334 + t480;
t435 = t245 * t333;
t434 = t245 * t335;
t272 = t278 + t327;
t433 = t272 * t333;
t432 = t272 * t335;
t431 = t282 * t337;
t430 = t282 * t340;
t428 = t328 * t333;
t427 = t328 * t335;
t421 = t339 * t310;
t311 = -t323 + t327;
t419 = t342 * t311;
t418 = -pkin(4) * t197 + pkin(10) * t122;
t417 = pkin(4) * t164 - pkin(10) * t140;
t416 = pkin(4) * t485 - pkin(10) * t146;
t404 = pkin(2) * t68 - pkin(3) * t120 + pkin(9) * t96;
t403 = pkin(2) * t82 - pkin(3) * t139 + pkin(9) * t106;
t402 = pkin(2) * t90 - pkin(3) * t145 + pkin(9) * t109;
t401 = t338 * t218;
t400 = t341 * t218;
t399 = t333 * t248;
t398 = t335 * t248;
t397 = t332 * t426;
t394 = t336 * t278;
t393 = -pkin(3) * t335 - pkin(2);
t392 = -pkin(2) * t95 + qJ(3) * t69;
t389 = -pkin(2) * t105 + qJ(3) * t83;
t388 = -pkin(2) * t108 + qJ(3) * t91;
t45 = t337 * t75 + t340 * t76;
t85 = t124 * t338 + t341 * t125;
t130 = t177 * t333 + t335 * t178;
t380 = -pkin(4) * t102 + pkin(10) * t45;
t379 = -t417 - t445;
t378 = -t416 + t446;
t44 = t337 * t76 - t340 * t75;
t84 = -t124 * t341 + t125 * t338;
t375 = qJD(1) * t328 - t336 * t343;
t374 = t45 + t418;
t366 = 0.2e1 * qJD(4) + t303;
t50 = -qJ(6) * t165 + (-t197 - t260) * pkin(5) + t372;
t257 = 0.2e1 * t410;
t52 = t257 - t362 + t493;
t365 = t337 * t52 + t340 * t50 + t418;
t127 = -pkin(5) * t485 - qJ(6) * t182;
t78 = -t190 * pkin(5) - t260 * qJ(6) + t230 * t264 + qJDD(6) + t102;
t77 = -qJ(6) * t215 + t78;
t364 = t127 * t340 + t337 * t77 - t416;
t65 = -pkin(5) * t164 + qJ(6) * t202 - t78;
t361 = -qJ(6) * t439 + t340 * t65 - t417;
t357 = t338 * t267;
t356 = t338 * t358;
t355 = t341 * t267;
t354 = t341 * t358;
t352 = t362 + t494;
t62 = -pkin(5) * t260 + t372;
t29 = t340 * t62 - t448;
t39 = -pkin(5) * t78 + qJ(6) * t62;
t351 = -pkin(4) * t78 + pkin(10) * t29 - qJ(6) * t448 + t340 * t39;
t207 = -t287 * t303 + t385;
t317 = t336 * t327;
t316 = t328 * t390;
t315 = t328 * t391;
t313 = (t332 - t478) * t426;
t312 = -t479 - t395;
t300 = -t479 - t397;
t292 = -t302 + t479;
t291 = t301 - t479;
t290 = -t315 + t359;
t289 = t315 + t359;
t288 = -t316 + t360;
t281 = -t302 - t479;
t277 = t302 - t301;
t270 = -t479 - t301;
t266 = -t284 + t363;
t265 = t283 - t363;
t259 = -t301 - t302;
t247 = t284 - t283;
t244 = -t284 - t363;
t240 = -t363 - t283;
t239 = -t281 * t333 - t432;
t238 = t281 * t335 - t433;
t233 = -t261 + t280;
t232 = t260 - t280;
t229 = t283 + t284;
t225 = t270 * t335 - t492;
t224 = t270 * t333 + t491;
t220 = -t355 + t356;
t219 = -t357 - t354;
t216 = t261 - t260;
t214 = t254 * t333 + t335 * t347;
t212 = t285 * t366 + t376;
t211 = t243 + t267;
t208 = -t287 * t366 - t385;
t206 = t341 * t243 - t356;
t205 = t338 * t243 + t354;
t204 = -t338 * t242 + t355;
t203 = t341 * t242 + t357;
t201 = t265 * t341 - t437;
t200 = -t266 * t338 + t489;
t199 = t265 * t338 + t436;
t198 = t266 * t341 + t490;
t195 = -t244 * t338 - t436;
t194 = t244 * t341 - t437;
t193 = (-t262 * t340 + t264 * t337) * t282;
t192 = (-t262 * t337 - t264 * t340) * t282;
t185 = t240 * t341 - t490;
t184 = t240 * t338 + t489;
t175 = -t207 * t341 + t211 * t338;
t174 = t208 * t341 - t210 * t338;
t173 = -t207 * t338 - t211 * t341;
t172 = t208 * t338 + t210 * t341;
t171 = t193 * t341 + t241 * t338;
t170 = t193 * t338 - t241 * t341;
t161 = t191 * t340 - t264 * t431;
t160 = t191 * t337 + t264 * t430;
t159 = -t190 * t337 + t262 * t430;
t158 = -t340 * t190 - t262 * t431;
t154 = t232 * t340 - t441;
t153 = -t233 * t337 + t438;
t152 = t232 * t337 + t440;
t151 = t233 * t340 + t439;
t150 = pkin(2) * t238 - t178;
t149 = pkin(2) * t224 - t177;
t148 = t195 * t335 - t212 * t333;
t147 = t195 * t333 + t212 * t335;
t142 = t185 * t335 - t208 * t333;
t141 = t185 * t333 + t208 * t335;
t136 = t175 * t335 - t229 * t333;
t135 = t175 * t333 + t229 * t335;
t134 = t161 * t341 + t401;
t133 = t159 * t341 - t401;
t132 = t161 * t338 - t400;
t131 = t159 * t338 + t400;
t128 = -pkin(9) * t194 + t442;
t126 = -pkin(9) * t184 + t443;
t121 = -t164 * t340 - t337 * t485;
t119 = -t164 * t337 + t340 * t485;
t113 = t154 * t341 - t165 * t338;
t112 = t153 * t341 + t168 * t338;
t111 = t154 * t338 + t165 * t341;
t110 = t153 * t338 - t168 * t341;
t100 = -pkin(3) * t194 + t125;
t99 = t121 * t341 + t216 * t338;
t98 = t121 * t338 - t216 * t341;
t97 = -pkin(3) * t184 + t124;
t93 = -pkin(4) * t120 + t462;
t92 = pkin(2) * t147 + pkin(3) * t212 + pkin(9) * t195 + t443;
t87 = pkin(2) * t141 + pkin(3) * t208 + pkin(9) * t185 - t442;
t86 = t445 - t456;
t79 = t446 - t457;
t73 = t336 * t170 + (t339 * (t171 * t335 + t192 * t333) + t342 * (t171 * t333 - t192 * t335)) * t334;
t72 = -pkin(9) * t173 - t84;
t71 = t156 * t333 + t335 * t85;
t70 = -t156 * t335 + t333 * t85;
t64 = t76 - t465;
t63 = t75 - t466;
t61 = pkin(2) * t135 + pkin(3) * t229 + pkin(9) * t175 + t85;
t59 = t336 * t132 + (t339 * (t134 * t335 + t160 * t333) + t342 * (t134 * t333 - t160 * t335)) * t334;
t58 = t336 * t131 + (t339 * (t133 * t335 - t158 * t333) + t342 * (t133 * t333 + t158 * t335)) * t334;
t57 = -t378 - t468;
t55 = -t379 - t469;
t54 = -t127 * t337 + t340 * t77 - t456;
t51 = -qJ(6) * t438 - t337 * t65 - t457;
t49 = -t465 - t488;
t48 = t257 - t352 - t466;
t47 = t336 * t111 + (t339 * (t113 * t335 + t152 * t333) + t342 * (t113 * t333 - t152 * t335)) * t334;
t46 = t336 * t110 + (t339 * (t112 * t335 + t151 * t333) + t342 * (t112 * t333 - t151 * t335)) * t334;
t42 = pkin(2) * t70 - pkin(3) * t156 + pkin(9) * t85;
t38 = -t364 - t468;
t37 = -t361 - t469;
t36 = -t44 - t458;
t35 = t102 * t338 + t341 * t45;
t34 = -t102 * t341 + t338 * t45;
t33 = -t338 * t64 + t341 * t86 - t459;
t32 = t336 * t98 + (t339 * (t119 * t333 + t335 * t99) + t342 * (-t119 * t335 + t333 * t99)) * t334;
t31 = -t338 * t63 + t341 * t79 - t460;
t28 = t337 * t62 + t447;
t27 = -t374 - t476;
t26 = t120 * t464 + t341 * t36 - t474;
t25 = t338 * t86 + t341 * t64 + t402;
t24 = t29 * t341 + t338 * t78;
t23 = t29 * t338 - t341 * t78;
t22 = t338 * t79 + t341 * t63 + t403;
t21 = -t337 * t50 + t340 * t52 - t458;
t20 = -t338 * t49 + t341 * t54 - t459;
t19 = -t338 * t48 + t341 * t51 - t460;
t18 = -pkin(4) * t28 - t475;
t17 = t333 * t44 + t335 * t35;
t16 = t333 * t35 - t335 * t44;
t15 = -t365 - t476;
t14 = t338 * t54 + t341 * t49 + t402;
t13 = -t120 * t463 + t338 * t36 + t404;
t12 = -pkin(3) * t34 - t380;
t11 = t21 * t341 - t338 * t93 - t474;
t10 = t338 * t51 + t341 * t48 + t403;
t9 = -pkin(10) * t28 - qJ(6) * t447 - t337 * t39;
t8 = t24 * t335 + t28 * t333;
t7 = t24 * t333 - t28 * t335;
t6 = -pkin(9) * t34 + (-pkin(10) * t341 + t464) * t44;
t5 = t21 * t338 + t341 * t93 + t404;
t4 = -pkin(3) * t23 - t351;
t3 = pkin(2) * t16 + pkin(9) * t35 + (-pkin(10) * t338 - pkin(3) - t463) * t44;
t2 = -pkin(9) * t23 - t18 * t338 + t341 * t9;
t1 = pkin(2) * t7 - pkin(3) * t28 + pkin(9) * t24 + t18 * t341 + t338 * t9;
t30 = [0, 0, 0, 0, 0, qJDD(1), t368, t369, 0, 0, (t342 * t375 + t371) * t331 * t339, t336 * t313 + (t339 * t290 + (t316 + t360) * t342) * t334, t336 * t288 + (t421 + t342 * (t479 - t397)) * t334, (-t339 * t375 + t370) * t331 * t342, t336 * t289 + (t339 * (-t479 + t395) + t419) * t334, t317, (-t268 + pkin(1) * (t310 * t342 + t312 * t339)) * t336 + (t342 * t294 + pkin(1) * t290 + pkin(8) * (t312 * t342 - t421)) * t334, -t294 * t425 - t336 * t269 + pkin(1) * ((t300 * t342 - t311 * t339) * t336 + (t342 * t383 - t406) * t331) + (-t300 * t339 - t419) * t461, pkin(1) * ((-t288 * t342 + t289 * t339) * t336 - (-t332 - t478) * t331 * t423) + (t288 * t339 + t289 * t342) * t461 + t481 * t334, pkin(1) * (t334 * t294 + (-t268 * t342 + t269 * t339) * t336) + t481 * t461, t394 + (t339 * (t279 * t335 - t305 * t428) + t342 * (t279 * t333 + t305 * t427)) * t334, t336 * t277 + (t339 * (-t249 * t335 - t333 * t482) + t342 * (-t249 * t333 + t335 * t482)) * t334, t336 * t254 + (t339 * (-t292 * t333 + t491) + t342 * (t292 * t335 + t492)) * t334, (t303 * t427 + t333 * t350) * t425 + (t303 * t428 - t335 * t350) * t424 - t394, t336 * t347 + (t339 * (t291 * t335 - t433) + t342 * (t291 * t333 + t432)) * t334, t317 + (t339 * (-t303 * t335 + t305 * t333) + t342 * (-t303 * t333 - t305 * t335)) * t334 * t328, (t149 + pkin(1) * (t224 * t342 + t225 * t339)) * t336 + (t339 * (-qJ(3) * t224 - t435) + t342 * (-pkin(2) * t249 + qJ(3) * t225 + t434) - pkin(1) * t249 + pkin(8) * (-t224 * t339 + t225 * t342)) * t334, (t150 + pkin(1) * (t238 * t342 + t239 * t339)) * t336 + (t339 * (-qJ(3) * t238 - t434) + t342 * (-pkin(2) * t482 + qJ(3) * t239 - t435) - pkin(1) * t482 + pkin(8) * (-t238 * t339 + t239 * t342)) * t334, (t470 + pkin(1) * (t213 * t342 + t214 * t339)) * t336 + (t339 * (-qJ(3) * t213 - t129) + t342 * (-pkin(2) * t259 + qJ(3) * t214 + t130) - pkin(1) * t259 + pkin(8) * (-t213 * t339 + t214 * t342)) * t334, (t471 + pkin(1) * (t129 * t342 + t130 * t339)) * t336 + (-qJ(3) * t444 + t342 * (pkin(2) * t245 + qJ(3) * t130) + pkin(1) * t245 + pkin(8) * (t130 * t342 - t444)) * t334, t336 * t205 + (t339 * (t206 * t335 + t399) + t342 * (t206 * t333 - t398)) * t334, t336 * t172 + (t339 * (t174 * t335 + t247 * t333) + t342 * (t174 * t333 - t247 * t335)) * t334, t336 * t198 + (t339 * (t200 * t335 + t211 * t333) + t342 * (t200 * t333 - t211 * t335)) * t334, t336 * t203 + (t339 * (t204 * t335 - t399) + t342 * (t204 * t333 + t398)) * t334, t336 * t199 + (t339 * (t201 * t335 - t207 * t333) + t342 * (t201 * t333 + t207 * t335)) * t334, (t335 * t220 - t333 * t349) * t425 + (t333 * t220 + t335 * t349) * t424 + t336 * t219, (t87 + pkin(1) * (t141 * t342 + t142 * t339)) * t336 + (t339 * (-qJ(3) * t141 + t126 * t335 - t333 * t97) + t342 * (-pkin(2) * t184 + qJ(3) * t142 + t126 * t333 + t335 * t97) - pkin(1) * t184 + pkin(8) * (-t141 * t339 + t142 * t342)) * t334, (t92 + pkin(1) * (t147 * t342 + t148 * t339)) * t336 + (t339 * (-qJ(3) * t147 - t100 * t333 + t128 * t335) + t342 * (-pkin(2) * t194 + qJ(3) * t148 + t100 * t335 + t128 * t333) - pkin(1) * t194 + pkin(8) * (-t147 * t339 + t148 * t342)) * t334, (t61 + pkin(1) * (t135 * t342 + t136 * t339)) * t336 + (t339 * (-qJ(3) * t135 + t335 * t72) + t342 * (qJ(3) * t136 + t333 * t72) + pkin(8) * (-t135 * t339 + t136 * t342) + (t339 * t467 + t342 * t393 - pkin(1)) * t173) * t334, (t42 + pkin(1) * (t339 * t71 + t342 * t70)) * t336 + ((t339 * (-pkin(9) * t335 + t467) + t342 * (-pkin(9) * t333 + t393) - pkin(1)) * t84 + (pkin(8) + qJ(3)) * (-t339 * t70 + t342 * t71)) * t334, t59, t32, t46, t58, t47, t73, t336 * t22 + (t339 * (t31 * t335 - t333 * t55 - t450) + t342 * (t31 * t333 + t335 * t55 + t389)) * t334 + t453, t336 * t25 + (t339 * (t33 * t335 - t333 * t57 - t449) + t342 * (t33 * t333 + t335 * t57 + t388)) * t334 + t452, t336 * t13 + (t339 * (t26 * t335 - t27 * t333 - t451) + t342 * (t26 * t333 + t27 * t335 + t392)) * t334 + t454, (t3 + pkin(1) * (t16 * t342 + t17 * t339)) * t336 + (t339 * (-qJ(3) * t16 - t12 * t333 + t335 * t6) + t342 * (-pkin(2) * t34 + qJ(3) * t17 + t12 * t335 + t333 * t6) - pkin(1) * t34 + pkin(8) * (-t16 * t339 + t17 * t342)) * t334, t59, t32, t46, t58, t47, t73, t336 * t10 + (t339 * (t19 * t335 - t333 * t37 - t450) + t342 * (t19 * t333 + t335 * t37 + t389)) * t334 + t453, t336 * t14 + (t339 * (t20 * t335 - t333 * t38 - t449) + t342 * (t20 * t333 + t335 * t38 + t388)) * t334 + t452, t336 * t5 + (t339 * (t11 * t335 - t15 * t333 - t451) + t342 * (t11 * t333 + t15 * t335 + t392)) * t334 + t454, (t1 + pkin(1) * (t339 * t8 + t342 * t7)) * t336 + (t339 * (-qJ(3) * t7 + t2 * t335 - t333 * t4) + t342 * (-pkin(2) * t23 + qJ(3) * t8 + t2 * t333 + t335 * t4) - pkin(1) * t23 + pkin(8) * (-t339 * t7 + t342 * t8)) * t334; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t323, t313, t288, t323, t289, t327, -t268, -t269, 0, 0, t278, t277, t254, -t278, t347, t327, t149, t150, t470, t471, t205, t172, t198, t203, t199, t219, t87, t92, t61, t42, t132, t98, t110, t131, t111, t170, t22, t25, t13, t3, t132, t98, t110, t131, t111, t170, t10, t14, t5, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t249, t482, t259, -t245, 0, 0, 0, 0, 0, 0, t184, t194, t173, t84, 0, 0, 0, 0, 0, 0, t105, t108, t95, t34, 0, 0, 0, 0, 0, 0, t105, t108, t95, t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t248, t247, t211, -t248, -t207, -t349, -t124, -t125, 0, 0, t160, t119, t151, -t158, t152, t192, t379, t378, t374, t380, t160, t119, t151, -t158, t152, t192, t361, t364, t365, t351; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t218, t216, t168, -t218, -t165, t241, -t75, -t76, 0, 0, t218, t216, t168, -t218, -t165, t241, t256 + t352, t488, -t462, t475; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t164, t485, t197, t78;];
tauJ_reg  = t30;
