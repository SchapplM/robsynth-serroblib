% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% tau_reg [6x35]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRRRR10_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR10_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR10_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR10_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR10_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR10_invdynJ_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:34:15
% EndTime: 2019-03-09 07:34:47
% DurationCPUTime: 14.21s
% Computational Cost: add. (24016->662), mult. (75919->942), div. (0->0), fcn. (67578->18), ass. (0->323)
t309 = sin(pkin(13));
t312 = cos(pkin(13));
t320 = cos(qJ(1));
t466 = cos(pkin(6));
t409 = t320 * t466;
t477 = sin(qJ(1));
t257 = t309 * t409 + t312 * t477;
t316 = sin(qJ(3));
t311 = sin(pkin(6));
t446 = t311 * t320;
t310 = sin(pkin(7));
t448 = t310 * t316;
t478 = cos(qJ(3));
t354 = t477 * t309 - t312 * t409;
t465 = cos(pkin(7));
t496 = t354 * t465;
t175 = -t257 * t478 + t316 * t496 + t446 * t448;
t411 = t311 * t465;
t231 = t310 * t354 - t320 * t411;
t308 = qJ(4) + qJ(5);
t303 = sin(t308);
t304 = cos(t308);
t133 = t175 * t304 - t231 * t303;
t422 = t310 * t478;
t172 = t257 * t316 + t422 * t446 + t478 * t496;
t313 = sin(qJ(6));
t317 = cos(qJ(6));
t532 = t133 * t313 + t172 * t317;
t531 = t133 * t317 - t172 * t313;
t410 = t316 * t465;
t332 = t311 * (-t309 * t410 + t312 * t478);
t245 = qJD(1) * t332;
t418 = qJD(3) * t478;
t393 = t310 * t418;
t530 = t245 - t393;
t431 = qJD(6) * t317;
t386 = t312 * t410;
t407 = t466 * t310;
t387 = t316 * t407;
t225 = t311 * (t309 * t478 + t386) + t387;
t217 = t225 * qJD(1);
t375 = t466 * t465;
t439 = qJD(1) * t311;
t420 = t312 * t439;
t252 = -qJD(1) * t375 + t310 * t420 - qJD(3);
t315 = sin(qJ(4));
t319 = cos(qJ(4));
t161 = -t315 * t217 - t252 * t319;
t162 = t217 * t319 - t252 * t315;
t314 = sin(qJ(5));
t318 = cos(qJ(5));
t104 = -t318 * t161 + t162 * t314;
t510 = t104 * t317;
t529 = t431 + t510;
t384 = t465 * t478;
t357 = t312 * t384;
t360 = t478 * t407;
t526 = t311 * t357 + t360;
t525 = t175 * t303 + t231 * t304;
t524 = t175 * t315 + t231 * t319;
t523 = t175 * t319 - t231 * t315;
t479 = pkin(10) + pkin(11);
t423 = qJD(4) * t479;
t404 = qJD(1) * t466;
t394 = pkin(1) * t404;
t255 = qJ(2) * t420 + t309 * t394;
t331 = (t312 * t411 + t407) * pkin(9);
t195 = qJD(1) * t331 + t255;
t451 = t309 * t310;
t248 = (-pkin(2) * t312 - pkin(9) * t451 - pkin(1)) * t311;
t237 = qJD(1) * t248 + qJD(2);
t290 = t312 * t394;
t450 = t309 * t311;
t325 = t466 * pkin(2) + (-pkin(9) * t465 - qJ(2)) * t450;
t208 = qJD(1) * t325 + t290;
t359 = t208 * t384;
t120 = -t316 * t195 + t237 * t422 + t359;
t421 = t309 * t439;
t504 = qJD(1) * t526 - t316 * t421;
t155 = pkin(3) * t217 - pkin(10) * t504;
t443 = t319 * t120 + t315 * t155;
t457 = t504 * t315;
t522 = -pkin(11) * t457 + t315 * t423 + t443;
t151 = t319 * t155;
t521 = pkin(4) * t217 - t120 * t315 + t151 + (-pkin(11) * t504 + t423) * t319;
t262 = t315 * t465 + t319 * t448;
t397 = t310 * t421;
t516 = qJD(4) * t262 - t315 * t530 + t319 * t397;
t261 = -t315 * t448 + t319 * t465;
t515 = -qJD(4) * t261 + t315 * t397 + t319 * t530;
t364 = t161 * t314 + t318 * t162;
t428 = qJDD(1) * t311;
t416 = t309 * t428;
t153 = qJD(3) * t504 + qJDD(1) * t387 + t386 * t428 + t478 * t416;
t427 = qJDD(1) * t312;
t415 = t311 * t427;
t251 = -qJDD(1) * t375 + t310 * t415 - qJDD(3);
t435 = qJD(4) * t319;
t436 = qJD(4) * t315;
t88 = t319 * t153 - t217 * t436 - t315 * t251 - t252 * t435;
t89 = qJD(4) * t162 + t153 * t315 + t319 * t251;
t46 = qJD(5) * t364 + t314 * t88 + t318 * t89;
t44 = qJDD(6) + t46;
t42 = t317 * t44;
t430 = -qJD(6) - t104;
t207 = qJD(4) - t504;
t203 = qJD(5) + t207;
t81 = -t317 * t203 + t313 * t364;
t514 = -t430 ^ 2 * t313 + t364 * t81 + t42;
t216 = t225 * qJD(3);
t154 = qJD(1) * t216 - qJDD(1) * t526 + t316 * t416;
t152 = qJDD(4) + t154;
t149 = qJDD(5) + t152;
t432 = qJD(6) * t313;
t433 = qJD(5) * t318;
t434 = qJD(5) * t314;
t45 = t161 * t433 - t162 * t434 - t314 * t89 + t318 * t88;
t27 = t313 * t149 + t203 * t431 + t317 * t45 - t364 * t432;
t25 = t27 * t313;
t83 = t203 * t313 + t317 * t364;
t513 = t529 * t83 + t25;
t512 = t313 * t44 - t364 * t83 - t430 * t529;
t385 = t466 * t477;
t258 = -t309 * t385 + t320 * t312;
t337 = t320 * t309 + t312 * t385;
t489 = -t477 * t311 * t310 + t337 * t465;
t177 = t258 * t478 - t316 * t489;
t232 = t310 * t337 + t411 * t477;
t134 = -t177 * t303 + t232 * t304;
t447 = t311 * t312;
t256 = t310 * t447 - t375;
t346 = -g(3) * (-t225 * t303 - t256 * t304) - g(2) * t525 - g(1) * t134;
t156 = -t208 * t310 + t465 * t237;
t97 = -pkin(3) * t504 - pkin(10) * t217 + t156;
t121 = t478 * t195 + t208 * t410 + t237 * t448;
t99 = -pkin(10) * t252 + t121;
t65 = t315 * t97 + t319 * t99;
t54 = pkin(11) * t161 + t65;
t469 = t318 * t54;
t64 = -t315 * t99 + t319 * t97;
t53 = -pkin(11) * t162 + t64;
t51 = pkin(4) * t207 + t53;
t23 = t314 * t51 + t469;
t424 = pkin(1) * t466;
t389 = qJDD(1) * t424;
t429 = qJD(1) * qJD(2);
t417 = t311 * t429;
t239 = qJ(2) * t415 + t309 * t389 + t312 * t417;
t182 = qJDD(1) * t331 + t239;
t288 = t312 * t389;
t183 = qJDD(1) * t325 - t309 * t417 + t288;
t233 = qJDD(1) * t248 + qJDD(2);
t437 = qJD(3) * t316;
t335 = qJD(3) * t359 + t478 * t182 + t183 * t410 - t195 * t437 + t233 * t448 + t237 * t393;
t67 = -pkin(10) * t251 + t335;
t143 = -t183 * t310 + t465 * t233;
t76 = pkin(3) * t154 - pkin(10) * t153 + t143;
t412 = -t315 * t67 + t319 * t76;
t329 = -qJD(4) * t65 + t412;
t14 = pkin(4) * t152 - pkin(11) * t88 + t329;
t352 = -t315 * t76 - t319 * t67 - t97 * t435 + t436 * t99;
t16 = -pkin(11) * t89 - t352;
t414 = -t318 * t14 + t16 * t314;
t481 = -qJD(5) * t23 - t414;
t5 = -pkin(5) * t149 - t481;
t339 = t346 - t5;
t470 = t314 * t54;
t22 = t318 * t51 - t470;
t20 = -pkin(5) * t203 - t22;
t511 = t104 * t20;
t509 = t364 * t104;
t268 = t314 * t315 - t318 * t319;
t145 = t268 * t504;
t487 = qJD(4) + qJD(5);
t228 = t487 * t268;
t442 = t145 - t228;
t269 = t314 * t319 + t315 * t318;
t441 = (-t504 + t487) * t269;
t506 = t309 * t384 + t312 * t316;
t505 = -t104 ^ 2 + t364 ^ 2;
t71 = pkin(5) * t364 + pkin(12) * t104;
t503 = t104 * t203 + t45;
t135 = t177 * t304 + t232 * t303;
t166 = t225 * t304 - t256 * t303;
t406 = -t314 * t14 - t318 * t16 - t51 * t433 + t54 * t434;
t98 = t252 * pkin(3) - t120;
t78 = -t161 * pkin(4) + t98;
t502 = g(1) * t135 - g(2) * t133 + g(3) * t166 + t104 * t78 + t406;
t28 = qJD(6) * t83 - t317 * t149 + t313 * t45;
t370 = t313 * t83 + t317 * t81;
t501 = -t104 * t370 + t27 * t317 - t313 * t28 - t81 * t431 - t432 * t83;
t392 = -t121 + (t436 - t457) * pkin(4);
t498 = t430 * t364;
t306 = t311 ^ 2;
t497 = t306 * (t309 ^ 2 + t312 ^ 2);
t296 = t312 * t424;
t226 = t296 + t325;
t164 = -t226 * t310 + t465 * t248;
t449 = t309 * t316;
t224 = t311 * t449 - t526;
t113 = pkin(3) * t224 - pkin(10) * t225 + t164;
t260 = qJ(2) * t447 + t309 * t424;
t220 = t331 + t260;
t425 = t478 * t220 + t226 * t410 + t248 * t448;
t119 = -pkin(10) * t256 + t425;
t444 = t315 * t113 + t319 * t119;
t362 = t261 * t318 - t262 * t314;
t495 = -qJD(5) * t362 + t314 * t516 + t318 * t515;
t197 = t261 * t314 + t262 * t318;
t494 = qJD(5) * t197 - t314 * t515 + t318 * t516;
t285 = t479 * t315;
t286 = t479 * t319;
t361 = -t285 * t318 - t286 * t314;
t493 = -qJD(5) * t361 + t521 * t314 + t318 * t522;
t243 = -t285 * t314 + t286 * t318;
t492 = -qJD(5) * t243 + t314 * t522 - t521 * t318;
t323 = -t316 * t220 + t226 * t384 + t248 * t422;
t170 = t225 * t315 + t256 * t319;
t171 = t225 * t319 - t256 * t315;
t123 = -t170 * t314 + t171 * t318;
t219 = t224 * t317;
t491 = -t123 * t313 + t219;
t244 = t506 * t439;
t419 = t310 * t437;
t490 = t244 - t419;
t21 = pkin(12) * t203 + t23;
t47 = t104 * pkin(5) - pkin(12) * t364 + t78;
t373 = t21 * t313 - t317 * t47;
t486 = t20 * t432 + t364 * t373;
t11 = t21 * t317 + t313 * t47;
t484 = t11 * t364 + t20 * t431 - t313 * t339;
t483 = -t364 * t78 + t346 + t481;
t480 = t203 * t364 - t46;
t4 = pkin(12) * t149 - t406;
t383 = qJD(3) * t410;
t353 = t316 * t182 - t183 * t384 + t195 * t418 + t208 * t383 - t233 * t422 + t237 * t419;
t68 = pkin(3) * t251 + t353;
t48 = pkin(4) * t89 + t68;
t9 = pkin(5) * t46 - pkin(12) * t45 + t48;
t1 = -t373 * qJD(6) + t313 * t9 + t317 * t4;
t476 = pkin(1) * t306;
t401 = t319 * t113 - t119 * t315;
t57 = pkin(4) * t224 - pkin(11) * t171 + t401;
t63 = -pkin(11) * t170 + t444;
t367 = t314 * t57 + t318 * t63;
t467 = pkin(5) * t217 - t492;
t459 = t161 * t207;
t458 = t162 * t207;
t456 = t224 * t313;
t455 = t269 * t313;
t454 = t269 * t317;
t453 = t304 * t313;
t452 = t304 * t317;
t438 = qJD(2) * t311;
t302 = -pkin(4) * t319 - pkin(3);
t321 = qJD(1) ^ 2;
t408 = t321 * t466;
t300 = pkin(4) * t314 + pkin(12);
t403 = pkin(4) * t162 + qJD(6) * t300 + t71;
t108 = qJD(2) * t332 + qJD(3) * t323;
t215 = (t360 + (t357 - t449) * t311) * qJD(3);
t395 = t438 * t451;
t142 = pkin(3) * t216 - pkin(10) * t215 + t395;
t402 = -t108 * t315 + t319 * t142;
t399 = t207 * t319;
t29 = t314 * t53 + t469;
t391 = pkin(4) * t434 - t29;
t30 = t318 * t53 - t470;
t390 = -pkin(4) * t433 + t30;
t380 = -pkin(5) * t441 + pkin(12) * t442 + qJD(6) * t243 - t392;
t221 = pkin(5) * t268 - pkin(12) * t269 + t302;
t379 = pkin(12) * t217 - qJD(6) * t221 + t493;
t378 = qJD(2) * t404;
t374 = -t300 * t44 + t511;
t32 = pkin(12) * t224 + t367;
t122 = t318 * t170 + t171 * t314;
t118 = t256 * pkin(3) - t323;
t84 = t170 * pkin(4) + t118;
t49 = t122 * pkin(5) - t123 * pkin(12) + t84;
t372 = t313 * t49 + t317 * t32;
t371 = -t313 * t32 + t317 * t49;
t130 = -qJD(4) * t170 + t215 * t319;
t34 = pkin(4) * t216 - pkin(11) * t130 - qJD(4) * t444 + t402;
t129 = qJD(4) * t171 + t215 * t315;
t347 = t319 * t108 + t113 * t435 - t119 * t436 + t315 * t142;
t36 = -pkin(11) * t129 + t347;
t369 = -t314 * t36 + t318 * t34;
t368 = -t314 * t63 + t318 * t57;
t109 = t220 * t418 + t226 * t383 + t248 * t419 + t438 * t506;
t94 = t123 * t317 + t456;
t363 = (-qJ(2) * t421 + t290) * t309 - t255 * t312;
t356 = g(1) * t477 - g(2) * t320;
t355 = -g(1) * t320 - g(2) * t477;
t351 = t314 * t34 + t318 * t36 + t57 * t433 - t434 * t63;
t349 = -pkin(10) * t152 + t207 * t98;
t176 = t258 * t316 + t478 * t489;
t345 = g(1) * t176 + g(2) * t172 + g(3) * t224;
t343 = -t313 * t197 - t317 * t422;
t342 = -t317 * t197 + t313 * t422;
t114 = -t145 * t313 - t317 * t217;
t341 = -t228 * t313 + t269 * t431 - t114;
t115 = -t145 * t317 + t217 * t313;
t340 = -t228 * t317 - t269 * t432 - t115;
t77 = pkin(4) * t129 + t109;
t2 = -qJD(6) * t11 - t313 * t4 + t317 * t9;
t326 = pkin(10) * qJD(4) * t207 - t345 + t68;
t301 = -pkin(4) * t318 - pkin(5);
t291 = -pkin(1) * t428 + qJDD(2);
t259 = -qJ(2) * t450 + t296;
t238 = t288 + (-qJ(2) * qJDD(1) - t429) * t450;
t139 = t177 * t319 + t232 * t315;
t138 = -t177 * t315 + t232 * t319;
t92 = t135 * t317 + t176 * t313;
t91 = -t135 * t313 + t176 * t317;
t61 = qJD(5) * t123 + t318 * t129 + t130 * t314;
t60 = -qJD(5) * t122 - t129 * t314 + t130 * t318;
t40 = qJD(6) * t94 - t216 * t317 + t313 * t60;
t39 = qJD(6) * t491 + t216 * t313 + t317 * t60;
t31 = -pkin(5) * t224 - t368;
t17 = pkin(5) * t61 - pkin(12) * t60 + t77;
t7 = -pkin(5) * t216 + qJD(5) * t367 - t369;
t6 = pkin(12) * t216 + t351;
t3 = [qJDD(1), t356, -t355, t238 * t466 + g(1) * t257 - g(2) * t258 + (-t291 * t312 - t309 * t378) * t311 + (t259 * t466 + t312 * t476) * qJDD(1), -t239 * t466 - g(1) * t354 + g(2) * t337 + (t291 * t309 - t312 * t378) * t311 + (-t260 * t466 - t309 * t476) * qJDD(1), t429 * t497 + (-t238 * t309 + t239 * t312 + (-t259 * t309 + t260 * t312) * qJDD(1) + t355) * t311, t238 * t259 + t239 * t260 + t356 * pkin(1) + (-t291 * pkin(1) + qJ(2) * t355 - qJD(2) * t363) * t311, t153 * t225 + t215 * t217, -t153 * t224 - t154 * t225 + t215 * t504 - t216 * t217, -t153 * t256 - t215 * t252 - t225 * t251, t154 * t256 + t216 * t252 + t224 * t251, t251 * t256, -g(1) * t175 - g(2) * t177 + t109 * t252 + t143 * t224 + t164 * t154 + t156 * t216 - t251 * t323 + t256 * t353 - t395 * t504, -g(1) * t172 + g(2) * t176 + t108 * t252 + t143 * t225 + t164 * t153 + t156 * t215 + t217 * t395 + t251 * t425 + t256 * t335, t130 * t162 + t171 * t88, -t129 * t162 + t130 * t161 - t170 * t88 - t171 * t89, t130 * t207 + t152 * t171 + t162 * t216 + t224 * t88, -t129 * t207 - t152 * t170 + t161 * t216 - t224 * t89, t152 * t224 + t207 * t216, t402 * t207 + t401 * t152 + t412 * t224 + t64 * t216 - t109 * t161 + t118 * t89 + t68 * t170 + t98 * t129 - g(1) * t523 - g(2) * t139 + (-t207 * t444 - t224 * t65) * qJD(4), g(1) * t524 - g(2) * t138 + t109 * t162 + t118 * t88 + t98 * t130 - t444 * t152 + t68 * t171 - t347 * t207 - t65 * t216 + t352 * t224, t123 * t45 + t364 * t60, -t104 * t60 - t122 * t45 - t123 * t46 - t364 * t61, t123 * t149 + t203 * t60 + t216 * t364 + t224 * t45, -t104 * t216 - t122 * t149 - t203 * t61 - t224 * t46, t149 * t224 + t203 * t216, t369 * t203 + t368 * t149 - t414 * t224 + t22 * t216 + t77 * t104 + t84 * t46 + t48 * t122 + t78 * t61 - g(1) * t133 - g(2) * t135 + (-t203 * t367 - t224 * t23) * qJD(5), g(1) * t525 - g(2) * t134 + t48 * t123 - t367 * t149 - t351 * t203 - t23 * t216 + t406 * t224 + t77 * t364 + t84 * t45 + t78 * t60, t27 * t94 + t39 * t83, t27 * t491 - t28 * t94 - t39 * t81 - t40 * t83, t122 * t27 - t39 * t430 + t44 * t94 + t61 * t83, -t122 * t28 + t40 * t430 + t44 * t491 - t61 * t81, t122 * t44 - t430 * t61 -(-qJD(6) * t372 + t17 * t317 - t313 * t6) * t430 + t371 * t44 + t2 * t122 - t373 * t61 + t7 * t81 + t31 * t28 - t5 * t491 + t20 * t40 - g(1) * t531 - g(2) * t92 (qJD(6) * t371 + t17 * t313 + t317 * t6) * t430 - t372 * t44 - t1 * t122 - t11 * t61 + t7 * t83 + t31 * t27 + t5 * t94 + t20 * t39 + g(1) * t532 - g(2) * t91; 0, 0, 0 (t309 * t408 - t427) * t311 (qJDD(1) * t309 + t312 * t408) * t311, -t321 * t497, -g(3) * t466 + qJDD(2) + (-pkin(1) * qJDD(1) + qJD(1) * t363 - t356) * t311, 0, 0, 0, 0, 0, t465 * t154 - t244 * t252 + (-t251 * t478 + t252 * t437 + t421 * t504) * t310, t465 * t153 - t245 * t252 + (-t217 * t421 + t251 * t316 + t252 * t418) * t310, 0, 0, 0, 0, 0, t261 * t152 + t490 * t161 - t207 * t516 - t89 * t422, -t262 * t152 - t490 * t162 + t207 * t515 - t88 * t422, 0, 0, 0, 0, 0, -t104 * t490 + t149 * t362 - t203 * t494 - t46 * t422, -t197 * t149 + t203 * t495 - t364 * t490 - t45 * t422, 0, 0, 0, 0, 0, -t362 * t28 + t343 * t44 + t494 * t81 - (qJD(6) * t342 + t313 * t495 - t317 * t490) * t430, -t362 * t27 + t342 * t44 + t494 * t83 - (-qJD(6) * t343 + t313 * t490 + t317 * t495) * t430; 0, 0, 0, 0, 0, 0, 0, -t217 * t504, t217 ^ 2 - t504 ^ 2, t252 * t504 + t153, -t217 * t252 - t154, -t251, -t121 * t252 - t156 * t217 + t345 - t353, g(1) * t177 - g(2) * t175 + g(3) * t225 - t120 * t252 - t156 * t504 - t335, t162 * t399 + t315 * t88 (t88 + t459) * t319 + (-t89 - t458) * t315, t315 * t152 - t162 * t217 + t207 * t399, -t207 ^ 2 * t315 + t152 * t319 - t161 * t217, -t207 * t217, -pkin(3) * t89 + t121 * t161 - t151 * t207 - t64 * t217 + (t120 * t207 + t349) * t315 - t326 * t319, -pkin(3) * t88 - t121 * t162 + t207 * t443 + t65 * t217 + t315 * t326 + t319 * t349, t269 * t45 + t364 * t442, -t104 * t442 - t268 * t45 - t269 * t46 - t364 * t441, t149 * t269 + t203 * t442 - t217 * t364, t104 * t217 - t149 * t268 - t203 * t441, -t203 * t217, t392 * t104 + t149 * t361 + t203 * t492 - t22 * t217 + t48 * t268 + t302 * t46 + t345 * t304 + t441 * t78, -t243 * t149 + t203 * t493 + t23 * t217 + t48 * t269 + t302 * t45 - t345 * t303 + t392 * t364 + t442 * t78, t27 * t454 + t340 * t83, t114 * t83 + t115 * t81 + t370 * t228 + (-t25 - t28 * t317 + (t313 * t81 - t317 * t83) * qJD(6)) * t269, t268 * t27 - t340 * t430 + t44 * t454 + t441 * t83, -t268 * t28 + t341 * t430 - t44 * t455 - t441 * t81, t268 * t44 - t430 * t441 (t221 * t317 - t243 * t313) * t44 + t2 * t268 - t361 * t28 + t5 * t455 - g(1) * (-t176 * t452 + t177 * t313) - g(2) * (-t172 * t452 - t175 * t313) - g(3) * (-t224 * t452 + t225 * t313) + t467 * t81 - (t313 * t379 - t317 * t380) * t430 - t441 * t373 + t341 * t20 -(t221 * t313 + t243 * t317) * t44 - t1 * t268 - t361 * t27 + t5 * t454 - g(1) * (t176 * t453 + t177 * t317) - g(2) * (t172 * t453 - t175 * t317) - g(3) * (t224 * t453 + t225 * t317) + t467 * t83 - t441 * t11 - (t313 * t380 + t317 * t379) * t430 + t340 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t162 * t161, -t161 ^ 2 + t162 ^ 2, t88 - t459, t458 - t89, t152, -g(1) * t138 - g(2) * t524 + g(3) * t170 - t98 * t162 + t65 * t207 + t329, g(1) * t139 - g(2) * t523 + g(3) * t171 - t98 * t161 + t64 * t207 + t352, t509, t505, t503, t480, t149, t203 * t29 + (-t104 * t162 + t149 * t318 - t203 * t434) * pkin(4) + t483, t203 * t30 + (-t149 * t314 - t162 * t364 - t203 * t433) * pkin(4) + t502, t513, t501, t512, t514, t498, t301 * t28 + t391 * t81 + (-t390 * t430 + t374) * t313 + (t403 * t430 + t339) * t317 + t486, t301 * t27 + t391 * t83 + t374 * t317 - (t313 * t403 + t317 * t390) * t430 + t484; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t509, t505, t503, t480, t149, t203 * t23 + t483, t203 * t22 + t502, t513, t501, t512, t514, t498, -pkin(5) * t28 - t23 * t81 + (-pkin(12) * t44 - t22 * t430 + t511) * t313 + (-(-pkin(12) * qJD(6) - t71) * t430 + t339) * t317 + t486, -pkin(5) * t27 - (t22 * t317 + t313 * t71) * t430 - t23 * t83 + t20 * t510 + (-t430 * t432 - t42) * pkin(12) + t484; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83 * t81, -t81 ^ 2 + t83 ^ 2, -t430 * t81 + t27, -t430 * t83 - t28, t44, -t11 * t430 - t20 * t83 - g(1) * t91 - g(2) * t532 - g(3) * (-t166 * t313 + t219) + t2, t373 * t430 + t20 * t81 + g(1) * t92 - g(2) * t531 - g(3) * (-t166 * t317 - t456) - t1;];
tau_reg  = t3;
