% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRRPRP10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% tau_reg [6x32]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRRPRP10_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP10_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP10_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP10_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP10_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP10_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:37:01
% EndTime: 2019-03-09 17:37:37
% DurationCPUTime: 18.50s
% Computational Cost: add. (18373->765), mult. (44864->1012), div. (0->0), fcn. (36090->14), ass. (0->337)
t323 = cos(qJ(2));
t316 = sin(pkin(6));
t447 = qJD(1) * t316;
t424 = t323 * t447;
t277 = -qJD(3) + t424;
t320 = sin(qJ(3));
t322 = cos(qJ(3));
t484 = cos(pkin(6));
t407 = t484 * qJD(1);
t380 = t407 + qJD(2);
t321 = sin(qJ(2));
t425 = t321 * t447;
t219 = t320 * t380 + t322 * t425;
t315 = sin(pkin(11));
t317 = cos(pkin(11));
t166 = t219 * t317 - t277 * t315;
t319 = sin(qJ(5));
t204 = t315 * t219;
t405 = t277 * t317 + t204;
t506 = cos(qJ(5));
t363 = t506 * t405;
t90 = t166 * t319 + t363;
t540 = t90 ^ 2;
t217 = t320 * t425 - t322 * t380;
t211 = qJD(5) + t217;
t539 = t211 * t90;
t289 = pkin(8) * t424;
t397 = pkin(1) * t407;
t242 = t321 * t397 + t289;
t538 = qJD(4) * t320 + t242 + t277 * (pkin(3) * t320 - qJ(4) * t322);
t383 = t319 * t405;
t534 = t166 * t506 - t383;
t508 = t534 ^ 2;
t537 = t211 * t534;
t239 = -pkin(8) * t425 + t323 * t397;
t390 = pkin(2) * t321 - pkin(9) * t323;
t240 = t390 * t447;
t453 = t322 * t239 + t320 * t240;
t133 = qJ(4) * t425 + t453;
t444 = qJD(3) * t320;
t434 = pkin(9) * t444;
t485 = -t538 * t317 + (t133 + t434) * t315;
t536 = t317 * t133 + t315 * t538;
t398 = t322 * t424;
t205 = t315 * t398 - t317 * t425;
t442 = qJD(3) * t322;
t531 = -t315 * t442 + t205;
t459 = t322 * t323;
t465 = t315 * t321;
t206 = (t317 * t459 + t465) * t447;
t535 = -t317 * t442 + t206;
t347 = qJD(3) * t380;
t403 = t484 * qJDD(1);
t375 = t403 + qJDD(2);
t445 = qJD(2) * t323;
t422 = t320 * t445;
t436 = t321 * qJDD(1);
t117 = t316 * (qJD(1) * (t321 * t442 + t422) + t320 * t436) + t320 * t347 - t322 * t375;
t111 = qJDD(5) + t117;
t399 = t320 * t424;
t460 = t317 * t322;
t533 = -pkin(4) * t399 + pkin(10) * t206 + (pkin(4) * t320 - pkin(10) * t460) * qJD(3) + t485;
t461 = t317 * t320;
t464 = t315 * t322;
t532 = -pkin(10) * t205 - (-pkin(9) * t461 - pkin(10) * t464) * qJD(3) + t536;
t262 = t315 * t506 + t319 * t317;
t251 = t262 * qJD(5);
t457 = -t251 * t320 + t319 * t531 - t506 * t535;
t415 = qJD(5) * t506;
t440 = qJD(5) * t319;
t517 = -t315 * t440 + t317 * t415;
t456 = -t205 * t506 - t206 * t319 + t262 * t442 + t320 * t517;
t358 = -t315 * t319 + t317 * t506;
t528 = -t358 * t217 - t517;
t451 = t262 * t217 + t251;
t220 = t320 * t239;
t135 = -pkin(3) * t425 - t240 * t322 + t220;
t527 = pkin(9) * t442 - t135;
t514 = t444 - t399;
t507 = cos(qJ(1));
t393 = t484 * t507;
t505 = sin(qJ(1));
t255 = t321 * t393 + t323 * t505;
t427 = t316 * t507;
t188 = t255 * t322 - t320 * t427;
t254 = t321 * t505 - t323 * t393;
t312 = pkin(11) + qJ(5);
t307 = sin(t312);
t308 = cos(t312);
t122 = t188 * t307 - t254 * t308;
t123 = t188 * t308 + t254 * t307;
t311 = t316 ^ 2;
t526 = 0.2e1 * t311;
t525 = -t317 * t434 - t536;
t524 = -pkin(4) * t531 + t527;
t200 = pkin(9) * t380 + t242;
t378 = -pkin(2) * t323 - pkin(9) * t321 - pkin(1);
t210 = t378 * t447;
t114 = -t320 * t200 + t322 * t210;
t148 = pkin(3) * t219 + qJ(4) * t217;
t77 = t317 * t114 + t315 * t148;
t523 = -qJD(4) * t317 + t77;
t496 = pkin(10) + qJ(4);
t272 = t496 * t315;
t273 = t496 * t317;
t359 = -t272 * t506 - t319 * t273;
t477 = t217 * t317;
t76 = -t114 * t315 + t317 * t148;
t61 = pkin(4) * t219 + pkin(10) * t477 + t76;
t478 = t217 * t315;
t69 = pkin(10) * t478 + t77;
t522 = -qJD(4) * t358 - qJD(5) * t359 + t319 * t61 + t506 * t69;
t197 = -t319 * t272 + t273 * t506;
t521 = -qJD(4) * t262 - qJD(5) * t197 + t319 * t69 - t506 * t61;
t504 = pkin(3) * t322;
t385 = qJ(4) * t320 + t504;
t269 = -pkin(2) - t385;
t260 = t317 * t269;
t175 = -pkin(10) * t461 + t260 + (-pkin(9) * t315 - pkin(4)) * t322;
t216 = pkin(9) * t460 + t315 * t269;
t466 = t315 * t320;
t186 = -pkin(10) * t466 + t216;
t520 = t319 * t175 + t506 * t186;
t428 = pkin(1) * t484;
t463 = t316 * t321;
t373 = -pkin(8) * t463 + t323 * t428;
t329 = t320 * t375 + t322 * t347;
t408 = t322 * t436;
t421 = t322 * t445;
t443 = qJD(3) * t321;
t327 = (t408 + (-t320 * t443 + t421) * qJD(1)) * t316 + t329;
t519 = -t175 * t415 + t186 * t440 - t319 * t533 + t532 * t506;
t518 = (qJDD(2) + 0.2e1 * t403) * t316;
t516 = t398 - t442;
t438 = qJDD(1) * t323;
t291 = t316 * t438;
t439 = qJD(1) * qJD(2);
t413 = t321 * t439;
t395 = t316 * t413;
t236 = qJDD(3) - t291 + t395;
t515 = -t236 * t322 - t277 * t444;
t498 = t236 * pkin(3);
t513 = qJDD(4) - t498;
t187 = t255 * t320 + t322 * t427;
t392 = t484 * t505;
t257 = -t321 * t392 + t323 * t507;
t426 = t316 * t505;
t191 = t257 * t320 - t322 * t426;
t252 = t320 * t463 - t322 * t484;
t349 = g(1) * t191 + g(2) * t187 + g(3) * t252;
t512 = -t197 * t111 - t349 * t307;
t253 = t320 * t484 + t322 * t463;
t184 = qJD(3) * t253 + t316 * t422;
t185 = -qJD(3) * t252 + t316 * t421;
t446 = qJD(2) * t321;
t423 = t316 * t446;
t345 = t185 * t317 + t315 * t423;
t462 = t316 * t323;
t449 = pkin(8) * t462 + t321 * t428;
t234 = pkin(9) * t484 + t449;
t450 = pkin(2) * t462 + pkin(9) * t463;
t235 = -pkin(1) * t316 - t450;
t364 = t390 * qJD(2);
t241 = t316 * t364;
t243 = t373 * qJD(2);
t353 = -t234 * t444 + t235 * t442 + t320 * t241 + t322 * t243;
t80 = (qJ(4) * t446 - qJD(4) * t323) * t316 + t353;
t244 = t449 * qJD(2);
t86 = t184 * pkin(3) - t185 * qJ(4) - t253 * qJD(4) + t244;
t40 = -t315 * t80 + t317 * t86;
t27 = t184 * pkin(4) - pkin(10) * t345 + t40;
t346 = -t185 * t315 + t317 * t423;
t41 = t315 * t86 + t317 * t80;
t35 = pkin(10) * t346 + t41;
t183 = t253 * t317 - t315 * t462;
t233 = -pkin(2) * t484 - t373;
t129 = t252 * pkin(3) - t253 * qJ(4) + t233;
t454 = t322 * t234 + t320 * t235;
t130 = -qJ(4) * t462 + t454;
t74 = t317 * t129 - t130 * t315;
t57 = pkin(4) * t252 - pkin(10) * t183 + t74;
t366 = t253 * t315 + t317 * t462;
t75 = t315 * t129 + t317 * t130;
t65 = -pkin(10) * t366 + t75;
t367 = t319 * t57 + t506 * t65;
t510 = -qJD(5) * t367 + t27 * t506 - t319 * t35;
t509 = -qJD(5) * t520 + t532 * t319 + t506 * t533;
t324 = qJD(1) ^ 2;
t503 = pkin(5) * t111;
t499 = g(3) * t316;
t309 = t320 * pkin(9);
t497 = t534 * t90;
t377 = qJD(2) * t397;
t394 = pkin(1) * t403;
t430 = pkin(8) * t291 + t321 * t394 + t323 * t377;
t341 = -pkin(8) * t395 + t430;
t146 = pkin(9) * t375 + t341;
t155 = (qJD(1) * t364 + qJDD(1) * t378) * t316;
t354 = t322 * t146 + t320 * t155 - t200 * t444 + t210 * t442;
t48 = qJ(4) * t236 - qJD(4) * t277 + t354;
t410 = t316 * t436;
t401 = pkin(8) * t410 + qJD(2) * t289 + t321 * t377 - t323 * t394;
t147 = -pkin(2) * t375 + t401;
t50 = t117 * pkin(3) - qJ(4) * t327 - t219 * qJD(4) + t147;
t20 = t315 * t50 + t317 * t48;
t495 = qJ(6) * t514 - qJD(6) * t322 - t519;
t494 = -pkin(5) * t514 - t509;
t238 = t358 * t320;
t491 = pkin(5) * t456 - qJ(6) * t457 - qJD(6) * t238 + t524;
t199 = -pkin(2) * t380 - t239;
t103 = t217 * pkin(3) - t219 * qJ(4) + t199;
t115 = t322 * t200 + t320 * t210;
t106 = -qJ(4) * t277 + t115;
t63 = t317 * t103 - t106 * t315;
t39 = pkin(4) * t217 - pkin(10) * t166 + t63;
t64 = t315 * t103 + t317 * t106;
t49 = -pkin(10) * t405 + t64;
t18 = t319 * t39 + t49 * t506;
t490 = t18 * t211;
t402 = -t320 * t146 + t322 * t155 - t200 * t442 - t210 * t444;
t51 = -t402 + t513;
t489 = t51 * t315;
t97 = -pkin(4) * t478 + t115;
t488 = pkin(5) * t451 + qJ(6) * t528 - qJD(6) * t262 - t97;
t487 = -qJ(6) * t219 - t522;
t486 = t219 * pkin(5) - t521;
t483 = qJ(4) * t315;
t482 = qJ(4) * t317;
t481 = qJ(6) * t111;
t480 = t115 * t277;
t476 = t219 * t277;
t472 = t254 * t320;
t256 = t321 * t507 + t323 * t392;
t471 = t256 * t320;
t356 = t277 * t320;
t469 = t307 * t322;
t468 = t308 * t322;
t467 = t311 * t324;
t17 = -t319 * t49 + t39 * t506;
t458 = qJD(6) - t17;
t263 = pkin(4) * t466 + t309;
t313 = t321 ^ 2;
t448 = -t323 ^ 2 + t313;
t435 = g(3) * t462;
t433 = t323 * t467;
t432 = t307 * t462;
t431 = t507 * pkin(1) + t257 * pkin(2) + pkin(8) * t426;
t305 = pkin(4) * t317 + pkin(3);
t429 = pkin(4) * t315 + pkin(9);
t105 = t277 * pkin(3) + qJD(4) - t114;
t416 = t105 * t442;
t414 = pkin(1) * t526;
t412 = qJD(1) * t443;
t411 = t323 * t439;
t19 = -t315 * t48 + t317 * t50;
t326 = t315 * t236 + t317 * t327;
t11 = t117 * pkin(4) - pkin(10) * t326 + t19;
t110 = t315 * t327;
t406 = t236 * t317 - t110;
t14 = pkin(10) * t406 + t20;
t4 = t506 * t11 - t319 * t14 - t39 * t440 - t49 * t415;
t404 = -t320 * t234 + t322 * t235;
t391 = t316 * t324 * t484;
t192 = t257 * t322 + t320 * t426;
t126 = t192 * t307 - t256 * t308;
t389 = -g(1) * t122 + g(2) * t126;
t127 = t192 * t308 + t256 * t307;
t388 = g(1) * t123 - g(2) * t127;
t387 = -g(1) * t187 + g(2) * t191;
t386 = g(1) * t257 + g(2) * t255;
t131 = pkin(3) * t462 - t404;
t382 = t305 * t322 + t320 * t496;
t379 = 0.2e1 * t407 + qJD(2);
t376 = -pkin(1) * t505 - t255 * pkin(2) + pkin(8) * t427;
t374 = -t234 * t442 - t235 * t444 + t322 * t241 - t320 * t243;
t371 = -t319 * t65 + t506 * t57;
t3 = t319 * t11 + t506 * t14 + t39 * t415 - t440 * t49;
t365 = t319 * t27 + t506 * t35 + t57 * t415 - t440 * t65;
t362 = t175 * t506 - t319 * t186;
t32 = qJD(5) * t363 + t166 * t440 - t319 * t406 - t506 * t326;
t355 = t411 + t436;
t151 = -t254 * t469 - t255 * t308;
t153 = -t256 * t469 - t257 * t308;
t202 = -t308 * t463 + t322 * t432;
t351 = -g(1) * t153 - g(2) * t151 - g(3) * t202;
t152 = -t254 * t468 + t255 * t307;
t154 = -t256 * t468 + t257 * t307;
t203 = (t307 * t321 + t308 * t459) * t316;
t350 = -g(1) * t154 - g(2) * t152 - g(3) * t203;
t348 = g(1) * t192 + g(2) * t188 + g(3) * t253;
t343 = t349 - t51;
t342 = g(1) * t256 + g(2) * t254 - t435;
t340 = t506 * t366;
t339 = t32 + t539;
t85 = -pkin(3) * t423 - t374;
t173 = t253 * t307 + t308 * t462;
t337 = g(1) * t126 + g(2) * t122 + g(3) * t173 + t4;
t336 = -g(1) * t471 - g(2) * t472 + t320 * t435;
t335 = t349 + t402;
t334 = t111 * t359 + t308 * t349;
t333 = -t19 * t315 + t20 * t317 - t348;
t108 = t183 * t506 - t319 * t366;
t79 = pkin(4) * t405 + t105;
t96 = pkin(4) * t366 + t131;
t174 = t253 * t308 - t432;
t332 = -g(1) * t127 - g(2) * t123 - g(3) * t174 + t3;
t31 = -pkin(4) * t406 + t51;
t28 = t90 * pkin(5) - qJ(6) * t534 + t79;
t331 = t28 * t534 + qJDD(6) - t337;
t66 = -pkin(4) * t346 + t85;
t33 = -qJD(5) * t383 + t166 * t415 + t319 * t326 - t506 * t406;
t7 = t33 * pkin(5) + t32 * qJ(6) - qJD(6) * t534 + t31;
t325 = t33 + t537;
t247 = t256 * pkin(2);
t245 = t254 * pkin(2);
t237 = t262 * t320;
t215 = -pkin(9) * t464 + t260;
t171 = -pkin(5) * t358 - qJ(6) * t262 - t305;
t132 = pkin(5) * t237 - qJ(6) * t238 + t263;
t107 = t183 * t319 + t340;
t102 = t322 * pkin(5) - t362;
t101 = -qJ(6) * t322 + t520;
t60 = qJD(5) * t108 + t319 * t345 - t346 * t506;
t59 = qJD(5) * t340 + t183 * t440 - t319 * t346 - t345 * t506;
t52 = pkin(5) * t534 + qJ(6) * t90;
t36 = t107 * pkin(5) - t108 * qJ(6) + t96;
t23 = -t252 * pkin(5) - t371;
t22 = qJ(6) * t252 + t367;
t21 = -t32 + t539;
t16 = t211 * qJ(6) + t18;
t15 = -t211 * pkin(5) + t458;
t8 = t60 * pkin(5) + t59 * qJ(6) - t108 * qJD(6) + t66;
t6 = -t184 * pkin(5) - t510;
t5 = qJ(6) * t184 + qJD(6) * t252 + t365;
t2 = qJDD(6) - t4 - t503;
t1 = qJD(6) * t211 + t3 + t481;
t9 = [qJDD(1), g(1) * t505 - g(2) * t507, g(1) * t507 + g(2) * t505 (qJDD(1) * t313 + 0.2e1 * t321 * t411) * t311 (t323 * t436 - t439 * t448) * t526, t316 * t379 * t445 + t321 * t518, t323 * t518 - t379 * t423, t375 * t484, -t244 * t380 + t373 * t375 - t401 * t484 + g(1) * t255 - g(2) * t257 + (-t413 + t438) * t414, -g(1) * t254 + g(2) * t256 - t243 * t380 - t341 * t484 - t355 * t414 - t375 * t449, t219 * t185 + t253 * t327, -t253 * t117 - t219 * t184 - t185 * t217 - t252 * t327, -t185 * t277 + t253 * t236 + (-t329 * t323 + t219 * t446 - (-t320 * t412 + t322 * t411 + t408) * t462) * t316, t184 * t277 - t236 * t252 + (t117 * t323 - t217 * t446) * t316 (-t236 * t323 - t277 * t446) * t316, -t374 * t277 + t404 * t236 + t244 * t217 + t233 * t117 + t147 * t252 + t199 * t184 + g(1) * t188 - g(2) * t192 + (t114 * t446 - t323 * t402) * t316, -t115 * t423 + t147 * t253 + t199 * t185 + t244 * t219 + t233 * t327 - t236 * t454 + t277 * t353 + t354 * t462 + t387, t40 * t217 + t74 * t117 + t19 * t252 + t63 * t184 + t85 * t405 - t131 * t406 + t51 * t366 - t105 * t346 - g(1) * (-t188 * t317 - t254 * t315) - g(2) * (t192 * t317 + t256 * t315) -t41 * t217 - t75 * t117 - t20 * t252 - t64 * t184 + t85 * t166 + t131 * t326 + t51 * t183 + t105 * t345 - g(1) * (t188 * t315 - t254 * t317) - g(2) * (-t192 * t315 + t256 * t317) -t40 * t166 - t19 * t183 - t20 * t366 - t326 * t74 - t345 * t63 + t346 * t64 - t405 * t41 + t406 * t75 - t387, t20 * t75 + t64 * t41 + t19 * t74 + t63 * t40 + t51 * t131 + t105 * t85 - g(1) * (-pkin(3) * t188 - t254 * pkin(9) - qJ(4) * t187 + t376) - g(2) * (pkin(3) * t192 + pkin(9) * t256 + qJ(4) * t191 + t431) -t108 * t32 - t534 * t59, t107 * t32 - t108 * t33 - t534 * t60 + t59 * t90, t108 * t111 + t184 * t534 - t211 * t59 - t252 * t32, -t107 * t111 - t184 * t90 - t211 * t60 - t252 * t33, t111 * t252 + t184 * t211, t31 * t107 + t371 * t111 + t17 * t184 + t211 * t510 + t4 * t252 + t96 * t33 + t79 * t60 + t66 * t90 + t388, t31 * t108 - t111 * t367 - t18 * t184 - t211 * t365 - t3 * t252 - t96 * t32 + t534 * t66 - t79 * t59 + t389, t107 * t7 - t111 * t23 - t15 * t184 - t2 * t252 - t211 * t6 + t28 * t60 + t33 * t36 + t8 * t90 + t388, -t1 * t107 + t108 * t2 - t15 * t59 - t16 * t60 - t22 * t33 - t23 * t32 - t5 * t90 + t534 * t6 - t387, t1 * t252 - t108 * t7 + t111 * t22 + t16 * t184 + t211 * t5 + t28 * t59 + t32 * t36 - t534 * t8 - t389, t1 * t22 + t16 * t5 + t7 * t36 + t28 * t8 + t2 * t23 + t15 * t6 - g(1) * (-pkin(5) * t123 - qJ(6) * t122 - t187 * t496 - t188 * t305 - t254 * t429 + t376) - g(2) * (pkin(5) * t127 + qJ(6) * t126 + t191 * t496 + t192 * t305 + t256 * t429 + t431); 0, 0, 0, -t321 * t433, t448 * t467, -t323 * t391 + t410, t321 * t391 + t291, t375, pkin(1) * t321 * t467 + t242 * t380 + t342 - t401, pkin(1) * t433 + t239 * t380 + (pkin(8) * t439 + g(3)) * t463 + t386 - t430 (-t316 * t412 + t375) * t320 ^ 2 + ((t316 * t355 + t347) * t320 - t476) * t322, -t320 * t117 + t217 * t516 - t219 * t514 + t322 * t327, -t277 * t442 + t320 * t236 + (-t219 * t321 + t277 * t459) * t447 (t217 * t321 - t323 * t356) * t447 - t515, t277 * t425, -t114 * t425 - pkin(2) * t117 - t242 * t217 - t220 * t277 + (-pkin(9) * t236 - t199 * t277) * t320 + (-t147 + (pkin(9) * qJD(3) + t240) * t277 + t342) * t322, -pkin(2) * t327 + pkin(9) * t515 + t115 * t425 + t147 * t320 - t199 * t516 - t242 * t219 - t277 * t453 + t336, t215 * t117 - t135 * t405 - t105 * t205 + t485 * t217 + (-g(3) * t463 - t386) * t315 + (-t19 + t342 * t317 + (pkin(9) * t405 + t105 * t315) * qJD(3)) * t322 + (-pkin(9) * t406 - t277 * t63 + t489) * t320, -t216 * t117 + t20 * t322 + t326 * t309 + t51 * t461 + t317 * t416 - t105 * t206 - g(1) * (t256 * t464 + t257 * t317) - g(2) * (t254 * t464 + t255 * t317) - (-t315 * t459 + t317 * t321) * t499 - t514 * t64 - t525 * t217 + t527 * t166, -t485 * t166 - t19 * t461 - t20 * t466 - t215 * t326 + t216 * t406 - t525 * t405 + t531 * t64 + t535 * t63 - t336, t20 * t216 + t19 * t215 - t105 * t135 - g(1) * (-qJ(4) * t471 - t256 * t504 - t247) - g(2) * (-qJ(4) * t472 - t254 * t504 - t245) - g(3) * (t385 * t462 + t450) + t525 * t64 + t485 * t63 + (t51 * t320 - t386 + t416) * pkin(9), -t238 * t32 + t457 * t534, t237 * t32 - t238 * t33 - t456 * t534 - t457 * t90, t111 * t238 + t211 * t457 + t32 * t322 - t356 * t534, -t111 * t237 - t211 * t456 + t322 * t33 + t356 * t90, -t111 * t322 - t211 * t356, t362 * t111 + t17 * t514 + t211 * t509 + t31 * t237 + t263 * t33 - t4 * t322 + t456 * t79 + t524 * t90 + t350, -t111 * t520 + t18 * t356 + t211 * t519 + t31 * t238 - t263 * t32 + t3 * t322 + t457 * t79 + t524 * t534 - t351, -t102 * t111 + t132 * t33 + t15 * t356 + t2 * t322 - t211 * t494 + t237 * t7 + t28 * t456 + t491 * t90 + t350, -t1 * t237 - t101 * t33 - t102 * t32 + t15 * t457 - t16 * t456 + t2 * t238 + t320 * t342 + t494 * t534 - t495 * t90, -t1 * t322 + t101 * t111 + t132 * t32 - t16 * t356 + t211 * t495 - t238 * t7 - t28 * t457 - t491 * t534 + t351, t1 * t101 + t7 * t132 + t2 * t102 - g(1) * (pkin(5) * t154 + qJ(6) * t153 - t256 * t382 + t257 * t429 - t247) - g(2) * (pkin(5) * t152 + qJ(6) * t151 - t254 * t382 + t255 * t429 - t245) - g(3) * (pkin(5) * t203 + qJ(6) * t202 + t450) - (pkin(4) * t465 + t323 * t382) * t499 + t491 * t28 + t495 * t16 + t494 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t219 * t217, -t217 ^ 2 + t219 ^ 2, -t217 * t277 + t327, -t117 - t476, t236, -t199 * t219 + t335 - t480, -t114 * t277 + t199 * t217 + t348 - t354, -t117 * t483 - pkin(3) * t110 - t115 * t204 - t63 * t219 + (-t76 + (-qJD(4) + t105) * t315) * t217 + (t343 - t480 + t498) * t317, -pkin(3) * t326 + t105 * t477 - t115 * t166 - t117 * t482 + t217 * t523 + t64 * t219 - t315 * t349 + t489, t326 * t483 + t406 * t482 - t63 * t477 - t64 * t478 + t333 + t523 * t405 + (qJD(4) * t315 + t76) * t166, -t105 * t115 - t63 * t76 - t64 * t77 + (-t315 * t63 + t317 * t64) * qJD(4) + t343 * pkin(3) + t333 * qJ(4), -t262 * t32 - t528 * t534, -t262 * t33 - t32 * t358 - t451 * t534 + t528 * t90, t111 * t262 - t211 * t528 - t219 * t534, t111 * t358 - t211 * t451 + t219 * t90, -t211 * t219, -t17 * t219 + t211 * t521 - t305 * t33 - t31 * t358 + t451 * t79 - t97 * t90 + t334, t18 * t219 + t211 * t522 + t31 * t262 + t305 * t32 - t528 * t79 - t534 * t97 + t512, t15 * t219 + t171 * t33 - t211 * t486 + t28 * t451 - t358 * t7 + t488 * t90 + t334, t1 * t358 - t15 * t528 - t16 * t451 - t197 * t33 + t2 * t262 + t32 * t359 + t486 * t534 - t487 * t90 - t348, -t16 * t219 + t171 * t32 + t211 * t487 - t262 * t7 + t28 * t528 - t488 * t534 - t512, t1 * t197 + t15 * t486 + t16 * t487 + t7 * t171 - t2 * t359 + t28 * t488 - t496 * t348 + t349 * (pkin(5) * t308 + qJ(6) * t307 + t305); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t166 * t217 - t406, -t217 * t405 + t326, -t166 ^ 2 - t405 ^ 2, t166 * t63 + t405 * t64 - t335 + t513, 0, 0, 0, 0, 0, t325, -t339, t325, -t508 - t540, t339, -t15 * t534 + t16 * t90 - t349 + t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t497, t508 - t540, t21, -t33 + t537, t111, -t534 * t79 + t337 + t490, t17 * t211 + t79 * t90 - t332, -t52 * t90 - t331 + t490 + 0.2e1 * t503, pkin(5) * t32 - qJ(6) * t33 + (t16 - t18) * t534 + (t15 - t458) * t90, 0.2e1 * t481 - t28 * t90 + t52 * t534 + (0.2e1 * qJD(6) - t17) * t211 + t332, t1 * qJ(6) - t2 * pkin(5) - t28 * t52 - t15 * t18 - g(1) * (-pkin(5) * t126 + qJ(6) * t127) - g(2) * (-pkin(5) * t122 + qJ(6) * t123) - g(3) * (-pkin(5) * t173 + qJ(6) * t174) + t458 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t497 - t111, t21, -t211 ^ 2 - t508, -t16 * t211 + t331 - t503;];
tau_reg  = t9;
