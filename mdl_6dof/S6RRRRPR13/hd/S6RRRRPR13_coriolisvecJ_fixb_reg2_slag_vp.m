% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RRRRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 00:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRRRPR13_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR13_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR13_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR13_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 00:04:04
% EndTime: 2019-03-10 00:04:41
% DurationCPUTime: 15.88s
% Computational Cost: add. (21522->790), mult. (54957->1031), div. (0->0), fcn. (42502->10), ass. (0->342)
t318 = sin(pkin(6));
t325 = cos(qJ(2));
t444 = qJD(1) * t325;
t301 = t318 * t444;
t321 = sin(qJ(3));
t394 = t321 * t301;
t441 = qJD(3) * t321;
t519 = t394 - t441;
t324 = cos(qJ(3));
t475 = cos(pkin(6));
t407 = t475 * qJD(1);
t364 = t407 + qJD(2);
t322 = sin(qJ(2));
t445 = qJD(1) * t322;
t422 = t318 * t445;
t395 = t321 * t422;
t225 = -t324 * t364 + t395;
t219 = qJD(4) + t225;
t425 = pkin(1) * t475;
t396 = t325 * t425;
t249 = -pkin(8) * t422 + qJD(1) * t396;
t360 = t318 * (pkin(2) * t322 - pkin(9) * t325);
t250 = qJD(1) * t360;
t170 = t324 * t249 + t321 * t250;
t159 = pkin(10) * t422 + t170;
t309 = t322 * t425;
t385 = pkin(3) * t321 - pkin(10) * t324;
t464 = t318 * t325;
t172 = (t309 + (pkin(8) + t385) * t464) * qJD(1);
t320 = sin(qJ(4));
t323 = cos(qJ(4));
t94 = t323 * t159 + t320 * t172;
t546 = -qJ(5) * t519 - t94;
t274 = t385 * qJD(3);
t287 = -pkin(3) * t324 - pkin(10) * t321 - pkin(2);
t457 = t323 * t324;
t311 = pkin(9) * t457;
t438 = qJD(4) * t320;
t527 = -qJD(4) * t311 + t320 * t159 - t287 * t438 + (-t172 + t274) * t323;
t198 = -pkin(2) * t364 - t249;
t227 = t321 * t364 + t324 * t422;
t125 = t225 * pkin(3) - t227 * pkin(10) + t198;
t448 = pkin(8) * t464 + t309;
t243 = t475 * pkin(9) + t448;
t199 = qJD(2) * pkin(9) + qJD(1) * t243;
t244 = (-pkin(2) * t325 - pkin(9) * t322 - pkin(1)) * t318;
t217 = qJD(1) * t244;
t144 = t324 * t199 + t321 * t217;
t378 = t301 - qJD(3);
t128 = -pkin(10) * t378 + t144;
t60 = t323 * t125 - t320 * t128;
t454 = qJD(5) - t60;
t319 = sin(qJ(6));
t504 = pkin(4) + pkin(5);
t177 = t323 * t227 - t320 * t378;
t540 = -t177 * pkin(11) + t454;
t40 = -t219 * t504 + t540;
t206 = t219 * qJ(5);
t175 = t227 * t320 + t323 * t378;
t61 = t125 * t320 + t128 * t323;
t50 = pkin(11) * t175 + t61;
t46 = t206 + t50;
t502 = cos(qJ(6));
t10 = -t319 * t46 + t40 * t502;
t432 = qJD(6) - t219;
t545 = t10 * t432;
t446 = qJD(1) * t318;
t456 = t323 * t325;
t212 = (t320 * t322 + t324 * t456) * t446;
t416 = t321 * t438;
t427 = -pkin(9) * t320 - pkin(4);
t544 = -t394 * t504 - (-pkin(5) + t427) * t441 + t527 + (qJD(3) * t457 - t212 - t416) * pkin(11);
t393 = t324 * t301;
t211 = t320 * t393 - t323 * t422;
t437 = qJD(4) * t323;
t449 = t320 * t274 + t287 * t437;
t459 = t321 * t323;
t543 = -pkin(11) * t211 + (-pkin(9) * qJD(3) + pkin(11) * qJD(4)) * t459 + (-qJD(5) + (-pkin(9) * qJD(4) + pkin(11) * qJD(3)) * t320) * t324 + t449 + t546;
t498 = pkin(11) * t225;
t503 = pkin(10) - pkin(11);
t143 = -t321 * t199 + t324 * t217;
t164 = pkin(3) * t227 + pkin(10) * t225;
t80 = t323 * t143 + t320 * t164;
t65 = t227 * qJ(5) + t80;
t542 = -t320 * t498 + t503 * t438 + t65;
t132 = t320 * t143;
t296 = t503 * t323;
t541 = qJD(4) * t296 - t132 - (-t164 + t498) * t323 + t504 * t227;
t110 = t319 * t175 + t177 * t502;
t357 = -t175 * t502 + t319 * t177;
t471 = t110 * t357;
t439 = qJD(3) * t324;
t539 = -t320 * t439 + t211;
t381 = t323 * t439 - t212;
t251 = qJD(2) * t360;
t238 = qJD(1) * t251;
t465 = t318 * t322;
t262 = -pkin(8) * t465 + t396;
t253 = t262 * qJD(2);
t239 = qJD(1) * t253;
t349 = t199 * t441 - t217 * t439 - t321 * t238 - t324 * t239;
t443 = qJD(2) * t322;
t421 = t318 * t443;
t389 = qJD(1) * t421;
t74 = pkin(10) * t389 - t349;
t431 = qJD(1) * qJD(2);
t408 = t325 * t431;
t388 = t318 * t408;
t365 = t321 * t388;
t181 = qJD(3) * t227 + t365;
t373 = pkin(8) * t388;
t383 = qJD(2) * t407;
t511 = qJD(3) * t395 - t324 * (qJD(3) * t364 + t388);
t98 = pkin(1) * t322 * t383 + t181 * pkin(3) + pkin(10) * t511 + t373;
t400 = t125 * t438 + t128 * t437 + t320 * t74 - t323 * t98;
t358 = t219 * t61 - t400;
t538 = qJD(4) * t378 + t511;
t537 = t110 ^ 2 - t357 ^ 2;
t11 = t319 * t40 + t46 * t502;
t102 = t227 * t438 - t320 * t389 + t323 * t538;
t7 = t102 * pkin(11) - t181 * t504 + t400;
t103 = t227 * t437 - t320 * t538 - t323 * t389;
t179 = t181 * qJ(5);
t205 = t219 * qJD(5);
t351 = -t125 * t437 + t128 * t438 - t320 * t98 - t323 * t74;
t14 = t179 + t205 - t351;
t9 = pkin(11) * t103 + t14;
t2 = -qJD(6) * t11 - t319 * t9 + t502 * t7;
t127 = pkin(3) * t378 - t143;
t341 = t177 * qJ(5) - t127;
t48 = -t175 * t504 + t341;
t536 = -t48 * t110 + t2;
t406 = -t319 * t102 - t502 * t103;
t37 = qJD(6) * t110 + t406;
t535 = t110 * t432 - t37;
t534 = t227 * t439 - t321 * t511;
t472 = t103 * t323;
t473 = t102 * t320;
t533 = t321 * (qJD(4) * (t175 * t320 - t177 * t323) - t472 + t473) - (t175 * t323 + t177 * t320) * t439 + t212 * t175 + t177 * t211;
t461 = t320 * t181;
t532 = t324 * t103 + t539 * t219 + t321 * (t175 * t378 - t219 * t437 - t461);
t402 = t177 * t219;
t403 = t175 * t219;
t531 = (t103 + t402) * t320 + (t102 + t403) * t323;
t315 = t318 ^ 2;
t529 = -0.2e1 * t315 * t431;
t501 = pkin(4) * t181;
t16 = t400 - t501;
t55 = t206 + t61;
t528 = -t219 * t55 + t16;
t526 = t378 * t321;
t523 = t320 * qJD(5) + t144;
t430 = pkin(9) * t441;
t522 = t170 + t430;
t259 = t321 * t465 - t324 * t475;
t420 = qJD(2) * t464;
t190 = -qJD(3) * t259 + t324 * t420;
t521 = -qJD(4) * t464 + t190;
t520 = t393 - t439;
t423 = t502 * t320;
t271 = -t319 * t323 + t423;
t410 = qJD(6) * t502;
t411 = qJD(4) * t502;
t435 = qJD(6) * t319;
t518 = -t319 * t437 + t323 * t435 + (-t410 + t411) * t320;
t356 = -t319 * t320 - t323 * t502;
t517 = -t356 * qJD(6) - t319 * t438 - t323 * t411;
t1 = t319 * t7 + t40 * t410 - t435 * t46 + t502 * t9;
t516 = t357 * t48 - t1;
t36 = t502 * t102 - t319 * t103 - t175 * t410 + t177 * t435;
t515 = t357 * t432 - t36;
t499 = pkin(10) * t181;
t59 = t175 * pkin(4) - t341;
t514 = t219 * t59 - t499;
t405 = t219 ^ 2;
t458 = t323 * t181;
t510 = -t227 * t175 + t320 * t405 - t458;
t474 = qJ(5) * t323;
t509 = t320 * t504 - t474;
t462 = t320 * qJ(5);
t508 = -t323 * t504 - t462;
t505 = t177 ^ 2;
t327 = qJD(1) ^ 2;
t500 = pkin(10) * t177;
t497 = t324 * pkin(9);
t310 = t320 * t497;
t314 = t324 * pkin(4);
t180 = t324 * pkin(5) + t310 + t314 + (-pkin(11) * t321 - t287) * t323;
t237 = t320 * t287 + t311;
t203 = -qJ(5) * t324 + t237;
t460 = t320 * t321;
t182 = pkin(11) * t460 + t203;
t119 = t319 * t180 + t182 * t502;
t496 = qJD(6) * t119 + t543 * t319 + t544 * t502;
t118 = t180 * t502 - t319 * t182;
t495 = -qJD(6) * t118 + t544 * t319 - t543 * t502;
t494 = t177 * t59;
t85 = -t199 * t439 - t217 * t441 + t324 * t238 - t321 * t239;
t337 = pkin(3) * t389 + t85;
t331 = -qJ(5) * t102 + qJD(5) * t177 + t337;
t27 = pkin(4) * t103 - t331;
t491 = t27 * t320;
t490 = t27 * t323;
t489 = t320 * t337;
t488 = t323 * t337;
t169 = -t321 * t249 + t324 * t250;
t352 = pkin(3) * t422 + t169;
t339 = t212 * qJ(5) + t352;
t348 = -pkin(9) - t509;
t436 = qJD(5) * t323;
t487 = t504 * t211 + (qJD(4) * t508 + t436) * t321 + t348 * t439 - t339;
t295 = t503 * t320;
t201 = t295 * t502 - t319 * t296;
t486 = qJD(6) * t201 + t541 * t319 - t502 * t542;
t202 = t319 * t295 + t296 * t502;
t485 = -qJD(6) * t202 + t319 * t542 + t541 * t502;
t440 = qJD(3) * t323;
t162 = (-t321 * t440 - t324 * t438) * pkin(9) + t449;
t484 = -t324 * qJD(5) + t162 + t546;
t483 = pkin(4) * t394 + t427 * t441 - t527;
t379 = pkin(4) * t320 - t474;
t362 = pkin(9) + t379;
t380 = pkin(4) * t323 + t462;
t482 = -pkin(4) * t211 + (qJD(4) * t380 - t436) * t321 + t362 * t439 + t339;
t481 = t162 - t94;
t480 = t320 * t430 + t527;
t280 = qJ(5) * t502 - t319 * t504;
t479 = qJD(6) * t280 + t319 * t540 + t50 * t502;
t279 = -t319 * qJ(5) - t502 * t504;
t478 = -qJD(6) * t279 + t319 * t50 - t502 * t540;
t477 = -t219 * t509 + t523;
t476 = t219 * t379 - t523;
t470 = t175 * qJ(5);
t469 = t177 * t175;
t157 = t181 * t259;
t174 = t181 * t324;
t468 = t219 * t227;
t467 = t227 * t225;
t466 = t315 * t327;
t392 = t502 * t439;
t453 = t212 * t502 + t319 * t539 + t321 * t518 - t323 * t392;
t452 = t211 * t502 + t319 * t381 - t320 * t392 + t321 * t517;
t242 = -pkin(2) * t475 - t262;
t260 = t321 * t475 + t324 * t465;
t153 = t259 * pkin(3) - t260 * pkin(10) + t242;
t166 = t324 * t243 + t321 * t244;
t155 = -pkin(10) * t464 + t166;
t77 = t320 * t153 + t323 * t155;
t451 = t225 * t271 + t518;
t450 = t356 * t225 + t517;
t165 = -t321 * t243 + t324 * t244;
t447 = t322 ^ 2 - t325 ^ 2;
t434 = t198 * qJD(3);
t429 = pkin(10) * t438;
t428 = pkin(10) * t437;
t68 = t259 * qJ(5) + t77;
t154 = pkin(3) * t464 - t165;
t412 = t175 ^ 2 - t505;
t79 = t164 * t323 - t132;
t76 = t323 * t153 - t320 * t155;
t236 = t287 * t323 - t310;
t404 = t325 * t378;
t401 = t378 * t318;
t399 = qJD(3) * t378;
t397 = t322 * t325 * t466;
t105 = -t243 * t439 - t244 * t441 + t324 * t251 - t321 * t253;
t387 = t318 * t327 * t475;
t386 = pkin(1) * t529;
t54 = -pkin(4) * t219 + t454;
t377 = -t320 * t55 + t323 * t54;
t376 = -t320 * t61 - t323 * t60;
t375 = -t85 * t321 - t324 * t349;
t374 = (qJD(4) * t175 - t102) * pkin(10);
t121 = t260 * t437 + t521 * t320 - t323 * t421;
t191 = t260 * t320 + t318 * t456;
t371 = t103 * t191 + t121 * t175;
t366 = t315 * t322 * t408;
t363 = 0.2e1 * t407 + qJD(2);
t189 = qJD(3) * t260 + t321 * t420;
t254 = t448 * qJD(2);
t113 = t189 * pkin(3) - t190 * pkin(10) + t254;
t104 = -t243 * t441 + t244 * t439 + t321 * t251 + t324 * t253;
t91 = pkin(10) * t421 + t104;
t32 = t323 * t113 - t153 * t438 - t155 * t437 - t320 * t91;
t192 = t260 * t323 - t320 * t464;
t361 = t192 * qJ(5) - t154;
t52 = -t192 * pkin(11) - t259 * t504 - t76;
t56 = pkin(11) * t191 + t68;
t23 = -t319 * t56 + t502 * t52;
t25 = t319 * t52 + t502 * t56;
t129 = -t191 * t502 + t192 * t319;
t130 = t191 * t319 + t192 * t502;
t350 = t127 * t219 - t499;
t31 = t320 * t113 + t153 * t437 - t155 * t438 + t323 * t91;
t347 = -t219 * t502 + t410;
t346 = t320 * t403 - t472;
t343 = pkin(3) * t421 + t105;
t22 = t189 * qJ(5) + t259 * qJD(5) + t31;
t340 = t219 * t60 + t351;
t338 = pkin(1) * (-t383 + t466);
t122 = -t260 * t438 + t320 * t421 + t521 * t323;
t336 = t102 * t191 - t103 * t192 - t121 * t177 - t122 * t175;
t335 = t103 * t259 + t121 * t219 + t175 * t189 + t181 * t191;
t333 = qJ(5) * t122 + qJD(5) * t192 + t343;
t332 = t103 * t460 + (t321 * t437 - t539) * t175;
t328 = t103 - t402;
t282 = -pkin(3) - t380;
t266 = pkin(3) - t508;
t256 = t356 * t321;
t255 = t319 * t459 - t321 * t423;
t252 = t448 * qJD(1);
t245 = t362 * t321;
t240 = qJD(1) * t254;
t204 = -t236 + t314;
t200 = t348 * t321;
t112 = pkin(4) * t177 + t470;
t97 = pkin(10) * t472;
t90 = -t219 * t526 - t174;
t86 = t189 * t219 + t157;
t78 = pkin(4) * t191 - t361;
t71 = -t177 * t504 - t470;
t69 = -pkin(4) * t259 - t76;
t67 = -pkin(4) * t227 - t79;
t63 = -t191 * t504 + t361;
t58 = -t102 + t403;
t53 = -t177 * t227 + t323 * t405 + t461;
t47 = t323 * t402 - t473;
t45 = -qJD(6) * t129 + t121 * t319 + t122 * t502;
t44 = qJD(6) * t130 - t121 * t502 + t122 * t319;
t41 = -t102 * t192 + t122 * t177;
t38 = -t102 * t459 + (t381 - t416) * t177;
t35 = pkin(4) * t121 - t333;
t30 = -pkin(4) * t189 - t32;
t29 = t102 * t324 + t381 * t219 + (-t177 * t378 - t219 * t438 + t458) * t321;
t28 = -t102 * t259 + t122 * t219 + t177 * t189 + t181 * t192;
t21 = -t121 * t504 + t333;
t15 = pkin(11) * t121 + t22;
t13 = -t122 * pkin(11) - t189 * t504 - t32;
t12 = -t103 * t504 + t331;
t4 = -qJD(6) * t25 + t13 * t502 - t319 * t15;
t3 = qJD(6) * t23 + t319 * t13 + t15 * t502;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t366, t447 * t529, t363 * t420, -0.2e1 * t366, -t363 * t421, 0, -t240 * t475 - t254 * t364 + t322 * t386, -t239 * t475 - t253 * t364 + t325 * t386 (t239 * t325 + t240 * t322 + (-t249 * t325 - t252 * t322) * qJD(2) + (t253 * t325 + t254 * t322 + (-t262 * t325 - t322 * t448) * qJD(2)) * qJD(1)) * t318, t239 * t448 - t240 * t262 - t249 * t254 + t252 * t253, t227 * t190 - t260 * t511, -t260 * t181 - t227 * t189 - t190 * t225 + t259 * t511, -t190 * t378 + t227 * t421 + t260 * t389 + t464 * t511, t189 * t225 + t157, t189 * t378 + (t181 * t325 + (-qJD(1) * t259 - t225) * t443) * t318 (-t315 * t444 - t401) * t443, -t105 * t378 + t254 * t225 + t242 * t181 + t240 * t259 + t198 * t189 + (-t85 * t325 + (qJD(1) * t165 + t143) * t443) * t318, t104 * t378 - t144 * t421 - t166 * t389 + t198 * t190 + t254 * t227 + t240 * t260 - t242 * t511 - t349 * t464, -t104 * t225 - t105 * t227 - t143 * t190 - t144 * t189 + t165 * t511 - t166 * t181 + t259 * t349 - t85 * t260, t104 * t144 + t105 * t143 + t165 * t85 - t166 * t349 + t198 * t254 + t240 * t242, t41, t336, t28, t371, -t335, t86, t103 * t154 + t121 * t127 - t175 * t343 + t181 * t76 + t189 * t60 - t191 * t337 + t219 * t32 - t259 * t400, -t102 * t154 + t122 * t127 - t177 * t343 - t181 * t77 - t189 * t61 - t192 * t337 - t219 * t31 + t259 * t351, t102 * t76 - t103 * t77 - t121 * t61 - t122 * t60 - t175 * t31 - t177 * t32 + t191 * t351 + t192 * t400, -t127 * t343 - t154 * t337 + t31 * t61 + t32 * t60 - t351 * t77 - t400 * t76, t41, t28, -t336, t86, t335, t371, t103 * t78 + t121 * t59 - t16 * t259 + t175 * t35 - t181 * t69 - t189 * t54 + t191 * t27 - t219 * t30, -t102 * t69 - t103 * t68 - t121 * t55 + t122 * t54 - t14 * t191 + t16 * t192 - t175 * t22 + t177 * t30, t102 * t78 - t122 * t59 + t14 * t259 - t177 * t35 + t181 * t68 + t189 * t55 - t192 * t27 + t219 * t22, t14 * t68 + t16 * t69 + t22 * t55 + t27 * t78 + t30 * t54 + t35 * t59, t110 * t45 - t130 * t36, -t110 * t44 + t129 * t36 - t130 * t37 - t357 * t45, -t110 * t189 - t130 * t181 + t259 * t36 + t432 * t45, t129 * t37 + t357 * t44, t129 * t181 + t189 * t357 + t259 * t37 - t432 * t44, -t189 * t432 + t157, -t10 * t189 + t12 * t129 - t181 * t23 - t2 * t259 + t21 * t357 + t37 * t63 + t4 * t432 + t44 * t48, t1 * t259 + t11 * t189 + t110 * t21 + t12 * t130 + t181 * t25 - t3 * t432 - t36 * t63 + t45 * t48, -t1 * t129 - t10 * t45 - t11 * t44 - t110 * t4 - t130 * t2 + t23 * t36 - t25 * t37 - t3 * t357, t1 * t25 + t10 * t4 + t11 * t3 + t12 * t63 + t2 * t23 + t21 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t397, t447 * t466, -t325 * t387, t397, t322 * t387, 0, t252 * t364 + t322 * t338 - t373, pkin(8) * t389 + t249 * t364 + t325 * t338, 0, 0, -t227 * t393 + t534, -t321 * t181 + t225 * t520 + t519 * t227 - t324 * t511, -t324 * t399 + (t324 * t404 + (qJD(2) * t321 - t227) * t322) * t446, -t225 * t526 - t174, t321 * t399 + (-t321 * t404 + (qJD(2) * t324 + t225) * t322) * t446, t401 * t445, -pkin(2) * t181 + t321 * t434 + t169 * t378 - t252 * t225 + (pkin(9) * t399 - t240) * t324 + (-t143 * t322 + (-pkin(9) * t443 - t198 * t325) * t321) * t446, pkin(2) * t511 + t144 * t422 - t198 * t393 - t252 * t227 + t240 * t321 + t324 * t434 - t522 * t378 - t389 * t497, t169 * t227 + t375 + t522 * t225 + t519 * t144 + t520 * t143 + (-t174 + t534) * pkin(9), -pkin(2) * t240 - t143 * t169 - t144 * t170 - t198 * t252 + ((-t143 * t324 - t144 * t321) * qJD(3) + t375) * pkin(9), t38, t533, t29, t332, t532, t90, -t127 * t211 + t352 * t175 + t181 * t236 + t480 * t219 + (t400 + (pkin(9) * t175 + t127 * t320) * qJD(3)) * t324 + (pkin(9) * t103 + t127 * t437 - t378 * t60 - t489) * t321, -t127 * t212 + t352 * t177 - t181 * t237 - t481 * t219 + (-t351 + (pkin(9) * t177 + t127 * t323) * qJD(3)) * t324 + (-pkin(9) * t102 - t127 * t438 + t378 * t61 - t488) * t321, t102 * t236 - t103 * t237 + t211 * t61 + t212 * t60 - t480 * t177 - t481 * t175 + t376 * t439 + (t351 * t320 + t400 * t323 + (t320 * t60 - t323 * t61) * qJD(4)) * t321, t127 * t352 - t351 * t237 - t400 * t236 + t481 * t61 + t480 * t60 + (t127 * t439 - t321 * t337) * pkin(9), t38, t29, -t533, t90, -t532, t332, t103 * t245 - t181 * t204 - t211 * t59 + (qJD(3) * t320 * t59 + t16) * t324 - t483 * t219 + t482 * t175 + (t378 * t54 + t437 * t59 + t491) * t321, -t102 * t204 - t103 * t203 + t211 * t55 - t212 * t54 + t483 * t177 - t484 * t175 + t377 * t439 + (-t14 * t320 + t16 * t323 + (-t320 * t54 - t323 * t55) * qJD(4)) * t321, t102 * t245 + t181 * t203 + t212 * t59 + (-t440 * t59 - t14) * t324 + t484 * t219 - t482 * t177 + (-t378 * t55 + t438 * t59 - t490) * t321, t14 * t203 + t16 * t204 + t245 * t27 + t482 * t59 + t483 * t54 + t484 * t55, -t110 * t453 + t36 * t256, -t110 * t452 + t36 * t255 + t256 * t37 + t357 * t453, t110 * t526 + t256 * t181 - t36 * t324 - t432 * t453, t37 * t255 + t357 * t452, t255 * t181 - t37 * t324 - t357 * t526 - t432 * t452, t432 * t526 - t174, t10 * t526 - t118 * t181 + t12 * t255 + t2 * t324 + t200 * t37 + t357 * t487 - t432 * t496 + t452 * t48, -t1 * t324 - t11 * t526 + t110 * t487 + t119 * t181 - t12 * t256 - t200 * t36 + t432 * t495 - t453 * t48, -t1 * t255 + t10 * t453 - t11 * t452 + t110 * t496 + t118 * t36 - t119 * t37 + t2 * t256 + t357 * t495, t1 * t119 - t10 * t496 - t11 * t495 + t118 * t2 + t12 * t200 + t48 * t487; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t467, -t225 ^ 2 + t227 ^ 2, -t225 * t378 - t511, -t467, -t227 * t301 - t365, t389, -t144 * t378 - t198 * t227 + t85, -t143 * t378 + t198 * t225 + t349, 0, 0, t47, -t531, t53, t346, -t510, -t468, -pkin(3) * t103 - t144 * t175 - t227 * t60 + t488 + (-t79 - t428) * t219 + t350 * t320, pkin(3) * t102 - t144 * t177 + t227 * t61 - t489 + (t80 + t429) * t219 + t350 * t323, t175 * t80 + t177 * t79 - t97 + (-t225 * t60 - t351 + (-t60 + t500) * qJD(4)) * t323 + (t374 - t358) * t320, pkin(3) * t337 - t127 * t144 - t60 * t79 - t61 * t80 + (qJD(4) * t376 + t320 * t400 - t323 * t351) * pkin(10), t47, t53, t531, -t468, t510, t346, t103 * t282 + t227 * t54 - t490 + (t67 - t428) * t219 + t476 * t175 + t514 * t320, t175 * t65 - t177 * t67 - t97 + (t225 * t54 + t14 + (t54 + t500) * qJD(4)) * t323 + (t374 + t528) * t320, t102 * t282 - t227 * t55 - t491 + (-t65 - t429) * t219 - t476 * t177 - t514 * t323, t27 * t282 - t54 * t67 - t55 * t65 + t476 * t59 + (qJD(4) * t377 + t14 * t323 + t16 * t320) * pkin(10), -t110 * t450 - t36 * t271, t110 * t451 - t271 * t37 - t356 * t36 + t357 * t450, t110 * t227 - t271 * t181 - t432 * t450, -t356 * t37 - t357 * t451, -t181 * t356 - t227 * t357 + t432 * t451, t432 * t227, t10 * t227 - t12 * t356 - t181 * t201 + t266 * t37 + t357 * t477 + t432 * t485 - t451 * t48, -t11 * t227 + t110 * t477 + t12 * t271 + t181 * t202 - t266 * t36 - t432 * t486 - t450 * t48, t1 * t356 + t10 * t450 + t11 * t451 - t110 * t485 - t2 * t271 + t201 * t36 - t202 * t37 - t357 * t486, t1 * t202 + t10 * t485 + t11 * t486 + t12 * t266 + t2 * t201 + t477 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t469, -t412, t58, -t469, -t328, t181, -t127 * t177 + t358, t127 * t175 + t340, 0, 0, t469, t58, t412, t181, t328, -t469, -t112 * t175 + t358 - t494 + 0.2e1 * t501, pkin(4) * t102 - t103 * qJ(5) + (t55 - t61) * t177 + (t54 - t454) * t175, t112 * t177 - t175 * t59 + 0.2e1 * t179 + 0.2e1 * t205 - t340, -t16 * pkin(4) + t14 * qJ(5) - t59 * t112 + t454 * t55 - t54 * t61, -t471, -t537, -t515, t471, -t535, t181, -t279 * t181 - t357 * t71 - t432 * t479 - t536, -t110 * t71 + t181 * t280 + t432 * t478 - t516, t279 * t36 - t280 * t37 + (t10 + t478) * t357 + (-t11 + t479) * t110, t1 * t280 - t10 * t479 - t11 * t478 + t2 * t279 - t48 * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t469 - t181, t58, -t405 - t505, t494 + t528, 0, 0, 0, 0, 0, 0, -t319 * t432 ^ 2 - t177 * t357 - t181 * t502, -t177 * t110 + t319 * t181 - t347 * t432, t319 * t535 - t347 * t357 + t502 * t36, t2 * t502 - t48 * t177 + t347 * t11 + (t1 - t545) * t319; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t471, t537, t515, -t471, -t406 + (-qJD(6) + t432) * t110, -t181, t11 * t432 + t536, t516 + t545, 0, 0;];
tauc_reg  = t5;
