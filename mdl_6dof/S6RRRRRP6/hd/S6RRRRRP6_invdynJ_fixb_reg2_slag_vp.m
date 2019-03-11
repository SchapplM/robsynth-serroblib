% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RRRRRP6
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRRRRP6_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP6_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP6_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP6_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP6_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP6_invdynJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:31:35
% EndTime: 2019-03-10 01:32:05
% DurationCPUTime: 15.58s
% Computational Cost: add. (25399->859), mult. (56682->1057), div. (0->0), fcn. (40717->14), ass. (0->389)
t336 = sin(qJ(4));
t337 = sin(qJ(3));
t340 = cos(qJ(4));
t341 = cos(qJ(3));
t262 = t336 * t341 + t337 * t340;
t342 = cos(qJ(2));
t472 = qJD(1) * t342;
t207 = t262 * t472;
t391 = t262 * qJD(4);
t392 = t262 * qJD(3);
t354 = -t392 - t391;
t595 = t207 + t354;
t410 = t336 * t337 - t340 * t341;
t208 = t410 * t472;
t390 = t410 * qJD(3);
t352 = -t410 * qJD(4) - t390;
t569 = t208 + t352;
t334 = qJ(3) + qJ(4);
t326 = qJ(5) + t334;
t312 = cos(t326);
t311 = sin(t326);
t343 = cos(qJ(1));
t483 = t343 * t311;
t339 = sin(qJ(1));
t486 = t339 * t342;
t201 = t312 * t486 - t483;
t484 = t342 * t343;
t203 = t311 * t339 + t312 * t484;
t322 = t342 * qJDD(1);
t338 = sin(qJ(2));
t461 = qJD(1) * qJD(2);
t566 = -t338 * t461 + t322;
t256 = qJDD(3) - t566;
t248 = qJDD(4) + t256;
t431 = pkin(2) * t338 - pkin(8) * t342;
t272 = t431 * qJD(2);
t432 = pkin(2) * t342 + pkin(8) * t338;
t279 = -pkin(1) - t432;
t188 = qJD(1) * t272 + qJDD(1) * t279;
t179 = t341 * t188;
t249 = t279 * qJD(1);
t318 = pkin(7) * t472;
t285 = qJD(2) * pkin(8) + t318;
t185 = t249 * t337 + t285 * t341;
t229 = pkin(7) * t566 + qJDD(2) * pkin(8);
t101 = -qJD(3) * t185 - t229 * t337 + t179;
t460 = t338 * qJDD(1);
t374 = qJD(2) * (qJD(3) + t472) + t460;
t467 = qJD(3) * t338;
t581 = qJD(1) * t467 - qJDD(2);
t158 = t337 * t581 - t341 * t374;
t72 = pkin(3) * t256 + pkin(9) * t158 + t101;
t466 = qJD(3) * t341;
t468 = qJD(3) * t337;
t100 = t337 * t188 + t341 * t229 + t249 * t466 - t285 * t468;
t159 = t337 * t374 + t581 * t341;
t77 = -pkin(9) * t159 + t100;
t184 = t341 * t249 - t285 * t337;
t473 = qJD(1) * t338;
t450 = t341 * t473;
t471 = qJD(2) * t337;
t261 = t450 + t471;
t142 = -pkin(9) * t261 + t184;
t304 = -qJD(3) + t472;
t133 = -pkin(3) * t304 + t142;
t451 = t337 * t473;
t462 = t341 * qJD(2);
t259 = t451 - t462;
t143 = -pkin(9) * t259 + t185;
t513 = t143 * t340;
t85 = t133 * t336 + t513;
t25 = -qJD(4) * t85 - t336 * t77 + t340 * t72;
t464 = qJD(4) * t340;
t465 = qJD(4) * t336;
t82 = t340 * t158 + t336 * t159 + t259 * t464 + t261 * t465;
t15 = pkin(4) * t248 + pkin(10) * t82 + t25;
t24 = t133 * t464 - t143 * t465 + t336 * t72 + t340 * t77;
t408 = -t336 * t158 + t159 * t340 - t259 * t465 + t261 * t464;
t18 = -pkin(10) * t408 + t24;
t335 = sin(qJ(5));
t550 = cos(qJ(5));
t445 = qJD(5) * t550;
t463 = qJD(5) * t335;
t295 = -qJD(4) + t304;
t416 = -t259 * t336 + t340 * t261;
t576 = pkin(10) * t416;
t136 = t336 * t143;
t84 = t340 * t133 - t136;
t70 = t84 - t576;
t63 = -pkin(4) * t295 + t70;
t417 = -t259 * t340 - t336 * t261;
t575 = pkin(10) * t417;
t71 = t85 + t575;
t3 = t335 * t15 + t550 * t18 + t63 * t445 - t71 * t463;
t498 = t312 * t338;
t372 = g(1) * t203 + g(2) * t201 + g(3) * t498 - t3;
t233 = qJDD(5) + t248;
t216 = t233 * qJ(6);
t282 = -qJD(5) + t295;
t274 = t282 * qJD(6);
t568 = t216 - t274;
t110 = t335 * t416 - t417 * t550;
t284 = -qJD(2) * pkin(2) + pkin(7) * t473;
t195 = pkin(3) * t259 + t284;
t130 = -pkin(4) * t417 + t195;
t554 = t335 * t417 + t416 * t550;
t51 = pkin(5) * t110 - qJ(6) * t554 + t130;
t590 = t110 * t51;
t599 = -t372 + t568 - t590;
t316 = pkin(7) * t460;
t443 = t342 * t461;
t230 = -qJDD(2) * pkin(2) + pkin(7) * t443 + t316;
t428 = g(1) * t343 + g(2) * t339;
t539 = g(3) * t342;
t377 = t338 * t428 - t539;
t368 = -t230 + t377;
t598 = pkin(8) * qJD(3) * t304 + t368;
t587 = t110 * t130;
t597 = t372 + t587;
t28 = qJD(5) * t554 - t335 * t82 + t550 * t408;
t588 = t282 * t554 + t28;
t596 = t110 * t554;
t344 = -pkin(9) - pkin(8);
t287 = t344 * t337;
t288 = t344 * t341;
t191 = t336 * t287 - t340 * t288;
t594 = t410 * pkin(10) - t191;
t551 = t554 ^ 2;
t446 = t110 ^ 2 - t551;
t27 = t335 * t408 + t416 * t463 - t417 * t445 + t550 * t82;
t593 = -t110 * t282 - t27;
t530 = t335 * t71;
t31 = t550 * t63 - t530;
t457 = t550 * t71;
t32 = t335 * t63 + t457;
t592 = -t110 * t31 + t32 * t554;
t64 = pkin(5) * t554 + qJ(6) * t110;
t495 = t336 * t288;
t190 = t340 * t287 + t495;
t547 = pkin(10) * t262;
t152 = t190 - t547;
t103 = t335 * t152 - t550 * t594;
t591 = -t103 * t233 - t311 * t377;
t181 = t262 * t550 - t335 * t410;
t528 = -qJD(5) * t181 - t569 * t335 + t550 * t595;
t382 = t550 * t410;
t527 = qJD(5) * t382 + t262 * t463 - t335 * t595 - t569 * t550;
t585 = t130 * t554;
t549 = pkin(3) * t337;
t251 = t472 * t549 + t318;
t583 = pkin(3) * t468 - t251;
t469 = qJD(2) * t342;
t449 = t337 * t469;
t385 = t338 * t466 + t449;
t218 = t233 * pkin(5);
t563 = t218 - qJDD(6);
t579 = t51 * t554 - t563;
t578 = pkin(4) * t416;
t453 = qJD(3) * t344;
t270 = t337 * t453;
t271 = t341 * t453;
t414 = -t336 * t270 + t340 * t271;
t567 = qJD(4) * t594 + t414;
t269 = t431 * qJD(1);
t192 = pkin(7) * t451 + t341 * t269;
t485 = t341 * t342;
t409 = pkin(3) * t338 - pkin(9) * t485;
t154 = qJD(1) * t409 + t192;
t242 = t337 * t269;
t490 = t338 * t341;
t493 = t337 * t342;
t182 = t242 + (-pkin(7) * t490 - pkin(9) * t493) * qJD(1);
t106 = t340 * t154 - t182 * t336;
t86 = pkin(4) * t473 + pkin(10) * t208 + t106;
t107 = t336 * t154 + t340 * t182;
t90 = -pkin(10) * t207 + t107;
t454 = t340 * t270 + t336 * t271 + t287 * t464;
t96 = -pkin(10) * t392 + (t495 - t547) * qJD(4) + t454;
t533 = t103 * qJD(5) + (-pkin(10) * t390 - t567 + t86) * t550 + (t96 - t90) * t335;
t314 = pkin(3) * t340 + pkin(4);
t91 = -t142 * t336 - t513;
t384 = t91 - t575;
t459 = t335 * t336 * pkin(3);
t92 = t340 * t142 - t136;
t75 = t92 - t576;
t525 = -t335 * t384 + t314 * t445 + (-qJD(4) - qJD(5)) * t459 + (pkin(3) * t464 - t75) * t550;
t573 = t262 * t408;
t510 = t417 * t416;
t570 = t184 * t304 + t100;
t479 = -pkin(4) * t595 + t583;
t258 = t341 * t279;
t548 = pkin(7) * t337;
t183 = -pkin(9) * t490 + t258 + (-pkin(3) - t548) * t342;
t307 = pkin(7) * t485;
t206 = t337 * t279 + t307;
t494 = t337 * t338;
t189 = -pkin(9) * t494 + t206;
t120 = t336 * t183 + t340 * t189;
t564 = qJD(3) + qJD(4);
t351 = t564 * t410;
t393 = qJD(2) * t262;
t381 = t342 * t393;
t348 = t351 * t338 - t381;
t237 = t337 * t486 + t341 * t343;
t239 = -t337 * t484 + t339 * t341;
t565 = -g(1) * t239 + g(2) * t237;
t562 = t416 ^ 2 - t417 ^ 2;
t561 = t295 * t417 - t82;
t200 = t311 * t486 + t312 * t343;
t487 = t339 * t312;
t202 = t342 * t483 - t487;
t4 = t550 * t15 - t335 * t18 - t71 * t445 - t63 * t463;
t499 = t311 * t338;
t371 = g(1) * t202 + g(2) * t200 + g(3) * t499 + t4;
t359 = -t371 + t579;
t560 = t371 - t585;
t324 = sin(t334);
t325 = cos(t334);
t496 = t325 * t343;
t212 = t324 * t486 + t496;
t497 = t325 * t339;
t214 = -t324 * t484 + t497;
t540 = g(3) * t338;
t559 = -g(1) * t214 + g(2) * t212 - t416 * t195 + t324 * t540 + t25;
t213 = t324 * t343 - t325 * t486;
t215 = t324 * t339 + t325 * t484;
t558 = g(1) * t215 - g(2) * t213 - t417 * t195 + t325 * t540 - t24;
t180 = t262 * t335 + t382;
t557 = -t110 * t473 + t180 * t233 + t282 * t528;
t556 = -t295 * t416 - t408;
t555 = t342 * t428 + t540;
t126 = -t338 * t354 + t410 * t469;
t470 = qJD(2) * t338;
t477 = t341 * t272 + t470 * t548;
t118 = t409 * qJD(2) + (-t307 + (pkin(9) * t338 - t279) * t337) * qJD(3) + t477;
t140 = t337 * t272 + t279 * t466 + (-t338 * t462 - t342 * t468) * pkin(7);
t125 = -pkin(9) * t385 + t140;
t50 = -qJD(4) * t120 + t340 * t118 - t125 * t336;
t40 = pkin(4) * t470 + pkin(10) * t126 + t50;
t386 = -t337 * t467 + t342 * t462;
t455 = t336 * t118 + t340 * t125 + t183 * t464;
t44 = (-t490 * t564 - t449) * pkin(10) * t340 + (-qJD(4) * t189 + (qJD(4) * t494 - t386) * pkin(10)) * t336 + t455;
t394 = t338 * t262;
t104 = -pkin(10) * t394 + t120;
t119 = t340 * t183 - t189 * t336;
t220 = t410 * t338;
t97 = -pkin(4) * t342 + pkin(10) * t220 + t119;
t532 = t550 * t104 + t335 * t97;
t9 = -qJD(5) * t532 - t335 * t44 + t40 * t550;
t553 = t110 * t527 + t180 * t27 - t181 * t28 + t528 * t554;
t552 = -0.2e1 * pkin(1);
t545 = g(1) * t339;
t541 = g(2) * t343;
t538 = t341 * pkin(3);
t537 = -pkin(5) * t528 + qJ(6) * t527 - t181 * qJD(6) + t479;
t402 = t152 * t550 + t335 * t594;
t38 = t402 * qJD(5) + t550 * t96 + (-(-t336 * t468 + t340 * t466) * pkin(10) + t567) * t335;
t48 = t335 * t86 + t550 * t90;
t45 = qJ(6) * t473 + t48;
t536 = t38 - t45;
t535 = t38 - t48;
t534 = pkin(5) * t473 + t533;
t531 = t282 * t32;
t529 = t82 * t262;
t526 = qJD(6) + t525;
t452 = t550 * t336;
t524 = -t335 * t75 + t384 * t550 + t314 * t463 + (t336 * t445 + (t335 * t340 + t452) * qJD(4)) * pkin(3);
t35 = t550 * t70 - t530;
t523 = pkin(4) * t445 + qJD(6) - t35;
t522 = pkin(7) * qJDD(1);
t132 = t159 * pkin(3) + t230;
t514 = t132 * t262;
t512 = t158 * t337;
t511 = t159 * t341;
t508 = t185 * t304;
t507 = t259 * t304;
t506 = t259 * t337;
t505 = t261 * t259;
t504 = t261 * t304;
t503 = t261 * t341;
t502 = t262 * t248;
t492 = t337 * t343;
t491 = t338 * t339;
t489 = t338 * t343;
t488 = t338 * t344;
t276 = pkin(4) * t325 + t538;
t268 = pkin(2) + t276;
t247 = t342 * t268;
t275 = pkin(4) * t324 + t549;
t254 = t343 * t275;
t482 = qJD(6) - t31;
t481 = t288 * t465 - t107 + t454;
t480 = -qJD(4) * t191 - t106 + t414;
t478 = -t342 * t254 + t339 * t276;
t232 = pkin(3) * t452 + t335 * t314;
t273 = pkin(3) * t494 + t338 * pkin(7);
t476 = t343 * pkin(1) + t339 * pkin(7);
t332 = t338 ^ 2;
t333 = t342 ^ 2;
t475 = t332 - t333;
t474 = t332 + t333;
t346 = qJD(1) ^ 2;
t456 = t338 * t346 * t342;
t198 = pkin(3) * t385 + pkin(7) * t469;
t315 = pkin(2) + t538;
t447 = t282 * t473;
t440 = t524 * t554;
t439 = -t200 * pkin(5) + qJ(6) * t201;
t438 = -t202 * pkin(5) + qJ(6) * t203;
t437 = -t275 * t486 - t276 * t343;
t435 = t338 * t443;
t34 = t335 * t70 + t457;
t434 = pkin(4) * t463 - t34;
t308 = g(1) * t491;
t433 = -g(2) * t489 + t308;
t135 = pkin(3) * t261 + t578;
t430 = g(1) * t200 - g(2) * t202;
t429 = g(1) * t201 - g(2) * t203;
t427 = pkin(5) * t312 + qJ(6) * t311;
t426 = pkin(7) * t259 + t284 * t337;
t425 = pkin(7) * t261 + t284 * t341;
t378 = t550 * t394;
t144 = -t220 * t335 + t378;
t145 = -t220 * t550 - t335 * t394;
t59 = qJD(5) * t145 - t335 * t126 - t348 * t550;
t424 = t110 * t59 + t144 * t28;
t422 = -pkin(8) * t256 + qJD(3) * t284;
t418 = -t184 * t341 - t185 * t337;
t412 = t315 * t342 - t488;
t56 = -t335 * t104 + t550 * t97;
t403 = -pkin(7) * qJDD(2) + t461 * t552;
t399 = t256 * t337 - t304 * t466;
t398 = t256 * t341 + t304 * t468;
t8 = -t104 * t463 + t335 * t40 + t550 * t44 + t97 * t445;
t396 = -t110 * t528 + t180 * t28;
t395 = t195 * t262;
t388 = pkin(1) * t346 + t428;
t231 = t314 * t550 - t459;
t345 = qJD(2) ^ 2;
t387 = pkin(7) * t345 + qJDD(1) * t552 + t541;
t331 = -pkin(10) + t344;
t383 = t268 * t484 + t339 * t275 - t331 * t489 + t476;
t380 = g(2) * t338 * t487 + t402 * t233 + (g(1) * t489 - t539) * t312;
t329 = t343 * pkin(7);
t379 = t254 + t331 * t491 + t329 + (-pkin(1) - t247) * t339;
t58 = qJD(5) * t378 + t126 * t550 - t220 * t463 - t335 * t348;
t369 = t110 * t58 + t144 * t27 - t145 * t28 - t554 * t59;
t209 = pkin(4) * t410 - t315;
t186 = pkin(4) * t394 + t273;
t366 = t110 * t470 + t144 * t233 - t28 * t342 - t282 * t59;
t364 = -t282 * t31 + t372;
t278 = qJ(6) * t498;
t363 = -g(2) * t439 - g(3) * (-pkin(5) * t499 + t278);
t361 = t282 * t524 + t371;
t360 = -t103 * t28 - t38 * t110 + t27 * t402 - t555;
t60 = pkin(4) * t408 + t132;
t356 = g(2) * t496 + t324 * t555;
t105 = -pkin(4) * t348 + t198;
t313 = -pkin(4) * t550 - pkin(5);
t309 = pkin(4) * t335 + qJ(6);
t292 = pkin(4) * t497;
t240 = t337 * t339 + t341 * t484;
t238 = -t339 * t485 + t492;
t222 = -pkin(5) - t231;
t219 = qJ(6) + t232;
t205 = -pkin(7) * t493 + t258;
t193 = -pkin(7) * t450 + t242;
t170 = -t233 * t342 - t282 * t470;
t141 = -qJD(3) * t206 + t477;
t99 = t180 * pkin(5) - t181 * qJ(6) + t209;
t78 = t144 * pkin(5) - t145 * qJ(6) + t186;
t55 = t578 + t64;
t54 = t342 * pkin(5) - t56;
t53 = -qJ(6) * t342 + t532;
t52 = t135 + t64;
t49 = -t189 * t465 + t455;
t33 = t181 * t233 + t282 * t527 - t473 * t554;
t30 = -t282 * qJ(6) + t32;
t29 = t282 * pkin(5) + t482;
t19 = t59 * pkin(5) + t58 * qJ(6) - t145 * qJD(6) + t105;
t12 = -t145 * t27 - t554 * t58;
t11 = t145 * t233 + t27 * t342 + t282 * t58 + t470 * t554;
t10 = -t181 * t27 - t527 * t554;
t7 = -pkin(5) * t470 - t9;
t6 = qJ(6) * t470 - qJD(6) * t342 + t8;
t5 = t28 * pkin(5) + t27 * qJ(6) - qJD(6) * t554 + t60;
t2 = -t4 - t563;
t1 = t3 + t568;
t13 = [0, 0, 0, 0, 0, qJDD(1), -t541 + t545, t428, 0, 0, qJDD(1) * t332 + 0.2e1 * t435, 0.2e1 * t322 * t338 - 0.2e1 * t461 * t475, qJDD(2) * t338 + t342 * t345, qJDD(1) * t333 - 0.2e1 * t435, qJDD(2) * t342 - t338 * t345, 0, t403 * t338 + (-t387 + t545) * t342, t338 * t387 + t342 * t403 - t308, 0.2e1 * t474 * t522 - t428, -g(1) * (-pkin(1) * t339 + t329) - g(2) * t476 + (pkin(7) ^ 2 * t474 + pkin(1) ^ 2) * qJDD(1), -t158 * t490 + t261 * t386 (-t259 * t341 - t261 * t337) * t469 + (t512 - t511 + (-t503 + t506) * qJD(3)) * t338 (-t304 * t462 + t158) * t342 + (qJD(2) * t261 + t398) * t338, t159 * t494 + t259 * t385 (t304 * t471 + t159) * t342 + (-qJD(2) * t259 - t399) * t338, -t256 * t342 - t304 * t470, -g(1) * t238 - g(2) * t240 - t141 * t304 + t205 * t256 + (qJD(2) * t426 - t101) * t342 + (pkin(7) * t159 + qJD(2) * t184 + t230 * t337 + t284 * t466) * t338, -g(1) * t237 - g(2) * t239 + t140 * t304 - t206 * t256 + (qJD(2) * t425 + t100) * t342 + (-pkin(7) * t158 - qJD(2) * t185 + t230 * t341 - t284 * t468) * t338, -t140 * t259 - t141 * t261 + t158 * t205 - t159 * t206 + t308 + t418 * t469 + (-t541 - t100 * t337 - t101 * t341 + (t184 * t337 - t185 * t341) * qJD(3)) * t338, t100 * t206 + t185 * t140 + t101 * t205 + t184 * t141 - g(1) * t329 - g(2) * (t343 * t432 + t476) - t279 * t545 + (t230 * t338 + t284 * t469) * pkin(7), -t126 * t416 + t220 * t82, -t126 * t417 + t220 * t408 - t416 * t381 + (t351 * t416 + t529) * t338, t126 * t295 - t220 * t248 + t342 * t82 + t416 * t470, -t417 * t381 + (t351 * t417 + t573) * t338 (t295 * t393 + t408) * t342 + (qJD(2) * t417 - t295 * t351 - t502) * t338, -t248 * t342 - t295 * t470, -t50 * t295 + t119 * t248 - t198 * t417 + t273 * t408 - g(1) * t213 - g(2) * t215 + (qJD(2) * t395 - t25) * t342 + (t84 * qJD(2) + t195 * t352 + t514) * t338, -g(1) * t212 - g(2) * t214 - t120 * t248 - t126 * t195 - t132 * t220 + t198 * t416 + t24 * t342 - t273 * t82 + t295 * t49 - t470 * t85, t49 * t417 - t120 * t408 - t50 * t416 + t119 * t82 + t25 * t220 + t84 * t126 + t308 - t85 * t381 + (-t24 * t262 + t351 * t85 - t541) * t338, t24 * t120 + t85 * t49 + t25 * t119 + t84 * t50 + t132 * t273 + t195 * t198 - g(1) * (pkin(3) * t492 + t329) - g(2) * (t315 * t484 - t343 * t488 + t476) + (-g(1) * (-pkin(1) - t412) - g(2) * t549) * t339, t12, t369, t11, t424, -t366, t170, t105 * t110 + t130 * t59 + t144 * t60 + t186 * t28 + t233 * t56 - t282 * t9 + t31 * t470 - t342 * t4 + t429, t105 * t554 - t130 * t58 + t145 * t60 - t186 * t27 - t233 * t532 + t282 * t8 + t3 * t342 - t32 * t470 - t430, -t110 * t8 - t144 * t3 - t145 * t4 + t27 * t56 - t28 * t532 + t31 * t58 - t32 * t59 - t554 * t9 + t433, -g(1) * t379 - g(2) * t383 + t130 * t105 + t60 * t186 + t3 * t532 + t31 * t9 + t32 * t8 + t4 * t56, t12, t11, -t369, t170, t366, t424, t110 * t19 + t144 * t5 + t2 * t342 - t233 * t54 + t28 * t78 + t282 * t7 - t29 * t470 + t51 * t59 + t429, -t1 * t144 - t110 * t6 + t145 * t2 - t27 * t54 - t28 * t53 - t29 * t58 - t30 * t59 + t554 * t7 + t433, -t1 * t342 - t145 * t5 - t19 * t554 + t233 * t53 + t27 * t78 - t282 * t6 + t30 * t470 + t51 * t58 + t430, t1 * t53 + t30 * t6 + t5 * t78 + t51 * t19 + t2 * t54 + t29 * t7 - g(1) * (-pkin(5) * t201 - qJ(6) * t200 + t379) - g(2) * (pkin(5) * t203 + qJ(6) * t202 + t383); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t456, t475 * t346, t460, t456, t322, qJDD(2), t338 * t388 - t316 - t539, t540 + (t388 - t522) * t342, 0, 0, -t304 * t503 - t512 (-t158 + t507) * t341 + (-t159 + t504) * t337 (-t261 * t338 + t304 * t485) * qJD(1) + t399, -t304 * t506 - t511 (t259 * t338 - t304 * t493) * qJD(1) + t398, t304 * t473, -pkin(2) * t159 + t192 * t304 + t422 * t337 + (-t184 * t338 - t342 * t426) * qJD(1) + t598 * t341, pkin(2) * t158 - t193 * t304 + t422 * t341 + (t185 * t338 - t342 * t425) * qJD(1) - t598 * t337, t192 * t261 + t193 * t259 + ((qJD(3) * t261 - t159) * pkin(8) + t570) * t341 + (-t101 + t508 + (qJD(3) * t259 - t158) * pkin(8)) * t337 - t555, -t284 * t318 - t184 * t192 - t185 * t193 + t368 * pkin(2) + (qJD(3) * t418 + t100 * t341 - t101 * t337 - t555) * pkin(8), t416 * t569 - t529, t208 * t417 + t416 * t207 - t573 + t410 * t82 + t564 * (-t262 * t416 - t410 * t417) -t416 * t473 + t502 + (-t208 + t351) * t295, t408 * t410 + t417 * t595, -t410 * t248 - t295 * t595 - t417 * t473, t295 * t473, t190 * t248 - t315 * t408 + t132 * t410 - t84 * t473 + t251 * t417 - t480 * t295 + (-t417 * t549 + t395) * qJD(3) + t377 * t325 + (-t207 + t391) * t195, -t191 * t248 + t481 * t295 + t315 * t82 - t377 * t324 + t85 * t473 + t514 + (qJD(3) * t549 - t251) * t416 + t569 * t195, t417 * t481 - t416 * t480 + t190 * t82 - t191 * t408 + t85 * t207 - t84 * t208 - t24 * t410 - t25 * t262 - t555 + t564 * (-t262 * t85 + t410 * t84) t24 * t191 + t25 * t190 - t132 * t315 - g(3) * t412 + t481 * t85 + t480 * t84 + t583 * t195 + t428 * (t315 * t338 + t342 * t344) t10, t553, t33, t396, -t557, t447, t110 * t479 - t130 * t528 + t180 * t60 + t209 * t28 + t282 * t533 - t31 * t473 + t380, -t130 * t527 + t181 * t60 - t209 * t27 + t282 * t535 + t32 * t473 + t479 * t554 + t591, t110 * t48 - t180 * t3 - t181 * t4 + t31 * t527 + t32 * t528 + t533 * t554 + t360, t3 * t103 + t4 * t402 + t60 * t209 - g(3) * (-t331 * t338 + t247) + t535 * t32 - t533 * t31 + t479 * t130 + t428 * (t268 * t338 + t331 * t342) t10, t33, -t553, t447, t557, t396, t110 * t537 + t180 * t5 + t28 * t99 + t282 * t534 + t29 * t473 - t51 * t528 + t380, -t1 * t180 + t110 * t45 + t181 * t2 - t29 * t527 + t30 * t528 + t534 * t554 + t360, -t181 * t5 + t27 * t99 - t282 * t536 - t30 * t473 + t51 * t527 - t537 * t554 - t591, -g(3) * t247 + t1 * t103 - t2 * t402 + t5 * t99 + t537 * t51 + t536 * t30 + t534 * t29 + (-g(3) * t427 + t331 * t428) * t342 + (g(3) * t331 + t428 * (t268 + t427)) * t338; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t505, -t259 ^ 2 + t261 ^ 2, -t158 - t507, -t505, -t159 - t504, t256, -t285 * t466 - t508 - t261 * t284 + t179 + (-qJD(3) * t249 - t229 + t540) * t337 + t565, g(1) * t240 - g(2) * t238 + g(3) * t490 + t259 * t284 - t570, 0, 0, -t510, t562, t561, t510, t556, t248, t295 * t91 + (t248 * t340 + t261 * t417 + t295 * t465) * pkin(3) + t559, -t295 * t92 + (-t248 * t336 - t261 * t416 + t295 * t464) * pkin(3) + t558, t85 * t416 - t92 * t417 + t84 * t417 + t91 * t416 + (-t408 * t336 + t82 * t340 + (t336 * t416 + t340 * t417) * qJD(4)) * pkin(3), -t84 * t91 - t85 * t92 + (t24 * t336 + t25 * t340 - t195 * t261 + g(3) * t494 + (-t336 * t84 + t340 * t85) * qJD(4) + t565) * pkin(3), t596, -t446, t593, -t596, -t588, t233, -t110 * t135 + t231 * t233 + t361 - t585, -t135 * t554 - t232 * t233 + t282 * t525 + t597, -t110 * t525 + t231 * t27 - t232 * t28 + t440 + t592, -g(1) * t478 - g(2) * t437 - t130 * t135 + t4 * t231 + t3 * t232 + t275 * t540 - t31 * t524 + t32 * t525, t596, t593, t446, t233, t588, -t596, -t110 * t52 - t222 * t233 + t361 - t579, -t219 * t28 - t222 * t27 + t30 * t554 + t440 + (t29 - t526) * t110, t219 * t233 - t282 * t526 + t52 * t554 + t599, t1 * t219 + t2 * t222 - t51 * t52 - g(1) * (t438 + t478) - g(2) * (t437 + t439) - g(3) * (t278 + (-pkin(5) * t311 - t275) * t338) + t526 * t30 + t524 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t510, t562, t561, t510, t556, t248, -t295 * t85 + t559, -t295 * t84 + t558, 0, 0, t596, -t446, t593, -t596, -t588, t233, -t34 * t282 + (-t110 * t416 + t233 * t550 + t282 * t463) * pkin(4) + t560, -t35 * t282 + (-t233 * t335 + t282 * t445 - t416 * t554) * pkin(4) + t597, t35 * t110 - t34 * t554 + (t550 * t27 - t28 * t335 + (-t110 * t550 + t335 * t554) * qJD(5)) * pkin(4) + t592, -g(1) * t292 + t31 * t34 - t32 * t35 + (t4 * t550 - t130 * t416 + t3 * t335 + (-t31 * t335 + t32 * t550) * qJD(5) + t356) * pkin(4), t596, t593, t446, t233, t588, -t596, -t110 * t55 - t233 * t313 + t282 * t434 - t359, -t27 * t313 - t28 * t309 + (t30 + t434) * t554 + (-t523 + t29) * t110, t233 * t309 - t282 * t523 + t55 * t554 + t599, t1 * t309 + t2 * t313 - t51 * t55 - t29 * t34 - g(1) * (t292 + t438) + t523 * t30 + (t29 * t463 + t356) * pkin(4) + t363; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t596, -t446, t593, -t596, -t588, t233, -t531 + t560, t364 + t587, 0, 0, t596, t593, t446, t233, t588, -t596, -t110 * t64 + t218 - t359 - t531, pkin(5) * t27 - qJ(6) * t28 + (t30 - t32) * t554 + (t29 - t482) * t110, t554 * t64 + 0.2e1 * t216 - 0.2e1 * t274 - t364 - t590, -t2 * pkin(5) - g(1) * t438 + t1 * qJ(6) - t29 * t32 + t30 * t482 - t51 * t64 + t363; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t233 + t596, t593, -t282 ^ 2 - t551, t282 * t30 + t359;];
tau_reg  = t13;
