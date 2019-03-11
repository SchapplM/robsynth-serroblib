% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 19:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRPRR12_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR12_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR12_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR12_coriolisvecJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR12_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR12_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR12_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:38:13
% EndTime: 2019-03-09 19:39:16
% DurationCPUTime: 32.59s
% Computational Cost: add. (31393->988), mult. (80557->1394), div. (0->0), fcn. (65208->12), ass. (0->403)
t395 = sin(qJ(2));
t389 = sin(pkin(6));
t439 = qJD(1) * t389;
t427 = t395 * t439;
t391 = cos(pkin(6));
t399 = cos(qJ(2));
t498 = pkin(1) * t399;
t431 = t391 * t498;
t342 = -pkin(8) * t427 + qJD(1) * t431;
t407 = (pkin(2) * t395 - pkin(9) * t399) * t389;
t343 = qJD(1) * t407;
t394 = sin(qJ(3));
t398 = cos(qJ(3));
t256 = t398 * t342 + t394 * t343;
t234 = qJ(4) * t427 + t256;
t499 = pkin(1) * t395;
t383 = t391 * t499;
t410 = pkin(3) * t394 - qJ(4) * t398;
t455 = t389 * t399;
t262 = (t383 + (pkin(8) + t410) * t455) * qJD(1);
t388 = sin(pkin(12));
t390 = cos(pkin(12));
t172 = -t388 * t234 + t390 * t262;
t348 = qJD(3) * t410 - qJD(4) * t394;
t435 = qJD(3) * t394;
t430 = pkin(9) * t435;
t285 = t390 * t348 + t388 * t430;
t576 = t285 - t172;
t451 = t398 * t399;
t305 = (t388 * t395 + t390 * t451) * t439;
t426 = t399 * t439;
t418 = t394 * t426;
t452 = t390 * t398;
t575 = -pkin(4) * t418 + t305 * pkin(10) + (pkin(4) * t394 - pkin(10) * t452) * qJD(3) + t576;
t173 = t390 * t234 + t388 * t262;
t304 = (-t388 * t451 + t390 * t395) * t439;
t333 = t388 * t348;
t454 = t390 * t394;
t457 = t388 * t398;
t574 = pkin(10) * t304 + t173 - t333 - (-pkin(9) * t454 - pkin(10) * t457) * qJD(3);
t393 = sin(qJ(5));
t397 = cos(qJ(5));
t223 = t304 * t397 - t305 * t393;
t453 = t390 * t397;
t408 = t388 * t393 - t453;
t349 = t408 * qJD(5);
t364 = t388 * t397 + t390 * t393;
t434 = qJD(3) * t398;
t260 = t349 * t394 - t364 * t434;
t446 = t223 - t260;
t224 = t304 * t393 + t305 * t397;
t350 = t364 * qJD(5);
t259 = -t350 * t394 - t408 * t434;
t445 = t224 - t259;
t367 = -pkin(3) * t398 - qJ(4) * t394 - pkin(2);
t358 = t390 * t367;
t281 = -pkin(10) * t454 + t358 + (-pkin(9) * t388 - pkin(4)) * t398;
t320 = pkin(9) * t452 + t388 * t367;
t458 = t388 * t394;
t293 = -pkin(10) * t458 + t320;
t432 = qJD(5) * t397;
t433 = qJD(5) * t393;
t561 = t281 * t432 - t293 * t433 + t393 * t575 - t574 * t397;
t205 = t393 * t281 + t397 * t293;
t560 = -qJD(5) * t205 + t574 * t393 + t397 * t575;
t419 = qJD(1) * t391 + qJD(2);
t540 = -t394 * t427 + t398 * t419;
t225 = t364 * t540;
t573 = t225 - t350;
t226 = t408 * t540;
t443 = -t226 + t349;
t572 = t560 + t445 * pkin(11) + (-t418 + t435) * pkin(5);
t571 = pkin(11) * t446 - t561;
t382 = pkin(8) * t455;
t302 = qJD(2) * pkin(9) + (t382 + (pkin(9) + t499) * t391) * qJD(1);
t337 = (-pkin(2) * t399 - pkin(9) * t395 - pkin(1)) * t389;
t313 = qJD(1) * t337;
t220 = -t394 * t302 + t313 * t398;
t322 = t394 * t419 + t398 * t427;
t242 = pkin(3) * t322 - qJ(4) * t540;
t159 = -t220 * t388 + t390 * t242;
t120 = -pkin(10) * t390 * t540 + pkin(4) * t322 + t159;
t160 = t390 * t220 + t388 * t242;
t459 = t540 * t388;
t138 = -pkin(10) * t459 + t160;
t495 = pkin(10) + qJ(4);
t368 = t495 * t388;
t369 = t495 * t390;
t545 = qJD(4) * t453 - t397 * t138 - t368 * t432 + (-qJD(4) * t388 - qJD(5) * t369 - t120) * t393;
t297 = -t393 * t368 + t397 * t369;
t544 = -t364 * qJD(4) - qJD(5) * t297 - t397 * t120 + t138 * t393;
t372 = qJD(3) - t426;
t265 = t322 * t390 + t372 * t388;
t420 = -t322 * t388 + t390 * t372;
t190 = t265 * t397 + t393 * t420;
t392 = sin(qJ(6));
t396 = cos(qJ(6));
t556 = -t265 * t393 + t397 * t420;
t570 = -t190 * t392 + t396 * t556;
t108 = t190 * t396 + t392 * t556;
t569 = pkin(11) * t573 + t545;
t568 = -pkin(5) * t322 + pkin(11) * t443 + t544;
t255 = -t394 * t342 + t343 * t398;
t235 = -pkin(3) * t427 - t255;
t558 = pkin(4) * t304 - t235 + (pkin(4) * t388 + pkin(9)) * t434;
t103 = Ifges(7,4) * t570;
t315 = qJD(5) - t540;
t301 = -pkin(2) * t419 - t342;
t201 = -pkin(3) * t540 - t322 * qJ(4) + t301;
t221 = t398 * t302 + t394 * t313;
t206 = qJ(4) * t372 + t221;
t126 = t388 * t201 + t390 * t206;
t100 = pkin(10) * t420 + t126;
t125 = t390 * t201 - t206 * t388;
t84 = -pkin(4) * t540 - pkin(10) * t265 + t125;
t45 = -t100 * t393 + t397 * t84;
t38 = -pkin(11) * t190 + t45;
t37 = pkin(5) * t315 + t38;
t46 = t100 * t397 + t393 * t84;
t39 = pkin(11) * t556 + t46;
t463 = t39 * t392;
t16 = t37 * t396 - t463;
t462 = t39 * t396;
t17 = t37 * t392 + t462;
t489 = Ifges(7,4) * t108;
t306 = qJD(6) + t315;
t511 = -t306 / 0.2e1;
t524 = -t108 / 0.2e1;
t526 = -t570 / 0.2e1;
t53 = Ifges(7,1) * t108 + Ifges(7,5) * t306 + t103;
t203 = -pkin(3) * t372 + qJD(4) - t220;
t162 = -pkin(4) * t420 + t203;
t99 = -pkin(5) * t556 + t162;
t565 = (Ifges(7,5) * t570 - Ifges(7,6) * t108) * t511 + (t108 * t17 + t16 * t570) * mrSges(7,3) + (-Ifges(7,2) * t108 + t103 + t53) * t526 - t99 * (mrSges(7,1) * t108 + mrSges(7,2) * t570) + (Ifges(7,1) * t570 - t489) * t524;
t438 = qJD(2) * t389;
t422 = qJD(1) * t438;
t415 = t399 * t422;
t275 = qJD(3) * t540 + t398 * t415;
t416 = t395 * t422;
t233 = t275 * t390 + t388 * t416;
t276 = qJD(3) * t322 + t394 * t415;
t344 = qJD(2) * t407;
t330 = qJD(1) * t344;
t456 = t389 * t395;
t379 = pkin(8) * t456;
t353 = -t379 + t431;
t346 = t353 * qJD(2);
t331 = qJD(1) * t346;
t164 = -t302 * t435 + t313 * t434 + t394 * t330 + t398 * t331;
t144 = qJ(4) * t416 + qJD(4) * t372 + t164;
t440 = t382 + t383;
t347 = t440 * qJD(2);
t332 = qJD(1) * t347;
t153 = t276 * pkin(3) - t275 * qJ(4) - t322 * qJD(4) + t332;
t76 = -t144 * t388 + t390 * t153;
t59 = pkin(4) * t276 - pkin(10) * t233 + t76;
t232 = -t275 * t388 + t390 * t416;
t77 = t390 * t144 + t388 * t153;
t65 = pkin(10) * t232 + t77;
t13 = -qJD(5) * t46 - t393 * t65 + t397 * t59;
t97 = qJD(5) * t556 + t232 * t393 + t233 * t397;
t6 = pkin(5) * t276 - pkin(11) * t97 + t13;
t12 = -t100 * t433 + t393 * t59 + t397 * t65 + t84 * t432;
t98 = -qJD(5) * t190 + t232 * t397 - t233 * t393;
t7 = pkin(11) * t98 + t12;
t2 = qJD(6) * t16 + t392 * t6 + t396 * t7;
t3 = -qJD(6) * t17 - t392 * t7 + t396 * t6;
t564 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t204 = t397 * t281 - t293 * t393;
t339 = t408 * t394;
t176 = -pkin(5) * t398 + pkin(11) * t339 + t204;
t338 = t364 * t394;
t183 = -pkin(11) * t338 + t205;
t102 = t176 * t392 + t183 * t396;
t563 = -qJD(6) * t102 + t571 * t392 + t396 * t572;
t101 = t176 * t396 - t183 * t392;
t562 = qJD(6) * t101 + t392 * t572 - t571 * t396;
t559 = pkin(5) * t446 + t558;
t314 = Ifges(4,4) * t540;
t466 = t372 * Ifges(4,5);
t468 = t322 * Ifges(4,1);
t216 = t314 + t466 + t468;
t557 = -t314 / 0.2e1 + t220 * mrSges(4,3) - t216 / 0.2e1 - t301 * mrSges(4,2) - t466 / 0.2e1;
t35 = qJD(6) * t570 + t392 * t98 + t396 * t97;
t535 = t35 / 0.2e1;
t36 = -qJD(6) * t108 - t392 * t97 + t396 * t98;
t534 = t36 / 0.2e1;
t52 = Ifges(7,2) * t570 + Ifges(7,6) * t306 + t489;
t554 = t52 / 0.2e1;
t531 = t97 / 0.2e1;
t530 = t98 / 0.2e1;
t513 = t276 / 0.2e1;
t423 = -Ifges(3,6) * qJD(2) / 0.2e1;
t92 = Ifges(6,6) * t98;
t93 = Ifges(6,5) * t97;
t42 = Ifges(6,3) * t276 + t92 + t93;
t33 = Ifges(7,6) * t36;
t34 = Ifges(7,5) * t35;
t8 = Ifges(7,3) * t276 + t33 + t34;
t553 = t42 + t8;
t552 = -t13 * mrSges(6,1) + t12 * mrSges(6,2) - t564;
t296 = -t397 * t368 - t369 * t393;
t257 = -pkin(11) * t364 + t296;
t258 = -pkin(11) * t408 + t297;
t179 = t257 * t396 - t258 * t392;
t551 = qJD(6) * t179 + t568 * t392 + t396 * t569;
t180 = t257 * t392 + t258 * t396;
t550 = -qJD(6) * t180 - t392 * t569 + t568 * t396;
t335 = t379 + (-pkin(2) - t498) * t391;
t351 = -t391 * t398 + t394 * t456;
t352 = t391 * t394 + t398 * t456;
t228 = t351 * pkin(3) - t352 * qJ(4) + t335;
t336 = pkin(9) * t391 + t440;
t245 = t398 * t336 + t394 * t337;
t230 = -qJ(4) * t455 + t245;
t155 = t390 * t228 - t230 * t388;
t290 = t352 * t390 - t388 * t455;
t115 = pkin(4) * t351 - pkin(10) * t290 + t155;
t156 = t388 * t228 + t390 * t230;
t289 = -t352 * t388 - t390 * t455;
t127 = pkin(10) * t289 + t156;
t63 = t393 * t115 + t397 * t127;
t193 = pkin(4) * t459 + t221;
t541 = -pkin(5) * t573 - t193;
t475 = t265 * Ifges(5,5);
t476 = t420 * Ifges(5,6);
t167 = -Ifges(5,3) * t540 + t475 + t476;
t465 = t372 * Ifges(4,6);
t493 = Ifges(4,4) * t322;
t215 = Ifges(4,2) * t540 + t465 + t493;
t470 = t315 * Ifges(6,3);
t471 = t306 * Ifges(7,3);
t479 = t190 * Ifges(6,5);
t480 = t556 * Ifges(6,6);
t483 = t108 * Ifges(7,5);
t484 = t570 * Ifges(7,6);
t51 = t471 + t483 + t484;
t94 = t470 + t479 + t480;
t400 = -t301 * mrSges(4,1) - t476 / 0.2e1 - t475 / 0.2e1 + t465 / 0.2e1 + t215 / 0.2e1 - t470 / 0.2e1 - t471 / 0.2e1 - t479 / 0.2e1 - t480 / 0.2e1 - t167 / 0.2e1 - t483 / 0.2e1 - t484 / 0.2e1 - t94 / 0.2e1 - t51 / 0.2e1 + t221 * mrSges(4,3) - t125 * mrSges(5,1) + t126 * mrSges(5,2) + t46 * mrSges(6,2) - t45 * mrSges(6,1) - t16 * mrSges(7,1) + t17 * mrSges(7,2) + t493 / 0.2e1;
t429 = Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1;
t539 = t429 * t540 + t400;
t168 = t265 * Ifges(5,4) + Ifges(5,2) * t420 - Ifges(5,6) * t540;
t169 = t265 * Ifges(5,1) + Ifges(5,4) * t420 - Ifges(5,5) * t540;
t491 = Ifges(5,4) * t390;
t412 = -Ifges(5,2) * t388 + t491;
t492 = Ifges(5,4) * t388;
t413 = Ifges(5,1) * t390 - t492;
t414 = mrSges(5,1) * t388 + mrSges(5,2) * t390;
t500 = t390 / 0.2e1;
t501 = -t388 / 0.2e1;
t514 = t265 / 0.2e1;
t515 = t420 / 0.2e1;
t538 = t203 * t414 + t412 * t515 + t413 * t514 + (-t125 * t390 - t126 * t388) * mrSges(5,3) + t168 * t501 + t169 * t500 - t557;
t537 = Ifges(7,4) * t535 + Ifges(7,2) * t534 + Ifges(7,6) * t513;
t536 = Ifges(7,1) * t535 + Ifges(7,4) * t534 + Ifges(7,5) * t513;
t533 = Ifges(6,4) * t531 + Ifges(6,2) * t530 + Ifges(6,6) * t513;
t532 = Ifges(6,1) * t531 + Ifges(6,4) * t530 + Ifges(6,5) * t513;
t529 = pkin(1) * mrSges(3,1);
t528 = pkin(1) * mrSges(3,2);
t525 = t570 / 0.2e1;
t523 = t108 / 0.2e1;
t130 = t233 * Ifges(5,1) + t232 * Ifges(5,4) + t276 * Ifges(5,5);
t522 = t130 / 0.2e1;
t521 = -t556 / 0.2e1;
t520 = t556 / 0.2e1;
t519 = -t190 / 0.2e1;
t518 = t190 / 0.2e1;
t517 = t232 / 0.2e1;
t516 = t233 / 0.2e1;
t510 = t306 / 0.2e1;
t508 = -t315 / 0.2e1;
t507 = t315 / 0.2e1;
t506 = t540 / 0.2e1;
t505 = -t540 / 0.2e1;
t504 = -t351 / 0.2e1;
t502 = t352 / 0.2e1;
t387 = t394 * pkin(9);
t496 = qJD(2) / 0.2e1;
t494 = Ifges(3,4) * t395;
t490 = Ifges(6,4) * t190;
t488 = Ifges(3,5) * t391;
t487 = Ifges(5,5) * t390;
t486 = Ifges(3,6) * t391;
t485 = Ifges(5,6) * t388;
t482 = t164 * mrSges(4,2);
t165 = -t302 * t434 - t313 * t435 + t330 * t398 - t394 * t331;
t481 = t165 * mrSges(4,1);
t478 = t232 * Ifges(5,6);
t477 = t233 * Ifges(5,5);
t474 = t275 * Ifges(4,1);
t473 = t275 * Ifges(4,4);
t472 = t276 * Ifges(4,4);
t248 = -t338 * t396 + t339 * t392;
t140 = qJD(6) * t248 + t259 * t396 + t260 * t392;
t150 = t223 * t392 + t224 * t396;
t450 = t140 - t150;
t249 = -t338 * t392 - t339 * t396;
t141 = -qJD(6) * t249 - t259 * t392 + t260 * t396;
t149 = t223 * t396 - t224 * t392;
t449 = t141 - t149;
t151 = -t225 * t396 + t226 * t392;
t280 = t364 * t396 - t392 * t408;
t200 = -qJD(6) * t280 + t349 * t392 - t350 * t396;
t448 = t151 - t200;
t152 = -t225 * t392 - t226 * t396;
t279 = -t364 * t392 - t396 * t408;
t199 = qJD(6) * t279 - t349 * t396 - t350 * t392;
t447 = t152 - t199;
t181 = -t336 * t435 + t337 * t434 + t394 * t344 + t398 * t346;
t437 = qJD(2) * t395;
t166 = (qJ(4) * t437 - qJD(4) * t399) * t389 + t181;
t436 = qJD(2) * t399;
t424 = t389 * t436;
t291 = qJD(3) * t352 + t394 * t424;
t292 = -qJD(3) * t351 + t398 * t424;
t175 = t291 * pkin(3) - t292 * qJ(4) - t352 * qJD(4) + t347;
t86 = t390 * t166 + t388 * t175;
t442 = -mrSges(3,1) * t419 - mrSges(4,1) * t540 + mrSges(4,2) * t322 + mrSges(3,3) * t427;
t194 = -mrSges(5,1) * t420 + mrSges(5,2) * t265;
t278 = mrSges(4,1) * t372 - mrSges(4,3) * t322;
t441 = t278 - t194;
t365 = pkin(4) * t458 + t387;
t428 = Ifges(4,5) * t275 - Ifges(4,6) * t276 + Ifges(4,3) * t416;
t385 = -pkin(4) * t390 - pkin(3);
t425 = t389 * t437;
t48 = -t98 * mrSges(6,1) + t97 * mrSges(6,2);
t11 = -t36 * mrSges(7,1) + t35 * mrSges(7,2);
t163 = -t232 * mrSges(5,1) + t233 * mrSges(5,2);
t62 = t397 * t115 - t127 * t393;
t85 = -t166 * t388 + t390 * t175;
t244 = -t394 * t336 + t337 * t398;
t231 = pkin(3) * t455 - t244;
t411 = -t485 + t487;
t409 = -t388 * t76 + t390 * t77;
t208 = t289 * t393 + t290 * t397;
t47 = pkin(5) * t351 - pkin(11) * t208 + t62;
t207 = t289 * t397 - t290 * t393;
t50 = pkin(11) * t207 + t63;
t22 = -t392 * t50 + t396 * t47;
t23 = t392 * t47 + t396 * t50;
t136 = t207 * t396 - t208 * t392;
t137 = t207 * t392 + t208 * t396;
t182 = -t336 * t434 - t337 * t435 + t344 * t398 - t394 * t346;
t254 = t292 * t390 + t388 * t425;
t71 = pkin(4) * t291 - pkin(10) * t254 + t85;
t253 = -t292 * t388 + t390 * t425;
t75 = pkin(10) * t253 + t86;
t20 = t115 * t432 - t127 * t433 + t393 * t71 + t397 * t75;
t192 = -pkin(4) * t289 + t231;
t373 = Ifges(3,4) * t426;
t404 = -t342 * mrSges(3,3) + Ifges(3,1) * t427 / 0.2e1 + t373 / 0.2e1 + (t419 / 0.2e1 + t496) * Ifges(3,5);
t174 = -pkin(3) * t425 - t182;
t21 = -qJD(5) * t63 - t393 * t75 + t397 * t71;
t154 = -pkin(3) * t416 - t165;
t131 = -pkin(4) * t253 + t174;
t111 = -pkin(4) * t232 + t154;
t345 = t440 * qJD(1);
t402 = t220 * mrSges(4,1) + t372 * Ifges(4,3) + t322 * Ifges(4,5) + t540 * Ifges(4,6) + t423 - (t486 + (Ifges(3,2) * t399 + t494) * t389) * qJD(1) / 0.2e1 - t221 * mrSges(4,2) - t345 * mrSges(3,3);
t371 = Ifges(3,5) * t415;
t341 = -mrSges(3,2) * t419 + mrSges(3,3) * t426;
t329 = pkin(5) * t408 + t385;
t319 = -pkin(9) * t457 + t358;
t286 = -t390 * t430 + t333;
t283 = pkin(5) * t338 + t365;
t277 = -mrSges(4,2) * t372 + mrSges(4,3) * t540;
t239 = -mrSges(4,2) * t416 - mrSges(4,3) * t276;
t238 = mrSges(4,1) * t416 - mrSges(4,3) * t275;
t210 = -mrSges(5,1) * t540 - mrSges(5,3) * t265;
t209 = mrSges(5,2) * t540 + mrSges(5,3) * t420;
t198 = mrSges(4,1) * t276 + mrSges(4,2) * t275;
t187 = Ifges(4,5) * t416 - t472 + t474;
t186 = -t276 * Ifges(4,2) + Ifges(4,6) * t416 + t473;
t185 = Ifges(6,4) * t556;
t178 = mrSges(5,1) * t276 - mrSges(5,3) * t233;
t177 = -mrSges(5,2) * t276 + mrSges(5,3) * t232;
t158 = mrSges(6,1) * t315 - mrSges(6,3) * t190;
t157 = -mrSges(6,2) * t315 + mrSges(6,3) * t556;
t129 = t233 * Ifges(5,4) + t232 * Ifges(5,2) + t276 * Ifges(5,6);
t128 = t276 * Ifges(5,3) + t477 + t478;
t123 = -pkin(5) * t207 + t192;
t119 = -qJD(5) * t208 + t253 * t397 - t254 * t393;
t118 = qJD(5) * t207 + t253 * t393 + t254 * t397;
t110 = -mrSges(6,1) * t556 + mrSges(6,2) * t190;
t96 = Ifges(6,1) * t190 + Ifges(6,5) * t315 + t185;
t95 = Ifges(6,2) * t556 + Ifges(6,6) * t315 + t490;
t89 = mrSges(7,1) * t306 - mrSges(7,3) * t108;
t88 = -mrSges(7,2) * t306 + mrSges(7,3) * t570;
t79 = -mrSges(6,2) * t276 + mrSges(6,3) * t98;
t78 = mrSges(6,1) * t276 - mrSges(6,3) * t97;
t66 = -pkin(5) * t119 + t131;
t56 = -mrSges(7,1) * t570 + mrSges(7,2) * t108;
t54 = -pkin(5) * t98 + t111;
t41 = -qJD(6) * t137 - t118 * t392 + t119 * t396;
t40 = qJD(6) * t136 + t118 * t396 + t119 * t392;
t29 = -mrSges(7,2) * t276 + mrSges(7,3) * t36;
t28 = mrSges(7,1) * t276 - mrSges(7,3) * t35;
t19 = t38 * t396 - t463;
t18 = -t38 * t392 - t462;
t15 = pkin(11) * t119 + t20;
t14 = pkin(5) * t291 - pkin(11) * t118 + t21;
t5 = -qJD(6) * t23 + t14 * t396 - t15 * t392;
t4 = qJD(6) * t22 + t14 * t392 + t15 * t396;
t1 = [(t167 + t94 + t51) * t291 / 0.2e1 + m(3) * (t331 * t440 - t332 * t353 - t342 * t347 + t345 * t346) + ((-t353 * mrSges(3,3) + t488 + (-0.2e1 * t528 + 0.3e1 / 0.2e1 * Ifges(3,4) * t399) * t389) * t436 + (-0.3e1 / 0.2e1 * t486 + Ifges(4,5) * t502 + Ifges(4,6) * t504 - t440 * mrSges(3,3) + (-0.2e1 * t529 - 0.3e1 / 0.2e1 * t494 + (-Ifges(4,3) / 0.2e1 + 0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t399) * t389) * t437) * t439 + t372 * (Ifges(4,5) * t292 - Ifges(4,6) * t291) / 0.2e1 + (Ifges(5,5) * t290 + Ifges(6,5) * t208 + Ifges(7,5) * t137 + Ifges(5,6) * t289 + Ifges(6,6) * t207 + Ifges(7,6) * t136 + (Ifges(5,3) + Ifges(6,3) + Ifges(7,3)) * t351) * t513 + t391 * t371 / 0.2e1 + m(5) * (t125 * t85 + t126 * t86 + t154 * t231 + t155 * t76 + t156 * t77 + t174 * t203) + m(6) * (t111 * t192 + t12 * t63 + t13 * t62 + t131 * t162 + t20 * t46 + t21 * t45) + m(7) * (t123 * t54 + t16 * t5 + t17 * t4 + t2 * t23 + t22 * t3 + t66 * t99) + m(4) * (t164 * t245 + t165 * t244 + t181 * t221 + t182 * t220 + t301 * t347 + t332 * t335) + t41 * t554 + (-t164 * t351 - t165 * t352 - t220 * t292 - t221 * t291) * mrSges(4,3) + t322 * (Ifges(4,1) * t292 - Ifges(4,4) * t291) / 0.2e1 + t13 * (mrSges(6,1) * t351 - mrSges(6,3) * t208) + t2 * (-mrSges(7,2) * t351 + mrSges(7,3) * t136) + t3 * (mrSges(7,1) * t351 - mrSges(7,3) * t137) + t77 * (-mrSges(5,2) * t351 + mrSges(5,3) * t289) + t76 * (mrSges(5,1) * t351 - mrSges(5,3) * t290) + t12 * (-mrSges(6,2) * t351 + mrSges(6,3) * t207) + t346 * t341 + t335 * t198 + t301 * (mrSges(4,1) * t291 + mrSges(4,2) * t292) - t291 * t215 / 0.2e1 + t292 * t216 / 0.2e1 + t289 * t129 / 0.2e1 + t154 * (-mrSges(5,1) * t289 + mrSges(5,2) * t290) + t126 * (-mrSges(5,2) * t291 + mrSges(5,3) * t253) + t125 * (mrSges(5,1) * t291 - mrSges(5,3) * t254) + t45 * (mrSges(6,1) * t291 - mrSges(6,3) * t118) + t46 * (-mrSges(6,2) * t291 + mrSges(6,3) * t119) + t17 * (-mrSges(7,2) * t291 + mrSges(7,3) * t41) + t16 * (mrSges(7,1) * t291 - mrSges(7,3) * t40) + t181 * t277 + t182 * t278 + t253 * t168 / 0.2e1 + t203 * (-mrSges(5,1) * t253 + mrSges(5,2) * t254) + t254 * t169 / 0.2e1 + t244 * t238 + t245 * t239 + t231 * t163 + t111 * (-mrSges(6,1) * t207 + mrSges(6,2) * t208) + t86 * t209 + t85 * t210 + t192 * t48 + t174 * t194 + t155 * t178 + t156 * t177 + t162 * (-mrSges(6,1) * t119 + mrSges(6,2) * t118) + t20 * t157 + t21 * t158 + t54 * (-mrSges(7,1) * t136 + mrSges(7,2) * t137) + t131 * t110 + t123 * t11 + t118 * t96 / 0.2e1 + t119 * t95 / 0.2e1 + t99 * (-mrSges(7,1) * t41 + mrSges(7,2) * t40) + t4 * t88 + t5 * t89 + t63 * t79 + t62 * t78 + t66 * t56 + t40 * t53 / 0.2e1 + t207 * t533 + (Ifges(7,4) * t137 + Ifges(7,2) * t136 + Ifges(7,6) * t351) * t534 + (Ifges(7,1) * t137 + Ifges(7,4) * t136 + Ifges(7,5) * t351) * t535 + t137 * t536 + t136 * t537 + t290 * t522 + (Ifges(7,1) * t40 + Ifges(7,4) * t41 + Ifges(7,5) * t291) * t523 + (Ifges(7,4) * t40 + Ifges(7,2) * t41 + Ifges(7,6) * t291) * t525 + (Ifges(6,4) * t208 + Ifges(6,2) * t207 + Ifges(6,6) * t351) * t530 + (Ifges(6,1) * t208 + Ifges(6,4) * t207 + Ifges(6,5) * t351) * t531 + t208 * t532 + (Ifges(5,1) * t254 + Ifges(5,4) * t253 + Ifges(5,5) * t291) * t514 + (Ifges(5,4) * t254 + Ifges(5,2) * t253 + Ifges(5,6) * t291) * t515 + (Ifges(5,1) * t290 + Ifges(5,4) * t289 + Ifges(5,5) * t351) * t516 + (Ifges(5,4) * t290 + Ifges(5,2) * t289 + Ifges(5,6) * t351) * t517 + (Ifges(6,1) * t118 + Ifges(6,4) * t119 + Ifges(6,5) * t291) * t518 + (Ifges(6,4) * t118 + Ifges(6,2) * t119 + Ifges(6,6) * t291) * t520 + t187 * t502 + t186 * t504 + (Ifges(5,5) * t254 + Ifges(5,6) * t253 + Ifges(5,3) * t291) * t505 + (Ifges(4,4) * t292 - Ifges(4,2) * t291) * t506 + (Ifges(6,5) * t118 + Ifges(6,6) * t119 + Ifges(6,3) * t291) * t507 + (Ifges(7,5) * t40 + Ifges(7,6) * t41 + Ifges(7,3) * t291) * t510 + t22 * t28 + t23 * t29 + t442 * t347 + (t128 + t553) * t351 / 0.2e1 - t428 * t455 / 0.2e1 - t276 * (Ifges(4,4) * t352 - Ifges(4,2) * t351 - Ifges(4,6) * t455) / 0.2e1 + t275 * (Ifges(4,1) * t352 - Ifges(4,4) * t351 - Ifges(4,5) * t455) / 0.2e1 + t331 * (-t391 * mrSges(3,2) + mrSges(3,3) * t455) + (-mrSges(3,1) * t391 + mrSges(4,1) * t351 + mrSges(4,2) * t352 + mrSges(3,3) * t456) * t332 + (t404 * t399 + (t423 + t402) * t395) * t438 + t455 * t482 - t455 * t481; t576 * t210 + (t125 * t305 - t126 * t304) * mrSges(5,3) - t420 * (Ifges(5,4) * t305 + Ifges(5,2) * t304) / 0.2e1 - t265 * (Ifges(5,1) * t305 + Ifges(5,4) * t304) / 0.2e1 + t371 + (-t150 / 0.2e1 + t140 / 0.2e1) * t53 + (t260 / 0.2e1 - t223 / 0.2e1) * t95 + (-t149 / 0.2e1 + t141 / 0.2e1) * t52 + (t286 - t173) * t209 + t111 * (mrSges(6,1) * t338 - mrSges(6,2) * t339) + (-Ifges(6,4) * t339 - Ifges(6,2) * t338) * t530 + (-Ifges(6,1) * t339 - Ifges(6,4) * t338) * t531 + (-t12 * t338 + t13 * t339 + t445 * t45 - t446 * t46) * mrSges(6,3) + (t259 / 0.2e1 - t224 / 0.2e1) * t96 - m(5) * (t125 * t172 + t126 * t173 + t203 * t235) - m(4) * (t220 * t255 + t221 * t256 + t301 * t345) + (t187 / 0.2e1 + t413 * t516 + t332 * mrSges(4,2) - t472 / 0.2e1 + t474 / 0.2e1 + t412 * t517 + t154 * t414 - t165 * mrSges(4,3) + t129 * t501 + t130 * t500 + (-t238 + t163) * pkin(9) + (-t388 * t77 - t390 * t76) * mrSges(5,3)) * t394 + t365 * t48 - t342 * t341 - t331 * mrSges(3,2) - t332 * mrSges(3,1) + t319 * t178 + t320 * t177 - t203 * (-mrSges(5,1) * t304 + mrSges(5,2) * t305) - t305 * t169 / 0.2e1 - t304 * t168 / 0.2e1 + t283 * t11 - t256 * t277 - t255 * t278 + t54 * (-mrSges(7,1) * t248 + mrSges(7,2) * t249) - t235 * t194 + t204 * t78 + t205 * t79 - pkin(2) * t198 + t102 * t29 + t101 * t28 + (Ifges(7,4) * t249 + Ifges(7,2) * t248) * t534 + (Ifges(7,1) * t249 + Ifges(7,4) * t248) * t535 + t249 * t536 + t248 * t537 + (Ifges(7,1) * t140 + Ifges(7,4) * t141) * t523 + (Ifges(7,1) * t150 + Ifges(7,4) * t149) * t524 + (Ifges(7,4) * t140 + Ifges(7,2) * t141) * t525 + (Ifges(7,4) * t150 + Ifges(7,2) * t149) * t526 - t339 * t532 - t338 * t533 + (Ifges(6,1) * t259 + Ifges(6,4) * t260) * t518 + (Ifges(6,1) * t224 + Ifges(6,4) * t223) * t519 + (Ifges(6,4) * t259 + Ifges(6,2) * t260) * t520 + (Ifges(6,4) * t224 + Ifges(6,2) * t223) * t521 + (Ifges(5,5) * t305 + Ifges(5,6) * t304) * t506 + (Ifges(6,5) * t259 + Ifges(6,6) * t260) * t507 + (Ifges(6,5) * t224 + Ifges(6,6) * t223) * t508 + (Ifges(7,5) * t140 + Ifges(7,6) * t141) * t510 + (Ifges(7,5) * t150 + Ifges(7,6) * t149) * t511 + (-Ifges(6,5) * t339 + Ifges(7,5) * t249 - Ifges(6,6) * t338 + Ifges(7,6) * t248 + t394 * t411) * t513 + ((t411 * t505 + t468 / 0.2e1 + (-m(4) * t220 + m(5) * t203 - t441) * pkin(9) + t538) * t398 + ((-m(4) * t221 - t277) * pkin(9) - t539) * t394) * qJD(3) - t442 * t345 + (mrSges(6,1) * t446 - mrSges(6,2) * t445) * t162 + (t164 * mrSges(4,3) + t473 / 0.2e1 - t332 * mrSges(4,1) + t77 * mrSges(5,2) - t478 / 0.2e1 - t477 / 0.2e1 - t76 * mrSges(5,1) - t128 / 0.2e1 - t8 / 0.2e1 - t42 / 0.2e1 + t186 / 0.2e1 + pkin(9) * t239 - t92 / 0.2e1 - t93 / 0.2e1 - t33 / 0.2e1 - t34 / 0.2e1 + (-Ifges(6,3) / 0.2e1 - Ifges(7,3) / 0.2e1 - t429) * t276 + t552) * t398 + (-mrSges(7,1) * t449 + mrSges(7,2) * t450) * t99 + (-t16 * t450 + t17 * t449 + t2 * t248 - t249 * t3) * mrSges(7,3) + ((t423 + (Ifges(4,5) * t394 + Ifges(4,6) * t398) * t496 + (t486 / 0.2e1 + (t529 + t494 / 0.2e1) * t389) * qJD(1) - t402) * t395 + (-t373 / 0.2e1 + (-t488 / 0.2e1 + (t528 + (-Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t395) * t389) * qJD(1) + (-t468 / 0.2e1 + t557) * t398 + t539 * t394 - t404) * t399) * t439 + t558 * t110 + t559 * t56 + t560 * t158 + t561 * t157 + (t111 * t365 + t12 * t205 + t13 * t204 + t162 * t558 + t45 * t560 + t46 * t561) * m(6) + m(4) * (pkin(9) * t164 * t398 - pkin(2) * t332 - t165 * t387) + m(5) * (t125 * t285 + t126 * t286 + t154 * t387 + t319 * t76 + t320 * t77) + t562 * t88 + t563 * t89 + (t101 * t3 + t102 * t2 + t16 * t563 + t17 * t562 + t283 * t54 + t559 * t99) * m(7); (t200 / 0.2e1 - t151 / 0.2e1) * t52 + (t177 * t390 - t178 * t388) * qJ(4) + (t209 * t390 - t210 * t388) * qJD(4) + (-Ifges(6,5) * t349 - Ifges(6,6) * t350) * t507 + t428 - (-(t487 / 0.2e1 - t485 / 0.2e1) * t540 + (Ifges(4,1) / 0.2e1 - t429) * t322 + t538) * t540 + (-t350 / 0.2e1 + t225 / 0.2e1) * t95 + (-Ifges(6,1) * t226 - Ifges(6,4) * t225) * t519 + (-Ifges(6,4) * t226 - Ifges(6,2) * t225) * t521 + (-Ifges(6,5) * t226 - Ifges(6,6) * t225) * t508 + (-t349 / 0.2e1 + t226 / 0.2e1) * t96 + (-pkin(3) * t154 + (-t125 * t388 + t126 * t390) * qJD(4) + t409 * qJ(4) - t125 * t159 - t126 * t160 - t203 * t221) * m(5) + (Ifges(5,5) * t388 + Ifges(6,5) * t364 + Ifges(7,5) * t280 + Ifges(5,6) * t390 - Ifges(6,6) * t408 + Ifges(7,6) * t279) * t513 + t111 * (mrSges(6,1) * t408 + mrSges(6,2) * t364) + (Ifges(6,4) * t364 - Ifges(6,2) * t408) * t530 + (Ifges(6,1) * t364 - Ifges(6,4) * t408) * t531 - t408 * t533 + t409 * mrSges(5,3) + t400 * t322 + (-Ifges(6,1) * t349 - Ifges(6,4) * t350) * t518 + (-Ifges(6,4) * t349 - Ifges(6,2) * t350) * t520 + t481 - t482 + t154 * (-mrSges(5,1) * t390 + mrSges(5,2) * t388) + t385 * t48 + t329 * t11 + t296 * t78 + t297 * t79 + t54 * (-mrSges(7,1) * t279 + mrSges(7,2) * t280) - t220 * t277 - t160 * t209 - t159 * t210 - t193 * t110 + t179 * t28 + t180 * t29 - pkin(3) * t163 + (Ifges(7,4) * t280 + Ifges(7,2) * t279) * t534 + (Ifges(7,1) * t280 + Ifges(7,4) * t279) * t535 + t280 * t536 + t279 * t537 + t388 * t522 + (Ifges(7,1) * t199 + Ifges(7,4) * t200) * t523 + (Ifges(7,1) * t152 + Ifges(7,4) * t151) * t524 + (Ifges(7,4) * t199 + Ifges(7,2) * t200) * t525 + (Ifges(7,4) * t152 + Ifges(7,2) * t151) * t526 + t364 * t532 + (Ifges(5,1) * t388 + t491) * t516 + (Ifges(5,2) * t390 + t492) * t517 + t129 * t500 + (Ifges(7,5) * t199 + Ifges(7,6) * t200) * t510 + (Ifges(7,5) * t152 + Ifges(7,6) * t151) * t511 + (t199 / 0.2e1 - t152 / 0.2e1) * t53 + (-t12 * t408 - t13 * t364 + t443 * t45 + t46 * t573) * mrSges(6,3) + (-mrSges(6,1) * t573 - mrSges(6,2) * t443) * t162 + t541 * t56 + t544 * t158 + t545 * t157 + (t111 * t385 + t12 * t297 + t13 * t296 - t162 * t193 + t45 * t544 + t46 * t545) * m(6) + t441 * t221 + (mrSges(7,1) * t448 - mrSges(7,2) * t447) * t99 + (t16 * t447 - t17 * t448 + t2 * t279 - t280 * t3) * mrSges(7,3) + t550 * t89 + t551 * t88 + (t16 * t550 + t17 * t551 + t179 * t3 + t180 * t2 + t329 * t54 + t541 * t99) * m(7); t108 * t89 - t570 * t88 - t556 * t157 + t190 * t158 - t420 * t209 + t265 * t210 + t11 + t163 + t48 + (t108 * t16 - t17 * t570 + t54) * m(7) + (t190 * t45 - t46 * t556 + t111) * m(6) + (t125 * t265 - t126 * t420 + t154) * m(5); -t552 + t108 * t554 + t553 + (-Ifges(6,2) * t190 + t185 + t96) * t521 - t162 * (mrSges(6,1) * t190 + mrSges(6,2) * t556) - t45 * t157 + t46 * t158 + (-t190 * t56 + t396 * t28 + t392 * t29 + (-t392 * t89 + t396 * t88) * qJD(6) + (-t190 * t99 + t2 * t392 + t3 * t396 + (-t16 * t392 + t17 * t396) * qJD(6)) * m(7)) * pkin(5) - t19 * t88 - t18 * t89 + (t190 * t46 + t45 * t556) * mrSges(6,3) + t95 * t518 + (Ifges(6,1) * t556 - t490) * t519 + (Ifges(6,5) * t556 - Ifges(6,6) * t190) * t508 - m(7) * (t16 * t18 + t17 * t19) + t565; -t16 * t88 + t17 * t89 + t52 * t523 + t564 + t565 + t8;];
tauc  = t1(:);
