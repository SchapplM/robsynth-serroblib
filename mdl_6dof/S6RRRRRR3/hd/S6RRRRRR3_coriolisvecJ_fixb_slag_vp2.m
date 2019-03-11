% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 03:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRRRR3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR3_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR3_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR3_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR3_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:39:14
% EndTime: 2019-03-10 03:39:55
% DurationCPUTime: 20.50s
% Computational Cost: add. (31747->800), mult. (74266->1088), div. (0->0), fcn. (55409->10), ass. (0->389)
t429 = cos(qJ(2));
t573 = -pkin(8) - pkin(7);
t407 = t573 * t429;
t394 = qJD(1) * t407;
t423 = sin(qJ(3));
t373 = t423 * t394;
t424 = sin(qJ(2));
t405 = t573 * t424;
t393 = qJD(1) * t405;
t379 = qJD(2) * pkin(2) + t393;
t428 = cos(qJ(3));
t321 = t379 * t428 + t373;
t419 = qJD(2) + qJD(3);
t291 = -pkin(3) * t419 - t321;
t388 = t423 * t429 + t428 * t424;
t371 = t388 * qJD(1);
t422 = sin(qJ(4));
t427 = cos(qJ(4));
t339 = -t371 * t422 + t419 * t427;
t242 = -pkin(4) * t339 + t291;
t340 = t371 * t427 + t419 * t422;
t421 = sin(qJ(5));
t426 = cos(qJ(5));
t453 = t426 * t339 - t340 * t421;
t145 = -pkin(5) * t453 + t242;
t253 = t339 * t421 + t340 * t426;
t420 = sin(qJ(6));
t425 = cos(qJ(6));
t152 = t253 * t425 + t420 * t453;
t614 = pkin(11) * t453;
t478 = qJD(1) * t429;
t479 = qJD(1) * t424;
t370 = -t423 * t479 + t428 * t478;
t414 = -pkin(2) * t429 - pkin(1);
t403 = qJD(1) * t414;
t280 = -t370 * pkin(3) - t371 * pkin(9) + t403;
t481 = t428 * t394;
t322 = t423 * t379 - t481;
t292 = pkin(9) * t419 + t322;
t200 = t427 * t280 - t292 * t422;
t163 = -pkin(10) * t340 + t200;
t364 = qJD(4) - t370;
t142 = pkin(4) * t364 + t163;
t201 = t280 * t422 + t292 * t427;
t164 = pkin(10) * t339 + t201;
t160 = t426 * t164;
t88 = t142 * t421 + t160;
t70 = t88 + t614;
t505 = t420 * t70;
t354 = qJD(5) + t364;
t634 = pkin(11) * t253;
t158 = t421 * t164;
t87 = t426 * t142 - t158;
t69 = t87 - t634;
t61 = pkin(5) * t354 + t69;
t23 = t425 * t61 - t505;
t503 = t425 * t70;
t24 = t420 * t61 + t503;
t332 = t419 * t388;
t312 = t332 * qJD(1);
t386 = t423 * t424 - t428 * t429;
t331 = t419 * t386;
t311 = t331 * qJD(1);
t239 = qJD(4) * t339 - t311 * t427;
t240 = -qJD(4) * t340 + t311 * t422;
t112 = qJD(5) * t453 + t239 * t426 + t240 * t421;
t113 = -qJD(5) * t253 - t239 * t421 + t240 * t426;
t619 = -t253 * t420 + t425 * t453;
t37 = qJD(6) * t619 + t112 * t425 + t113 * t420;
t38 = -qJD(6) * t152 - t112 * t420 + t113 * t425;
t469 = Ifges(7,5) * t37 + Ifges(7,6) * t38 + Ifges(7,3) * t312;
t522 = Ifges(7,4) * t152;
t351 = qJD(6) + t354;
t548 = -t351 / 0.2e1;
t560 = -t152 / 0.2e1;
t473 = qJD(5) * t426;
t474 = qJD(5) * t421;
t477 = qJD(2) * t424;
t468 = pkin(2) * t477;
t622 = qJD(1) * t468;
t214 = pkin(3) * t312 + pkin(9) * t311 + t622;
t464 = qJD(2) * t573;
t452 = qJD(1) * t464;
t380 = t424 * t452;
t439 = t429 * t452;
t232 = qJD(3) * t321 + t428 * t380 + t423 * t439;
t95 = -qJD(4) * t201 + t427 * t214 - t232 * t422;
t59 = pkin(4) * t312 - pkin(10) * t239 + t95;
t475 = qJD(4) * t427;
t476 = qJD(4) * t422;
t94 = t422 * t214 + t427 * t232 + t280 * t475 - t292 * t476;
t68 = pkin(10) * t240 + t94;
t16 = t142 * t473 - t164 * t474 + t421 * t59 + t426 * t68;
t10 = pkin(11) * t113 + t16;
t17 = -qJD(5) * t88 - t421 * t68 + t426 * t59;
t7 = pkin(5) * t312 - pkin(11) * t112 + t17;
t3 = qJD(6) * t23 + t10 * t425 + t420 * t7;
t4 = -qJD(6) * t24 - t10 * t420 + t425 * t7;
t613 = t4 * mrSges(7,1) - t3 * mrSges(7,2);
t659 = t469 + t613 + (Ifges(7,1) * t619 - t522) * t560 + (Ifges(7,5) * t619 - Ifges(7,6) * t152) * t548 + (t152 * t24 + t23 * t619) * mrSges(7,3) - t145 * (mrSges(7,1) * t152 + mrSges(7,2) * t619);
t317 = pkin(3) * t371 - pkin(9) * t370;
t285 = pkin(2) * t479 + t317;
t327 = t393 * t428 + t373;
t222 = t427 * t285 - t327 * t422;
t489 = t370 * t427;
t451 = t371 * pkin(4) - pkin(10) * t489;
t410 = pkin(2) * t423 + pkin(9);
t529 = -pkin(10) - t410;
t456 = qJD(4) * t529;
t521 = pkin(2) * qJD(3);
t467 = t428 * t521;
t658 = -t422 * t467 + t427 * t456 - t222 - t451;
t226 = t427 * t317 - t321 * t422;
t572 = -pkin(10) - pkin(9);
t463 = qJD(4) * t572;
t657 = t427 * t463 - t226 - t451;
t223 = t422 * t285 + t427 * t327;
t490 = t370 * t422;
t470 = pkin(10) * t490;
t656 = -t422 * t456 - t427 * t467 + t223 - t470;
t227 = t422 * t317 + t427 * t321;
t655 = -t422 * t463 + t227 - t470;
t387 = t421 * t427 + t422 * t426;
t278 = t387 * t370;
t585 = qJD(4) + qJD(5);
t330 = t585 * t387;
t628 = t278 - t330;
t440 = t421 * t422 - t426 * t427;
t279 = t440 * t370;
t329 = t585 * t440;
t627 = t279 - t329;
t466 = Ifges(6,5) * t112 + Ifges(6,6) * t113 + Ifges(6,3) * t312;
t523 = Ifges(6,4) * t253;
t546 = -t354 / 0.2e1;
t554 = -t253 / 0.2e1;
t562 = -t619 / 0.2e1;
t144 = Ifges(7,4) * t619;
t81 = t152 * Ifges(7,1) + t351 * Ifges(7,5) + t144;
t575 = -t81 / 0.2e1;
t80 = Ifges(7,2) * t619 + Ifges(7,6) * t351 + t522;
t577 = -t80 / 0.2e1;
t599 = t17 * mrSges(6,1) - t16 * mrSges(6,2);
t620 = -Ifges(7,2) * t152 + t144;
t654 = -t152 * t577 + t466 + t599 + (Ifges(6,5) * t453 - Ifges(6,6) * t253) * t546 + (t253 * t88 + t453 * t87) * mrSges(6,3) - t242 * (mrSges(6,1) * t253 + mrSges(6,2) * t453) + t619 * t575 + (Ifges(6,1) * t453 - t523) * t554 + t620 * t562 + t659;
t540 = -t422 / 0.2e1;
t382 = t529 * t422;
t418 = t427 * pkin(10);
t383 = t410 * t427 + t418;
t314 = t421 * t382 + t426 * t383;
t597 = -qJD(5) * t314 + t656 * t421 + t426 * t658;
t596 = t382 * t473 - t383 * t474 + t421 * t658 - t656 * t426;
t404 = t572 * t422;
t406 = pkin(9) * t427 + t418;
t343 = t421 * t404 + t426 * t406;
t595 = -qJD(5) * t343 + t421 * t655 + t426 * t657;
t594 = t404 * t473 - t406 * t474 + t421 * t657 - t426 * t655;
t651 = t628 * pkin(11);
t650 = -t371 * pkin(5) - pkin(11) * t627;
t647 = t419 * Ifges(4,6) / 0.2e1;
t643 = t651 + t596;
t642 = t650 + t597;
t641 = t651 + t594;
t640 = t650 + t595;
t526 = Ifges(5,4) * t340;
t229 = t339 * Ifges(5,2) + t364 * Ifges(5,6) + t526;
t449 = mrSges(5,1) * t422 + mrSges(5,2) * t427;
t639 = t229 * t540 + t291 * t449;
t507 = t419 * Ifges(4,5);
t638 = t403 * mrSges(4,2) + t507 / 0.2e1;
t508 = t371 * Ifges(4,4);
t637 = t647 + t508 / 0.2e1 + t370 * Ifges(4,2) / 0.2e1;
t411 = pkin(4) * t426 + pkin(5);
t471 = qJD(6) * t425;
t472 = qJD(6) * t420;
t484 = t421 * t425;
t99 = -t163 * t421 - t160;
t71 = t99 - t614;
t100 = t426 * t163 - t158;
t72 = t100 - t634;
t632 = t420 * t72 - t425 * t71 - t411 * t472 + (-t421 * t471 + (-t420 * t426 - t484) * qJD(5)) * pkin(4);
t485 = t420 * t421;
t631 = -t420 * t71 - t425 * t72 + t411 * t471 + (-t421 * t472 + (t425 * t426 - t485) * qJD(5)) * pkin(4);
t323 = -t387 * t420 - t425 * t440;
t178 = qJD(6) * t323 - t329 * t425 - t330 * t420;
t194 = -t278 * t420 - t279 * t425;
t630 = t178 - t194;
t324 = t387 * t425 - t420 * t440;
t179 = -qJD(6) * t324 + t329 * t420 - t330 * t425;
t193 = -t278 * t425 + t279 * t420;
t629 = t179 - t193;
t416 = pkin(4) * t476;
t626 = -pkin(5) * t628 + t416;
t326 = t393 * t423 - t481;
t352 = pkin(4) * t490;
t625 = t423 * t521 - t326 - t352;
t624 = -t331 * t422 + t388 * t475;
t623 = -t422 * t95 + t427 * t94;
t249 = Ifges(6,4) * t453;
t621 = -Ifges(6,2) * t253 + t249;
t444 = Ifges(5,5) * t427 - Ifges(5,6) * t422;
t524 = Ifges(5,4) * t427;
t446 = -Ifges(5,2) * t422 + t524;
t525 = Ifges(5,4) * t422;
t448 = Ifges(5,1) * t427 - t525;
t338 = Ifges(5,4) * t339;
t230 = t340 * Ifges(5,1) + t364 * Ifges(5,5) + t338;
t482 = t427 * t230;
t549 = t340 / 0.2e1;
t618 = t364 * t444 / 0.2e1 + t448 * t549 + t339 * t446 / 0.2e1 + t482 / 0.2e1 + t639;
t615 = -t403 * mrSges(4,1) - t200 * mrSges(5,1) - t87 * mrSges(6,1) - t23 * mrSges(7,1) + t201 * mrSges(5,2) + t88 * mrSges(6,2) + t24 * mrSges(7,2) + t637;
t581 = t37 / 0.2e1;
t580 = t38 / 0.2e1;
t568 = t112 / 0.2e1;
t567 = t113 / 0.2e1;
t552 = t312 / 0.2e1;
t313 = t426 * t382 - t383 * t421;
t533 = pkin(11) * t387;
t267 = t313 - t533;
t381 = t440 * pkin(11);
t268 = -t381 + t314;
t175 = t267 * t425 - t268 * t420;
t612 = qJD(6) * t175 + t420 * t642 + t425 * t643;
t176 = t267 * t420 + t268 * t425;
t611 = -qJD(6) * t176 - t420 * t643 + t425 * t642;
t341 = t426 * t404 - t406 * t421;
t289 = t341 - t533;
t290 = -t381 + t343;
t203 = t289 * t425 - t290 * t420;
t610 = qJD(6) * t203 + t420 * t640 + t425 * t641;
t204 = t289 * t420 + t290 * t425;
t609 = -qJD(6) * t204 - t420 * t641 + t425 * t640;
t608 = mrSges(4,2) * t371;
t298 = t440 * t388;
t261 = t352 + t322;
t593 = -t261 + t626;
t592 = t625 + t626;
t320 = t386 * pkin(3) - t388 * pkin(9) + t414;
t344 = t405 * t423 - t407 * t428;
t243 = t427 * t320 - t344 * t422;
t186 = pkin(4) * t386 - t388 * t418 + t243;
t334 = t427 * t344;
t244 = t422 * t320 + t334;
t487 = t388 * t422;
t210 = -pkin(10) * t487 + t244;
t120 = t421 * t186 + t426 * t210;
t591 = Ifges(5,5) * t239 + Ifges(5,6) * t240;
t590 = t428 * t405 + t407 * t423;
t589 = t416 + t625;
t588 = -t200 * t422 + t201 * t427;
t587 = t95 * mrSges(5,1) - t94 * mrSges(5,2);
t140 = -mrSges(5,1) * t240 + mrSges(5,2) * t239;
t233 = qJD(3) * t322 + t380 * t423 - t428 * t439;
t586 = m(5) * t233 + t140;
t584 = m(4) / 0.2e1;
t583 = Ifges(7,4) * t581 + Ifges(7,2) * t580 + Ifges(7,6) * t552;
t582 = Ifges(7,1) * t581 + Ifges(7,4) * t580 + Ifges(7,5) * t552;
t579 = Ifges(6,4) * t568 + Ifges(6,2) * t567 + Ifges(6,6) * t552;
t578 = Ifges(6,1) * t568 + Ifges(6,4) * t567 + Ifges(6,5) * t552;
t576 = t80 / 0.2e1;
t574 = t81 / 0.2e1;
t571 = pkin(1) * mrSges(3,1);
t570 = pkin(1) * mrSges(3,2);
t136 = Ifges(6,2) * t453 + Ifges(6,6) * t354 + t523;
t566 = -t136 / 0.2e1;
t565 = t136 / 0.2e1;
t137 = Ifges(6,1) * t253 + Ifges(6,5) * t354 + t249;
t564 = -t137 / 0.2e1;
t563 = t137 / 0.2e1;
t561 = t619 / 0.2e1;
t559 = t152 / 0.2e1;
t558 = t239 / 0.2e1;
t557 = t240 / 0.2e1;
t556 = -t453 / 0.2e1;
t555 = t453 / 0.2e1;
t553 = t253 / 0.2e1;
t551 = -t339 / 0.2e1;
t550 = -t340 / 0.2e1;
t547 = t351 / 0.2e1;
t545 = t354 / 0.2e1;
t544 = -t364 / 0.2e1;
t543 = -t370 / 0.2e1;
t539 = t427 / 0.2e1;
t538 = m(4) * t403;
t536 = pkin(2) * t428;
t528 = mrSges(4,3) * t370;
t527 = Ifges(3,4) * t424;
t361 = Ifges(4,4) * t370;
t520 = t619 * Ifges(7,6);
t519 = t152 * Ifges(7,5);
t518 = t453 * Ifges(6,6);
t517 = t253 * Ifges(6,5);
t516 = t339 * Ifges(5,6);
t515 = t340 * Ifges(5,5);
t514 = t351 * Ifges(7,3);
t513 = t354 * Ifges(6,3);
t512 = t364 * Ifges(5,3);
t510 = t371 * mrSges(4,3);
t509 = t371 * Ifges(4,1);
t501 = Ifges(3,5) * qJD(2);
t500 = Ifges(3,6) * qJD(2);
t499 = qJD(2) * mrSges(3,1);
t498 = qJD(2) * mrSges(3,2);
t495 = t233 * t590;
t480 = -mrSges(4,1) * t419 - mrSges(5,1) * t339 + mrSges(5,2) * t340 + t510;
t465 = Ifges(5,3) * t312 + t591;
t413 = -pkin(4) * t427 - pkin(3);
t461 = t501 / 0.2e1;
t460 = -t500 / 0.2e1;
t119 = t426 * t186 - t210 * t421;
t238 = pkin(3) * t332 + pkin(9) * t331 + t468;
t398 = t424 * t464;
t399 = t429 * t464;
t256 = qJD(3) * t590 + t398 * t428 + t399 * t423;
t455 = t427 * t238 - t256 * t422;
t281 = pkin(4) * t487 - t590;
t450 = mrSges(5,1) * t427 - mrSges(5,2) * t422;
t447 = Ifges(5,1) * t422 + t524;
t445 = Ifges(5,2) * t427 + t525;
t443 = Ifges(5,5) * t422 + Ifges(5,6) * t427;
t297 = t387 * t388;
t105 = -pkin(11) * t297 + t120;
t98 = pkin(5) * t386 + pkin(11) * t298 + t119;
t42 = t105 * t425 + t420 * t98;
t41 = -t105 * t420 + t425 * t98;
t172 = mrSges(5,1) * t312 - mrSges(5,3) * t239;
t173 = -mrSges(5,2) * t312 + mrSges(5,3) * t240;
t442 = -t172 * t422 + t173 * t427;
t441 = -t200 * t427 - t201 * t422;
t212 = -t297 * t425 + t298 * t420;
t213 = -t297 * t420 - t298 * t425;
t350 = pkin(5) * t440 + t413;
t77 = t331 * t418 + pkin(4) * t332 + (-t334 + (pkin(10) * t388 - t320) * t422) * qJD(4) + t455;
t114 = t422 * t238 + t427 * t256 + t320 * t475 - t344 * t476;
t93 = -pkin(10) * t624 + t114;
t21 = t186 * t473 - t210 * t474 + t421 * t77 + t426 * t93;
t436 = t441 * mrSges(5,3);
t22 = -qJD(5) * t120 - t421 * t93 + t426 * t77;
t257 = qJD(3) * t344 + t398 * t423 - t428 * t399;
t171 = pkin(4) * t624 + t257;
t138 = -pkin(4) * t240 + t233;
t431 = m(5) * (qJD(4) * t441 + t623);
t123 = t239 * Ifges(5,4) + t240 * Ifges(5,2) + t312 * Ifges(5,6);
t124 = Ifges(5,1) * t239 + Ifges(5,4) * t240 + Ifges(5,5) * t312;
t135 = t513 + t517 + t518;
t228 = t512 + t515 + t516;
t287 = t361 + t507 + t509;
t60 = -pkin(5) * t113 + t138;
t79 = t514 + t519 + t520;
t430 = (-Ifges(6,4) * t279 - Ifges(6,2) * t278) * t556 + (Ifges(7,4) * t194 + Ifges(7,2) * t193) * t562 + (-t450 - mrSges(4,1)) * t233 + (-Ifges(6,5) * t279 - Ifges(6,6) * t278) * t546 + (-Ifges(6,1) * t279 - Ifges(6,4) * t278) * t554 + (t482 + t287 + t361) * t543 + (Ifges(7,1) * t194 + Ifges(7,4) * t193) * t560 - (Ifges(4,1) * t370 + t135 + t228 - t508 + t79) * t371 / 0.2e1 + (Ifges(6,5) * t387 + Ifges(7,5) * t324 - Ifges(6,6) * t440 + Ifges(7,6) * t323 + t443) * t552 + t138 * (mrSges(6,1) * t440 + mrSges(6,2) * t387) + (Ifges(6,4) * t387 - Ifges(6,2) * t440) * t567 + (Ifges(6,1) * t387 - Ifges(6,4) * t440) * t568 - t440 * t579 + (-Ifges(6,5) * t329 - Ifges(6,6) * t330) * t545 + (-Ifges(6,1) * t329 - Ifges(6,4) * t330) * t553 + (-Ifges(6,4) * t329 - Ifges(6,2) * t330) * t555 + (-t16 * t440 - t17 * t387 - t627 * t87 + t628 * t88) * mrSges(6,3) + (-mrSges(6,1) * t628 + mrSges(6,2) * t627) * t242 + (t444 * t544 + t446 * t551 + t448 * t550 - t638 - t639) * t370 + (t200 * t489 + t201 * t490 + t623) * mrSges(5,3) + t321 * t528 + t618 * qJD(4) + (-mrSges(7,1) * t629 + mrSges(7,2) * t630) * t145 + (-t23 * t630 + t24 * t629 + t3 * t323 - t324 * t4) * mrSges(7,3) + t123 * t539 + (Ifges(7,5) * t178 + Ifges(7,6) * t179) * t547 + t445 * t557 + t447 * t558 + (Ifges(7,1) * t178 + Ifges(7,4) * t179) * t559 + (Ifges(7,4) * t178 + Ifges(7,2) * t179) * t561 - t329 * t563 + t422 * t124 / 0.2e1 + (Ifges(5,5) * t550 + Ifges(6,5) * t554 + Ifges(7,5) * t560 - Ifges(4,2) * t543 + Ifges(5,6) * t551 + Ifges(6,6) * t556 + Ifges(7,6) * t562 + Ifges(5,3) * t544 + Ifges(6,3) * t546 + Ifges(7,3) * t548 + t615 + t647) * t371 + (Ifges(7,5) * t194 + Ifges(7,6) * t193) * t548 - t279 * t564 - t330 * t565 - t278 * t566 + t178 * t574 + t194 * t575 + t179 * t576 + t193 * t577 + t387 * t578 + (Ifges(7,4) * t324 + Ifges(7,2) * t323) * t580 + (Ifges(7,1) * t324 + Ifges(7,4) * t323) * t581 + t324 * t582 + t323 * t583 - t232 * mrSges(4,2) - Ifges(4,5) * t311 - Ifges(4,6) * t312 + t60 * (-mrSges(7,1) * t323 + mrSges(7,2) * t324);
t415 = Ifges(3,4) * t478;
t402 = mrSges(3,3) * t478 - t498;
t401 = -mrSges(3,3) * t479 + t499;
t400 = t413 - t536;
t369 = Ifges(3,1) * t479 + t415 + t501;
t368 = t500 + (Ifges(3,2) * t429 + t527) * qJD(1);
t367 = pkin(4) * t484 + t411 * t420;
t366 = -pkin(4) * t485 + t411 * t425;
t348 = -mrSges(4,2) * t419 + t528;
t347 = t350 - t536;
t316 = -mrSges(4,1) * t370 + t608;
t271 = mrSges(5,1) * t364 - mrSges(5,3) * t340;
t270 = -mrSges(5,2) * t364 + mrSges(5,3) * t339;
t217 = pkin(5) * t297 + t281;
t216 = mrSges(6,1) * t354 - mrSges(6,3) * t253;
t215 = -mrSges(6,2) * t354 + mrSges(6,3) * t453;
t199 = pkin(4) * t340 + pkin(5) * t253;
t157 = -mrSges(6,1) * t453 + mrSges(6,2) * t253;
t134 = mrSges(7,1) * t351 - mrSges(7,3) * t152;
t133 = -mrSges(7,2) * t351 + mrSges(7,3) * t619;
t132 = t298 * t585 + t387 * t331;
t131 = -t330 * t388 + t331 * t440;
t115 = -qJD(4) * t244 + t455;
t97 = -mrSges(6,2) * t312 + mrSges(6,3) * t113;
t96 = mrSges(6,1) * t312 - mrSges(6,3) * t112;
t92 = -mrSges(7,1) * t619 + mrSges(7,2) * t152;
t86 = -pkin(5) * t132 + t171;
t50 = -qJD(6) * t213 - t131 * t420 + t132 * t425;
t49 = qJD(6) * t212 + t131 * t425 + t132 * t420;
t48 = -mrSges(6,1) * t113 + mrSges(6,2) * t112;
t32 = -mrSges(7,2) * t312 + mrSges(7,3) * t38;
t31 = mrSges(7,1) * t312 - mrSges(7,3) * t37;
t26 = t425 * t69 - t505;
t25 = -t420 * t69 - t503;
t19 = pkin(11) * t132 + t21;
t18 = pkin(5) * t332 - pkin(11) * t131 + t22;
t13 = -mrSges(7,1) * t38 + mrSges(7,2) * t37;
t6 = -qJD(6) * t42 + t18 * t425 - t19 * t420;
t5 = qJD(6) * t41 + t18 * t420 + t19 * t425;
t1 = [(Ifges(7,1) * t213 + Ifges(7,4) * t212) * t581 + (t321 * t331 - t322 * t332) * mrSges(4,3) + m(4) * (t232 * t344 + t256 * t322 - t257 * t321 - t495) + m(5) * (t114 * t201 + t115 * t200 + t243 * t95 + t244 * t94 + t257 * t291 - t495) - (t509 / 0.2e1 + t287 / 0.2e1 + t361 / 0.2e1 + t436 + t618 + t638) * t331 + (t212 * t3 - t213 * t4 - t23 * t49 + t24 * t50) * mrSges(7,3) + (t124 * t539 + t123 * t540 + t448 * t558 + t446 * t557 - Ifges(4,1) * t311 + (mrSges(4,3) + t449) * t233 + (-t422 * t94 - t427 * t95) * mrSges(5,3) + (-t427 * t229 / 0.2e1 + t230 * t540 + t291 * t450 + t445 * t551 + t447 * t550 + t443 * t544 - t588 * mrSges(5,3)) * qJD(4)) * t388 + (mrSges(4,1) * t622 - mrSges(4,3) * t232 + Ifges(4,4) * t311 + Ifges(6,5) * t568 + Ifges(7,5) * t581 + Ifges(6,6) * t567 + Ifges(7,6) * t580 + (Ifges(6,3) + Ifges(7,3)) * t552 + t587 + t599 + t613) * t386 + (-Ifges(6,4) * t298 - Ifges(6,2) * t297) * t567 + (-Ifges(6,1) * t298 - Ifges(6,4) * t297) * t568 + (-t131 * t87 + t132 * t88 - t16 * t297 + t17 * t298) * mrSges(6,3) + (-Ifges(6,5) * t298 + Ifges(7,5) * t213 - Ifges(6,6) * t297 + Ifges(7,6) * t212 + t388 * t444) * t552 + t138 * (mrSges(6,1) * t297 - mrSges(6,2) * t298) + t480 * t257 + (Ifges(7,4) * t213 + Ifges(7,2) * t212) * t580 + (t515 / 0.2e1 + t516 / 0.2e1 + t519 / 0.2e1 + t517 / 0.2e1 + t518 / 0.2e1 + t520 / 0.2e1 + t512 / 0.2e1 + t513 / 0.2e1 + t514 / 0.2e1 - t615 + t79 / 0.2e1 + t135 / 0.2e1 + t228 / 0.2e1 - t637) * t332 + (t469 + t466 + t465 + t591) * t386 / 0.2e1 + t86 * t92 + (-pkin(7) * t402 - t368 / 0.2e1 + t460 + (-0.2e1 * t571 - 0.3e1 / 0.2e1 * t527 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t429) * qJD(1) + (t316 + 0.2e1 * t538 + t608) * pkin(2)) * t477 + (Ifges(6,5) * t131 + Ifges(6,6) * t132) * t545 + (Ifges(7,5) * t49 + Ifges(7,6) * t50) * t547 + ((Ifges(5,3) / 0.2e1 + Ifges(4,2)) * t386 - Ifges(4,4) * t388 - t344 * mrSges(4,3) + t414 * mrSges(4,1)) * t312 + (Ifges(6,1) * t131 + Ifges(6,4) * t132) * t553 + (Ifges(6,4) * t131 + Ifges(6,2) * t132) * t555 + (Ifges(7,1) * t49 + Ifges(7,4) * t50) * t559 + (Ifges(7,4) * t49 + Ifges(7,2) * t50) * t561 + (-pkin(7) * t401 + t369 / 0.2e1 + t461 + (-0.2e1 * t570 + 0.3e1 / 0.2e1 * Ifges(3,4) * t429) * qJD(1)) * t429 * qJD(2) + t41 * t31 + t42 * t32 - (mrSges(4,2) * t414 - mrSges(4,3) * t590) * t311 - t590 * t140 + t119 * t96 + t120 * t97 + t5 * t133 + t6 * t134 + t145 * (-mrSges(7,1) * t50 + mrSges(7,2) * t49) + t171 * t157 + m(6) * (t119 * t17 + t120 * t16 + t138 * t281 + t171 * t242 + t21 * t88 + t22 * t87) + m(7) * (t145 * t86 + t217 * t60 + t23 * t6 + t24 * t5 + t3 * t42 + t4 * t41) + t131 * t563 + t132 * t565 + t49 * t574 + t50 * t576 - t298 * t578 - t297 * t579 + t213 * t582 + t212 * t583 + t60 * (-mrSges(7,1) * t212 + mrSges(7,2) * t213) + t21 * t215 + t22 * t216 + t217 * t13 + t242 * (-mrSges(6,1) * t132 + mrSges(6,2) * t131) + t243 * t172 + t244 * t173 + t114 * t270 + t115 * t271 + t281 * t48 + t256 * t348; ((t311 * t428 - t312 * t423) * mrSges(4,3) + (t480 * t423 + (t270 * t427 - t271 * t422 + t348) * t428) * qJD(3)) * pkin(2) - t480 * t326 - m(5) * (t200 * t222 + t201 * t223 + t291 * t326) + t586 * (-pkin(3) - t536) + 0.2e1 * ((t232 * t423 - t233 * t428) * t584 + (m(5) * (t291 * t423 + t428 * t588) / 0.2e1 + (-t321 * t423 + t322 * t428) * t584) * qJD(3)) * pkin(2) + t589 * t157 + t322 * t510 + t596 * t215 + (t138 * t400 + t16 * t314 + t17 * t313 + t242 * t589 + t596 * t88 + t597 * t87) * m(6) + t597 * t216 + t592 * t92 + ((t461 - t415 / 0.2e1 - t369 / 0.2e1 + qJD(1) * t570 + (t401 - t499) * pkin(7)) * t429 + (t460 + t368 / 0.2e1 + (t571 + t527 / 0.2e1 + (-Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t429) * qJD(1) + (t402 + t498) * pkin(7) + (-t316 - t538) * pkin(2)) * t424) * qJD(1) + t430 + t400 * t48 + t175 * t31 + t176 * t32 - m(4) * (-t321 * t326 + t322 * t327) + t611 * t134 + t612 * t133 + (t145 * t592 + t175 * t4 + t176 * t3 + t23 * t611 + t24 * t612 + t347 * t60) * m(7) - t223 * t270 - t222 * t271 + t313 * t96 + t314 * t97 + t347 * t13 - t327 * t348 + t436 * qJD(4) + ((-t270 * t422 - t271 * t427) * qJD(4) + t431 + t442) * t410; ((-t200 * mrSges(5,3) - pkin(9) * t271) * t427 + (-mrSges(5,3) * t201 + pkin(4) * t157 - pkin(9) * t270) * t422) * qJD(4) + t593 * t92 - m(5) * (t200 * t226 + t201 * t227 + t291 * t322) + pkin(9) * t431 + (-t480 + t510) * t322 + t430 + t413 * t48 + t442 * pkin(9) + t595 * t216 + t594 * t215 + t203 * t31 + t204 * t32 - t261 * t157 - t227 * t270 + t610 * t133 - t226 * t271 + t609 * t134 + t341 * t96 + t343 * t97 - t321 * t348 + t350 * t13 - t586 * pkin(3) + (t145 * t593 + t203 * t4 + t204 * t3 + t23 * t609 + t24 * t610 + t350 * t60) * m(7) + (t138 * t413 + t16 * t343 + t17 * t341 + t594 * t88 + t595 * t87 + (-t261 + t416) * t242) * m(6); -m(6) * (t100 * t88 + t87 * t99) + t453 * t564 + (t200 * t339 + t201 * t340) * mrSges(5,3) + t587 + (-t340 * t157 + t421 * t97 + t426 * t96 + (t215 * t426 - t216 * t421) * qJD(5) + (t16 * t421 + t17 * t426 - t242 * t340 + t473 * t88 - t474 * t87) * m(6)) * pkin(4) + t465 + t621 * t556 - t253 * t566 + t631 * t133 + t632 * t134 + (-t145 * t199 + t23 * t632 + t24 * t631 + t3 * t367 + t366 * t4) * m(7) + (Ifges(5,5) * t339 - Ifges(5,6) * t340) * t544 + t229 * t549 + (Ifges(5,1) * t339 - t526) * t550 + (-Ifges(5,2) * t340 + t230 + t338) * t551 + t654 - t199 * t92 - t100 * t215 - t99 * t216 - t200 * t270 + t201 * t271 - t291 * (mrSges(5,1) * t340 + mrSges(5,2) * t339) + t366 * t31 + t367 * t32; -m(7) * (t23 * t25 + t24 * t26) + (t137 + t621) * t556 + (-t253 * t92 + t425 * t31 + t420 * t32 + (t133 * t425 - t134 * t420) * qJD(6) + (-t145 * t253 - t23 * t472 + t24 * t471 + t3 * t420 + t4 * t425) * m(7)) * pkin(5) + t136 * t553 - t26 * t133 - t25 * t134 - t87 * t215 + t88 * t216 + t654; t80 * t559 - t23 * t133 + t24 * t134 + (t620 + t81) * t562 + t659;];
tauc  = t1(:);
