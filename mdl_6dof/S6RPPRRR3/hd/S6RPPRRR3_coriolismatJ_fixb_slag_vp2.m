% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2]';
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
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPPRRR3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR3_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR3_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR3_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR3_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR3_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR3_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:22:58
% EndTime: 2019-03-09 02:23:10
% DurationCPUTime: 7.94s
% Computational Cost: add. (13125->536), mult. (26324->752), div. (0->0), fcn. (25076->8), ass. (0->293)
t353 = sin(qJ(4));
t531 = t353 * pkin(8);
t355 = cos(qJ(4));
t535 = pkin(4) * t355;
t324 = t531 + t535;
t354 = cos(qJ(5));
t309 = t354 * t324;
t336 = -cos(pkin(10)) * pkin(1) - pkin(2) - pkin(7);
t352 = sin(qJ(5));
t417 = -t336 * t352 + pkin(5);
t474 = t353 * t354;
t183 = pkin(9) * t474 + t355 * t417 + t309;
t477 = t352 * t355;
t217 = -t336 * t477 + t309;
t472 = t354 * t355;
t218 = t352 * t324 + t336 * t472;
t539 = sin(qJ(6));
t540 = cos(qJ(6));
t383 = t352 * t540 + t354 * t539;
t275 = t383 * t353;
t401 = -mrSges(7,2) * t355 + t275 * mrSges(7,3);
t594 = t539 * t352 - t540 * t354;
t277 = t594 * t353;
t402 = mrSges(7,1) * t355 - t277 * mrSges(7,3);
t478 = t352 * t353;
t424 = t478 / 0.2e1;
t456 = pkin(5) * t539;
t457 = pkin(5) * t540;
t572 = m(7) * pkin(5);
t461 = t572 / 0.2e1;
t541 = t355 / 0.2e1;
t571 = -mrSges(6,2) / 0.2e1;
t192 = pkin(9) * t478 + t218;
t382 = -t183 * t540 + t192 * t539;
t384 = t183 * t539 + t192 * t540;
t553 = -t277 / 0.2e1;
t558 = -t275 / 0.2e1;
t394 = Ifges(7,5) * t553 + Ifges(7,6) * t558;
t542 = -t355 / 0.2e1;
t585 = -t382 * mrSges(7,1) / 0.2e1 - t384 * mrSges(7,2) / 0.2e1 - Ifges(7,3) * t542 - t394;
t619 = Ifges(6,3) * t541 + t217 * mrSges(6,1) / 0.2e1 + t218 * t571 - Ifges(6,5) * t474 / 0.2e1 + Ifges(6,6) * t424 + t402 * t457 / 0.2e1 + t401 * t456 / 0.2e1 + (t539 ^ 2 + t540 ^ 2) * t183 * t461 + t585;
t497 = t354 * mrSges(6,1);
t498 = t352 * mrSges(6,2);
t500 = t383 * mrSges(7,2);
t501 = t594 * mrSges(7,1);
t598 = -t501 / 0.2e1 - t500 / 0.2e1;
t618 = -t498 / 0.2e1 + t497 / 0.2e1 + (t383 * t539 - t540 * t594) * t461 + t598;
t338 = sin(pkin(10)) * pkin(1) + qJ(3);
t532 = t353 * pkin(4);
t290 = -pkin(8) * t355 + t338 + t532;
t273 = t354 * t290;
t458 = pkin(9) * t472;
t167 = t353 * t417 + t273 - t458;
t194 = t290 * t352 + t336 * t474;
t179 = -pkin(9) * t477 + t194;
t435 = t539 * t179;
t83 = t167 * t540 - t435;
t193 = -t336 * t478 + t273;
t178 = t193 - t458;
t88 = t178 * t540 - t435;
t614 = t83 - t88;
t569 = -pkin(9) - pkin(8);
t323 = t569 * t352;
t325 = t569 * t354;
t219 = -t323 * t539 + t325 * t540;
t406 = t540 * t323 + t325 * t539;
t465 = -Ifges(7,5) * t594 - Ifges(7,6) * t383;
t41 = t219 * mrSges(7,1) - t406 * mrSges(7,2) + t465;
t617 = t41 * qJD(6);
t341 = t355 * t353;
t350 = t352 ^ 2;
t351 = t354 ^ 2;
t464 = t350 + t351;
t413 = t464 * t355;
t274 = t594 * t355;
t482 = t277 * t274;
t276 = t383 * t355;
t483 = t275 * t276;
t66 = m(7) * (t341 - t482 - t483) - m(6) * (t353 * t413 - t341);
t616 = t66 * qJD(4);
t168 = -t274 * mrSges(7,1) - t276 * mrSges(7,2);
t503 = t276 * mrSges(7,3);
t224 = -t353 * mrSges(7,2) - t503;
t225 = t353 * mrSges(7,1) + t274 * mrSges(7,3);
t554 = t277 / 0.2e1;
t359 = t168 * t542 + t224 * t558 + t225 * t554 + (-t482 / 0.2e1 - t483 / 0.2e1) * mrSges(7,3);
t575 = -m(7) / 0.2e1;
t615 = pkin(5) * t575;
t438 = t540 * t179;
t84 = t167 * t539 + t438;
t87 = -t178 * t539 - t438;
t603 = t84 + t87;
t349 = Ifges(6,4) * t354;
t319 = Ifges(6,1) * t352 + t349;
t516 = Ifges(6,2) * t352;
t597 = t349 - t516;
t612 = t319 + t597;
t591 = t498 - t497;
t291 = t591 * t355;
t533 = pkin(5) * t354;
t574 = m(7) / 0.2e1;
t579 = t355 ^ 2;
t611 = (-t275 * t603 + t277 * t614 - t533 * t579) * t574 - t291 * t542 + t359;
t499 = t383 * mrSges(7,3);
t442 = -t499 / 0.2e1;
t610 = t442 + t499 / 0.2e1;
t519 = Ifges(7,4) * t383;
t202 = -Ifges(7,2) * t594 + t519;
t203 = -Ifges(7,1) * t594 - t519;
t421 = -t203 / 0.4e1 + t202 / 0.4e1;
t562 = t224 / 0.2e1;
t430 = t406 * t562;
t552 = t291 / 0.2e1;
t608 = t430 + t219 * t225 / 0.2e1 + pkin(4) * t552 + t421 * t274;
t607 = qJD(6) * t168;
t415 = t277 * mrSges(7,1) + t275 * mrSges(7,2);
t606 = qJD(6) * t415;
t570 = -mrSges(7,3) / 0.2e1;
t604 = -t319 / 0.4e1;
t534 = pkin(5) * t352;
t316 = mrSges(6,1) * t352 + mrSges(6,2) * t354;
t601 = t336 * t316;
t249 = Ifges(7,4) * t276;
t162 = -t274 * Ifges(7,1) + Ifges(7,5) * t353 - t249;
t169 = Ifges(7,2) * t274 - t249;
t599 = t162 + t169;
t311 = t353 * mrSges(6,1) - mrSges(6,3) * t472;
t473 = t354 * t311;
t455 = mrSges(6,3) * t477;
t310 = -t353 * mrSges(6,2) - t455;
t480 = t352 * t310;
t596 = -t473 / 0.2e1 - t480 / 0.2e1;
t484 = t219 * t274;
t485 = t406 * t276;
t595 = -t484 / 0.2e1 + t485 / 0.2e1;
t592 = -t217 * t352 + t218 * t354;
t581 = t276 ^ 2;
t582 = t274 ^ 2;
t590 = t581 + t582;
t504 = t276 * mrSges(7,1);
t507 = t274 * mrSges(7,2);
t467 = -t504 / 0.2e1 + t507 / 0.2e1;
t502 = t277 * mrSges(7,2);
t505 = t275 * mrSges(7,1);
t468 = t505 / 0.2e1 - t502 / 0.2e1;
t588 = mrSges(6,3) * t464 - mrSges(5,2);
t496 = t355 * mrSges(6,3);
t586 = t596 - t464 * t496 / 0.2e1;
t506 = t274 * Ifges(7,4);
t161 = -t276 * Ifges(7,2) + t353 * Ifges(7,6) - t506;
t170 = -Ifges(7,1) * t276 + t506;
t199 = mrSges(7,1) * t383 - mrSges(7,2) * t594;
t416 = -t336 + t534;
t279 = t416 * t355;
t544 = t353 / 0.4e1;
t342 = -pkin(4) - t533;
t548 = t342 / 0.2e1;
t584 = t279 * t199 / 0.2e1 + t168 * t548 + t465 * t544 - (-t170 / 0.4e1 + t161 / 0.4e1) * t383;
t580 = t353 ^ 2;
t578 = 0.2e1 * qJD(4);
t577 = m(6) / 0.2e1;
t576 = m(6) / 0.4e1;
t573 = -pkin(8) / 0.2e1;
t200 = t500 + t501;
t566 = t200 / 0.2e1;
t563 = -t406 / 0.2e1;
t561 = -t225 / 0.2e1;
t560 = -t274 / 0.2e1;
t559 = t274 / 0.2e1;
t557 = t275 / 0.2e1;
t556 = -t276 / 0.2e1;
t551 = -t594 / 0.2e1;
t520 = Ifges(6,4) * t352;
t320 = Ifges(6,1) * t354 - t520;
t549 = t320 / 0.4e1;
t547 = -t352 / 0.2e1;
t546 = -t353 / 0.2e1;
t545 = t353 / 0.2e1;
t108 = -t274 * t383 + t276 * t594;
t537 = m(7) * t108;
t536 = m(7) * t279;
t530 = t83 * mrSges(7,2);
t529 = t84 * mrSges(7,1);
t528 = t87 * mrSges(7,1);
t527 = t88 * mrSges(7,2);
t518 = Ifges(6,5) * t353;
t348 = Ifges(6,5) * t354;
t515 = Ifges(6,6) * t352;
t514 = Ifges(6,6) * t353;
t495 = t83 * t594;
t494 = t84 * t274;
t493 = t591 - mrSges(5,1);
t492 = mrSges(7,3) * qJD(4);
t13 = t353 * t586 + t611 - t618;
t491 = t13 * qJD(1);
t390 = t224 * t556 + t225 * t559;
t419 = t351 / 0.2e1 + t350 / 0.2e1;
t14 = (t274 * t614 - t603 * t276) * t575 + (-t168 / 0.2e1 + t552) * t353 + (t582 / 0.2e1 + t581 / 0.2e1) * mrSges(7,3) + (t480 / 0.2e1 + (t311 / 0.2e1 + t546 * t572) * t354 + t419 * t496) * t355 - t390;
t490 = t14 * qJD(1);
t16 = t359 - t598;
t489 = t16 * qJD(1);
t19 = t168 * t545 + t570 * t590 + t390;
t488 = t19 * qJD(1);
t29 = t353 * mrSges(5,1) + t355 * mrSges(5,2) + t383 * t224 - t594 * t225 + t480 + t473 + mrSges(4,3) + m(7) * (t383 * t84 - t495) + m(6) * (t193 * t354 + t194 * t352) + 0.4e1 * (m(5) / 0.4e1 + m(4) / 0.4e1) * t338;
t481 = t29 * qJD(1);
t317 = Ifges(6,2) * t354 + t520;
t479 = t352 * t317;
t476 = t353 * t199;
t475 = t353 * t316;
t471 = t355 * t199;
t466 = -Ifges(7,5) * t276 + Ifges(7,6) * t274;
t463 = qJD(4) * t353;
t462 = qJD(4) * t355;
t460 = -0.1e1 + t464;
t52 = (-t275 ^ 2 - t277 ^ 2 + t590) * t574 + 0.2e1 * (t460 * t576 - m(7) / 0.4e1) * t579 + 0.2e1 * (-m(6) * t460 / 0.4e1 + m(7) / 0.4e1) * t580;
t459 = t52 * qJD(4);
t454 = -t537 / 0.2e1;
t453 = t537 / 0.2e1;
t451 = -Ifges(6,2) / 0.4e1 + Ifges(6,1) / 0.4e1;
t450 = -t83 / 0.2e1 + t88 / 0.2e1;
t449 = t84 / 0.2e1 + t87 / 0.2e1;
t440 = t476 / 0.2e1 + t610 * t274;
t439 = -t471 / 0.2e1 - t610 * t277;
t437 = t540 * t276;
t434 = t539 * t274;
t426 = -t601 / 0.2e1;
t295 = Ifges(7,4) * t594;
t201 = -Ifges(7,2) * t383 - t295;
t204 = Ifges(7,1) * t383 - t295;
t420 = t204 / 0.4e1 + t201 / 0.4e1;
t418 = Ifges(5,4) - t348;
t414 = t348 - t515;
t412 = mrSges(7,3) * t457;
t411 = mrSges(7,3) * t456;
t408 = mrSges(6,3) * t419;
t407 = qJD(4) * (t200 + t493);
t404 = -t504 + t507;
t403 = t502 - t505;
t400 = Ifges(7,1) * t277 + Ifges(7,4) * t275;
t398 = Ifges(7,4) * t277 + Ifges(7,2) * t275;
t278 = t416 * t353;
t378 = t475 / 0.2e1 + t311 * t547 + t354 * t310 / 0.2e1;
t389 = m(6) * t592;
t11 = t224 * t553 + t225 * t558 + (t274 * t557 - t276 * t554) * mrSges(7,3) + (t378 + t468) * t353 - m(6) * (t193 * t478 - t194 * t474 + t336 * t580) / 0.2e1 + (-t274 * t384 + t83 * t275 + t276 * t382 + t84 * t277 - t278 * t353) * t575 + 0.2e1 * (-t389 / 0.4e1 - t536 / 0.4e1 + t336 * t355 * t576) * t355;
t4 = -t84 * t401 - t83 * t402 - t279 * t403 - t278 * t404 - t193 * (mrSges(6,1) * t355 + mrSges(6,3) * t474) - t338 * (t355 * mrSges(5,1) - t353 * mrSges(5,2)) - t218 * t310 - t217 * t311 - t194 * (-mrSges(6,2) * t355 + mrSges(6,3) * t478) + t400 * t559 + t276 * t398 / 0.2e1 + t162 * t553 + t161 * t558 - t384 * t224 + t382 * t225 - t341 * t601 - m(6) * (t193 * t217 + t194 * t218) - m(7) * (-t279 * t278 - t382 * t83 + t384 * t84) + ((-t418 - t515) * t353 + t394) * t353 + (-t336 * t475 + Ifges(7,5) * t274 + Ifges(7,6) * t276 + t418 * t355 + (m(6) * t336 ^ 2 + t351 * Ifges(6,1) + Ifges(5,1) - Ifges(5,2) - Ifges(6,3) - Ifges(7,3)) * t353 + (Ifges(6,6) * t355 + (t516 - 0.2e1 * t349) * t353) * t352) * t355;
t397 = -t4 * qJD(1) - t11 * qJD(2);
t386 = mrSges(7,3) * t494 + t279 * t168 + t466 * t545;
t392 = pkin(5) * t404;
t5 = t392 * t472 - m(7) * (t83 * t87 + t84 * t88) - t88 * t224 - t87 * t225 + t194 * t311 + (t170 / 0.2e1 - t161 / 0.2e1) * t274 - (-t162 / 0.2e1 - t169 / 0.2e1 + t83 * mrSges(7,3)) * t276 + (-t336 * t291 + (-Ifges(6,4) * t477 + t518) * t352 + (t514 - pkin(5) * t536 + t194 * mrSges(6,3) + Ifges(6,4) * t472 + (Ifges(6,1) - Ifges(6,2)) * t477) * t354) * t355 - t386 + (-t310 - t455) * t193;
t396 = -t5 * qJD(1) - t14 * qJD(2);
t8 = -t225 * t84 + t161 * t559 + t170 * t560 + (t224 + t503) * t83 + t386 + t599 * t556;
t395 = t8 * qJD(1) + t19 * qJD(2);
t388 = m(7) * (-t219 * t383 - t406 * t594);
t385 = -t392 / 0.2e1;
t356 = (-t278 * t575 + (-t193 * t352 + t194 * t354 - 0.2e1 * t336 * t353) * t577 + t378) * t355 + (-t217 * t478 + t218 * t474) * t577 + (t275 * t382 - t83 * t276 - t277 * t384 + t279 * t353 - t494) * t574 + t224 * t560 + t225 * t556 + t404 * t546;
t10 = -t388 / 0.2e1 + t356;
t380 = t10 * qJD(1) + t52 * qJD(2) - qJD(3) * t66;
t379 = -t11 * qJD(1) + t66 * qJD(2) + t52 * qJD(3);
t364 = (t540 * t562 + t539 * t561 + (t434 / 0.2e1 + t437 / 0.2e1) * mrSges(7,3)) * pkin(5);
t18 = -mrSges(7,1) * t449 + mrSges(7,2) * t450 + t364;
t312 = (mrSges(7,1) * t539 + mrSges(7,2) * t540) * pkin(5);
t42 = (t563 + t406 / 0.2e1) * mrSges(7,2);
t377 = -t18 * qJD(1) - t42 * qJD(4) + t312 * qJD(5);
t373 = t614 * t219 + t603 * t406;
t1 = t584 + t586 * pkin(8) + t354 * (t355 * t320 + t518) / 0.4e1 - t352 * (t355 * t597 + t514) / 0.4e1 - t599 * t594 / 0.4e1 - (t204 + t201) * t276 / 0.4e1 + (-t317 / 0.2e1 + t549 + (t342 * t574 + t566) * pkin(5)) * t472 + (t279 * t534 + t373) * t574 + t603 * t442 + t414 * t544 + t352 * t385 + t355 * t426 - t495 * t570 + (t551 * t88 + t595) * mrSges(7,3) + t608 + (t604 - t612 / 0.4e1) * t477 - t619;
t34 = t342 * t199 - (-t203 / 0.2e1 + t202 / 0.2e1) * t383 - (t204 / 0.2e1 + t201 / 0.2e1) * t594;
t20 = -t479 / 0.2e1 - pkin(4) * t316 + t352 * t320 / 0.2e1 + (t597 / 0.2e1 + t319 / 0.2e1) * t354 + t34 + (m(7) * t342 + t200) * t534;
t363 = (t275 * t540 + t277 * t539) * t461 + mrSges(6,1) * t424 + mrSges(6,2) * t474 / 0.2e1 + t468;
t366 = t478 * t615 - t475 / 0.2e1;
t30 = -t476 / 0.2e1 + t363 + t366;
t362 = (-t437 - t434) * t461 - mrSges(6,1) * t477 / 0.2e1 + t472 * t571 + t467;
t365 = t316 * t541 - t477 * t615;
t31 = t471 / 0.2e1 + t362 + t365;
t370 = t1 * qJD(1) - t30 * qJD(2) - t31 * qJD(3) + t20 * qJD(4);
t37 = -t468 + t440;
t39 = t439 - t467;
t361 = -(t162 / 0.4e1 + t169 / 0.4e1) * t594 + t584;
t358 = -(mrSges(7,3) * t563 + t420) * t276 + (t219 * t570 + t421) * t274 + t430 - t219 * t561 + t361;
t7 = t358 - t585;
t369 = t7 * qJD(1) + t37 * qJD(2) + t39 * qJD(3) + t34 * qJD(4);
t303 = t312 * qJD(6);
t40 = t439 + t467;
t38 = t440 + t468;
t33 = t362 - t365 + t439;
t32 = t363 - t366 + t440;
t17 = t359 + t598;
t15 = -t530 / 0.2e1 - t529 / 0.2e1 - t527 / 0.2e1 + t528 / 0.2e1 + t364 + t466;
t12 = (-t355 * t408 + t596) * t353 + t611 + t618;
t9 = t388 / 0.2e1 + t356;
t6 = t358 + t585;
t3 = qJD(3) * t453 - t11 * qJD(4) - t14 * qJD(5) + t19 * qJD(6);
t2 = t348 * t544 - t420 * t276 + t361 + (-pkin(8) * t408 + t426) * t355 + (-t383 * t449 - t450 * t594 + t595) * mrSges(7,3) + (t279 * t461 + t310 * t573 + t385 - t514 / 0.2e1 + (t604 - t597 / 0.4e1 - t349 - t451 * t352) * t355) * t352 + (t518 / 0.4e1 + t311 * t573 + (t549 - t317 / 0.4e1 + t451 * t354 + (m(7) * t548 + t566) * pkin(5)) * t355) * t354 + t373 * t574 + t608 + t619;
t21 = [qJD(3) * t29 - qJD(4) * t4 - qJD(5) * t5 + qJD(6) * t8, t3, t481 + t9 * qJD(4) + t12 * qJD(5) + t17 * qJD(6) + 0.2e1 * (t108 * qJD(2) / 0.4e1 + (t275 * t594 - t277 * t383) * qJD(3) / 0.2e1) * m(7), t9 * qJD(3) + (t383 * t400 / 0.2e1 + t398 * t551 + t342 * t403 + pkin(4) * t475 + t204 * t554 - t278 * t200 + t202 * t557 + t592 * mrSges(6,3)) * qJD(4) + t2 * qJD(5) + t6 * qJD(6) + ((-t219 * t384 - t342 * t278 - t382 * t406) * t574 + pkin(8) * t389 / 0.2e1) * t578 + (-t219 * t275 - t277 * t406 + t382 * t383 - t384 * t594) * t492 + (t479 / 0.2e1 + t320 * t547 - Ifges(5,5) + (-m(6) * pkin(4) + t493) * t336 - t612 * t354 / 0.2e1) * t463 + (mrSges(7,1) * t406 - mrSges(5,2) * t336 + mrSges(7,2) * t219 + Ifges(6,5) * t352 + Ifges(7,5) * t383 + Ifges(6,6) * t354 - Ifges(7,6) * t594 - pkin(8) * t316 - Ifges(5,6)) * t462 + t397, t12 * qJD(3) + t2 * qJD(4) + (-Ifges(6,5) * t477 - Ifges(6,6) * t472 + (t539 * t88 + t540 * t87) * t572 + t274 * t411 + t276 * t412 - t527 + t528 - t193 * mrSges(6,2) - t194 * mrSges(6,1) + t466) * qJD(5) + t15 * qJD(6) + t396, t17 * qJD(3) + t6 * qJD(4) + t15 * qJD(5) + (t466 - t529 - t530) * qJD(6) + t395; t3, t616, qJD(1) * t453 + t459, t32 * qJD(5) + t38 * qJD(6) + (-t275 * t383 - t277 * t594) * t492 + t355 * t407 - t588 * t463 + ((-t464 * t531 - t535) * t577 + (-t219 * t277 + t275 * t406 + t342 * t355) * t574) * t578 + t379, -t490 + t32 * qJD(4) + (m(7) * (t274 * t457 - t276 * t456) + t291 - t168) * qJD(5) - t607, qJD(4) * t38 - qJD(5) * t168 + t488 - t607; qJD(2) * t454 + t10 * qJD(4) + t13 * qJD(5) + t16 * qJD(6) - t481, qJD(1) * t454 + t459, -t616, t33 * qJD(5) + t40 * qJD(6) + (t274 * t594 + t276 * t383) * t492 + t588 * t462 + t353 * t407 + ((pkin(8) * t413 - t532) * t577 + (t342 * t353 + t484 - t485) * t574) * t578 + t380, t491 + t33 * qJD(4) + (m(7) * (-t275 * t456 + t277 * t457) + t591 * t353 + t415) * qJD(5) + t606, t40 * qJD(4) + qJD(5) * t415 + t489 + t606; -qJD(3) * t10 + qJD(5) * t1 + qJD(6) * t7 - t397, -qJD(5) * t30 + qJD(6) * t37 - t379, -qJD(5) * t31 + qJD(6) * t39 - t380, qJD(5) * t20 + qJD(6) * t34 (-t383 * t411 + t594 * t412 + (t219 * t540 + t406 * t539) * t572 + t414 + t591 * pkin(8) + t41) * qJD(5) + t617 + t370, t41 * qJD(5) + t369 + t617; -qJD(3) * t13 - qJD(4) * t1 + qJD(6) * t18 - t396, t30 * qJD(4) + t490, t31 * qJD(4) - t491, qJD(6) * t42 - t370, -t303, -t303 - t377; -qJD(3) * t16 - qJD(4) * t7 - qJD(5) * t18 - t395, -t37 * qJD(4) - t488, -t39 * qJD(4) - t489, -qJD(5) * t42 - t369, t377, 0;];
Cq  = t21;
