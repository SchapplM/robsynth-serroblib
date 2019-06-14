% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRRPR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
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
% tauJ [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-08 02:41
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRRPR14_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR14_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR14_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR14_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR14_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR14_invdynJ_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR14_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR14_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR14_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 02:13:57
% EndTime: 2019-05-08 02:14:53
% DurationCPUTime: 40.25s
% Computational Cost: add. (639119->379), mult. (1584605->511), div. (0->0), fcn. (1344086->16), ass. (0->167)
t563 = cos(pkin(6));
t556 = qJD(1) * t563 + qJD(2);
t559 = sin(pkin(7));
t562 = cos(pkin(7));
t560 = sin(pkin(6));
t571 = cos(qJ(2));
t593 = qJD(1) * t571;
t589 = t560 * t593;
t543 = (t556 * t559 + t562 * t589) * pkin(10);
t567 = sin(qJ(2));
t595 = qJD(1) * t560;
t606 = pkin(10) * t559;
t547 = (-pkin(2) * t571 - t567 * t606) * t595;
t592 = qJD(1) * qJD(2);
t553 = (qJDD(1) * t567 + t571 * t592) * t560;
t555 = qJDD(1) * t563 + qJDD(2);
t573 = qJD(1) ^ 2;
t568 = sin(qJ(1));
t572 = cos(qJ(1));
t584 = -g(1) * t572 - g(2) * t568;
t607 = pkin(9) * t560;
t551 = -pkin(1) * t573 + qJDD(1) * t607 + t584;
t588 = t568 * g(1) - g(2) * t572;
t550 = qJDD(1) * pkin(1) + t573 * t607 + t588;
t603 = t550 * t563;
t585 = -t567 * t551 + t571 * t603;
t594 = qJD(1) * t567;
t605 = pkin(10) * t562;
t499 = -t553 * t605 + t555 * pkin(2) + t556 * t543 + (-g(3) * t571 - t547 * t594) * t560 + t585;
t590 = t560 * t594;
t546 = pkin(2) * t556 - t590 * t605;
t554 = (qJDD(1) * t571 - t567 * t592) * t560;
t581 = t554 * t562 + t555 * t559;
t596 = t571 * t551 + t567 * t603;
t500 = -t556 * t546 + (-g(3) * t567 + t547 * t593) * t560 + t581 * pkin(10) + t596;
t604 = t563 * g(3);
t505 = -t553 * t606 - t554 * pkin(2) - t604 + (-t550 + (-t543 * t571 + t546 * t567) * qJD(1)) * t560;
t566 = sin(qJ(3));
t570 = cos(qJ(3));
t470 = -t566 * t500 + (t499 * t562 + t505 * t559) * t570;
t597 = t562 * t571;
t602 = t559 * t566;
t533 = t556 * t602 + (t566 * t597 + t567 * t570) * t595;
t518 = -t533 * qJD(3) - t566 * t553 + t570 * t581;
t601 = t559 * t570;
t532 = (-t566 * t567 + t570 * t597) * t595 + t556 * t601;
t608 = cos(qJ(4));
t600 = t560 * t567;
t599 = t560 * t571;
t598 = t562 * t566;
t471 = t499 * t598 + t570 * t500 + t505 * t602;
t521 = -pkin(3) * t532 - pkin(11) * t533;
t534 = -t554 * t559 + t555 * t562 + qJDD(3);
t544 = t556 * t562 - t559 * t589 + qJD(3);
t542 = t544 ^ 2;
t462 = -pkin(3) * t542 + pkin(11) * t534 + t521 * t532 + t471;
t481 = -t559 * t499 + t562 * t505;
t519 = t532 * qJD(3) + t570 * t553 + t566 * t581;
t464 = (-t532 * t544 - t519) * pkin(11) + (t533 * t544 - t518) * pkin(3) + t481;
t565 = sin(qJ(4));
t455 = t462 * t608 + t565 * t464;
t523 = t533 * t565 - t544 * t608;
t524 = t533 * t608 + t565 * t544;
t501 = pkin(4) * t523 - qJ(5) * t524;
t517 = qJDD(4) - t518;
t530 = qJD(4) - t532;
t529 = t530 ^ 2;
t450 = -pkin(4) * t529 + qJ(5) * t517 - t501 * t523 + t455;
t461 = -t534 * pkin(3) - t542 * pkin(11) + t533 * t521 - t470;
t487 = qJD(4) * t524 + t519 * t565 - t534 * t608;
t488 = -t523 * qJD(4) + t519 * t608 + t565 * t534;
t453 = (t523 * t530 - t488) * qJ(5) + (t524 * t530 + t487) * pkin(4) + t461;
t558 = sin(pkin(13));
t561 = cos(pkin(13));
t511 = t524 * t561 + t530 * t558;
t445 = -0.2e1 * qJD(5) * t511 - t558 * t450 + t561 * t453;
t474 = t488 * t561 + t517 * t558;
t510 = -t524 * t558 + t530 * t561;
t443 = (t510 * t523 - t474) * pkin(12) + (t510 * t511 + t487) * pkin(5) + t445;
t446 = 0.2e1 * qJD(5) * t510 + t561 * t450 + t558 * t453;
t473 = -t488 * t558 + t517 * t561;
t491 = pkin(5) * t523 - pkin(12) * t511;
t509 = t510 ^ 2;
t444 = -pkin(5) * t509 + pkin(12) * t473 - t491 * t523 + t446;
t564 = sin(qJ(6));
t569 = cos(qJ(6));
t441 = t443 * t569 - t444 * t564;
t482 = t510 * t569 - t511 * t564;
t458 = qJD(6) * t482 + t473 * t564 + t474 * t569;
t483 = t510 * t564 + t511 * t569;
t469 = -mrSges(7,1) * t482 + mrSges(7,2) * t483;
t522 = qJD(6) + t523;
t475 = -mrSges(7,2) * t522 + mrSges(7,3) * t482;
t486 = qJDD(6) + t487;
t438 = m(7) * t441 + mrSges(7,1) * t486 - mrSges(7,3) * t458 - t469 * t483 + t475 * t522;
t442 = t443 * t564 + t444 * t569;
t457 = -qJD(6) * t483 + t473 * t569 - t474 * t564;
t476 = mrSges(7,1) * t522 - mrSges(7,3) * t483;
t439 = m(7) * t442 - mrSges(7,2) * t486 + mrSges(7,3) * t457 + t469 * t482 - t476 * t522;
t430 = t569 * t438 + t564 * t439;
t484 = -mrSges(6,1) * t510 + mrSges(6,2) * t511;
t583 = -mrSges(6,2) * t523 + mrSges(6,3) * t510;
t428 = m(6) * t445 + t487 * mrSges(6,1) - t474 * mrSges(6,3) - t511 * t484 + t523 * t583 + t430;
t490 = mrSges(6,1) * t523 - mrSges(6,3) * t511;
t586 = -t438 * t564 + t569 * t439;
t429 = m(6) * t446 - mrSges(6,2) * t487 + mrSges(6,3) * t473 + t484 * t510 - t490 * t523 + t586;
t426 = -t428 * t558 + t561 * t429;
t502 = mrSges(5,1) * t523 + mrSges(5,2) * t524;
t513 = mrSges(5,1) * t530 - mrSges(5,3) * t524;
t424 = m(5) * t455 - mrSges(5,2) * t517 - mrSges(5,3) * t487 - t502 * t523 - t513 * t530 + t426;
t454 = -t565 * t462 + t464 * t608;
t449 = -t517 * pkin(4) - t529 * qJ(5) + t524 * t501 + qJDD(5) - t454;
t447 = -t473 * pkin(5) - t509 * pkin(12) + t511 * t491 + t449;
t577 = m(7) * t447 - t457 * mrSges(7,1) + mrSges(7,2) * t458 - t482 * t475 + t476 * t483;
t440 = m(6) * t449 - t473 * mrSges(6,1) + mrSges(6,2) * t474 + t490 * t511 - t510 * t583 + t577;
t512 = -mrSges(5,2) * t530 - mrSges(5,3) * t523;
t434 = m(5) * t454 + mrSges(5,1) * t517 - mrSges(5,3) * t488 - t502 * t524 + t512 * t530 - t440;
t416 = t565 * t424 + t434 * t608;
t520 = -mrSges(4,1) * t532 + mrSges(4,2) * t533;
t526 = mrSges(4,1) * t544 - mrSges(4,3) * t533;
t587 = t424 * t608 - t434 * t565;
t413 = m(4) * t471 - mrSges(4,2) * t534 + mrSges(4,3) * t518 + t520 * t532 - t526 * t544 + t587;
t525 = -mrSges(4,2) * t544 + mrSges(4,3) * t532;
t415 = m(4) * t481 - mrSges(4,1) * t518 + mrSges(4,2) * t519 - t525 * t532 + t526 * t533 + t416;
t425 = t428 * t561 + t429 * t558;
t576 = -m(5) * t461 - t487 * mrSges(5,1) - mrSges(5,2) * t488 - t523 * t512 - t513 * t524 - t425;
t421 = m(4) * t470 + mrSges(4,1) * t534 - mrSges(4,3) * t519 - t520 * t533 + t525 * t544 + t576;
t404 = t413 * t602 + t562 * t415 + t421 * t601;
t408 = t570 * t413 - t421 * t566;
t405 = t562 * t570 * t421 + t413 * t598 - t415 * t559;
t465 = Ifges(7,5) * t483 + Ifges(7,6) * t482 + Ifges(7,3) * t522;
t467 = Ifges(7,1) * t483 + Ifges(7,4) * t482 + Ifges(7,5) * t522;
t431 = -mrSges(7,1) * t447 + mrSges(7,3) * t442 + Ifges(7,4) * t458 + Ifges(7,2) * t457 + Ifges(7,6) * t486 - t465 * t483 + t467 * t522;
t466 = Ifges(7,4) * t483 + Ifges(7,2) * t482 + Ifges(7,6) * t522;
t432 = mrSges(7,2) * t447 - mrSges(7,3) * t441 + Ifges(7,1) * t458 + Ifges(7,4) * t457 + Ifges(7,5) * t486 + t465 * t482 - t466 * t522;
t477 = Ifges(6,5) * t511 + Ifges(6,6) * t510 + Ifges(6,3) * t523;
t479 = Ifges(6,1) * t511 + Ifges(6,4) * t510 + Ifges(6,5) * t523;
t417 = -mrSges(6,1) * t449 + mrSges(6,3) * t446 + Ifges(6,4) * t474 + Ifges(6,2) * t473 + Ifges(6,6) * t487 - pkin(5) * t577 + pkin(12) * t586 + t569 * t431 + t564 * t432 - t511 * t477 + t523 * t479;
t478 = Ifges(6,4) * t511 + Ifges(6,2) * t510 + Ifges(6,6) * t523;
t418 = mrSges(6,2) * t449 - mrSges(6,3) * t445 + Ifges(6,1) * t474 + Ifges(6,4) * t473 + Ifges(6,5) * t487 - pkin(12) * t430 - t431 * t564 + t432 * t569 + t477 * t510 - t478 * t523;
t492 = Ifges(5,5) * t524 - Ifges(5,6) * t523 + Ifges(5,3) * t530;
t493 = Ifges(5,4) * t524 - Ifges(5,2) * t523 + Ifges(5,6) * t530;
t406 = mrSges(5,2) * t461 - mrSges(5,3) * t454 + Ifges(5,1) * t488 - Ifges(5,4) * t487 + Ifges(5,5) * t517 - qJ(5) * t425 - t417 * t558 + t418 * t561 - t492 * t523 - t493 * t530;
t494 = Ifges(5,1) * t524 - Ifges(5,4) * t523 + Ifges(5,5) * t530;
t575 = mrSges(7,1) * t441 - mrSges(7,2) * t442 + Ifges(7,5) * t458 + Ifges(7,6) * t457 + Ifges(7,3) * t486 + t483 * t466 - t482 * t467;
t409 = -t524 * t492 + t530 * t494 + Ifges(5,6) * t517 + t510 * t479 - t511 * t478 + Ifges(5,4) * t488 - Ifges(6,6) * t473 - Ifges(6,5) * t474 - t575 - mrSges(5,1) * t461 + mrSges(5,3) * t455 + mrSges(6,2) * t446 - mrSges(6,1) * t445 + (-Ifges(5,2) - Ifges(6,3)) * t487 - pkin(5) * t430 - pkin(4) * t425;
t514 = Ifges(4,5) * t533 + Ifges(4,6) * t532 + Ifges(4,3) * t544;
t515 = Ifges(4,4) * t533 + Ifges(4,2) * t532 + Ifges(4,6) * t544;
t401 = mrSges(4,2) * t481 - mrSges(4,3) * t470 + Ifges(4,1) * t519 + Ifges(4,4) * t518 + Ifges(4,5) * t534 - pkin(11) * t416 + t406 * t608 - t565 * t409 + t532 * t514 - t544 * t515;
t516 = Ifges(4,1) * t533 + Ifges(4,4) * t532 + Ifges(4,5) * t544;
t574 = mrSges(5,1) * t454 - mrSges(5,2) * t455 + Ifges(5,5) * t488 - Ifges(5,6) * t487 + Ifges(5,3) * t517 - pkin(4) * t440 + qJ(5) * t426 + t561 * t417 + t558 * t418 + t524 * t493 + t523 * t494;
t402 = -mrSges(4,1) * t481 + mrSges(4,3) * t471 + Ifges(4,4) * t519 + Ifges(4,2) * t518 + Ifges(4,6) * t534 - pkin(3) * t416 - t533 * t514 + t544 * t516 - t574;
t578 = pkin(10) * t408 + t401 * t566 + t402 * t570;
t552 = (-mrSges(3,1) * t571 + mrSges(3,2) * t567) * t595;
t549 = -mrSges(3,2) * t556 + mrSges(3,3) * t589;
t548 = mrSges(3,1) * t556 - mrSges(3,3) * t590;
t538 = -t560 * t550 - t604;
t537 = Ifges(3,5) * t556 + (Ifges(3,1) * t567 + Ifges(3,4) * t571) * t595;
t536 = Ifges(3,6) * t556 + (Ifges(3,4) * t567 + Ifges(3,2) * t571) * t595;
t535 = Ifges(3,3) * t556 + (Ifges(3,5) * t567 + Ifges(3,6) * t571) * t595;
t528 = -g(3) * t600 + t596;
t527 = -g(3) * t599 + t585;
t407 = m(3) * t528 - mrSges(3,2) * t555 + mrSges(3,3) * t554 - t548 * t556 + t552 * t589 + t408;
t403 = m(3) * t527 + mrSges(3,1) * t555 - mrSges(3,3) * t553 + t549 * t556 - t552 * t590 + t405;
t400 = mrSges(4,1) * t470 - mrSges(4,2) * t471 + Ifges(4,5) * t519 + Ifges(4,6) * t518 + Ifges(4,3) * t534 + pkin(3) * t576 + pkin(11) * t587 + t565 * t406 + t409 * t608 + t533 * t515 - t532 * t516;
t399 = mrSges(3,1) * t527 - mrSges(3,2) * t528 + Ifges(3,5) * t553 + Ifges(3,6) * t554 + Ifges(3,3) * t555 + pkin(2) * t405 + t562 * t400 + (t536 * t567 - t537 * t571) * t595 + t578 * t559;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t588 - mrSges(2,2) * t584 + (t535 * t589 + mrSges(3,2) * t538 - mrSges(3,3) * t527 + Ifges(3,1) * t553 + Ifges(3,4) * t554 + Ifges(3,5) * t555 + t570 * t401 - t566 * t402 - t556 * t536 + (-t404 * t559 - t405 * t562) * pkin(10)) * t600 + (-mrSges(3,1) * t538 + mrSges(3,3) * t528 + Ifges(3,4) * t553 + Ifges(3,2) * t554 + Ifges(3,6) * t555 - pkin(2) * t404 - t559 * t400 - t535 * t590 + t556 * t537 + t562 * t578) * t599 + t563 * t399 + pkin(1) * ((t403 * t571 + t407 * t567) * t563 + (-m(3) * t538 + t554 * mrSges(3,1) - t553 * mrSges(3,2) + (-t548 * t567 + t549 * t571) * t595 - t404) * t560) + (-t403 * t567 + t407 * t571) * t607; t399; t400; t574; t440; t575;];
tauJ  = t1;
