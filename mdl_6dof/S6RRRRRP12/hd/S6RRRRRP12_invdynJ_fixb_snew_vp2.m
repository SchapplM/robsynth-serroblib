% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRRRP12
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
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
% Datum: 2019-05-08 07:55
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRRRP12_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP12_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP12_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP12_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP12_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP12_invdynJ_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP12_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP12_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP12_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 07:22:02
% EndTime: 2019-05-08 07:22:24
% DurationCPUTime: 18.00s
% Computational Cost: add. (285264->356), mult. (703090->468), div. (0->0), fcn. (591135->14), ass. (0->160)
t609 = Ifges(6,1) + Ifges(7,1);
t599 = Ifges(6,4) - Ifges(7,5);
t598 = -Ifges(6,5) - Ifges(7,4);
t608 = Ifges(6,2) + Ifges(7,3);
t597 = Ifges(6,6) - Ifges(7,6);
t607 = -Ifges(6,3) - Ifges(7,2);
t552 = cos(pkin(6));
t547 = qJD(1) * t552 + qJD(2);
t549 = sin(pkin(7));
t551 = cos(pkin(7));
t550 = sin(pkin(6));
t560 = cos(qJ(2));
t581 = qJD(1) * t560;
t576 = t550 * t581;
t534 = (t547 * t549 + t551 * t576) * pkin(10);
t556 = sin(qJ(2));
t583 = qJD(1) * t550;
t603 = pkin(10) * t549;
t538 = (-pkin(2) * t560 - t556 * t603) * t583;
t580 = qJD(1) * qJD(2);
t544 = (qJDD(1) * t556 + t560 * t580) * t550;
t546 = qJDD(1) * t552 + qJDD(2);
t562 = qJD(1) ^ 2;
t557 = sin(qJ(1));
t561 = cos(qJ(1));
t571 = -g(1) * t561 - g(2) * t557;
t604 = pkin(9) * t550;
t542 = -pkin(1) * t562 + qJDD(1) * t604 + t571;
t575 = t557 * g(1) - g(2) * t561;
t541 = qJDD(1) * pkin(1) + t562 * t604 + t575;
t595 = t541 * t552;
t573 = -t556 * t542 + t560 * t595;
t582 = qJD(1) * t556;
t602 = pkin(10) * t551;
t492 = -t544 * t602 + t546 * pkin(2) + t547 * t534 + (-g(3) * t560 - t538 * t582) * t550 + t573;
t577 = t550 * t582;
t537 = pkin(2) * t547 - t577 * t602;
t545 = (qJDD(1) * t560 - t556 * t580) * t550;
t569 = t545 * t551 + t546 * t549;
t584 = t560 * t542 + t556 * t595;
t493 = -t547 * t537 + (-g(3) * t556 + t538 * t581) * t550 + t569 * pkin(10) + t584;
t601 = t552 * g(3);
t498 = -t544 * t603 - t545 * pkin(2) - t601 + (-t541 + (-t534 * t560 + t537 * t556) * qJD(1)) * t550;
t555 = sin(qJ(3));
t559 = cos(qJ(3));
t458 = -t555 * t493 + (t492 * t551 + t498 * t549) * t559;
t589 = t551 * t560;
t594 = t549 * t555;
t525 = t547 * t594 + (t555 * t589 + t556 * t559) * t583;
t509 = -t525 * qJD(3) - t555 * t544 + t559 * t569;
t593 = t549 * t559;
t524 = (-t555 * t556 + t559 * t589) * t583 + t547 * t593;
t510 = t524 * qJD(3) + t559 * t544 + t555 * t569;
t535 = t547 * t551 - t549 * t576 + qJD(3);
t554 = sin(qJ(4));
t558 = cos(qJ(4));
t515 = -t525 * t554 + t535 * t558;
t526 = -t545 * t549 + t546 * t551 + qJDD(3);
t479 = qJD(4) * t515 + t510 * t558 + t526 * t554;
t516 = t525 * t558 + t535 * t554;
t523 = qJD(4) - t524;
t553 = sin(qJ(5));
t605 = cos(qJ(5));
t500 = t553 * t516 - t523 * t605;
t508 = qJDD(4) - t509;
t457 = -t500 * qJD(5) + t479 * t605 + t553 * t508;
t501 = t516 * t605 + t553 * t523;
t472 = mrSges(7,1) * t500 - mrSges(7,3) * t501;
t590 = t551 * t555;
t459 = t492 * t590 + t559 * t493 + t498 * t594;
t512 = -pkin(3) * t524 - pkin(11) * t525;
t533 = t535 ^ 2;
t451 = -pkin(3) * t533 + pkin(11) * t526 + t512 * t524 + t459;
t468 = -t549 * t492 + t551 * t498;
t453 = (-t524 * t535 - t510) * pkin(11) + (t525 * t535 - t509) * pkin(3) + t468;
t447 = t558 * t451 + t554 * t453;
t495 = -pkin(4) * t515 - pkin(12) * t516;
t522 = t523 ^ 2;
t443 = -pkin(4) * t522 + pkin(12) * t508 + t495 * t515 + t447;
t450 = -t526 * pkin(3) - t533 * pkin(11) + t525 * t512 - t458;
t478 = -qJD(4) * t516 - t510 * t554 + t526 * t558;
t445 = (-t515 * t523 - t479) * pkin(12) + (t516 * t523 - t478) * pkin(4) + t450;
t439 = -t553 * t443 + t445 * t605;
t471 = pkin(5) * t500 - qJ(6) * t501;
t477 = qJDD(5) - t478;
t514 = qJD(5) - t515;
t513 = t514 ^ 2;
t437 = -t477 * pkin(5) - t513 * qJ(6) + t501 * t471 + qJDD(6) - t439;
t481 = -mrSges(7,2) * t500 + mrSges(7,3) * t514;
t572 = -m(7) * t437 + t477 * mrSges(7,1) + t514 * t481;
t433 = t457 * mrSges(7,2) + t501 * t472 - t572;
t440 = t443 * t605 + t553 * t445;
t436 = -pkin(5) * t513 + qJ(6) * t477 + 0.2e1 * qJD(6) * t514 - t471 * t500 + t440;
t456 = t501 * qJD(5) + t553 * t479 - t508 * t605;
t484 = -mrSges(7,1) * t514 + mrSges(7,2) * t501;
t578 = m(7) * t436 + t477 * mrSges(7,3) + t514 * t484;
t586 = t599 * t500 - t609 * t501 + t598 * t514;
t587 = t608 * t500 - t599 * t501 - t597 * t514;
t606 = -t456 * t597 - t457 * t598 - t607 * t477 - t500 * t586 - t501 * t587 + mrSges(6,1) * t439 - mrSges(7,1) * t437 - mrSges(6,2) * t440 + mrSges(7,3) * t436 - pkin(5) * t433 + qJ(6) * (-t456 * mrSges(7,2) - t500 * t472 + t578);
t600 = -mrSges(6,3) - mrSges(7,2);
t592 = t550 * t556;
t591 = t550 * t560;
t483 = mrSges(6,1) * t514 - mrSges(6,3) * t501;
t585 = -mrSges(6,1) * t500 - mrSges(6,2) * t501 - t472;
t429 = m(6) * t440 - t477 * mrSges(6,2) + t456 * t600 - t514 * t483 + t500 * t585 + t578;
t482 = -mrSges(6,2) * t514 - mrSges(6,3) * t500;
t430 = m(6) * t439 + t477 * mrSges(6,1) + t457 * t600 + t514 * t482 + t501 * t585 + t572;
t425 = t429 * t605 - t430 * t553;
t494 = -mrSges(5,1) * t515 + mrSges(5,2) * t516;
t503 = mrSges(5,1) * t523 - mrSges(5,3) * t516;
t421 = m(5) * t447 - mrSges(5,2) * t508 + mrSges(5,3) * t478 + t494 * t515 - t503 * t523 + t425;
t446 = -t554 * t451 + t558 * t453;
t442 = -t508 * pkin(4) - t522 * pkin(12) + t516 * t495 - t446;
t438 = -0.2e1 * qJD(6) * t501 + (t500 * t514 - t457) * qJ(6) + (t501 * t514 + t456) * pkin(5) + t442;
t434 = m(7) * t438 + mrSges(7,1) * t456 - t457 * mrSges(7,3) + t481 * t500 - t501 * t484;
t431 = -m(6) * t442 - t456 * mrSges(6,1) - mrSges(6,2) * t457 - t500 * t482 - t483 * t501 - t434;
t502 = -mrSges(5,2) * t523 + mrSges(5,3) * t515;
t427 = m(5) * t446 + mrSges(5,1) * t508 - mrSges(5,3) * t479 - t494 * t516 + t502 * t523 + t431;
t415 = t554 * t421 + t558 * t427;
t588 = t597 * t500 + t598 * t501 + t607 * t514;
t511 = -mrSges(4,1) * t524 + mrSges(4,2) * t525;
t518 = mrSges(4,1) * t535 - mrSges(4,3) * t525;
t574 = t558 * t421 - t427 * t554;
t412 = m(4) * t459 - mrSges(4,2) * t526 + mrSges(4,3) * t509 + t511 * t524 - t518 * t535 + t574;
t517 = -mrSges(4,2) * t535 + mrSges(4,3) * t524;
t414 = m(4) * t468 - mrSges(4,1) * t509 + mrSges(4,2) * t510 - t517 * t524 + t518 * t525 + t415;
t424 = t553 * t429 + t430 * t605;
t564 = -m(5) * t450 + t478 * mrSges(5,1) - t479 * mrSges(5,2) + t515 * t502 - t516 * t503 - t424;
t418 = m(4) * t458 + t526 * mrSges(4,1) - t510 * mrSges(4,3) - t525 * t511 + t535 * t517 + t564;
t403 = t412 * t594 + t551 * t414 + t418 * t593;
t407 = t559 * t412 - t418 * t555;
t404 = t551 * t559 * t418 + t412 * t590 - t414 * t549;
t422 = -mrSges(6,1) * t442 - mrSges(7,1) * t438 + mrSges(7,2) * t436 + mrSges(6,3) * t440 - pkin(5) * t434 - t608 * t456 + t599 * t457 + t597 * t477 + t588 * t501 - t586 * t514;
t423 = mrSges(6,2) * t442 + mrSges(7,2) * t437 - mrSges(6,3) * t439 - mrSges(7,3) * t438 - qJ(6) * t434 - t599 * t456 + t609 * t457 - t598 * t477 + t588 * t500 + t587 * t514;
t485 = Ifges(5,5) * t516 + Ifges(5,6) * t515 + Ifges(5,3) * t523;
t486 = Ifges(5,4) * t516 + Ifges(5,2) * t515 + Ifges(5,6) * t523;
t405 = mrSges(5,2) * t450 - mrSges(5,3) * t446 + Ifges(5,1) * t479 + Ifges(5,4) * t478 + Ifges(5,5) * t508 - pkin(12) * t424 - t553 * t422 + t423 * t605 + t515 * t485 - t523 * t486;
t487 = Ifges(5,1) * t516 + Ifges(5,4) * t515 + Ifges(5,5) * t523;
t408 = -mrSges(5,1) * t450 + mrSges(5,3) * t447 + Ifges(5,4) * t479 + Ifges(5,2) * t478 + Ifges(5,6) * t508 - pkin(4) * t424 - t516 * t485 + t523 * t487 - t606;
t504 = Ifges(4,5) * t525 + Ifges(4,6) * t524 + Ifges(4,3) * t535;
t505 = Ifges(4,4) * t525 + Ifges(4,2) * t524 + Ifges(4,6) * t535;
t400 = mrSges(4,2) * t468 - mrSges(4,3) * t458 + Ifges(4,1) * t510 + Ifges(4,4) * t509 + Ifges(4,5) * t526 - pkin(11) * t415 + t405 * t558 - t408 * t554 + t504 * t524 - t505 * t535;
t506 = Ifges(4,1) * t525 + Ifges(4,4) * t524 + Ifges(4,5) * t535;
t563 = mrSges(5,1) * t446 - mrSges(5,2) * t447 + Ifges(5,5) * t479 + Ifges(5,6) * t478 + Ifges(5,3) * t508 + pkin(4) * t431 + pkin(12) * t425 + t422 * t605 + t553 * t423 + t516 * t486 - t515 * t487;
t401 = -mrSges(4,1) * t468 + mrSges(4,3) * t459 + Ifges(4,4) * t510 + Ifges(4,2) * t509 + Ifges(4,6) * t526 - pkin(3) * t415 - t525 * t504 + t535 * t506 - t563;
t566 = pkin(10) * t407 + t400 * t555 + t401 * t559;
t543 = (-mrSges(3,1) * t560 + mrSges(3,2) * t556) * t583;
t540 = -mrSges(3,2) * t547 + mrSges(3,3) * t576;
t539 = mrSges(3,1) * t547 - mrSges(3,3) * t577;
t530 = -t550 * t541 - t601;
t529 = Ifges(3,5) * t547 + (Ifges(3,1) * t556 + Ifges(3,4) * t560) * t583;
t528 = Ifges(3,6) * t547 + (Ifges(3,4) * t556 + Ifges(3,2) * t560) * t583;
t527 = Ifges(3,3) * t547 + (Ifges(3,5) * t556 + Ifges(3,6) * t560) * t583;
t520 = -g(3) * t592 + t584;
t519 = -g(3) * t591 + t573;
t406 = m(3) * t520 - mrSges(3,2) * t546 + mrSges(3,3) * t545 - t539 * t547 + t543 * t576 + t407;
t402 = m(3) * t519 + mrSges(3,1) * t546 - mrSges(3,3) * t544 + t540 * t547 - t543 * t577 + t404;
t399 = mrSges(4,1) * t458 - mrSges(4,2) * t459 + Ifges(4,5) * t510 + Ifges(4,6) * t509 + Ifges(4,3) * t526 + pkin(3) * t564 + pkin(11) * t574 + t554 * t405 + t558 * t408 + t525 * t505 - t524 * t506;
t398 = mrSges(3,1) * t519 - mrSges(3,2) * t520 + Ifges(3,5) * t544 + Ifges(3,6) * t545 + Ifges(3,3) * t546 + pkin(2) * t404 + t551 * t399 + (t528 * t556 - t529 * t560) * t583 + t566 * t549;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t575 - mrSges(2,2) * t571 + (t527 * t576 + mrSges(3,2) * t530 - mrSges(3,3) * t519 + Ifges(3,1) * t544 + Ifges(3,4) * t545 + Ifges(3,5) * t546 + t559 * t400 - t555 * t401 - t547 * t528 + (-t403 * t549 - t404 * t551) * pkin(10)) * t592 + (-mrSges(3,1) * t530 + mrSges(3,3) * t520 + Ifges(3,4) * t544 + Ifges(3,2) * t545 + Ifges(3,6) * t546 - pkin(2) * t403 - t549 * t399 - t527 * t577 + t547 * t529 + t551 * t566) * t591 + t552 * t398 + pkin(1) * ((t402 * t560 + t406 * t556) * t552 + (-m(3) * t530 + t545 * mrSges(3,1) - t544 * mrSges(3,2) + (-t539 * t556 + t540 * t560) * t583 - t403) * t550) + (-t402 * t556 + t406 * t560) * t604; t398; t399; t563; t606; t433;];
tauJ  = t1;
