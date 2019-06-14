% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRRRP12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-05-06 02:30
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRRRP12_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP12_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP12_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP12_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP12_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP12_invdynJ_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP12_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP12_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP12_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 02:17:11
% EndTime: 2019-05-06 02:17:24
% DurationCPUTime: 13.06s
% Computational Cost: add. (174510->332), mult. (541419->449), div. (0->0), fcn. (460597->14), ass. (0->158)
t556 = sin(pkin(12));
t558 = sin(pkin(6));
t559 = cos(pkin(12));
t561 = cos(pkin(6));
t564 = sin(qJ(3));
t560 = cos(pkin(7));
t567 = cos(qJ(3));
t601 = t560 * t567;
t557 = sin(pkin(7));
t605 = t557 * t567;
t573 = t558 * (-t556 * t564 + t559 * t601) + t561 * t605;
t533 = t573 * qJD(1);
t602 = t560 * t564;
t606 = t557 * t564;
t575 = t561 * t606 + (t556 * t567 + t559 * t602) * t558;
t534 = t575 * qJD(1);
t520 = -qJD(3) * t534 + qJDD(1) * t573;
t623 = Ifges(6,1) + Ifges(7,1);
t614 = Ifges(6,4) - Ifges(7,5);
t613 = -Ifges(6,5) - Ifges(7,4);
t622 = Ifges(6,2) + Ifges(7,3);
t612 = Ifges(6,6) - Ifges(7,6);
t621 = -Ifges(6,3) - Ifges(7,2);
t603 = t558 * t560;
t545 = (t557 * t561 + t559 * t603) * qJD(1) * pkin(9);
t569 = qJD(1) ^ 2;
t565 = sin(qJ(1));
t568 = cos(qJ(1));
t586 = -g(1) * t568 - g(2) * t565;
t609 = qJ(2) * t558;
t549 = -pkin(1) * t569 + qJDD(1) * t609 + t586;
t616 = pkin(9) * t557;
t582 = -pkin(2) * t559 - t556 * t616;
t596 = qJD(1) * t558;
t610 = pkin(9) * qJDD(1);
t577 = qJD(1) * t582 * t596 + t560 * t610;
t591 = qJD(2) * t596;
t604 = t558 * t559;
t590 = t565 * g(1) - g(2) * t568;
t548 = qJDD(1) * pkin(1) + t569 * t609 + t590;
t608 = t548 * t561;
t583 = -g(3) * t604 - 0.2e1 * t556 * t591 + t559 * t608;
t501 = (pkin(2) * qJDD(1) + qJD(1) * t545) * t561 + (-t558 * t577 - t549) * t556 + t583;
t617 = pkin(9) * t556;
t550 = (pkin(2) * t561 - t603 * t617) * qJD(1);
t619 = 0.2e1 * t559;
t592 = t559 * t549 + t556 * t608 + t591 * t619;
t502 = (-qJD(1) * t550 + t557 * t610) * t561 + (-g(3) * t556 + t559 * t577) * t558 + t592;
t589 = -g(3) * t561 + qJDD(2);
t510 = (-t548 + t582 * qJDD(1) + (-t545 * t559 + t550 * t556) * qJD(1)) * t558 + t589;
t463 = -t564 * t502 + (t501 * t560 + t510 * t557) * t567;
t521 = qJD(3) * t533 + qJDD(1) * t575;
t578 = -t557 * t604 + t560 * t561;
t546 = qJD(1) * t578 + qJD(3);
t563 = sin(qJ(4));
t566 = cos(qJ(4));
t526 = -t534 * t563 + t546 * t566;
t543 = qJDD(1) * t578 + qJDD(3);
t497 = qJD(4) * t526 + t521 * t566 + t543 * t563;
t527 = t534 * t566 + t546 * t563;
t532 = qJD(4) - t533;
t562 = sin(qJ(5));
t618 = cos(qJ(5));
t508 = t527 * t562 - t532 * t618;
t517 = qJDD(4) - t520;
t468 = -t508 * qJD(5) + t497 * t618 + t562 * t517;
t509 = t527 * t618 + t562 * t532;
t481 = mrSges(7,1) * t508 - mrSges(7,3) * t509;
t464 = t501 * t602 + t567 * t502 + t510 * t606;
t519 = -pkin(3) * t533 - pkin(10) * t534;
t542 = t546 ^ 2;
t460 = -pkin(3) * t542 + pkin(10) * t543 + t519 * t533 + t464;
t477 = -t501 * t557 + t560 * t510;
t462 = (-t533 * t546 - t521) * pkin(10) + (t534 * t546 - t520) * pkin(3) + t477;
t456 = t566 * t460 + t563 * t462;
t504 = -pkin(4) * t526 - pkin(11) * t527;
t531 = t532 ^ 2;
t452 = -pkin(4) * t531 + pkin(11) * t517 + t504 * t526 + t456;
t459 = -pkin(3) * t543 - pkin(10) * t542 + t534 * t519 - t463;
t496 = -qJD(4) * t527 - t521 * t563 + t543 * t566;
t454 = (-t526 * t532 - t497) * pkin(11) + (t527 * t532 - t496) * pkin(4) + t459;
t448 = -t562 * t452 + t454 * t618;
t480 = pkin(5) * t508 - qJ(6) * t509;
t494 = qJDD(5) - t496;
t523 = qJD(5) - t526;
t522 = t523 ^ 2;
t446 = -t494 * pkin(5) - t522 * qJ(6) + t509 * t480 + qJDD(6) - t448;
t484 = -mrSges(7,2) * t508 + mrSges(7,3) * t523;
t587 = -m(7) * t446 + t494 * mrSges(7,1) + t523 * t484;
t442 = mrSges(7,2) * t468 + t481 * t509 - t587;
t449 = t618 * t452 + t562 * t454;
t445 = -pkin(5) * t522 + qJ(6) * t494 + 0.2e1 * qJD(6) * t523 - t480 * t508 + t449;
t467 = qJD(5) * t509 + t497 * t562 - t517 * t618;
t487 = -mrSges(7,1) * t523 + mrSges(7,2) * t509;
t593 = m(7) * t445 + t494 * mrSges(7,3) + t523 * t487;
t598 = t614 * t508 - t623 * t509 + t613 * t523;
t599 = t622 * t508 - t614 * t509 - t612 * t523;
t620 = -t467 * t612 - t468 * t613 - t621 * t494 - t508 * t598 - t509 * t599 + mrSges(6,1) * t448 - mrSges(7,1) * t446 - mrSges(6,2) * t449 + mrSges(7,3) * t445 - pkin(5) * t442 + qJ(6) * (-mrSges(7,2) * t467 - t481 * t508 + t593);
t615 = -mrSges(6,3) - mrSges(7,2);
t607 = t556 * t558;
t486 = mrSges(6,1) * t523 - mrSges(6,3) * t509;
t597 = -mrSges(6,1) * t508 - mrSges(6,2) * t509 - t481;
t438 = m(6) * t449 - mrSges(6,2) * t494 + t467 * t615 - t486 * t523 + t508 * t597 + t593;
t485 = -mrSges(6,2) * t523 - mrSges(6,3) * t508;
t439 = m(6) * t448 + mrSges(6,1) * t494 + t468 * t615 + t485 * t523 + t509 * t597 + t587;
t434 = t618 * t438 - t439 * t562;
t503 = -mrSges(5,1) * t526 + mrSges(5,2) * t527;
t512 = mrSges(5,1) * t532 - mrSges(5,3) * t527;
t430 = m(5) * t456 - mrSges(5,2) * t517 + mrSges(5,3) * t496 + t503 * t526 - t512 * t532 + t434;
t455 = -t563 * t460 + t462 * t566;
t451 = -pkin(4) * t517 - pkin(11) * t531 + t527 * t504 - t455;
t447 = -0.2e1 * qJD(6) * t509 + (t508 * t523 - t468) * qJ(6) + (t509 * t523 + t467) * pkin(5) + t451;
t443 = m(7) * t447 + mrSges(7,1) * t467 - t468 * mrSges(7,3) + t484 * t508 - t509 * t487;
t440 = -m(6) * t451 - t467 * mrSges(6,1) - mrSges(6,2) * t468 - t508 * t485 - t486 * t509 - t443;
t511 = -mrSges(5,2) * t532 + mrSges(5,3) * t526;
t436 = m(5) * t455 + mrSges(5,1) * t517 - mrSges(5,3) * t497 - t503 * t527 + t511 * t532 + t440;
t424 = t563 * t430 + t566 * t436;
t600 = t612 * t508 + t613 * t509 + t621 * t523;
t518 = -mrSges(4,1) * t533 + mrSges(4,2) * t534;
t529 = mrSges(4,1) * t546 - mrSges(4,3) * t534;
t588 = t566 * t430 - t436 * t563;
t421 = m(4) * t464 - mrSges(4,2) * t543 + mrSges(4,3) * t520 + t518 * t533 - t529 * t546 + t588;
t528 = -mrSges(4,2) * t546 + mrSges(4,3) * t533;
t423 = m(4) * t477 - mrSges(4,1) * t520 + mrSges(4,2) * t521 - t528 * t533 + t529 * t534 + t424;
t433 = t562 * t438 + t439 * t618;
t571 = -m(5) * t459 + t496 * mrSges(5,1) - t497 * mrSges(5,2) + t526 * t511 - t527 * t512 - t433;
t427 = m(4) * t463 + t543 * mrSges(4,1) - t521 * mrSges(4,3) - t534 * t518 + t546 * t528 + t571;
t412 = t421 * t606 + t560 * t423 + t427 * t605;
t416 = t567 * t421 - t427 * t564;
t413 = t421 * t602 - t423 * t557 + t427 * t601;
t585 = -mrSges(3,1) * t559 + mrSges(3,2) * t556;
t581 = mrSges(3,1) * t561 - mrSges(3,3) * t607;
t580 = -mrSges(3,2) * t561 + mrSges(3,3) * t604;
t431 = -mrSges(6,1) * t451 - mrSges(7,1) * t447 + mrSges(7,2) * t445 + mrSges(6,3) * t449 - pkin(5) * t443 - t622 * t467 + t614 * t468 + t612 * t494 + t600 * t509 - t598 * t523;
t432 = mrSges(6,2) * t451 + mrSges(7,2) * t446 - mrSges(6,3) * t448 - mrSges(7,3) * t447 - qJ(6) * t443 - t614 * t467 + t623 * t468 - t613 * t494 + t600 * t508 + t599 * t523;
t490 = Ifges(5,5) * t527 + Ifges(5,6) * t526 + Ifges(5,3) * t532;
t491 = Ifges(5,4) * t527 + Ifges(5,2) * t526 + Ifges(5,6) * t532;
t414 = mrSges(5,2) * t459 - mrSges(5,3) * t455 + Ifges(5,1) * t497 + Ifges(5,4) * t496 + Ifges(5,5) * t517 - pkin(11) * t433 - t562 * t431 + t432 * t618 + t526 * t490 - t532 * t491;
t492 = Ifges(5,1) * t527 + Ifges(5,4) * t526 + Ifges(5,5) * t532;
t417 = -mrSges(5,1) * t459 + mrSges(5,3) * t456 + Ifges(5,4) * t497 + Ifges(5,2) * t496 + Ifges(5,6) * t517 - pkin(4) * t433 - t527 * t490 + t532 * t492 - t620;
t513 = Ifges(4,5) * t534 + Ifges(4,6) * t533 + Ifges(4,3) * t546;
t514 = Ifges(4,4) * t534 + Ifges(4,2) * t533 + Ifges(4,6) * t546;
t408 = mrSges(4,2) * t477 - mrSges(4,3) * t463 + Ifges(4,1) * t521 + Ifges(4,4) * t520 + Ifges(4,5) * t543 - pkin(10) * t424 + t414 * t566 - t417 * t563 + t513 * t533 - t514 * t546;
t515 = Ifges(4,1) * t534 + Ifges(4,4) * t533 + Ifges(4,5) * t546;
t570 = mrSges(5,1) * t455 - mrSges(5,2) * t456 + Ifges(5,5) * t497 + Ifges(5,6) * t496 + Ifges(5,3) * t517 + pkin(4) * t440 + pkin(11) * t434 + t431 * t618 + t562 * t432 + t527 * t491 - t526 * t492;
t409 = -mrSges(4,1) * t477 + mrSges(4,3) * t464 + Ifges(4,4) * t521 + Ifges(4,2) * t520 + Ifges(4,6) * t543 - pkin(3) * t424 - t534 * t513 + t546 * t515 - t570;
t576 = pkin(9) * t416 + t408 * t564 + t409 * t567;
t552 = t580 * qJD(1);
t551 = t581 * qJD(1);
t547 = t585 * t596;
t535 = -t548 * t558 + t589;
t525 = -g(3) * t607 + t592;
t524 = -t549 * t556 + t583;
t415 = m(3) * t525 + t580 * qJDD(1) + (t547 * t604 - t551 * t561) * qJD(1) + t416;
t411 = m(3) * t535 + (t585 * qJDD(1) + (t551 * t556 - t552 * t559) * qJD(1)) * t558 + t412;
t410 = m(3) * t524 + t581 * qJDD(1) + (-t547 * t607 + t552 * t561) * qJD(1) + t413;
t407 = mrSges(4,1) * t463 - mrSges(4,2) * t464 + Ifges(4,5) * t521 + Ifges(4,6) * t520 + Ifges(4,3) * t543 + pkin(3) * t571 + pkin(10) * t588 + t563 * t414 + t566 * t417 + t534 * t514 - t533 * t515;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t590 - mrSges(2,2) * t586 + (mrSges(3,1) * t524 - mrSges(3,2) * t525 + pkin(2) * t413 + t560 * t407 + pkin(1) * (t410 * t559 + t415 * t556) + Ifges(3,3) * qJDD(1) * t561 + t576 * t557) * t561 + (t556 * (mrSges(3,2) * t535 - mrSges(3,3) * t524 + t408 * t567 - t409 * t564 - t412 * t616) + t559 * (-mrSges(3,1) * t535 + mrSges(3,3) * t525 - pkin(2) * t412 - t557 * t407) - pkin(1) * t411 + qJ(2) * (-t410 * t556 + t415 * t559) + (-t413 * t617 + t559 * t576) * t560 + ((Ifges(3,2) * t559 ^ 2 + (Ifges(3,1) * t556 + Ifges(3,4) * t619) * t556) * t558 + 0.2e1 * t561 * (Ifges(3,5) * t556 + Ifges(3,6) * t559)) * qJDD(1)) * t558; t411; t407; t570; t620; t442;];
tauJ  = t1;
