% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRRRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-05-08 16:27
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRRRR9_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR9_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR9_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR9_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR9_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR9_invdynJ_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR9_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR9_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR9_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 15:49:45
% EndTime: 2019-05-08 15:50:50
% DurationCPUTime: 40.79s
% Computational Cost: add. (664766->381), mult. (1635018->510), div. (0->0), fcn. (1392501->16), ass. (0->170)
t571 = cos(pkin(6));
t566 = qJD(1) * t571 + qJD(2);
t568 = sin(pkin(7));
t570 = cos(pkin(7));
t569 = sin(pkin(6));
t582 = cos(qJ(2));
t604 = qJD(1) * t582;
t600 = t569 * t604;
t553 = (t566 * t568 + t570 * t600) * pkin(10);
t576 = sin(qJ(2));
t606 = qJD(1) * t569;
t617 = pkin(10) * t568;
t557 = (-pkin(2) * t582 - t576 * t617) * t606;
t603 = qJD(1) * qJD(2);
t563 = (qJDD(1) * t576 + t582 * t603) * t569;
t565 = qJDD(1) * t571 + qJDD(2);
t584 = qJD(1) ^ 2;
t577 = sin(qJ(1));
t583 = cos(qJ(1));
t595 = -g(1) * t583 - g(2) * t577;
t618 = pkin(9) * t569;
t561 = -pkin(1) * t584 + qJDD(1) * t618 + t595;
t599 = t577 * g(1) - g(2) * t583;
t560 = qJDD(1) * pkin(1) + t584 * t618 + t599;
t614 = t560 * t571;
t596 = -t576 * t561 + t582 * t614;
t605 = qJD(1) * t576;
t616 = pkin(10) * t570;
t509 = -t563 * t616 + t565 * pkin(2) + t566 * t553 + (-g(3) * t582 - t557 * t605) * t569 + t596;
t601 = t569 * t605;
t556 = pkin(2) * t566 - t601 * t616;
t564 = (qJDD(1) * t582 - t576 * t603) * t569;
t593 = t564 * t570 + t565 * t568;
t607 = t582 * t561 + t576 * t614;
t510 = -t566 * t556 + (-g(3) * t576 + t557 * t604) * t569 + t593 * pkin(10) + t607;
t615 = t571 * g(3);
t515 = -t563 * t617 - t564 * pkin(2) - t615 + (-t560 + (-t553 * t582 + t556 * t576) * qJD(1)) * t569;
t575 = sin(qJ(3));
t581 = cos(qJ(3));
t480 = -t575 * t510 + (t509 * t570 + t515 * t568) * t581;
t608 = t570 * t582;
t613 = t568 * t575;
t544 = t566 * t613 + (t575 * t608 + t576 * t581) * t606;
t527 = -t544 * qJD(3) - t575 * t563 + t581 * t593;
t612 = t568 * t581;
t543 = (-t575 * t576 + t581 * t608) * t606 + t566 * t612;
t611 = t569 * t576;
t610 = t569 * t582;
t609 = t570 * t575;
t481 = t509 * t609 + t581 * t510 + t515 * t613;
t530 = -pkin(3) * t543 - pkin(11) * t544;
t545 = -t564 * t568 + t565 * t570 + qJDD(3);
t554 = t566 * t570 - t568 * t600 + qJD(3);
t552 = t554 ^ 2;
t469 = -pkin(3) * t552 + pkin(11) * t545 + t530 * t543 + t481;
t488 = -t568 * t509 + t570 * t515;
t528 = t543 * qJD(3) + t581 * t563 + t575 * t593;
t471 = (-t543 * t554 - t528) * pkin(11) + (t544 * t554 - t527) * pkin(3) + t488;
t574 = sin(qJ(4));
t580 = cos(qJ(4));
t459 = t580 * t469 + t574 * t471;
t534 = -t574 * t544 + t554 * t580;
t535 = t544 * t580 + t554 * t574;
t512 = -pkin(4) * t534 - pkin(12) * t535;
t526 = qJDD(4) - t527;
t542 = qJD(4) - t543;
t541 = t542 ^ 2;
t454 = -pkin(4) * t541 + pkin(12) * t526 + t512 * t534 + t459;
t468 = -t545 * pkin(3) - t552 * pkin(11) + t544 * t530 - t480;
t496 = -t535 * qJD(4) - t574 * t528 + t545 * t580;
t497 = qJD(4) * t534 + t528 * t580 + t545 * t574;
t457 = (-t534 * t542 - t497) * pkin(12) + (t535 * t542 - t496) * pkin(4) + t468;
t573 = sin(qJ(5));
t579 = cos(qJ(5));
t449 = -t573 * t454 + t579 * t457;
t518 = -t535 * t573 + t542 * t579;
t479 = qJD(5) * t518 + t497 * t579 + t526 * t573;
t495 = qJDD(5) - t496;
t519 = t535 * t579 + t542 * t573;
t533 = qJD(5) - t534;
t447 = (t518 * t533 - t479) * pkin(13) + (t518 * t519 + t495) * pkin(5) + t449;
t450 = t579 * t454 + t573 * t457;
t478 = -qJD(5) * t519 - t497 * t573 + t526 * t579;
t501 = pkin(5) * t533 - pkin(13) * t519;
t517 = t518 ^ 2;
t448 = -pkin(5) * t517 + pkin(13) * t478 - t501 * t533 + t450;
t572 = sin(qJ(6));
t578 = cos(qJ(6));
t445 = t447 * t578 - t448 * t572;
t489 = t518 * t578 - t519 * t572;
t464 = qJD(6) * t489 + t478 * t572 + t479 * t578;
t490 = t518 * t572 + t519 * t578;
t476 = -mrSges(7,1) * t489 + mrSges(7,2) * t490;
t531 = qJD(6) + t533;
t482 = -mrSges(7,2) * t531 + mrSges(7,3) * t489;
t493 = qJDD(6) + t495;
t441 = m(7) * t445 + mrSges(7,1) * t493 - mrSges(7,3) * t464 - t476 * t490 + t482 * t531;
t446 = t447 * t572 + t448 * t578;
t463 = -qJD(6) * t490 + t478 * t578 - t479 * t572;
t483 = mrSges(7,1) * t531 - mrSges(7,3) * t490;
t442 = m(7) * t446 - mrSges(7,2) * t493 + mrSges(7,3) * t463 + t476 * t489 - t483 * t531;
t433 = t578 * t441 + t572 * t442;
t491 = -mrSges(6,1) * t518 + mrSges(6,2) * t519;
t499 = -mrSges(6,2) * t533 + mrSges(6,3) * t518;
t431 = m(6) * t449 + mrSges(6,1) * t495 - mrSges(6,3) * t479 - t491 * t519 + t499 * t533 + t433;
t500 = mrSges(6,1) * t533 - mrSges(6,3) * t519;
t597 = -t441 * t572 + t578 * t442;
t432 = m(6) * t450 - mrSges(6,2) * t495 + mrSges(6,3) * t478 + t491 * t518 - t500 * t533 + t597;
t429 = -t431 * t573 + t579 * t432;
t511 = -mrSges(5,1) * t534 + mrSges(5,2) * t535;
t521 = mrSges(5,1) * t542 - mrSges(5,3) * t535;
t427 = m(5) * t459 - mrSges(5,2) * t526 + mrSges(5,3) * t496 + t511 * t534 - t521 * t542 + t429;
t458 = -t574 * t469 + t471 * t580;
t453 = -pkin(4) * t526 - pkin(12) * t541 + t535 * t512 - t458;
t451 = -pkin(5) * t478 - pkin(13) * t517 + t501 * t519 + t453;
t589 = m(7) * t451 - t463 * mrSges(7,1) + mrSges(7,2) * t464 - t489 * t482 + t483 * t490;
t443 = -m(6) * t453 + t478 * mrSges(6,1) - mrSges(6,2) * t479 + t518 * t499 - t500 * t519 - t589;
t520 = -mrSges(5,2) * t542 + mrSges(5,3) * t534;
t437 = m(5) * t458 + mrSges(5,1) * t526 - mrSges(5,3) * t497 - t511 * t535 + t520 * t542 + t443;
t419 = t574 * t427 + t580 * t437;
t529 = -mrSges(4,1) * t543 + mrSges(4,2) * t544;
t537 = mrSges(4,1) * t554 - mrSges(4,3) * t544;
t598 = t580 * t427 - t437 * t574;
t416 = m(4) * t481 - mrSges(4,2) * t545 + mrSges(4,3) * t527 + t529 * t543 - t537 * t554 + t598;
t536 = -mrSges(4,2) * t554 + mrSges(4,3) * t543;
t418 = m(4) * t488 - mrSges(4,1) * t527 + mrSges(4,2) * t528 - t536 * t543 + t537 * t544 + t419;
t428 = t431 * t579 + t432 * t573;
t587 = -m(5) * t468 + t496 * mrSges(5,1) - mrSges(5,2) * t497 + t534 * t520 - t521 * t535 - t428;
t424 = m(4) * t480 + mrSges(4,1) * t545 - mrSges(4,3) * t528 - t529 * t544 + t536 * t554 + t587;
t407 = t416 * t613 + t570 * t418 + t424 * t612;
t411 = t581 * t416 - t424 * t575;
t408 = t570 * t581 * t424 + t416 * t609 - t418 * t568;
t472 = Ifges(7,5) * t490 + Ifges(7,6) * t489 + Ifges(7,3) * t531;
t474 = Ifges(7,1) * t490 + Ifges(7,4) * t489 + Ifges(7,5) * t531;
t434 = -mrSges(7,1) * t451 + mrSges(7,3) * t446 + Ifges(7,4) * t464 + Ifges(7,2) * t463 + Ifges(7,6) * t493 - t472 * t490 + t474 * t531;
t473 = Ifges(7,4) * t490 + Ifges(7,2) * t489 + Ifges(7,6) * t531;
t435 = mrSges(7,2) * t451 - mrSges(7,3) * t445 + Ifges(7,1) * t464 + Ifges(7,4) * t463 + Ifges(7,5) * t493 + t472 * t489 - t473 * t531;
t484 = Ifges(6,5) * t519 + Ifges(6,6) * t518 + Ifges(6,3) * t533;
t486 = Ifges(6,1) * t519 + Ifges(6,4) * t518 + Ifges(6,5) * t533;
t420 = -mrSges(6,1) * t453 + mrSges(6,3) * t450 + Ifges(6,4) * t479 + Ifges(6,2) * t478 + Ifges(6,6) * t495 - pkin(5) * t589 + pkin(13) * t597 + t578 * t434 + t572 * t435 - t519 * t484 + t533 * t486;
t485 = Ifges(6,4) * t519 + Ifges(6,2) * t518 + Ifges(6,6) * t533;
t421 = mrSges(6,2) * t453 - mrSges(6,3) * t449 + Ifges(6,1) * t479 + Ifges(6,4) * t478 + Ifges(6,5) * t495 - pkin(13) * t433 - t434 * t572 + t435 * t578 + t484 * t518 - t485 * t533;
t502 = Ifges(5,5) * t535 + Ifges(5,6) * t534 + Ifges(5,3) * t542;
t503 = Ifges(5,4) * t535 + Ifges(5,2) * t534 + Ifges(5,6) * t542;
t409 = mrSges(5,2) * t468 - mrSges(5,3) * t458 + Ifges(5,1) * t497 + Ifges(5,4) * t496 + Ifges(5,5) * t526 - pkin(12) * t428 - t420 * t573 + t421 * t579 + t502 * t534 - t503 * t542;
t504 = Ifges(5,1) * t535 + Ifges(5,4) * t534 + Ifges(5,5) * t542;
t588 = -mrSges(7,1) * t445 + mrSges(7,2) * t446 - Ifges(7,5) * t464 - Ifges(7,6) * t463 - Ifges(7,3) * t493 - t490 * t473 + t489 * t474;
t585 = mrSges(6,1) * t449 - mrSges(6,2) * t450 + Ifges(6,5) * t479 + Ifges(6,6) * t478 + Ifges(6,3) * t495 + pkin(5) * t433 + t519 * t485 - t518 * t486 - t588;
t412 = -mrSges(5,1) * t468 + mrSges(5,3) * t459 + Ifges(5,4) * t497 + Ifges(5,2) * t496 + Ifges(5,6) * t526 - pkin(4) * t428 - t535 * t502 + t542 * t504 - t585;
t522 = Ifges(4,5) * t544 + Ifges(4,6) * t543 + Ifges(4,3) * t554;
t523 = Ifges(4,4) * t544 + Ifges(4,2) * t543 + Ifges(4,6) * t554;
t404 = mrSges(4,2) * t488 - mrSges(4,3) * t480 + Ifges(4,1) * t528 + Ifges(4,4) * t527 + Ifges(4,5) * t545 - pkin(11) * t419 + t409 * t580 - t412 * t574 + t522 * t543 - t523 * t554;
t524 = Ifges(4,1) * t544 + Ifges(4,4) * t543 + Ifges(4,5) * t554;
t586 = mrSges(5,1) * t458 - mrSges(5,2) * t459 + Ifges(5,5) * t497 + Ifges(5,6) * t496 + Ifges(5,3) * t526 + pkin(4) * t443 + pkin(12) * t429 + t579 * t420 + t573 * t421 + t535 * t503 - t534 * t504;
t405 = -mrSges(4,1) * t488 + mrSges(4,3) * t481 + Ifges(4,4) * t528 + Ifges(4,2) * t527 + Ifges(4,6) * t545 - pkin(3) * t419 - t544 * t522 + t554 * t524 - t586;
t590 = pkin(10) * t411 + t404 * t575 + t405 * t581;
t562 = (-mrSges(3,1) * t582 + mrSges(3,2) * t576) * t606;
t559 = -mrSges(3,2) * t566 + mrSges(3,3) * t600;
t558 = mrSges(3,1) * t566 - mrSges(3,3) * t601;
t549 = -t569 * t560 - t615;
t548 = Ifges(3,5) * t566 + (Ifges(3,1) * t576 + Ifges(3,4) * t582) * t606;
t547 = Ifges(3,6) * t566 + (Ifges(3,4) * t576 + Ifges(3,2) * t582) * t606;
t546 = Ifges(3,3) * t566 + (Ifges(3,5) * t576 + Ifges(3,6) * t582) * t606;
t539 = -g(3) * t611 + t607;
t538 = -g(3) * t610 + t596;
t410 = m(3) * t539 - mrSges(3,2) * t565 + mrSges(3,3) * t564 - t558 * t566 + t562 * t600 + t411;
t406 = m(3) * t538 + mrSges(3,1) * t565 - mrSges(3,3) * t563 + t559 * t566 - t562 * t601 + t408;
t403 = mrSges(4,1) * t480 - mrSges(4,2) * t481 + Ifges(4,5) * t528 + Ifges(4,6) * t527 + Ifges(4,3) * t545 + pkin(3) * t587 + pkin(11) * t598 + t574 * t409 + t580 * t412 + t544 * t523 - t543 * t524;
t402 = mrSges(3,1) * t538 - mrSges(3,2) * t539 + Ifges(3,5) * t563 + Ifges(3,6) * t564 + Ifges(3,3) * t565 + pkin(2) * t408 + t570 * t403 + (t547 * t576 - t548 * t582) * t606 + t590 * t568;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t599 - mrSges(2,2) * t595 + (t546 * t600 + mrSges(3,2) * t549 - mrSges(3,3) * t538 + Ifges(3,1) * t563 + Ifges(3,4) * t564 + Ifges(3,5) * t565 + t581 * t404 - t575 * t405 - t566 * t547 + (-t407 * t568 - t408 * t570) * pkin(10)) * t611 + (-mrSges(3,1) * t549 + mrSges(3,3) * t539 + Ifges(3,4) * t563 + Ifges(3,2) * t564 + Ifges(3,6) * t565 - pkin(2) * t407 - t568 * t403 - t546 * t601 + t566 * t548 + t570 * t590) * t610 + t571 * t402 + pkin(1) * ((t406 * t582 + t410 * t576) * t571 + (-m(3) * t549 + t564 * mrSges(3,1) - t563 * mrSges(3,2) + (-t558 * t576 + t559 * t582) * t606 - t407) * t569) + (-t406 * t576 + t410 * t582) * t618; t402; t403; t586; t585; -t588;];
tauJ  = t1;
