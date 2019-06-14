% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRPPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta5]';
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
% Datum: 2019-05-07 07:09
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRPPR10_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR10_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR10_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR10_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR10_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR10_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR10_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR10_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR10_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 06:59:36
% EndTime: 2019-05-07 06:59:45
% DurationCPUTime: 6.47s
% Computational Cost: add. (88913->342), mult. (194122->429), div. (0->0), fcn. (150632->12), ass. (0->144)
t598 = Ifges(4,1) + Ifges(5,2);
t591 = Ifges(4,5) - Ifges(5,4);
t597 = -Ifges(4,2) - Ifges(5,3);
t590 = Ifges(4,6) - Ifges(5,5);
t589 = -Ifges(5,6) - Ifges(4,4);
t596 = Ifges(4,3) + Ifges(5,1);
t548 = sin(pkin(6));
t553 = sin(qJ(2));
t556 = cos(qJ(2));
t575 = qJD(1) * qJD(2);
t535 = (-qJDD(1) * t556 + t553 * t575) * t548;
t578 = qJD(1) * t548;
t533 = (-pkin(2) * t556 - pkin(9) * t553) * t578;
t550 = cos(pkin(6));
t544 = qJD(1) * t550 + qJD(2);
t542 = t544 ^ 2;
t543 = qJDD(1) * t550 + qJDD(2);
t577 = qJD(1) * t556;
t558 = qJD(1) ^ 2;
t554 = sin(qJ(1));
t557 = cos(qJ(1));
t569 = -g(1) * t557 - g(2) * t554;
t593 = pkin(8) * t548;
t531 = -pkin(1) * t558 + qJDD(1) * t593 + t569;
t572 = t554 * g(1) - g(2) * t557;
t530 = qJDD(1) * pkin(1) + t558 * t593 + t572;
t586 = t530 * t550;
t579 = t556 * t531 + t553 * t586;
t465 = -pkin(2) * t542 + pkin(9) * t543 + (-g(3) * t553 + t533 * t577) * t548 + t579;
t534 = (qJDD(1) * t553 + t556 * t575) * t548;
t592 = g(3) * t550;
t466 = pkin(2) * t535 - pkin(9) * t534 - t592 + (-t530 + (pkin(2) * t553 - pkin(9) * t556) * t544 * qJD(1)) * t548;
t552 = sin(qJ(3));
t594 = cos(qJ(3));
t449 = t594 * t465 + t552 * t466;
t574 = t553 * t578;
t521 = -t544 * t594 + t552 * t574;
t522 = t552 * t544 + t574 * t594;
t496 = pkin(3) * t521 - qJ(4) * t522;
t527 = qJDD(3) + t535;
t573 = t548 * t577;
t540 = -qJD(3) + t573;
t539 = t540 ^ 2;
t439 = pkin(3) * t539 - t527 * qJ(4) + 0.2e1 * qJD(4) * t540 + t521 * t496 - t449;
t448 = -t552 * t465 + t466 * t594;
t440 = -t527 * pkin(3) - t539 * qJ(4) + t522 * t496 + qJDD(4) - t448;
t493 = -t521 * qJD(3) + t534 * t594 + t552 * t543;
t587 = t521 * t540;
t434 = (t521 * t522 - t527) * qJ(5) + (t493 - t587) * pkin(4) + t440;
t492 = qJD(3) * t522 + t534 * t552 - t543 * t594;
t508 = pkin(4) * t522 + qJ(5) * t540;
t520 = t521 ^ 2;
t584 = t548 * t556;
t494 = -g(3) * t584 - t553 * t531 + t556 * t586;
t464 = -pkin(2) * t543 - pkin(9) * t542 + t533 * t574 - t494;
t561 = (-t493 - t587) * qJ(4) + t464 + (-t540 * pkin(3) - 0.2e1 * qJD(4)) * t522;
t438 = -pkin(4) * t520 - t508 * t522 + (pkin(3) + qJ(5)) * t492 + t561;
t547 = sin(pkin(11));
t549 = cos(pkin(11));
t505 = t521 * t547 - t540 * t549;
t428 = -0.2e1 * qJD(5) * t505 + t549 * t434 - t438 * t547;
t471 = t492 * t547 + t527 * t549;
t504 = t521 * t549 + t540 * t547;
t426 = (t504 * t522 - t471) * pkin(10) + (t504 * t505 + t493) * pkin(5) + t428;
t429 = 0.2e1 * qJD(5) * t504 + t547 * t434 + t549 * t438;
t470 = t492 * t549 - t527 * t547;
t478 = pkin(5) * t522 - pkin(10) * t505;
t503 = t504 ^ 2;
t427 = -pkin(5) * t503 + pkin(10) * t470 - t478 * t522 + t429;
t551 = sin(qJ(6));
t555 = cos(qJ(6));
t425 = t426 * t551 + t427 * t555;
t437 = -pkin(4) * t492 - qJ(5) * t520 - t540 * t508 + qJDD(5) - t439;
t431 = -pkin(5) * t470 - pkin(10) * t503 + t478 * t505 + t437;
t468 = t504 * t551 + t505 * t555;
t446 = -qJD(6) * t468 + t470 * t555 - t471 * t551;
t467 = t504 * t555 - t505 * t551;
t447 = qJD(6) * t467 + t470 * t551 + t471 * t555;
t519 = qJD(6) + t522;
t450 = Ifges(7,5) * t468 + Ifges(7,6) * t467 + Ifges(7,3) * t519;
t452 = Ifges(7,1) * t468 + Ifges(7,4) * t467 + Ifges(7,5) * t519;
t489 = qJDD(6) + t493;
t411 = -mrSges(7,1) * t431 + mrSges(7,3) * t425 + Ifges(7,4) * t447 + Ifges(7,2) * t446 + Ifges(7,6) * t489 - t450 * t468 + t452 * t519;
t424 = t426 * t555 - t427 * t551;
t451 = Ifges(7,4) * t468 + Ifges(7,2) * t467 + Ifges(7,6) * t519;
t412 = mrSges(7,2) * t431 - mrSges(7,3) * t424 + Ifges(7,1) * t447 + Ifges(7,4) * t446 + Ifges(7,5) * t489 + t450 * t467 - t451 * t519;
t457 = Ifges(6,5) * t505 + Ifges(6,6) * t504 + Ifges(6,3) * t522;
t459 = Ifges(6,1) * t505 + Ifges(6,4) * t504 + Ifges(6,5) * t522;
t455 = -mrSges(7,2) * t519 + mrSges(7,3) * t467;
t456 = mrSges(7,1) * t519 - mrSges(7,3) * t468;
t565 = m(7) * t431 - mrSges(7,1) * t446 + t447 * mrSges(7,2) - t455 * t467 + t468 * t456;
t454 = -mrSges(7,1) * t467 + mrSges(7,2) * t468;
t419 = m(7) * t424 + mrSges(7,1) * t489 - mrSges(7,3) * t447 - t454 * t468 + t455 * t519;
t420 = m(7) * t425 - mrSges(7,2) * t489 + mrSges(7,3) * t446 + t454 * t467 - t456 * t519;
t570 = -t419 * t551 + t555 * t420;
t398 = -mrSges(6,1) * t437 + mrSges(6,3) * t429 + Ifges(6,4) * t471 + Ifges(6,2) * t470 + Ifges(6,6) * t493 - pkin(5) * t565 + pkin(10) * t570 + t555 * t411 + t551 * t412 - t505 * t457 + t522 * t459;
t410 = t555 * t419 + t551 * t420;
t458 = Ifges(6,4) * t505 + Ifges(6,2) * t504 + Ifges(6,6) * t522;
t399 = mrSges(6,2) * t437 - mrSges(6,3) * t428 + Ifges(6,1) * t471 + Ifges(6,4) * t470 + Ifges(6,5) * t493 - pkin(10) * t410 - t411 * t551 + t412 * t555 + t457 * t504 - t458 * t522;
t509 = mrSges(5,1) * t521 + mrSges(5,3) * t540;
t472 = -mrSges(6,1) * t504 + mrSges(6,2) * t505;
t476 = -mrSges(6,2) * t522 + mrSges(6,3) * t504;
t408 = m(6) * t428 + mrSges(6,1) * t493 - mrSges(6,3) * t471 - t472 * t505 + t476 * t522 + t410;
t477 = mrSges(6,1) * t522 - mrSges(6,3) * t505;
t409 = m(6) * t429 - mrSges(6,2) * t493 + mrSges(6,3) * t470 + t472 * t504 - t477 * t522 + t570;
t405 = t408 * t549 + t409 * t547;
t498 = -mrSges(5,2) * t521 - mrSges(5,3) * t522;
t564 = -m(5) * t440 - t493 * mrSges(5,1) - t522 * t498 - t405;
t404 = mrSges(5,2) * t527 - t509 * t540 - t564;
t422 = m(6) * t437 - mrSges(6,1) * t470 + t471 * mrSges(6,2) - t476 * t504 + t505 * t477 + t565;
t510 = mrSges(5,1) * t522 - mrSges(5,2) * t540;
t560 = -m(5) * t439 + t527 * mrSges(5,3) - t540 * t510 + t422;
t580 = t589 * t521 + t598 * t522 - t591 * t540;
t581 = t597 * t521 - t589 * t522 - t590 * t540;
t595 = -t492 * t590 + t493 * t591 + t521 * t580 + t522 * t581 + t596 * t527 + mrSges(4,1) * t448 - mrSges(4,2) * t449 + mrSges(5,2) * t440 - mrSges(5,3) * t439 - pkin(3) * t404 + qJ(4) * (-mrSges(5,1) * t492 - t498 * t521 + t560) - qJ(5) * t405 - t547 * t398 + t549 * t399;
t585 = t548 * t553;
t583 = -t547 * t408 + t549 * t409;
t497 = mrSges(4,1) * t521 + mrSges(4,2) * t522;
t506 = mrSges(4,2) * t540 - mrSges(4,3) * t521;
t402 = m(4) * t448 - mrSges(4,3) * t493 - t497 * t522 + (-t506 + t509) * t540 + (mrSges(4,1) - mrSges(5,2)) * t527 + t564;
t507 = -mrSges(4,1) * t540 - mrSges(4,3) * t522;
t415 = t560 + (-t497 - t498) * t521 + t507 * t540 + (-mrSges(4,3) - mrSges(5,1)) * t492 + m(4) * t449 - mrSges(4,2) * t527;
t397 = t594 * t402 + t552 * t415;
t582 = t590 * t521 - t591 * t522 + t596 * t540;
t571 = -t402 * t552 + t594 * t415;
t441 = pkin(3) * t492 + t561;
t568 = -m(5) * t441 + t492 * mrSges(5,2) + t521 * t509 - t583;
t563 = mrSges(7,1) * t424 - mrSges(7,2) * t425 + Ifges(7,5) * t447 + Ifges(7,6) * t446 + Ifges(7,3) * t489 + t468 * t451 - t467 * t452;
t562 = -m(4) * t464 - t492 * mrSges(4,1) - t521 * t506 + (-t507 + t510) * t522 + (-mrSges(4,2) + mrSges(5,3)) * t493 + t568;
t532 = (-mrSges(3,1) * t556 + mrSges(3,2) * t553) * t578;
t529 = -mrSges(3,2) * t544 + mrSges(3,3) * t573;
t528 = mrSges(3,1) * t544 - mrSges(3,3) * t574;
t515 = -t530 * t548 - t592;
t514 = Ifges(3,5) * t544 + (Ifges(3,1) * t553 + Ifges(3,4) * t556) * t578;
t513 = Ifges(3,6) * t544 + (Ifges(3,4) * t553 + Ifges(3,2) * t556) * t578;
t512 = Ifges(3,3) * t544 + (Ifges(3,5) * t553 + Ifges(3,6) * t556) * t578;
t495 = -g(3) * t585 + t579;
t403 = -mrSges(5,3) * t493 - t510 * t522 - t568;
t400 = m(3) * t494 + mrSges(3,1) * t543 - mrSges(3,3) * t534 + t529 * t544 - t532 * t574 + t562;
t396 = m(3) * t495 - mrSges(3,2) * t543 - mrSges(3,3) * t535 - t528 * t544 + t532 * t573 + t571;
t395 = pkin(4) * t405 + t563 - qJ(4) * t403 - t504 * t459 + t505 * t458 + Ifges(6,6) * t470 + Ifges(6,5) * t471 + mrSges(4,2) * t464 - mrSges(4,3) * t448 + mrSges(5,1) * t440 - mrSges(5,3) * t441 + mrSges(6,1) * t428 - mrSges(6,2) * t429 + t589 * t492 + (Ifges(6,3) + t598) * t493 + t582 * t521 + t591 * t527 + t581 * t540 + pkin(5) * t410;
t394 = -mrSges(4,1) * t464 - mrSges(5,1) * t439 + mrSges(5,2) * t441 + mrSges(4,3) * t449 - pkin(3) * t403 + pkin(4) * t422 - qJ(5) * t583 - t549 * t398 - t547 * t399 + t597 * t492 - t589 * t493 + t582 * t522 + t590 * t527 - t580 * t540;
t393 = Ifges(3,5) * t534 - Ifges(3,6) * t535 + Ifges(3,3) * t543 + mrSges(3,1) * t494 - mrSges(3,2) * t495 + t552 * t395 + t594 * t394 + pkin(2) * t562 + pkin(9) * t571 + (t513 * t553 - t514 * t556) * t578;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t572 - mrSges(2,2) * t569 + (mrSges(3,2) * t515 - mrSges(3,3) * t494 + Ifges(3,1) * t534 - Ifges(3,4) * t535 + Ifges(3,5) * t543 - pkin(9) * t397 - t552 * t394 + t395 * t594 + t512 * t573 - t544 * t513) * t585 + (-mrSges(3,1) * t515 + mrSges(3,3) * t495 + Ifges(3,4) * t534 - Ifges(3,2) * t535 + Ifges(3,6) * t543 - pkin(2) * t397 - t512 * t574 + t544 * t514 - t595) * t584 + t550 * t393 + pkin(1) * ((t396 * t553 + t400 * t556) * t550 + (-m(3) * t515 - t535 * mrSges(3,1) - t534 * mrSges(3,2) + (-t528 * t553 + t529 * t556) * t578 - t397) * t548) + (t396 * t556 - t400 * t553) * t593; t393; t595; t404; t422; t563;];
tauJ  = t1;
