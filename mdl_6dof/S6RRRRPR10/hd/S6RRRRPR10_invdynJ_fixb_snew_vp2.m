% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRRPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
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
% Datum: 2019-05-07 23:04
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRRPR10_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR10_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR10_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR10_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR10_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR10_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR10_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR10_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR10_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 22:51:05
% EndTime: 2019-05-07 22:51:22
% DurationCPUTime: 7.86s
% Computational Cost: add. (99858->342), mult. (214030->427), div. (0->0), fcn. (171740->12), ass. (0->148)
t599 = Ifges(5,4) + Ifges(6,6);
t611 = -Ifges(5,2) - Ifges(6,3);
t607 = Ifges(5,6) - Ifges(6,5);
t560 = cos(pkin(6));
t556 = qJD(1) * t560 + qJD(2);
t563 = sin(qJ(3));
t567 = cos(qJ(3));
t564 = sin(qJ(2));
t559 = sin(pkin(6));
t590 = qJD(1) * t559;
t587 = t564 * t590;
t535 = t556 * t567 - t563 * t587;
t536 = t556 * t563 + t567 * t587;
t562 = sin(qJ(4));
t602 = cos(qJ(4));
t518 = -t535 * t602 + t536 * t562;
t568 = cos(qJ(2));
t589 = qJD(1) * t568;
t586 = t559 * t589;
t552 = qJD(3) - t586;
t550 = -qJD(4) - t552;
t498 = mrSges(6,1) * t518 + mrSges(6,3) * t550;
t588 = qJDD(1) * t559;
t548 = -qJD(2) * t587 + t568 * t588;
t540 = qJDD(3) - t548;
t539 = qJDD(4) + t540;
t546 = (-pkin(2) * t568 - pkin(9) * t564) * t590;
t554 = t556 ^ 2;
t555 = qJDD(1) * t560 + qJDD(2);
t570 = qJD(1) ^ 2;
t565 = sin(qJ(1));
t569 = cos(qJ(1));
t581 = -g(1) * t569 - g(2) * t565;
t544 = -pkin(1) * t570 + pkin(8) * t588 + t581;
t585 = t565 * g(1) - g(2) * t569;
t601 = pkin(8) * t559;
t543 = qJDD(1) * pkin(1) + t570 * t601 + t585;
t597 = t543 * t560;
t591 = t568 * t544 + t564 * t597;
t493 = -t554 * pkin(2) + t555 * pkin(9) + (-g(3) * t564 + t546 * t589) * t559 + t591;
t547 = (qJD(2) * t589 + qJDD(1) * t564) * t559;
t600 = t560 * g(3);
t494 = -t548 * pkin(2) - t547 * pkin(9) - t600 + (-t543 + (pkin(2) * t564 - pkin(9) * t568) * t556 * qJD(1)) * t559;
t452 = -t563 * t493 + t567 * t494;
t512 = qJD(3) * t535 + t547 * t567 + t555 * t563;
t441 = (t535 * t552 - t512) * pkin(10) + (t535 * t536 + t540) * pkin(3) + t452;
t453 = t567 * t493 + t563 * t494;
t511 = -qJD(3) * t536 - t547 * t563 + t555 * t567;
t523 = pkin(3) * t552 - pkin(10) * t536;
t534 = t535 ^ 2;
t444 = -pkin(3) * t534 + pkin(10) * t511 - t523 * t552 + t453;
t438 = t441 * t602 - t562 * t444;
t519 = t562 * t535 + t536 * t602;
t486 = pkin(4) * t518 - qJ(5) * t519;
t549 = t550 ^ 2;
t434 = -t539 * pkin(4) - t549 * qJ(5) + t519 * t486 + qJDD(5) - t438;
t469 = -t518 * qJD(4) + t562 * t511 + t512 * t602;
t598 = t518 * t550;
t428 = (t518 * t519 - t539) * pkin(11) + (t469 - t598) * pkin(5) + t434;
t468 = t519 * qJD(4) - t511 * t602 + t562 * t512;
t502 = pkin(5) * t519 + pkin(11) * t550;
t517 = t518 ^ 2;
t595 = t559 * t568;
t513 = -g(3) * t595 - t564 * t544 + t568 * t597;
t492 = -t555 * pkin(2) - t554 * pkin(9) + t546 * t587 - t513;
t451 = -t511 * pkin(3) - t534 * pkin(10) + t536 * t523 + t492;
t603 = -2 * qJD(5);
t573 = (-t469 - t598) * qJ(5) + t451 + (-t550 * pkin(4) + t603) * t519;
t431 = -t519 * t502 - t517 * pkin(5) + t573 + (pkin(4) + pkin(11)) * t468;
t561 = sin(qJ(6));
t566 = cos(qJ(6));
t426 = t428 * t566 - t431 * t561;
t496 = t518 * t566 + t550 * t561;
t449 = qJD(6) * t496 + t468 * t561 + t539 * t566;
t467 = qJDD(6) + t469;
t497 = t518 * t561 - t550 * t566;
t474 = -mrSges(7,1) * t496 + mrSges(7,2) * t497;
t516 = qJD(6) + t519;
t475 = -mrSges(7,2) * t516 + mrSges(7,3) * t496;
t423 = m(7) * t426 + mrSges(7,1) * t467 - mrSges(7,3) * t449 - t474 * t497 + t475 * t516;
t427 = t428 * t561 + t431 * t566;
t448 = -qJD(6) * t497 + t468 * t566 - t539 * t561;
t476 = mrSges(7,1) * t516 - mrSges(7,3) * t497;
t424 = m(7) * t427 - mrSges(7,2) * t467 + mrSges(7,3) * t448 + t474 * t496 - t476 * t516;
t412 = t566 * t423 + t561 * t424;
t488 = -mrSges(6,2) * t518 - mrSges(6,3) * t519;
t579 = -m(6) * t434 - t469 * mrSges(6,1) - t519 * t488 - t412;
t411 = t539 * mrSges(6,2) - t550 * t498 - t579;
t439 = t562 * t441 + t602 * t444;
t578 = -t549 * pkin(4) + t539 * qJ(5) - t518 * t486 + t439;
t430 = -t468 * pkin(5) - t517 * pkin(11) + (t603 - t502) * t550 + t578;
t454 = Ifges(7,5) * t497 + Ifges(7,6) * t496 + Ifges(7,3) * t516;
t456 = Ifges(7,1) * t497 + Ifges(7,4) * t496 + Ifges(7,5) * t516;
t414 = -mrSges(7,1) * t430 + mrSges(7,3) * t427 + Ifges(7,4) * t449 + Ifges(7,2) * t448 + Ifges(7,6) * t467 - t454 * t497 + t456 * t516;
t455 = Ifges(7,4) * t497 + Ifges(7,2) * t496 + Ifges(7,6) * t516;
t415 = mrSges(7,2) * t430 - mrSges(7,3) * t426 + Ifges(7,1) * t449 + Ifges(7,4) * t448 + Ifges(7,5) * t467 + t454 * t496 - t455 * t516;
t432 = 0.2e1 * qJD(5) * t550 - t578;
t499 = mrSges(6,1) * t519 - mrSges(6,2) * t550;
t580 = -m(7) * t430 + t448 * mrSges(7,1) - t449 * mrSges(7,2) + t496 * t475 - t497 * t476;
t576 = -m(6) * t432 + t539 * mrSges(6,3) - t550 * t499 - t580;
t608 = Ifges(5,5) - Ifges(6,4);
t609 = Ifges(5,1) + Ifges(6,2);
t592 = t599 * t518 - t519 * t609 + t608 * t550;
t605 = t611 * t518 + t599 * t519 - t607 * t550;
t606 = Ifges(5,3) + Ifges(6,1);
t610 = -mrSges(5,2) * t439 - mrSges(6,3) * t432 - pkin(4) * t411 - pkin(11) * t412 - t561 * t414 + t566 * t415 + qJ(5) * (-t518 * t488 + t576) + mrSges(6,2) * t434 + mrSges(5,1) * t438 + t606 * t539 + t605 * t519 + t608 * t469 + (-qJ(5) * mrSges(6,1) - t607) * t468 - t592 * t518;
t487 = mrSges(5,1) * t518 + mrSges(5,2) * t519;
t500 = mrSges(5,2) * t550 - mrSges(5,3) * t518;
t408 = m(5) * t438 - t469 * mrSges(5,3) - t519 * t487 + (t498 - t500) * t550 + (mrSges(5,1) - mrSges(6,2)) * t539 + t579;
t501 = -mrSges(5,1) * t550 - mrSges(5,3) * t519;
t418 = m(5) * t439 - t539 * mrSges(5,2) + t550 * t501 + (-t487 - t488) * t518 + (-mrSges(5,3) - mrSges(6,1)) * t468 + t576;
t405 = t602 * t408 + t562 * t418;
t506 = Ifges(4,4) * t536 + Ifges(4,2) * t535 + Ifges(4,6) * t552;
t507 = Ifges(4,1) * t536 + Ifges(4,4) * t535 + Ifges(4,5) * t552;
t604 = mrSges(4,1) * t452 - mrSges(4,2) * t453 + Ifges(4,5) * t512 + Ifges(4,6) * t511 + Ifges(4,3) * t540 + pkin(3) * t405 + t536 * t506 - t535 * t507 + t610;
t596 = t559 * t564;
t520 = -mrSges(4,1) * t535 + mrSges(4,2) * t536;
t521 = -mrSges(4,2) * t552 + mrSges(4,3) * t535;
t403 = m(4) * t452 + mrSges(4,1) * t540 - mrSges(4,3) * t512 - t520 * t536 + t521 * t552 + t405;
t522 = mrSges(4,1) * t552 - mrSges(4,3) * t536;
t582 = -t408 * t562 + t602 * t418;
t404 = m(4) * t453 - mrSges(4,2) * t540 + mrSges(4,3) * t511 + t520 * t535 - t522 * t552 + t582;
t398 = t567 * t403 + t563 * t404;
t594 = -t561 * t423 + t566 * t424;
t593 = t518 * t607 - t519 * t608 + t550 * t606;
t583 = -t403 * t563 + t567 * t404;
t436 = t468 * pkin(4) + t573;
t409 = m(6) * t436 - t468 * mrSges(6,2) - t469 * mrSges(6,3) - t518 * t498 - t519 * t499 + t594;
t577 = mrSges(7,1) * t426 - mrSges(7,2) * t427 + Ifges(7,5) * t449 + Ifges(7,6) * t448 + Ifges(7,3) * t467 + t497 * t455 - t496 * t456;
t575 = m(5) * t451 + t468 * mrSges(5,1) + t469 * mrSges(5,2) + t518 * t500 + t519 * t501 + t409;
t572 = -m(4) * t492 + t511 * mrSges(4,1) - t512 * mrSges(4,2) + t535 * t521 - t536 * t522 - t575;
t545 = (-mrSges(3,1) * t568 + mrSges(3,2) * t564) * t590;
t542 = -mrSges(3,2) * t556 + mrSges(3,3) * t586;
t541 = mrSges(3,1) * t556 - mrSges(3,3) * t587;
t527 = -t559 * t543 - t600;
t526 = Ifges(3,5) * t556 + (Ifges(3,1) * t564 + Ifges(3,4) * t568) * t590;
t525 = Ifges(3,6) * t556 + (Ifges(3,4) * t564 + Ifges(3,2) * t568) * t590;
t524 = Ifges(3,3) * t556 + (Ifges(3,5) * t564 + Ifges(3,6) * t568) * t590;
t514 = -g(3) * t596 + t591;
t505 = Ifges(4,5) * t536 + Ifges(4,6) * t535 + Ifges(4,3) * t552;
t406 = m(3) * t513 + t555 * mrSges(3,1) - t547 * mrSges(3,3) + t556 * t542 - t545 * t587 + t572;
t399 = mrSges(6,1) * t434 + mrSges(5,2) * t451 - mrSges(5,3) * t438 - mrSges(6,3) * t436 + pkin(5) * t412 - qJ(5) * t409 - t599 * t468 + t469 * t609 + t593 * t518 + t608 * t539 + t605 * t550 + t577;
t397 = m(3) * t514 - mrSges(3,2) * t555 + mrSges(3,3) * t548 - t541 * t556 + t545 * t586 + t583;
t396 = -mrSges(5,1) * t451 - mrSges(6,1) * t432 + mrSges(6,2) * t436 + mrSges(5,3) * t439 - pkin(4) * t409 - pkin(5) * t580 - pkin(11) * t594 - t566 * t414 - t561 * t415 + t611 * t468 + t599 * t469 + t593 * t519 + t607 * t539 + t592 * t550;
t395 = mrSges(4,2) * t492 - mrSges(4,3) * t452 + Ifges(4,1) * t512 + Ifges(4,4) * t511 + Ifges(4,5) * t540 - pkin(10) * t405 - t562 * t396 + t399 * t602 + t535 * t505 - t552 * t506;
t394 = -mrSges(4,1) * t492 + mrSges(4,3) * t453 + Ifges(4,4) * t512 + Ifges(4,2) * t511 + Ifges(4,6) * t540 - pkin(3) * t575 + pkin(10) * t582 + t396 * t602 + t562 * t399 - t536 * t505 + t552 * t507;
t393 = Ifges(3,5) * t547 + Ifges(3,6) * t548 + Ifges(3,3) * t555 + mrSges(3,1) * t513 - mrSges(3,2) * t514 + t563 * t395 + t567 * t394 + pkin(2) * t572 + pkin(9) * t583 + (t525 * t564 - t526 * t568) * t590;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t585 - mrSges(2,2) * t581 + (mrSges(3,2) * t527 - mrSges(3,3) * t513 + Ifges(3,1) * t547 + Ifges(3,4) * t548 + Ifges(3,5) * t555 - pkin(9) * t398 - t394 * t563 + t395 * t567 + t524 * t586 - t525 * t556) * t596 + (-mrSges(3,1) * t527 + mrSges(3,3) * t514 + Ifges(3,4) * t547 + Ifges(3,2) * t548 + Ifges(3,6) * t555 - pkin(2) * t398 - t524 * t587 + t556 * t526 - t604) * t595 + t560 * t393 + pkin(1) * ((t397 * t564 + t406 * t568) * t560 + (-m(3) * t527 + t548 * mrSges(3,1) - t547 * mrSges(3,2) + (-t541 * t564 + t542 * t568) * t590 - t398) * t559) + (t397 * t568 - t406 * t564) * t601; t393; t604; t610; t411; t577;];
tauJ  = t1;
