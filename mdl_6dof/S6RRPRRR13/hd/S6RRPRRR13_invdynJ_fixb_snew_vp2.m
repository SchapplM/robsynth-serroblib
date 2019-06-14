% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPRRR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
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
% Datum: 2019-05-07 01:21
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPRRR13_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR13_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR13_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR13_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR13_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR13_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR13_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR13_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR13_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 01:09:39
% EndTime: 2019-05-07 01:09:53
% DurationCPUTime: 7.77s
% Computational Cost: add. (81731->347), mult. (181164->432), div. (0->0), fcn. (134670->12), ass. (0->152)
t628 = -2 * qJD(3);
t627 = Ifges(3,1) + Ifges(4,2);
t619 = Ifges(3,4) + Ifges(4,6);
t618 = Ifges(3,5) - Ifges(4,4);
t626 = Ifges(3,2) + Ifges(4,3);
t617 = Ifges(3,6) - Ifges(4,5);
t625 = Ifges(3,3) + Ifges(4,1);
t573 = cos(pkin(6));
t567 = qJD(1) * t573 + qJD(2);
t577 = sin(qJ(2));
t572 = sin(pkin(6));
t608 = qJD(1) * t572;
t602 = t577 * t608;
t624 = (pkin(2) * t567 + t628) * t602;
t584 = qJD(1) ^ 2;
t578 = sin(qJ(1));
t583 = cos(qJ(1));
t596 = -g(1) * t583 - g(2) * t578;
t605 = qJDD(1) * t572;
t547 = -pkin(1) * t584 + pkin(8) * t605 + t596;
t582 = cos(qJ(2));
t614 = t572 * t577;
t600 = t578 * g(1) - g(2) * t583;
t622 = pkin(8) * t572;
t546 = qJDD(1) * pkin(1) + t584 * t622 + t600;
t616 = t546 * t573;
t509 = -g(3) * t614 + t582 * t547 + t577 * t616;
t548 = (-pkin(2) * t582 - qJ(3) * t577) * t608;
t565 = t567 ^ 2;
t566 = qJDD(1) * t573 + qJDD(2);
t607 = qJD(1) * t582;
t601 = t572 * t607;
t481 = pkin(2) * t565 - t566 * qJ(3) - t548 * t601 + t567 * t628 - t509;
t623 = -pkin(2) - pkin(9);
t621 = g(3) * t573;
t620 = mrSges(3,1) - mrSges(4,2);
t615 = t572 ^ 2 * t584;
t613 = t572 * t582;
t551 = pkin(3) * t602 - pkin(9) * t567;
t552 = (qJD(2) * t607 + qJDD(1) * t577) * t572;
t553 = -qJD(2) * t602 + t582 * t605;
t603 = t582 ^ 2 * t615;
t472 = -pkin(3) * t603 - t621 - qJ(3) * t552 + t623 * t553 + (-t546 + (-qJ(3) * t567 * t582 - t551 * t577) * qJD(1)) * t572 + t624;
t609 = g(3) * t613 + t577 * t547;
t594 = -qJ(3) * t565 + t548 * t602 + qJDD(3) + t609;
t474 = pkin(3) * t552 + t623 * t566 + (-pkin(3) * t567 * t608 - pkin(9) * t577 * t615 - t616) * t582 + t594;
t576 = sin(qJ(4));
t581 = cos(qJ(4));
t461 = t581 * t472 + t576 * t474;
t534 = -t576 * t567 - t581 * t601;
t535 = t567 * t581 - t576 * t601;
t511 = -pkin(4) * t534 - pkin(10) * t535;
t541 = qJDD(4) + t552;
t560 = qJD(4) + t602;
t557 = t560 ^ 2;
t450 = -pkin(4) * t557 + pkin(10) * t541 + t511 * t534 + t461;
t471 = pkin(3) * t553 - pkin(9) * t603 + t567 * t551 - t481;
t506 = -t535 * qJD(4) - t553 * t581 - t576 * t566;
t507 = qJD(4) * t534 - t553 * t576 + t566 * t581;
t456 = (-t534 * t560 - t507) * pkin(10) + (t535 * t560 - t506) * pkin(4) + t471;
t575 = sin(qJ(5));
t580 = cos(qJ(5));
t445 = -t450 * t575 + t580 * t456;
t513 = -t535 * t575 + t560 * t580;
t477 = qJD(5) * t513 + t507 * t580 + t541 * t575;
t504 = qJDD(5) - t506;
t514 = t535 * t580 + t560 * t575;
t533 = qJD(5) - t534;
t443 = (t513 * t533 - t477) * pkin(11) + (t513 * t514 + t504) * pkin(5) + t445;
t446 = t580 * t450 + t575 * t456;
t476 = -qJD(5) * t514 - t507 * t575 + t541 * t580;
t495 = pkin(5) * t533 - pkin(11) * t514;
t512 = t513 ^ 2;
t444 = -pkin(5) * t512 + pkin(11) * t476 - t495 * t533 + t446;
t574 = sin(qJ(6));
t579 = cos(qJ(6));
t441 = t443 * t579 - t444 * t574;
t488 = t513 * t579 - t514 * t574;
t458 = qJD(6) * t488 + t476 * t574 + t477 * t579;
t489 = t513 * t574 + t514 * t579;
t467 = -mrSges(7,1) * t488 + mrSges(7,2) * t489;
t531 = qJD(6) + t533;
t478 = -mrSges(7,2) * t531 + mrSges(7,3) * t488;
t497 = qJDD(6) + t504;
t437 = m(7) * t441 + mrSges(7,1) * t497 - mrSges(7,3) * t458 - t467 * t489 + t478 * t531;
t442 = t443 * t574 + t444 * t579;
t457 = -qJD(6) * t489 + t476 * t579 - t477 * t574;
t479 = mrSges(7,1) * t531 - mrSges(7,3) * t489;
t438 = m(7) * t442 - mrSges(7,2) * t497 + mrSges(7,3) * t457 + t467 * t488 - t479 * t531;
t430 = t579 * t437 + t574 * t438;
t490 = -mrSges(6,1) * t513 + mrSges(6,2) * t514;
t493 = -mrSges(6,2) * t533 + mrSges(6,3) * t513;
t428 = m(6) * t445 + mrSges(6,1) * t504 - mrSges(6,3) * t477 - t490 * t514 + t493 * t533 + t430;
t494 = mrSges(6,1) * t533 - mrSges(6,3) * t514;
t597 = -t437 * t574 + t579 * t438;
t429 = m(6) * t446 - mrSges(6,2) * t504 + mrSges(6,3) * t476 + t490 * t513 - t494 * t533 + t597;
t424 = t580 * t428 + t575 * t429;
t612 = (t577 * t618 + t582 * t617) * t608 + t625 * t567;
t611 = (t577 * t619 + t582 * t626) * t608 + t617 * t567;
t610 = (t577 * t627 + t582 * t619) * t608 + t618 * t567;
t604 = t582 * t616;
t510 = -mrSges(5,1) * t534 + mrSges(5,2) * t535;
t516 = mrSges(5,1) * t560 - mrSges(5,3) * t535;
t598 = -t428 * t575 + t580 * t429;
t422 = m(5) * t461 - mrSges(5,2) * t541 + mrSges(5,3) * t506 + t510 * t534 - t516 * t560 + t598;
t460 = -t576 * t472 + t474 * t581;
t515 = -mrSges(5,2) * t560 + mrSges(5,3) * t534;
t449 = -pkin(4) * t541 - pkin(10) * t557 + t535 * t511 - t460;
t447 = -pkin(5) * t476 - pkin(11) * t512 + t495 * t514 + t449;
t592 = m(7) * t447 - t457 * mrSges(7,1) + mrSges(7,2) * t458 - t488 * t478 + t479 * t489;
t586 = -m(6) * t449 + t476 * mrSges(6,1) - mrSges(6,2) * t477 + t513 * t493 - t494 * t514 - t592;
t433 = m(5) * t460 + mrSges(5,1) * t541 - mrSges(5,3) * t507 - t510 * t535 + t515 * t560 + t586;
t599 = t581 * t422 - t576 * t433;
t523 = -t546 * t572 - t621;
t415 = t422 * t576 + t433 * t581;
t482 = -pkin(2) * t553 + (-t567 * t601 - t552) * qJ(3) + t523 + t624;
t544 = -mrSges(4,1) * t601 - mrSges(4,3) * t567;
t595 = -m(4) * t482 + t552 * mrSges(4,3) - t544 * t601 - t599;
t487 = -pkin(2) * t566 + t594 - t604;
t593 = -m(4) * t487 - t552 * mrSges(4,1) - t415;
t590 = -m(5) * t471 + mrSges(5,1) * t506 - t507 * mrSges(5,2) + t515 * t534 - t535 * t516 - t424;
t464 = Ifges(7,4) * t489 + Ifges(7,2) * t488 + Ifges(7,6) * t531;
t465 = Ifges(7,1) * t489 + Ifges(7,4) * t488 + Ifges(7,5) * t531;
t589 = -mrSges(7,1) * t441 + mrSges(7,2) * t442 - Ifges(7,5) * t458 - Ifges(7,6) * t457 - Ifges(7,3) * t497 - t489 * t464 + t488 * t465;
t463 = Ifges(7,5) * t489 + Ifges(7,6) * t488 + Ifges(7,3) * t531;
t431 = -mrSges(7,1) * t447 + mrSges(7,3) * t442 + Ifges(7,4) * t458 + Ifges(7,2) * t457 + Ifges(7,6) * t497 - t463 * t489 + t465 * t531;
t432 = mrSges(7,2) * t447 - mrSges(7,3) * t441 + Ifges(7,1) * t458 + Ifges(7,4) * t457 + Ifges(7,5) * t497 + t463 * t488 - t464 * t531;
t483 = Ifges(6,5) * t514 + Ifges(6,6) * t513 + Ifges(6,3) * t533;
t485 = Ifges(6,1) * t514 + Ifges(6,4) * t513 + Ifges(6,5) * t533;
t417 = -mrSges(6,1) * t449 + mrSges(6,3) * t446 + Ifges(6,4) * t477 + Ifges(6,2) * t476 + Ifges(6,6) * t504 - pkin(5) * t592 + pkin(11) * t597 + t579 * t431 + t574 * t432 - t514 * t483 + t533 * t485;
t484 = Ifges(6,4) * t514 + Ifges(6,2) * t513 + Ifges(6,6) * t533;
t419 = mrSges(6,2) * t449 - mrSges(6,3) * t445 + Ifges(6,1) * t477 + Ifges(6,4) * t476 + Ifges(6,5) * t504 - pkin(11) * t430 - t431 * t574 + t432 * t579 + t483 * t513 - t484 * t533;
t499 = Ifges(5,4) * t535 + Ifges(5,2) * t534 + Ifges(5,6) * t560;
t500 = Ifges(5,1) * t535 + Ifges(5,4) * t534 + Ifges(5,5) * t560;
t588 = mrSges(5,1) * t460 - mrSges(5,2) * t461 + Ifges(5,5) * t507 + Ifges(5,6) * t506 + Ifges(5,3) * t541 + pkin(4) * t586 + pkin(10) * t598 + t580 * t417 + t575 * t419 + t535 * t499 - t534 * t500;
t545 = mrSges(4,1) * t602 + mrSges(4,2) * t567;
t549 = (mrSges(4,2) * t582 - mrSges(4,3) * t577) * t608;
t587 = -m(4) * t481 + t566 * mrSges(4,3) + t567 * t545 + t549 * t601 - t590;
t585 = mrSges(6,1) * t445 - mrSges(6,2) * t446 + Ifges(6,5) * t477 + Ifges(6,6) * t476 + Ifges(6,3) * t504 + pkin(5) * t430 + t514 * t484 - t513 * t485 - t589;
t550 = (-mrSges(3,1) * t582 + mrSges(3,2) * t577) * t608;
t543 = -mrSges(3,2) * t567 + mrSges(3,3) * t601;
t542 = mrSges(3,1) * t567 - mrSges(3,3) * t602;
t508 = t604 - t609;
t498 = Ifges(5,5) * t535 + Ifges(5,6) * t534 + Ifges(5,3) * t560;
t420 = t587 + (mrSges(3,3) + mrSges(4,1)) * t553 + t550 * t601 + m(3) * t509 - mrSges(3,2) * t566 - t542 * t567;
t414 = mrSges(4,2) * t566 + t544 * t567 + t549 * t602 - t593;
t413 = t553 * mrSges(4,2) - t545 * t602 - t595;
t412 = m(3) * t508 - mrSges(3,3) * t552 + (t543 - t544) * t567 + t620 * t566 + (-t549 - t550) * t602 + t593;
t411 = -mrSges(5,1) * t471 + mrSges(5,3) * t461 + Ifges(5,4) * t507 + Ifges(5,2) * t506 + Ifges(5,6) * t541 - pkin(4) * t424 - t535 * t498 + t560 * t500 - t585;
t410 = mrSges(5,2) * t471 - mrSges(5,3) * t460 + Ifges(5,1) * t507 + Ifges(5,4) * t506 + Ifges(5,5) * t541 - pkin(10) * t424 - t417 * t575 + t419 * t580 + t498 * t534 - t499 * t560;
t409 = mrSges(3,1) * t508 - mrSges(3,2) * t509 + mrSges(4,2) * t487 - mrSges(4,3) * t481 + t581 * t410 - t576 * t411 - pkin(9) * t415 - pkin(2) * t414 + qJ(3) * t587 + t625 * t566 + (mrSges(4,1) * qJ(3) + t617) * t553 + t618 * t552 + (t577 * t611 - t582 * t610) * t608;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t600 - mrSges(2,2) * t596 + (mrSges(4,1) * t487 + mrSges(3,2) * t523 - mrSges(3,3) * t508 - mrSges(4,3) * t482 + pkin(3) * t415 - qJ(3) * t413 + t627 * t552 + t619 * t553 + t618 * t566 - t611 * t567 + t612 * t601 + t588) * t614 + (-mrSges(3,1) * t523 - mrSges(4,1) * t481 + mrSges(4,2) * t482 + mrSges(3,3) * t509 - pkin(2) * t413 - pkin(3) * t590 - pkin(9) * t599 - t576 * t410 - t581 * t411 + t619 * t552 + t626 * t553 + t617 * t566 + t610 * t567 - t612 * t602) * t613 + t573 * t409 + pkin(1) * ((t412 * t582 + t420 * t577) * t573 + (-m(3) * t523 - t552 * mrSges(3,2) + t620 * t553 + (t543 * t582 + (-t542 + t545) * t577) * t608 + t595) * t572) + (-t412 * t577 + t420 * t582) * t622; t409; t414; t588; t585; -t589;];
tauJ  = t1;
