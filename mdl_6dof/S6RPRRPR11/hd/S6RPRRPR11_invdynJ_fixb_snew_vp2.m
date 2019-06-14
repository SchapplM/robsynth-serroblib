% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRRPR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-05-06 00:22
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRRPR11_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR11_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR11_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR11_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR11_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR11_invdynJ_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR11_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR11_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR11_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 00:09:03
% EndTime: 2019-05-06 00:09:29
% DurationCPUTime: 26.00s
% Computational Cost: add. (386997->355), mult. (1211269->492), div. (0->0), fcn. (1040689->16), ass. (0->165)
t565 = sin(pkin(12));
t567 = sin(pkin(6));
t569 = cos(pkin(12));
t571 = cos(pkin(6));
t574 = sin(qJ(3));
t570 = cos(pkin(7));
t577 = cos(qJ(3));
t608 = t570 * t577;
t566 = sin(pkin(7));
t612 = t566 * t577;
t583 = t567 * (-t565 * t574 + t569 * t608) + t571 * t612;
t539 = t583 * qJD(1);
t609 = t570 * t574;
t613 = t566 * t574;
t585 = t571 * t613 + (t565 * t577 + t569 * t609) * t567;
t540 = t585 * qJD(1);
t528 = -t540 * qJD(3) + qJDD(1) * t583;
t610 = t567 * t570;
t553 = (t566 * t571 + t569 * t610) * qJD(1) * pkin(9);
t579 = qJD(1) ^ 2;
t575 = sin(qJ(1));
t578 = cos(qJ(1));
t598 = -g(1) * t578 - g(2) * t575;
t616 = qJ(2) * t567;
t557 = -pkin(1) * t579 + qJDD(1) * t616 + t598;
t618 = pkin(9) * t566;
t593 = -pkin(2) * t569 - t565 * t618;
t607 = qJD(1) * t567;
t617 = pkin(9) * qJDD(1);
t588 = qJD(1) * t593 * t607 + t570 * t617;
t603 = qJD(2) * t607;
t611 = t567 * t569;
t602 = t575 * g(1) - g(2) * t578;
t556 = qJDD(1) * pkin(1) + t579 * t616 + t602;
t615 = t556 * t571;
t594 = -g(3) * t611 - 0.2e1 * t565 * t603 + t569 * t615;
t507 = (pkin(2) * qJDD(1) + qJD(1) * t553) * t571 + (-t567 * t588 - t557) * t565 + t594;
t619 = pkin(9) * t565;
t558 = (pkin(2) * t571 - t610 * t619) * qJD(1);
t621 = 0.2e1 * t569;
t604 = t569 * t557 + t565 * t615 + t603 * t621;
t508 = (-qJD(1) * t558 + t566 * t617) * t571 + (-g(3) * t565 + t569 * t588) * t567 + t604;
t601 = -t571 * g(3) + qJDD(2);
t519 = (-t556 + t593 * qJDD(1) + (-t553 * t569 + t558 * t565) * qJD(1)) * t567 + t601;
t478 = -t574 * t508 + (t507 * t570 + t519 * t566) * t577;
t620 = cos(qJ(4));
t614 = t565 * t567;
t479 = t507 * t609 + t577 * t508 + t519 * t613;
t527 = -pkin(3) * t539 - pkin(10) * t540;
t589 = -t566 * t611 + t570 * t571;
t554 = qJD(1) * t589 + qJD(3);
t550 = t554 ^ 2;
t551 = qJDD(1) * t589 + qJDD(3);
t470 = -pkin(3) * t550 + pkin(10) * t551 + t527 * t539 + t479;
t487 = -t566 * t507 + t570 * t519;
t529 = t539 * qJD(3) + qJDD(1) * t585;
t476 = (-t539 * t554 - t529) * pkin(10) + (t540 * t554 - t528) * pkin(3) + t487;
t573 = sin(qJ(4));
t463 = t620 * t470 + t573 * t476;
t533 = t540 * t573 - t620 * t554;
t534 = t540 * t620 + t573 * t554;
t509 = pkin(4) * t533 - qJ(5) * t534;
t525 = qJDD(4) - t528;
t538 = qJD(4) - t539;
t537 = t538 ^ 2;
t458 = -pkin(4) * t537 + qJ(5) * t525 - t509 * t533 + t463;
t469 = -t551 * pkin(3) - t550 * pkin(10) + t540 * t527 - t478;
t502 = qJD(4) * t534 + t529 * t573 - t620 * t551;
t503 = -t533 * qJD(4) + t529 * t620 + t573 * t551;
t461 = (t533 * t538 - t503) * qJ(5) + (t534 * t538 + t502) * pkin(4) + t469;
t564 = sin(pkin(13));
t568 = cos(pkin(13));
t518 = t534 * t568 + t538 * t564;
t453 = -0.2e1 * qJD(5) * t518 - t564 * t458 + t568 * t461;
t489 = t503 * t568 + t525 * t564;
t517 = -t534 * t564 + t538 * t568;
t451 = (t517 * t533 - t489) * pkin(11) + (t517 * t518 + t502) * pkin(5) + t453;
t454 = 0.2e1 * qJD(5) * t517 + t568 * t458 + t564 * t461;
t488 = -t503 * t564 + t525 * t568;
t495 = pkin(5) * t533 - pkin(11) * t518;
t516 = t517 ^ 2;
t452 = -pkin(5) * t516 + pkin(11) * t488 - t495 * t533 + t454;
t572 = sin(qJ(6));
t576 = cos(qJ(6));
t449 = t451 * t576 - t452 * t572;
t490 = t517 * t576 - t518 * t572;
t466 = qJD(6) * t490 + t488 * t572 + t489 * t576;
t491 = t517 * t572 + t518 * t576;
t477 = -mrSges(7,1) * t490 + mrSges(7,2) * t491;
t530 = qJD(6) + t533;
t480 = -mrSges(7,2) * t530 + mrSges(7,3) * t490;
t500 = qJDD(6) + t502;
t446 = m(7) * t449 + mrSges(7,1) * t500 - mrSges(7,3) * t466 - t477 * t491 + t480 * t530;
t450 = t451 * t572 + t452 * t576;
t465 = -qJD(6) * t491 + t488 * t576 - t489 * t572;
t481 = mrSges(7,1) * t530 - mrSges(7,3) * t491;
t447 = m(7) * t450 - mrSges(7,2) * t500 + mrSges(7,3) * t465 + t477 * t490 - t481 * t530;
t438 = t576 * t446 + t572 * t447;
t492 = -mrSges(6,1) * t517 + mrSges(6,2) * t518;
t596 = -mrSges(6,2) * t533 + mrSges(6,3) * t517;
t436 = m(6) * t453 + t502 * mrSges(6,1) - t489 * mrSges(6,3) - t518 * t492 + t533 * t596 + t438;
t494 = mrSges(6,1) * t533 - mrSges(6,3) * t518;
t599 = -t446 * t572 + t576 * t447;
t437 = m(6) * t454 - mrSges(6,2) * t502 + mrSges(6,3) * t488 + t492 * t517 - t494 * t533 + t599;
t434 = -t436 * t564 + t568 * t437;
t510 = mrSges(5,1) * t533 + mrSges(5,2) * t534;
t521 = mrSges(5,1) * t538 - mrSges(5,3) * t534;
t432 = m(5) * t463 - mrSges(5,2) * t525 - mrSges(5,3) * t502 - t510 * t533 - t521 * t538 + t434;
t462 = -t573 * t470 + t476 * t620;
t457 = -t525 * pkin(4) - t537 * qJ(5) + t534 * t509 + qJDD(5) - t462;
t455 = -t488 * pkin(5) - t516 * pkin(11) + t518 * t495 + t457;
t586 = m(7) * t455 - t465 * mrSges(7,1) + mrSges(7,2) * t466 - t490 * t480 + t481 * t491;
t448 = m(6) * t457 - t488 * mrSges(6,1) + mrSges(6,2) * t489 + t494 * t518 - t517 * t596 + t586;
t520 = -mrSges(5,2) * t538 - mrSges(5,3) * t533;
t442 = m(5) * t462 + mrSges(5,1) * t525 - mrSges(5,3) * t503 - t510 * t534 + t520 * t538 - t448;
t424 = t573 * t432 + t620 * t442;
t526 = -mrSges(4,1) * t539 + mrSges(4,2) * t540;
t536 = mrSges(4,1) * t554 - mrSges(4,3) * t540;
t600 = t620 * t432 - t442 * t573;
t421 = m(4) * t479 - mrSges(4,2) * t551 + mrSges(4,3) * t528 + t526 * t539 - t536 * t554 + t600;
t535 = -mrSges(4,2) * t554 + mrSges(4,3) * t539;
t423 = m(4) * t487 - mrSges(4,1) * t528 + mrSges(4,2) * t529 - t535 * t539 + t536 * t540 + t424;
t433 = t436 * t568 + t437 * t564;
t582 = -m(5) * t469 - t502 * mrSges(5,1) - mrSges(5,2) * t503 - t533 * t520 - t534 * t521 - t433;
t429 = m(4) * t478 + mrSges(4,1) * t551 - mrSges(4,3) * t529 - t526 * t540 + t535 * t554 + t582;
t412 = t421 * t613 + t570 * t423 + t429 * t612;
t416 = t577 * t421 - t574 * t429;
t413 = t421 * t609 - t423 * t566 + t429 * t608;
t597 = -mrSges(3,1) * t569 + mrSges(3,2) * t565;
t592 = mrSges(3,1) * t571 - mrSges(3,3) * t614;
t591 = -mrSges(3,2) * t571 + mrSges(3,3) * t611;
t472 = Ifges(7,5) * t491 + Ifges(7,6) * t490 + Ifges(7,3) * t530;
t474 = Ifges(7,1) * t491 + Ifges(7,4) * t490 + Ifges(7,5) * t530;
t439 = -mrSges(7,1) * t455 + mrSges(7,3) * t450 + Ifges(7,4) * t466 + Ifges(7,2) * t465 + Ifges(7,6) * t500 - t472 * t491 + t474 * t530;
t473 = Ifges(7,4) * t491 + Ifges(7,2) * t490 + Ifges(7,6) * t530;
t440 = mrSges(7,2) * t455 - mrSges(7,3) * t449 + Ifges(7,1) * t466 + Ifges(7,4) * t465 + Ifges(7,5) * t500 + t472 * t490 - t473 * t530;
t482 = Ifges(6,5) * t518 + Ifges(6,6) * t517 + Ifges(6,3) * t533;
t484 = Ifges(6,1) * t518 + Ifges(6,4) * t517 + Ifges(6,5) * t533;
t425 = -mrSges(6,1) * t457 + mrSges(6,3) * t454 + Ifges(6,4) * t489 + Ifges(6,2) * t488 + Ifges(6,6) * t502 - pkin(5) * t586 + pkin(11) * t599 + t576 * t439 + t572 * t440 - t518 * t482 + t533 * t484;
t483 = Ifges(6,4) * t518 + Ifges(6,2) * t517 + Ifges(6,6) * t533;
t426 = mrSges(6,2) * t457 - mrSges(6,3) * t453 + Ifges(6,1) * t489 + Ifges(6,4) * t488 + Ifges(6,5) * t502 - pkin(11) * t438 - t439 * t572 + t440 * t576 + t482 * t517 - t483 * t533;
t496 = Ifges(5,5) * t534 - Ifges(5,6) * t533 + Ifges(5,3) * t538;
t497 = Ifges(5,4) * t534 - Ifges(5,2) * t533 + Ifges(5,6) * t538;
t414 = mrSges(5,2) * t469 - mrSges(5,3) * t462 + Ifges(5,1) * t503 - Ifges(5,4) * t502 + Ifges(5,5) * t525 - qJ(5) * t433 - t425 * t564 + t426 * t568 - t496 * t533 - t497 * t538;
t498 = Ifges(5,1) * t534 - Ifges(5,4) * t533 + Ifges(5,5) * t538;
t581 = mrSges(7,1) * t449 - mrSges(7,2) * t450 + Ifges(7,5) * t466 + Ifges(7,6) * t465 + Ifges(7,3) * t500 + t491 * t473 - t490 * t474;
t417 = mrSges(5,3) * t463 - t581 - pkin(5) * t438 - mrSges(5,1) * t469 + Ifges(5,6) * t525 + t538 * t498 + Ifges(5,4) * t503 + (-Ifges(5,2) - Ifges(6,3)) * t502 - mrSges(6,1) * t453 - pkin(4) * t433 - t534 * t496 - Ifges(6,6) * t488 - Ifges(6,5) * t489 + mrSges(6,2) * t454 + t517 * t484 - t518 * t483;
t522 = Ifges(4,5) * t540 + Ifges(4,6) * t539 + Ifges(4,3) * t554;
t523 = Ifges(4,4) * t540 + Ifges(4,2) * t539 + Ifges(4,6) * t554;
t408 = mrSges(4,2) * t487 - mrSges(4,3) * t478 + Ifges(4,1) * t529 + Ifges(4,4) * t528 + Ifges(4,5) * t551 - pkin(10) * t424 + t414 * t620 - t573 * t417 + t539 * t522 - t554 * t523;
t524 = Ifges(4,1) * t540 + Ifges(4,4) * t539 + Ifges(4,5) * t554;
t580 = mrSges(5,1) * t462 - mrSges(5,2) * t463 + Ifges(5,5) * t503 - Ifges(5,6) * t502 + Ifges(5,3) * t525 - pkin(4) * t448 + qJ(5) * t434 + t568 * t425 + t564 * t426 + t534 * t497 + t533 * t498;
t409 = -mrSges(4,1) * t487 + mrSges(4,3) * t479 + Ifges(4,4) * t529 + Ifges(4,2) * t528 + Ifges(4,6) * t551 - pkin(3) * t424 - t540 * t522 + t554 * t524 - t580;
t587 = pkin(9) * t416 + t408 * t574 + t409 * t577;
t560 = t591 * qJD(1);
t559 = t592 * qJD(1);
t555 = t597 * t607;
t541 = -t567 * t556 + t601;
t532 = -g(3) * t614 + t604;
t531 = -t557 * t565 + t594;
t415 = m(3) * t532 + t591 * qJDD(1) + (t555 * t611 - t559 * t571) * qJD(1) + t416;
t411 = m(3) * t541 + (t597 * qJDD(1) + (t559 * t565 - t560 * t569) * qJD(1)) * t567 + t412;
t410 = m(3) * t531 + t592 * qJDD(1) + (-t555 * t614 + t560 * t571) * qJD(1) + t413;
t407 = mrSges(4,1) * t478 - mrSges(4,2) * t479 + Ifges(4,5) * t529 + Ifges(4,6) * t528 + Ifges(4,3) * t551 + pkin(3) * t582 + pkin(10) * t600 + t573 * t414 + t417 * t620 + t540 * t523 - t539 * t524;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t602 - mrSges(2,2) * t598 + (mrSges(3,1) * t531 - mrSges(3,2) * t532 + pkin(2) * t413 + t570 * t407 + pkin(1) * (t410 * t569 + t415 * t565) + Ifges(3,3) * t571 * qJDD(1) + t587 * t566) * t571 + (t565 * (mrSges(3,2) * t541 - mrSges(3,3) * t531 + t577 * t408 - t574 * t409 - t412 * t618) + t569 * (-mrSges(3,1) * t541 + mrSges(3,3) * t532 - pkin(2) * t412 - t566 * t407) - pkin(1) * t411 + qJ(2) * (-t410 * t565 + t415 * t569) + (-t413 * t619 + t569 * t587) * t570 + ((Ifges(3,2) * t569 ^ 2 + (Ifges(3,1) * t565 + Ifges(3,4) * t621) * t565) * t567 + 0.2e1 * t571 * (Ifges(3,5) * t565 + Ifges(3,6) * t569)) * qJDD(1)) * t567; t411; t407; t580; t448; t581;];
tauJ  = t1;
