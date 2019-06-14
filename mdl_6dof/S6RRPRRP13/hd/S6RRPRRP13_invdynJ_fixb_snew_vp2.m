% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPRRP13
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
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
% Datum: 2019-05-06 19:15
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPRRP13_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP13_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP13_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP13_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP13_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP13_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP13_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP13_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP13_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 19:07:01
% EndTime: 2019-05-06 19:07:13
% DurationCPUTime: 4.19s
% Computational Cost: add. (35703->324), mult. (79549->390), div. (0->0), fcn. (57456->10), ass. (0->143)
t619 = -2 * qJD(3);
t618 = Ifges(3,1) + Ifges(4,2);
t617 = Ifges(6,1) + Ifges(7,1);
t605 = Ifges(3,4) + Ifges(4,6);
t604 = Ifges(6,4) + Ifges(7,4);
t603 = Ifges(3,5) - Ifges(4,4);
t602 = Ifges(6,5) + Ifges(7,5);
t616 = Ifges(3,2) + Ifges(4,3);
t615 = Ifges(6,2) + Ifges(7,2);
t601 = Ifges(3,6) - Ifges(4,5);
t600 = Ifges(6,6) + Ifges(7,6);
t614 = Ifges(3,3) + Ifges(4,1);
t613 = Ifges(6,3) + Ifges(7,3);
t553 = cos(pkin(6));
t547 = qJD(1) * t553 + qJD(2);
t556 = sin(qJ(2));
t552 = sin(pkin(6));
t586 = qJD(1) * t552;
t578 = t556 * t586;
t612 = (pkin(2) * t547 + t619) * t578;
t562 = qJD(1) ^ 2;
t557 = sin(qJ(1));
t561 = cos(qJ(1));
t572 = -g(1) * t561 - g(2) * t557;
t583 = qJDD(1) * t552;
t530 = -pkin(1) * t562 + pkin(8) * t583 + t572;
t560 = cos(qJ(2));
t596 = t552 * t556;
t576 = t557 * g(1) - g(2) * t561;
t609 = pkin(8) * t552;
t529 = qJDD(1) * pkin(1) + t562 * t609 + t576;
t598 = t529 * t553;
t494 = -g(3) * t596 + t560 * t530 + t556 * t598;
t531 = (-pkin(2) * t560 - qJ(3) * t556) * t586;
t545 = t547 ^ 2;
t546 = qJDD(1) * t553 + qJDD(2);
t585 = qJD(1) * t560;
t577 = t552 * t585;
t461 = pkin(2) * t545 - t546 * qJ(3) - t531 * t577 + t547 * t619 - t494;
t555 = sin(qJ(4));
t559 = cos(qJ(4));
t517 = -t547 * t555 - t559 * t577;
t536 = -qJD(2) * t578 + t560 * t583;
t492 = qJD(4) * t517 - t536 * t555 + t546 * t559;
t518 = t547 * t559 - t555 * t577;
t540 = qJD(4) + t578;
t554 = sin(qJ(5));
t558 = cos(qJ(5));
t498 = -t518 * t554 + t540 * t558;
t535 = (qJD(2) * t585 + qJDD(1) * t556) * t552;
t524 = qJDD(4) + t535;
t458 = qJD(5) * t498 + t492 * t558 + t524 * t554;
t499 = t518 * t558 + t540 * t554;
t473 = -mrSges(7,1) * t498 + mrSges(7,2) * t499;
t534 = pkin(3) * t578 - pkin(9) * t547;
t597 = t552 ^ 2 * t562;
t581 = t560 ^ 2 * t597;
t608 = g(3) * t553;
t610 = -pkin(2) - pkin(9);
t451 = -pkin(3) * t581 - t608 - qJ(3) * t535 + t610 * t536 + (-t529 + (-qJ(3) * t547 * t560 - t534 * t556) * qJD(1)) * t552 + t612;
t595 = t552 * t560;
t587 = g(3) * t595 + t556 * t530;
t570 = -qJ(3) * t545 + t531 * t578 + qJDD(3) + t587;
t453 = pkin(3) * t535 + t610 * t546 + (-pkin(3) * t547 * t586 - pkin(9) * t556 * t597 - t598) * t560 + t570;
t446 = t559 * t451 + t555 * t453;
t496 = -pkin(4) * t517 - pkin(10) * t518;
t538 = t540 ^ 2;
t440 = -pkin(4) * t538 + pkin(10) * t524 + t496 * t517 + t446;
t450 = pkin(3) * t536 - pkin(9) * t581 + t547 * t534 - t461;
t491 = -qJD(4) * t518 - t536 * t559 - t546 * t555;
t443 = (-t517 * t540 - t492) * pkin(10) + (t518 * t540 - t491) * pkin(4) + t450;
t435 = -t440 * t554 + t558 * t443;
t489 = qJDD(5) - t491;
t516 = qJD(5) - t517;
t432 = -0.2e1 * qJD(6) * t499 + (t498 * t516 - t458) * qJ(6) + (t498 * t499 + t489) * pkin(5) + t435;
t477 = -mrSges(7,2) * t516 + mrSges(7,3) * t498;
t580 = m(7) * t432 + t489 * mrSges(7,1) + t516 * t477;
t429 = -mrSges(7,3) * t458 - t473 * t499 + t580;
t436 = t558 * t440 + t554 * t443;
t457 = -qJD(5) * t499 - t492 * t554 + t524 * t558;
t479 = pkin(5) * t516 - qJ(6) * t499;
t497 = t498 ^ 2;
t434 = -pkin(5) * t497 + qJ(6) * t457 + 0.2e1 * qJD(6) * t498 - t479 * t516 + t436;
t592 = t498 * t604 + t499 * t617 + t516 * t602;
t593 = -t498 * t615 - t499 * t604 - t516 * t600;
t611 = mrSges(6,1) * t435 + mrSges(7,1) * t432 - mrSges(6,2) * t436 - mrSges(7,2) * t434 + pkin(5) * t429 + t600 * t457 + t602 * t458 + t489 * t613 - t592 * t498 - t593 * t499;
t607 = mrSges(3,1) - mrSges(4,2);
t606 = -mrSges(6,2) - mrSges(7,2);
t474 = -mrSges(6,1) * t498 + mrSges(6,2) * t499;
t478 = -mrSges(6,2) * t516 + mrSges(6,3) * t498;
t423 = m(6) * t435 + mrSges(6,1) * t489 + t478 * t516 + (-t473 - t474) * t499 + (-mrSges(6,3) - mrSges(7,3)) * t458 + t580;
t579 = m(7) * t434 + t457 * mrSges(7,3) + t498 * t473;
t480 = mrSges(7,1) * t516 - mrSges(7,3) * t499;
t591 = -mrSges(6,1) * t516 + mrSges(6,3) * t499 - t480;
t426 = m(6) * t436 + mrSges(6,3) * t457 + t474 * t498 + t606 * t489 + t591 * t516 + t579;
t421 = t558 * t423 + t554 * t426;
t594 = -t498 * t600 - t499 * t602 - t516 * t613;
t590 = (t556 * t603 + t560 * t601) * t586 + t614 * t547;
t589 = (t556 * t605 + t560 * t616) * t586 + t601 * t547;
t588 = (t556 * t618 + t560 * t605) * t586 + t603 * t547;
t582 = t560 * t598;
t495 = -mrSges(5,1) * t517 + mrSges(5,2) * t518;
t501 = mrSges(5,1) * t540 - mrSges(5,3) * t518;
t574 = -t423 * t554 + t558 * t426;
t417 = m(5) * t446 - mrSges(5,2) * t524 + mrSges(5,3) * t491 + t495 * t517 - t501 * t540 + t574;
t445 = -t555 * t451 + t453 * t559;
t500 = -mrSges(5,2) * t540 + mrSges(5,3) * t517;
t439 = -pkin(4) * t524 - pkin(10) * t538 + t518 * t496 - t445;
t437 = -pkin(5) * t457 - qJ(6) * t497 + t479 * t499 + qJDD(6) + t439;
t573 = -m(7) * t437 + t457 * mrSges(7,1) + t498 * t477;
t563 = -m(6) * t439 + t457 * mrSges(6,1) + t606 * t458 + t498 * t478 + t591 * t499 + t573;
t427 = m(5) * t445 + mrSges(5,1) * t524 - mrSges(5,3) * t492 - t495 * t518 + t500 * t540 + t563;
t575 = t559 * t417 - t555 * t427;
t508 = -t529 * t552 - t608;
t412 = t417 * t555 + t427 * t559;
t462 = -pkin(2) * t536 + (-t547 * t577 - t535) * qJ(3) + t508 + t612;
t527 = -mrSges(4,1) * t577 - mrSges(4,3) * t547;
t571 = -m(4) * t462 + t535 * mrSges(4,3) - t527 * t577 - t575;
t471 = -pkin(2) * t546 + t570 - t582;
t569 = -m(4) * t471 - t535 * mrSges(4,1) - t412;
t567 = -m(5) * t450 + mrSges(5,1) * t491 - t492 * mrSges(5,2) + t500 * t517 - t518 * t501 - t421;
t430 = mrSges(7,2) * t458 + t480 * t499 - t573;
t414 = -mrSges(6,1) * t439 + mrSges(6,3) * t436 - mrSges(7,1) * t437 + mrSges(7,3) * t434 - pkin(5) * t430 + qJ(6) * t579 + (-qJ(6) * t480 + t592) * t516 + t594 * t499 + (-mrSges(7,2) * qJ(6) + t600) * t489 + t604 * t458 + t615 * t457;
t419 = mrSges(6,2) * t439 + mrSges(7,2) * t437 - mrSges(6,3) * t435 - mrSges(7,3) * t432 - qJ(6) * t429 + t604 * t457 + t458 * t617 + t602 * t489 - t594 * t498 + t593 * t516;
t484 = Ifges(5,4) * t518 + Ifges(5,2) * t517 + Ifges(5,6) * t540;
t485 = Ifges(5,1) * t518 + Ifges(5,4) * t517 + Ifges(5,5) * t540;
t565 = mrSges(5,1) * t445 - mrSges(5,2) * t446 + Ifges(5,5) * t492 + Ifges(5,6) * t491 + Ifges(5,3) * t524 + pkin(4) * t563 + pkin(10) * t574 + t558 * t414 + t554 * t419 + t518 * t484 - t517 * t485;
t528 = mrSges(4,1) * t578 + mrSges(4,2) * t547;
t532 = (mrSges(4,2) * t560 - mrSges(4,3) * t556) * t586;
t564 = -m(4) * t461 + t546 * mrSges(4,3) + t547 * t528 + t532 * t577 - t567;
t533 = (-mrSges(3,1) * t560 + mrSges(3,2) * t556) * t586;
t526 = -mrSges(3,2) * t547 + mrSges(3,3) * t577;
t525 = mrSges(3,1) * t547 - mrSges(3,3) * t578;
t493 = t582 - t587;
t483 = Ifges(5,5) * t518 + Ifges(5,6) * t517 + Ifges(5,3) * t540;
t415 = t564 + m(3) * t494 - mrSges(3,2) * t546 - t525 * t547 + (mrSges(3,3) + mrSges(4,1)) * t536 + t533 * t577;
t411 = mrSges(4,2) * t546 + t527 * t547 + t532 * t578 - t569;
t410 = t536 * mrSges(4,2) - t528 * t578 - t571;
t409 = m(3) * t493 - mrSges(3,3) * t535 + (t526 - t527) * t547 + t607 * t546 + (-t532 - t533) * t578 + t569;
t408 = -mrSges(5,1) * t450 + mrSges(5,3) * t446 + Ifges(5,4) * t492 + Ifges(5,2) * t491 + Ifges(5,6) * t524 - pkin(4) * t421 - t518 * t483 + t540 * t485 - t611;
t407 = mrSges(5,2) * t450 - mrSges(5,3) * t445 + Ifges(5,1) * t492 + Ifges(5,4) * t491 + Ifges(5,5) * t524 - pkin(10) * t421 - t414 * t554 + t419 * t558 + t483 * t517 - t484 * t540;
t406 = mrSges(3,1) * t493 - mrSges(3,2) * t494 + mrSges(4,2) * t471 - mrSges(4,3) * t461 + t559 * t407 - t555 * t408 - pkin(9) * t412 - pkin(2) * t411 + qJ(3) * t564 + t614 * t546 + (mrSges(4,1) * qJ(3) + t601) * t536 + t603 * t535 + (t589 * t556 - t588 * t560) * t586;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t576 - mrSges(2,2) * t572 + (mrSges(4,1) * t471 + mrSges(3,2) * t508 - mrSges(3,3) * t493 - mrSges(4,3) * t462 + pkin(3) * t412 - qJ(3) * t410 + t535 * t618 + t605 * t536 + t603 * t546 - t589 * t547 + t590 * t577 + t565) * t596 + (-mrSges(3,1) * t508 - mrSges(4,1) * t461 + mrSges(4,2) * t462 + mrSges(3,3) * t494 - pkin(2) * t410 - pkin(3) * t567 - pkin(9) * t575 - t555 * t407 - t559 * t408 + t605 * t535 + t536 * t616 + t601 * t546 + t588 * t547 - t590 * t578) * t595 + t553 * t406 + pkin(1) * ((t409 * t560 + t415 * t556) * t553 + (-m(3) * t508 - t535 * mrSges(3,2) + t607 * t536 + (t526 * t560 + (-t525 + t528) * t556) * t586 + t571) * t552) + (-t409 * t556 + t415 * t560) * t609; t406; t411; t565; t611; t430;];
tauJ  = t1;
