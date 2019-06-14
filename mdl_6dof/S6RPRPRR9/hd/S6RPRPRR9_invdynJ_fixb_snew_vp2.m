% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRPRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-05-05 19:44
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRPRR9_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR9_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR9_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR9_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_invdynJ_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR9_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR9_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR9_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 19:32:38
% EndTime: 2019-05-05 19:32:59
% DurationCPUTime: 22.06s
% Computational Cost: add. (322472->355), mult. (1075007->491), div. (0->0), fcn. (925075->16), ass. (0->165)
t566 = sin(pkin(7));
t569 = cos(pkin(12));
t571 = cos(pkin(6));
t567 = sin(pkin(6));
t570 = cos(pkin(7));
t606 = t567 * t570;
t553 = (t566 * t571 + t569 * t606) * qJD(1) * pkin(9);
t580 = qJD(1) ^ 2;
t575 = sin(qJ(1));
t579 = cos(qJ(1));
t594 = -g(1) * t579 - g(2) * t575;
t612 = qJ(2) * t567;
t557 = -pkin(1) * t580 + qJDD(1) * t612 + t594;
t565 = sin(pkin(12));
t614 = pkin(9) * t566;
t591 = -t569 * pkin(2) - t565 * t614;
t603 = qJD(1) * t567;
t613 = pkin(9) * qJDD(1);
t587 = qJD(1) * t591 * t603 + t570 * t613;
t600 = qJD(2) * t603;
t607 = t567 * t569;
t599 = t575 * g(1) - g(2) * t579;
t556 = qJDD(1) * pkin(1) + t580 * t612 + t599;
t611 = t556 * t571;
t592 = -g(3) * t607 - 0.2e1 * t565 * t600 + t569 * t611;
t511 = (pkin(2) * qJDD(1) + qJD(1) * t553) * t571 + (-t567 * t587 - t557) * t565 + t592;
t615 = pkin(9) * t565;
t558 = (pkin(2) * t571 - t606 * t615) * qJD(1);
t616 = 0.2e1 * t569;
t601 = t569 * t557 + t565 * t611 + t600 * t616;
t512 = (-qJD(1) * t558 + t566 * t613) * t571 + (-g(3) * t565 + t569 * t587) * t567 + t601;
t598 = -g(3) * t571 + qJDD(2);
t519 = (-t556 + t591 * qJDD(1) + (-t553 * t569 + t558 * t565) * qJD(1)) * t567 + t598;
t574 = sin(qJ(3));
t578 = cos(qJ(3));
t604 = t570 * t578;
t608 = t566 * t578;
t477 = t511 * t604 - t512 * t574 + t519 * t608;
t584 = t571 * t608 + (-t565 * t574 + t569 * t604) * t567;
t541 = t584 * qJD(1);
t605 = t570 * t574;
t609 = t566 * t574;
t585 = t571 * t609 + (t565 * t578 + t569 * t605) * t567;
t534 = qJD(3) * t541 + qJDD(1) * t585;
t542 = t585 * qJD(1);
t588 = -t566 * t607 + t570 * t571;
t551 = qJDD(1) * t588 + qJDD(3);
t554 = qJD(1) * t588 + qJD(3);
t467 = (t541 * t554 - t534) * qJ(4) + (t541 * t542 + t551) * pkin(3) + t477;
t478 = t511 * t605 + t578 * t512 + t519 * t609;
t533 = -qJD(3) * t542 + qJDD(1) * t584;
t538 = pkin(3) * t554 - qJ(4) * t542;
t540 = t541 ^ 2;
t470 = -pkin(3) * t540 + qJ(4) * t533 - t538 * t554 + t478;
t564 = sin(pkin(13));
t568 = cos(pkin(13));
t531 = t541 * t564 + t542 * t568;
t459 = -0.2e1 * qJD(4) * t531 + t467 * t568 - t564 * t470;
t610 = t565 * t567;
t530 = t541 * t568 - t542 * t564;
t460 = 0.2e1 * qJD(4) * t530 + t564 * t467 + t568 * t470;
t502 = -mrSges(5,1) * t530 + mrSges(5,2) * t531;
t506 = t533 * t568 - t534 * t564;
t521 = mrSges(5,1) * t554 - mrSges(5,3) * t531;
t503 = -pkin(4) * t530 - pkin(10) * t531;
t550 = t554 ^ 2;
t458 = -pkin(4) * t550 + pkin(10) * t551 + t503 * t530 + t460;
t490 = -t511 * t566 + t570 * t519;
t476 = -pkin(3) * t533 - qJ(4) * t540 + t542 * t538 + qJDD(4) + t490;
t507 = t533 * t564 + t534 * t568;
t462 = (-t530 * t554 - t507) * pkin(10) + (t531 * t554 - t506) * pkin(4) + t476;
t573 = sin(qJ(5));
t577 = cos(qJ(5));
t454 = t577 * t458 + t573 * t462;
t517 = -t531 * t573 + t554 * t577;
t518 = t531 * t577 + t554 * t573;
t493 = -pkin(5) * t517 - pkin(11) * t518;
t505 = qJDD(5) - t506;
t529 = qJD(5) - t530;
t528 = t529 ^ 2;
t452 = -pkin(5) * t528 + pkin(11) * t505 + t493 * t517 + t454;
t457 = -pkin(4) * t551 - pkin(10) * t550 + t531 * t503 - t459;
t488 = -qJD(5) * t518 - t507 * t573 + t551 * t577;
t489 = qJD(5) * t517 + t507 * t577 + t551 * t573;
t455 = (-t517 * t529 - t489) * pkin(11) + (t518 * t529 - t488) * pkin(5) + t457;
t572 = sin(qJ(6));
t576 = cos(qJ(6));
t448 = -t452 * t572 + t455 * t576;
t494 = -t518 * t572 + t529 * t576;
t465 = qJD(6) * t494 + t489 * t576 + t505 * t572;
t495 = t518 * t576 + t529 * t572;
t479 = -mrSges(7,1) * t494 + mrSges(7,2) * t495;
t516 = qJD(6) - t517;
t480 = -mrSges(7,2) * t516 + mrSges(7,3) * t494;
t487 = qJDD(6) - t488;
t446 = m(7) * t448 + mrSges(7,1) * t487 - mrSges(7,3) * t465 - t479 * t495 + t480 * t516;
t449 = t452 * t576 + t455 * t572;
t464 = -qJD(6) * t495 - t489 * t572 + t505 * t576;
t481 = mrSges(7,1) * t516 - mrSges(7,3) * t495;
t447 = m(7) * t449 - mrSges(7,2) * t487 + mrSges(7,3) * t464 + t479 * t494 - t481 * t516;
t440 = -t446 * t572 + t576 * t447;
t492 = -mrSges(6,1) * t517 + mrSges(6,2) * t518;
t497 = mrSges(6,1) * t529 - mrSges(6,3) * t518;
t438 = m(6) * t454 - mrSges(6,2) * t505 + mrSges(6,3) * t488 + t492 * t517 - t497 * t529 + t440;
t453 = -t458 * t573 + t462 * t577;
t451 = -pkin(5) * t505 - pkin(11) * t528 + t493 * t518 - t453;
t450 = -m(7) * t451 + t464 * mrSges(7,1) - mrSges(7,2) * t465 + t494 * t480 - t481 * t495;
t496 = -mrSges(6,2) * t529 + mrSges(6,3) * t517;
t444 = m(6) * t453 + mrSges(6,1) * t505 - mrSges(6,3) * t489 - t492 * t518 + t496 * t529 + t450;
t596 = t577 * t438 - t444 * t573;
t429 = m(5) * t460 - mrSges(5,2) * t551 + mrSges(5,3) * t506 + t502 * t530 - t521 * t554 + t596;
t520 = -mrSges(5,2) * t554 + mrSges(5,3) * t530;
t439 = t446 * t576 + t447 * t572;
t583 = -m(6) * t457 + t488 * mrSges(6,1) - mrSges(6,2) * t489 + t517 * t496 - t497 * t518 - t439;
t435 = m(5) * t459 + mrSges(5,1) * t551 - mrSges(5,3) * t507 - t502 * t531 + t520 * t554 + t583;
t424 = t564 * t429 + t568 * t435;
t433 = t573 * t438 + t577 * t444;
t532 = -mrSges(4,1) * t541 + mrSges(4,2) * t542;
t537 = -mrSges(4,2) * t554 + mrSges(4,3) * t541;
t422 = m(4) * t477 + mrSges(4,1) * t551 - mrSges(4,3) * t534 - t532 * t542 + t537 * t554 + t424;
t539 = mrSges(4,1) * t554 - mrSges(4,3) * t542;
t597 = t568 * t429 - t435 * t564;
t423 = m(4) * t478 - mrSges(4,2) * t551 + mrSges(4,3) * t533 + t532 * t541 - t539 * t554 + t597;
t432 = m(5) * t476 - mrSges(5,1) * t506 + t507 * mrSges(5,2) - t520 * t530 + t531 * t521 + t433;
t431 = m(4) * t490 - mrSges(4,1) * t533 + mrSges(4,2) * t534 - t537 * t541 + t539 * t542 + t432;
t411 = t422 * t608 + t423 * t609 + t570 * t431;
t415 = -t422 * t574 + t578 * t423;
t412 = t422 * t604 + t423 * t605 - t431 * t566;
t593 = -t569 * mrSges(3,1) + t565 * mrSges(3,2);
t590 = mrSges(3,1) * t571 - mrSges(3,3) * t610;
t589 = -mrSges(3,2) * t571 + mrSges(3,3) * t607;
t471 = Ifges(7,5) * t495 + Ifges(7,6) * t494 + Ifges(7,3) * t516;
t473 = Ifges(7,1) * t495 + Ifges(7,4) * t494 + Ifges(7,5) * t516;
t441 = -mrSges(7,1) * t451 + mrSges(7,3) * t449 + Ifges(7,4) * t465 + Ifges(7,2) * t464 + Ifges(7,6) * t487 - t471 * t495 + t473 * t516;
t472 = Ifges(7,4) * t495 + Ifges(7,2) * t494 + Ifges(7,6) * t516;
t442 = mrSges(7,2) * t451 - mrSges(7,3) * t448 + Ifges(7,1) * t465 + Ifges(7,4) * t464 + Ifges(7,5) * t487 + t471 * t494 - t472 * t516;
t482 = Ifges(6,5) * t518 + Ifges(6,6) * t517 + Ifges(6,3) * t529;
t483 = Ifges(6,4) * t518 + Ifges(6,2) * t517 + Ifges(6,6) * t529;
t425 = mrSges(6,2) * t457 - mrSges(6,3) * t453 + Ifges(6,1) * t489 + Ifges(6,4) * t488 + Ifges(6,5) * t505 - pkin(11) * t439 - t441 * t572 + t442 * t576 + t482 * t517 - t483 * t529;
t484 = Ifges(6,1) * t518 + Ifges(6,4) * t517 + Ifges(6,5) * t529;
t582 = mrSges(7,1) * t448 - mrSges(7,2) * t449 + Ifges(7,5) * t465 + Ifges(7,6) * t464 + Ifges(7,3) * t487 + t472 * t495 - t473 * t494;
t426 = -mrSges(6,1) * t457 + mrSges(6,3) * t454 + Ifges(6,4) * t489 + Ifges(6,2) * t488 + Ifges(6,6) * t505 - pkin(5) * t439 - t482 * t518 + t484 * t529 - t582;
t498 = Ifges(5,5) * t531 + Ifges(5,6) * t530 + Ifges(5,3) * t554;
t499 = Ifges(5,4) * t531 + Ifges(5,2) * t530 + Ifges(5,6) * t554;
t413 = mrSges(5,2) * t476 - mrSges(5,3) * t459 + Ifges(5,1) * t507 + Ifges(5,4) * t506 + Ifges(5,5) * t551 - pkin(10) * t433 + t425 * t577 - t426 * t573 + t498 * t530 - t499 * t554;
t500 = Ifges(5,1) * t531 + Ifges(5,4) * t530 + Ifges(5,5) * t554;
t581 = mrSges(6,1) * t453 - mrSges(6,2) * t454 + Ifges(6,5) * t489 + Ifges(6,6) * t488 + Ifges(6,3) * t505 + pkin(5) * t450 + pkin(11) * t440 + t576 * t441 + t572 * t442 + t518 * t483 - t517 * t484;
t416 = -mrSges(5,1) * t476 + mrSges(5,3) * t460 + Ifges(5,4) * t507 + Ifges(5,2) * t506 + Ifges(5,6) * t551 - pkin(4) * t433 - t531 * t498 + t554 * t500 - t581;
t522 = Ifges(4,5) * t542 + Ifges(4,6) * t541 + Ifges(4,3) * t554;
t524 = Ifges(4,1) * t542 + Ifges(4,4) * t541 + Ifges(4,5) * t554;
t406 = -mrSges(4,1) * t490 + mrSges(4,3) * t478 + Ifges(4,4) * t534 + Ifges(4,2) * t533 + Ifges(4,6) * t551 - pkin(3) * t432 + qJ(4) * t597 + t564 * t413 + t568 * t416 - t542 * t522 + t554 * t524;
t523 = Ifges(4,4) * t542 + Ifges(4,2) * t541 + Ifges(4,6) * t554;
t407 = mrSges(4,2) * t490 - mrSges(4,3) * t477 + Ifges(4,1) * t534 + Ifges(4,4) * t533 + Ifges(4,5) * t551 - qJ(4) * t424 + t413 * t568 - t416 * t564 + t522 * t541 - t523 * t554;
t586 = pkin(9) * t415 + t406 * t578 + t407 * t574;
t560 = t589 * qJD(1);
t559 = t590 * qJD(1);
t555 = t593 * t603;
t543 = -t556 * t567 + t598;
t536 = -g(3) * t610 + t601;
t535 = -t557 * t565 + t592;
t414 = m(3) * t536 + t589 * qJDD(1) + (t555 * t607 - t559 * t571) * qJD(1) + t415;
t410 = m(3) * t543 + (t593 * qJDD(1) + (t559 * t565 - t560 * t569) * qJD(1)) * t567 + t411;
t409 = m(3) * t535 + t590 * qJDD(1) + (-t555 * t610 + t560 * t571) * qJD(1) + t412;
t408 = Ifges(4,5) * t534 + Ifges(4,6) * t533 + t542 * t523 - t541 * t524 + mrSges(4,1) * t477 - mrSges(4,2) * t478 + Ifges(5,5) * t507 + Ifges(5,6) * t506 + t531 * t499 - t530 * t500 + mrSges(5,1) * t459 - mrSges(5,2) * t460 + t573 * t425 + t577 * t426 + pkin(4) * t583 + pkin(10) * t596 + pkin(3) * t424 + (Ifges(4,3) + Ifges(5,3)) * t551;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t599 - mrSges(2,2) * t594 + (mrSges(3,1) * t535 - mrSges(3,2) * t536 + pkin(2) * t412 + t570 * t408 + pkin(1) * (t409 * t569 + t414 * t565) + Ifges(3,3) * t571 * qJDD(1) + t586 * t566) * t571 + (t565 * (mrSges(3,2) * t543 - mrSges(3,3) * t535 - t574 * t406 + t578 * t407 - t411 * t614) + t569 * (-mrSges(3,1) * t543 + mrSges(3,3) * t536 - pkin(2) * t411 - t566 * t408) - pkin(1) * t410 + qJ(2) * (-t409 * t565 + t414 * t569) + (-t412 * t615 + t569 * t586) * t570 + ((Ifges(3,2) * t569 ^ 2 + (Ifges(3,1) * t565 + Ifges(3,4) * t616) * t565) * t567 + 0.2e1 * t571 * (Ifges(3,5) * t565 + Ifges(3,6) * t569)) * qJDD(1)) * t567; t410; t408; t432; t581; t582;];
tauJ  = t1;
