% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPRRR12
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
% Datum: 2019-05-07 00:47
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPRRR12_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR12_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR12_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR12_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR12_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR12_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR12_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR12_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR12_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 00:38:18
% EndTime: 2019-05-07 00:38:28
% DurationCPUTime: 7.48s
% Computational Cost: add. (80325->346), mult. (179303->432), div. (0->0), fcn. (134435->12), ass. (0->152)
t618 = -2 * qJD(3);
t617 = Ifges(3,1) + Ifges(4,2);
t609 = Ifges(3,4) + Ifges(4,6);
t608 = Ifges(3,5) - Ifges(4,4);
t616 = Ifges(3,2) + Ifges(4,3);
t607 = Ifges(3,6) - Ifges(4,5);
t615 = Ifges(3,3) + Ifges(4,1);
t568 = sin(qJ(2));
t563 = sin(pkin(6));
t598 = qJD(1) * t563;
t556 = t568 * t598;
t564 = cos(pkin(6));
t558 = qJD(1) * t564 + qJD(2);
t614 = (pkin(2) * t558 + t618) * t556;
t575 = qJD(1) ^ 2;
t569 = sin(qJ(1));
t574 = cos(qJ(1));
t587 = -g(1) * t574 - g(2) * t569;
t595 = qJDD(1) * t563;
t539 = -pkin(1) * t575 + pkin(8) * t595 + t587;
t573 = cos(qJ(2));
t604 = t563 * t568;
t591 = t569 * g(1) - g(2) * t574;
t612 = pkin(8) * t563;
t538 = qJDD(1) * pkin(1) + t575 * t612 + t591;
t606 = t538 * t564;
t501 = -g(3) * t604 + t573 * t539 + t568 * t606;
t540 = (-t573 * pkin(2) - t568 * qJ(3)) * t598;
t555 = t558 ^ 2;
t557 = qJDD(1) * t564 + qJDD(2);
t597 = qJD(1) * t573;
t592 = t563 * t597;
t480 = pkin(2) * t555 - t557 * qJ(3) - t540 * t592 + t558 * t618 - t501;
t613 = -pkin(2) - pkin(9);
t611 = g(3) * t564;
t610 = mrSges(3,1) - mrSges(4,2);
t605 = t563 ^ 2 * t575;
t603 = t563 * t573;
t543 = pkin(3) * t556 - pkin(9) * t558;
t544 = (qJD(2) * t597 + qJDD(1) * t568) * t563;
t545 = -qJD(2) * t556 + t573 * t595;
t593 = t573 ^ 2 * t605;
t468 = -pkin(3) * t593 - t611 - qJ(3) * t544 + t613 * t545 + (-t538 + (-qJ(3) * t558 * t573 - t543 * t568) * qJD(1)) * t563 + t614;
t599 = g(3) * t603 + t568 * t539;
t585 = -qJ(3) * t555 + t540 * t556 + qJDD(3) + t599;
t471 = pkin(3) * t544 + t613 * t557 + (-pkin(3) * t558 * t598 - pkin(9) * t568 * t605 - t606) * t573 + t585;
t567 = sin(qJ(4));
t572 = cos(qJ(4));
t447 = -t468 * t567 + t572 * t471;
t525 = -t558 * t567 - t572 * t592;
t499 = qJD(4) * t525 - t545 * t567 + t557 * t572;
t526 = t558 * t572 - t567 * t592;
t533 = qJDD(4) + t544;
t550 = t556 + qJD(4);
t443 = (t525 * t550 - t499) * pkin(10) + (t525 * t526 + t533) * pkin(4) + t447;
t448 = t572 * t468 + t567 * t471;
t498 = -qJD(4) * t526 - t545 * t572 - t557 * t567;
t508 = pkin(4) * t550 - pkin(10) * t526;
t524 = t525 ^ 2;
t445 = -pkin(4) * t524 + pkin(10) * t498 - t508 * t550 + t448;
t566 = sin(qJ(5));
t571 = cos(qJ(5));
t440 = t566 * t443 + t571 * t445;
t504 = t525 * t566 + t526 * t571;
t462 = -qJD(5) * t504 + t498 * t571 - t499 * t566;
t503 = t525 * t571 - t526 * t566;
t482 = -mrSges(6,1) * t503 + mrSges(6,2) * t504;
t548 = qJD(5) + t550;
t489 = mrSges(6,1) * t548 - mrSges(6,3) * t504;
t530 = qJDD(5) + t533;
t483 = -pkin(5) * t503 - pkin(11) * t504;
t547 = t548 ^ 2;
t437 = -pkin(5) * t547 + pkin(11) * t530 + t483 * t503 + t440;
t467 = pkin(3) * t545 - pkin(9) * t593 + t558 * t543 - t480;
t450 = -pkin(4) * t498 - pkin(10) * t524 + t526 * t508 + t467;
t463 = qJD(5) * t503 + t498 * t566 + t499 * t571;
t441 = (-t503 * t548 - t463) * pkin(11) + (t504 * t548 - t462) * pkin(5) + t450;
t565 = sin(qJ(6));
t570 = cos(qJ(6));
t434 = -t437 * t565 + t441 * t570;
t486 = -t504 * t565 + t548 * t570;
t453 = qJD(6) * t486 + t463 * t570 + t530 * t565;
t461 = qJDD(6) - t462;
t487 = t504 * t570 + t548 * t565;
t472 = -mrSges(7,1) * t486 + mrSges(7,2) * t487;
t502 = qJD(6) - t503;
t473 = -mrSges(7,2) * t502 + mrSges(7,3) * t486;
t431 = m(7) * t434 + mrSges(7,1) * t461 - mrSges(7,3) * t453 - t472 * t487 + t473 * t502;
t435 = t437 * t570 + t441 * t565;
t452 = -qJD(6) * t487 - t463 * t565 + t530 * t570;
t474 = mrSges(7,1) * t502 - mrSges(7,3) * t487;
t432 = m(7) * t435 - mrSges(7,2) * t461 + mrSges(7,3) * t452 + t472 * t486 - t474 * t502;
t588 = -t431 * t565 + t570 * t432;
t419 = m(6) * t440 - mrSges(6,2) * t530 + mrSges(6,3) * t462 + t482 * t503 - t489 * t548 + t588;
t439 = t443 * t571 - t445 * t566;
t488 = -mrSges(6,2) * t548 + mrSges(6,3) * t503;
t436 = -pkin(5) * t530 - pkin(11) * t547 + t483 * t504 - t439;
t583 = -m(7) * t436 + t452 * mrSges(7,1) - mrSges(7,2) * t453 + t486 * t473 - t474 * t487;
t427 = m(6) * t439 + mrSges(6,1) * t530 - mrSges(6,3) * t463 - t482 * t504 + t488 * t548 + t583;
t415 = t566 * t419 + t571 * t427;
t421 = t570 * t431 + t565 * t432;
t602 = (t568 * t608 + t573 * t607) * t598 + t615 * t558;
t601 = (t568 * t609 + t573 * t616) * t598 + t607 * t558;
t600 = (t568 * t617 + t573 * t609) * t598 + t608 * t558;
t594 = t573 * t606;
t505 = -mrSges(5,1) * t525 + mrSges(5,2) * t526;
t506 = -mrSges(5,2) * t550 + mrSges(5,3) * t525;
t412 = m(5) * t447 + mrSges(5,1) * t533 - mrSges(5,3) * t499 - t505 * t526 + t506 * t550 + t415;
t507 = mrSges(5,1) * t550 - mrSges(5,3) * t526;
t589 = t571 * t419 - t427 * t566;
t413 = m(5) * t448 - mrSges(5,2) * t533 + mrSges(5,3) * t498 + t505 * t525 - t507 * t550 + t589;
t590 = -t567 * t412 + t572 * t413;
t515 = -t538 * t563 - t611;
t408 = t412 * t572 + t413 * t567;
t481 = -pkin(2) * t545 + (-t558 * t592 - t544) * qJ(3) + t515 + t614;
t536 = -mrSges(4,1) * t592 - mrSges(4,3) * t558;
t586 = -m(4) * t481 + t544 * mrSges(4,3) - t536 * t592 - t590;
t484 = -pkin(2) * t557 + t585 - t594;
t584 = -m(4) * t484 - t544 * mrSges(4,1) - t408;
t581 = m(6) * t450 - t462 * mrSges(6,1) + t463 * mrSges(6,2) - t503 * t488 + t504 * t489 + t421;
t454 = Ifges(7,5) * t487 + Ifges(7,6) * t486 + Ifges(7,3) * t502;
t456 = Ifges(7,1) * t487 + Ifges(7,4) * t486 + Ifges(7,5) * t502;
t424 = -mrSges(7,1) * t436 + mrSges(7,3) * t435 + Ifges(7,4) * t453 + Ifges(7,2) * t452 + Ifges(7,6) * t461 - t454 * t487 + t456 * t502;
t455 = Ifges(7,4) * t487 + Ifges(7,2) * t486 + Ifges(7,6) * t502;
t425 = mrSges(7,2) * t436 - mrSges(7,3) * t434 + Ifges(7,1) * t453 + Ifges(7,4) * t452 + Ifges(7,5) * t461 + t454 * t486 - t455 * t502;
t476 = Ifges(6,4) * t504 + Ifges(6,2) * t503 + Ifges(6,6) * t548;
t477 = Ifges(6,1) * t504 + Ifges(6,4) * t503 + Ifges(6,5) * t548;
t580 = mrSges(6,1) * t439 - mrSges(6,2) * t440 + Ifges(6,5) * t463 + Ifges(6,6) * t462 + Ifges(6,3) * t530 + pkin(5) * t583 + pkin(11) * t588 + t570 * t424 + t565 * t425 + t504 * t476 - t503 * t477;
t579 = mrSges(7,1) * t434 - mrSges(7,2) * t435 + Ifges(7,5) * t453 + Ifges(7,6) * t452 + Ifges(7,3) * t461 + t455 * t487 - t456 * t486;
t578 = -m(5) * t467 + t498 * mrSges(5,1) - t499 * mrSges(5,2) + t525 * t506 - t526 * t507 - t581;
t537 = mrSges(4,1) * t556 + mrSges(4,2) * t558;
t541 = (mrSges(4,2) * t573 - mrSges(4,3) * t568) * t598;
t577 = -m(4) * t480 + t557 * mrSges(4,3) + t558 * t537 + t541 * t592 - t578;
t491 = Ifges(5,4) * t526 + Ifges(5,2) * t525 + Ifges(5,6) * t550;
t492 = Ifges(5,1) * t526 + Ifges(5,4) * t525 + Ifges(5,5) * t550;
t576 = mrSges(5,1) * t447 - mrSges(5,2) * t448 + Ifges(5,5) * t499 + Ifges(5,6) * t498 + Ifges(5,3) * t533 + pkin(4) * t415 + t526 * t491 - t525 * t492 + t580;
t542 = (-t573 * mrSges(3,1) + t568 * mrSges(3,2)) * t598;
t535 = -mrSges(3,2) * t558 + mrSges(3,3) * t592;
t534 = mrSges(3,1) * t558 - mrSges(3,3) * t556;
t500 = t594 - t599;
t490 = Ifges(5,5) * t526 + Ifges(5,6) * t525 + Ifges(5,3) * t550;
t475 = Ifges(6,5) * t504 + Ifges(6,6) * t503 + Ifges(6,3) * t548;
t416 = t542 * t592 + t577 - t557 * mrSges(3,2) - t558 * t534 + m(3) * t501 + (mrSges(3,3) + mrSges(4,1)) * t545;
t410 = -mrSges(6,1) * t450 + mrSges(6,3) * t440 + Ifges(6,4) * t463 + Ifges(6,2) * t462 + Ifges(6,6) * t530 - pkin(5) * t421 - t475 * t504 + t477 * t548 - t579;
t409 = mrSges(6,2) * t450 - mrSges(6,3) * t439 + Ifges(6,1) * t463 + Ifges(6,4) * t462 + Ifges(6,5) * t530 - pkin(11) * t421 - t424 * t565 + t425 * t570 + t475 * t503 - t476 * t548;
t407 = mrSges(4,2) * t557 + t536 * t558 + t541 * t556 - t584;
t406 = t545 * mrSges(4,2) - t537 * t556 - t586;
t405 = m(3) * t500 - mrSges(3,3) * t544 + (t535 - t536) * t558 + t610 * t557 + (-t541 - t542) * t556 + t584;
t404 = mrSges(5,2) * t467 - mrSges(5,3) * t447 + Ifges(5,1) * t499 + Ifges(5,4) * t498 + Ifges(5,5) * t533 - pkin(10) * t415 + t409 * t571 - t410 * t566 + t490 * t525 - t491 * t550;
t403 = -mrSges(5,1) * t467 + mrSges(5,3) * t448 + Ifges(5,4) * t499 + Ifges(5,2) * t498 + Ifges(5,6) * t533 - pkin(4) * t581 + pkin(10) * t589 + t566 * t409 + t571 * t410 - t526 * t490 + t550 * t492;
t402 = mrSges(3,1) * t500 - mrSges(3,2) * t501 + mrSges(4,2) * t484 - mrSges(4,3) * t480 + t572 * t404 - t567 * t403 - pkin(9) * t408 - pkin(2) * t407 + qJ(3) * t577 + t615 * t557 + (qJ(3) * mrSges(4,1) + t607) * t545 + t608 * t544 + (t601 * t568 - t600 * t573) * t598;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t591 - mrSges(2,2) * t587 + (mrSges(4,1) * t484 + mrSges(3,2) * t515 - mrSges(3,3) * t500 - mrSges(4,3) * t481 + pkin(3) * t408 - qJ(3) * t406 + t617 * t544 + t609 * t545 + t608 * t557 - t601 * t558 + t602 * t592 + t576) * t604 + (-mrSges(3,1) * t515 - mrSges(4,1) * t480 + mrSges(4,2) * t481 + mrSges(3,3) * t501 - pkin(2) * t406 - pkin(3) * t578 - pkin(9) * t590 - t572 * t403 - t567 * t404 + t609 * t544 + t616 * t545 - t602 * t556 + t607 * t557 + t600 * t558) * t603 + t564 * t402 + pkin(1) * ((t405 * t573 + t416 * t568) * t564 + (-m(3) * t515 - t544 * mrSges(3,2) + t610 * t545 + (t535 * t573 + (-t534 + t537) * t568) * t598 + t586) * t563) + (-t405 * t568 + t416 * t573) * t612; t402; t407; t576; t580; t579;];
tauJ  = t1;
