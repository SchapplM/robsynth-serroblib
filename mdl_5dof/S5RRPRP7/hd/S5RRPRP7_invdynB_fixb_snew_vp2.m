% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRPRP7_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP7_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP7_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP7_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP7_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP7_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP7_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP7_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP7_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:00:14
% EndTime: 2019-12-31 20:00:19
% DurationCPUTime: 3.91s
% Computational Cost: add. (35035->289), mult. (78918->357), div. (0->0), fcn. (51766->8), ass. (0->112)
t600 = -2 * qJD(3);
t599 = Ifges(5,1) + Ifges(6,1);
t592 = Ifges(5,4) - Ifges(6,5);
t598 = -Ifges(5,5) - Ifges(6,4);
t597 = Ifges(5,2) + Ifges(6,3);
t590 = Ifges(5,6) - Ifges(6,6);
t596 = -Ifges(5,3) - Ifges(6,2);
t562 = sin(qJ(2));
t564 = cos(qJ(2));
t580 = qJD(1) * qJD(2);
t548 = qJDD(1) * t562 + t564 * t580;
t563 = sin(qJ(1));
t565 = cos(qJ(1));
t554 = -g(1) * t565 - g(2) * t563;
t567 = qJD(1) ^ 2;
t543 = -pkin(1) * t567 + qJDD(1) * pkin(6) + t554;
t589 = t562 * t543;
t594 = pkin(2) * t567;
t503 = qJDD(2) * pkin(2) - t548 * qJ(3) - t589 + (qJ(3) * t580 + t562 * t594 - g(3)) * t564;
t529 = -g(3) * t562 + t564 * t543;
t549 = qJDD(1) * t564 - t562 * t580;
t583 = qJD(1) * t562;
t550 = qJD(2) * pkin(2) - qJ(3) * t583;
t558 = t564 ^ 2;
t504 = qJ(3) * t549 - qJD(2) * t550 - t558 * t594 + t529;
t559 = sin(pkin(8));
t560 = cos(pkin(8));
t538 = (t559 * t564 + t560 * t562) * qJD(1);
t484 = t560 * t503 - t559 * t504 + t538 * t600;
t537 = (t559 * t562 - t560 * t564) * qJD(1);
t595 = cos(qJ(4));
t593 = -mrSges(5,3) - mrSges(6,2);
t485 = t559 * t503 + t560 * t504 + t537 * t600;
t518 = mrSges(4,1) * t537 + mrSges(4,2) * t538;
t524 = -t548 * t559 + t549 * t560;
t531 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t538;
t519 = pkin(3) * t537 - pkin(7) * t538;
t566 = qJD(2) ^ 2;
t481 = -pkin(3) * t566 + qJDD(2) * pkin(7) - t519 * t537 + t485;
t553 = t563 * g(1) - t565 * g(2);
t571 = -qJDD(1) * pkin(1) - t553;
t509 = -t549 * pkin(2) + qJDD(3) + t550 * t583 + (-qJ(3) * t558 - pkin(6)) * t567 + t571;
t525 = t548 * t560 + t549 * t559;
t483 = (qJD(2) * t537 - t525) * pkin(7) + (qJD(2) * t538 - t524) * pkin(3) + t509;
t561 = sin(qJ(4));
t478 = t595 * t481 + t561 * t483;
t527 = t561 * qJD(2) + t595 * t538;
t496 = qJD(4) * t527 - t595 * qJDD(2) + t525 * t561;
t536 = qJD(4) + t537;
t512 = mrSges(5,1) * t536 - mrSges(5,3) * t527;
t523 = qJDD(4) - t524;
t526 = -t595 * qJD(2) + t538 * t561;
t505 = pkin(4) * t526 - qJ(5) * t527;
t535 = t536 ^ 2;
t474 = -pkin(4) * t535 + qJ(5) * t523 + 0.2e1 * qJD(5) * t536 - t505 * t526 + t478;
t513 = -mrSges(6,1) * t536 + mrSges(6,2) * t527;
t579 = m(6) * t474 + t523 * mrSges(6,3) + t536 * t513;
t506 = mrSges(6,1) * t526 - mrSges(6,3) * t527;
t584 = -mrSges(5,1) * t526 - mrSges(5,2) * t527 - t506;
t469 = m(5) * t478 - t523 * mrSges(5,2) + t593 * t496 - t536 * t512 + t584 * t526 + t579;
t477 = -t561 * t481 + t595 * t483;
t497 = -t526 * qJD(4) + t561 * qJDD(2) + t595 * t525;
t511 = -mrSges(5,2) * t536 - mrSges(5,3) * t526;
t475 = -t523 * pkin(4) - t535 * qJ(5) + t527 * t505 + qJDD(5) - t477;
t510 = -mrSges(6,2) * t526 + mrSges(6,3) * t536;
t573 = -m(6) * t475 + t523 * mrSges(6,1) + t536 * t510;
t471 = m(5) * t477 + t523 * mrSges(5,1) + t593 * t497 + t536 * t511 + t584 * t527 + t573;
t575 = t595 * t469 - t471 * t561;
t462 = m(4) * t485 - qJDD(2) * mrSges(4,2) + mrSges(4,3) * t524 - qJD(2) * t531 - t518 * t537 + t575;
t530 = -qJD(2) * mrSges(4,2) - mrSges(4,3) * t537;
t480 = -qJDD(2) * pkin(3) - t566 * pkin(7) + t538 * t519 - t484;
t476 = -0.2e1 * qJD(5) * t527 + (t526 * t536 - t497) * qJ(5) + (t527 * t536 + t496) * pkin(4) + t480;
t472 = m(6) * t476 + mrSges(6,1) * t496 - t497 * mrSges(6,3) + t510 * t526 - t527 * t513;
t569 = -m(5) * t480 - t496 * mrSges(5,1) - mrSges(5,2) * t497 - t526 * t511 - t512 * t527 - t472;
t466 = m(4) * t484 + qJDD(2) * mrSges(4,1) - mrSges(4,3) * t525 + qJD(2) * t530 - t518 * t538 + t569;
t456 = t559 * t462 + t560 * t466;
t528 = -t564 * g(3) - t589;
t547 = (-mrSges(3,1) * t564 + mrSges(3,2) * t562) * qJD(1);
t582 = qJD(1) * t564;
t552 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t582;
t454 = m(3) * t528 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t548 + qJD(2) * t552 - t547 * t583 + t456;
t551 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t583;
t576 = t560 * t462 - t466 * t559;
t455 = m(3) * t529 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t549 - qJD(2) * t551 + t547 * t582 + t576;
t577 = -t454 * t562 + t564 * t455;
t447 = m(2) * t554 - mrSges(2,1) * t567 - qJDD(1) * mrSges(2,2) + t577;
t542 = -t567 * pkin(6) + t571;
t464 = t561 * t469 + t595 * t471;
t570 = m(4) * t509 - t524 * mrSges(4,1) + mrSges(4,2) * t525 + t537 * t530 + t531 * t538 + t464;
t568 = -m(3) * t542 + t549 * mrSges(3,1) - mrSges(3,2) * t548 - t551 * t583 + t552 * t582 - t570;
t458 = m(2) * t553 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t567 + t568;
t588 = t563 * t447 + t565 * t458;
t448 = t564 * t454 + t562 * t455;
t587 = t597 * t526 - t592 * t527 - t590 * t536;
t586 = t590 * t526 + t598 * t527 + t596 * t536;
t585 = -t592 * t526 + t599 * t527 - t598 * t536;
t578 = t565 * t447 - t458 * t563;
t541 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t562 + Ifges(3,4) * t564) * qJD(1);
t540 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t562 + Ifges(3,2) * t564) * qJD(1);
t539 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t562 + Ifges(3,6) * t564) * qJD(1);
t516 = Ifges(4,1) * t538 - Ifges(4,4) * t537 + Ifges(4,5) * qJD(2);
t515 = Ifges(4,4) * t538 - Ifges(4,2) * t537 + Ifges(4,6) * qJD(2);
t514 = Ifges(4,5) * t538 - Ifges(4,6) * t537 + Ifges(4,3) * qJD(2);
t463 = mrSges(5,2) * t480 + mrSges(6,2) * t475 - mrSges(5,3) * t477 - mrSges(6,3) * t476 - qJ(5) * t472 - t592 * t496 + t599 * t497 - t523 * t598 + t586 * t526 + t587 * t536;
t459 = -mrSges(5,1) * t480 - mrSges(6,1) * t476 + mrSges(6,2) * t474 + mrSges(5,3) * t478 - pkin(4) * t472 - t597 * t496 + t592 * t497 + t590 * t523 + t586 * t527 + t585 * t536;
t450 = Ifges(4,4) * t525 + Ifges(4,2) * t524 + Ifges(4,6) * qJDD(2) - t538 * t514 + qJD(2) * t516 - mrSges(4,1) * t509 + mrSges(4,3) * t485 - mrSges(5,1) * t477 + mrSges(5,2) * t478 + mrSges(6,1) * t475 - mrSges(6,3) * t474 - pkin(4) * t573 - qJ(5) * t579 - pkin(3) * t464 + (pkin(4) * t506 + t587) * t527 + (qJ(5) * t506 - t585) * t526 + t596 * t523 + (mrSges(6,2) * pkin(4) + t598) * t497 + (mrSges(6,2) * qJ(5) + t590) * t496;
t449 = mrSges(4,2) * t509 - mrSges(4,3) * t484 + Ifges(4,1) * t525 + Ifges(4,4) * t524 + Ifges(4,5) * qJDD(2) - pkin(7) * t464 - qJD(2) * t515 - t561 * t459 + t595 * t463 - t537 * t514;
t444 = mrSges(3,2) * t542 - mrSges(3,3) * t528 + Ifges(3,1) * t548 + Ifges(3,4) * t549 + Ifges(3,5) * qJDD(2) - qJ(3) * t456 - qJD(2) * t540 + t449 * t560 - t450 * t559 + t539 * t582;
t443 = -mrSges(3,1) * t542 + mrSges(3,3) * t529 + Ifges(3,4) * t548 + Ifges(3,2) * t549 + Ifges(3,6) * qJDD(2) - pkin(2) * t570 + qJ(3) * t576 + qJD(2) * t541 + t559 * t449 + t560 * t450 - t539 * t583;
t442 = -pkin(1) * t448 + mrSges(2,3) * t554 - pkin(2) * t456 - Ifges(3,5) * t548 - Ifges(3,6) * t549 - mrSges(3,1) * t528 + mrSges(3,2) * t529 - Ifges(4,5) * t525 - Ifges(4,6) * t524 - mrSges(4,1) * t484 + mrSges(4,2) * t485 - t561 * t463 - t595 * t459 - pkin(3) * t569 - pkin(7) * t575 + mrSges(2,1) * g(3) + t567 * Ifges(2,5) - t538 * t515 - t537 * t516 + Ifges(2,6) * qJDD(1) + (-Ifges(3,3) - Ifges(4,3)) * qJDD(2) + (-t540 * t562 + t541 * t564) * qJD(1);
t441 = -mrSges(2,2) * g(3) - mrSges(2,3) * t553 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t567 - pkin(6) * t448 - t443 * t562 + t444 * t564;
t1 = [-m(1) * g(1) + t578; -m(1) * g(2) + t588; (-m(1) - m(2)) * g(3) + t448; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t588 + t565 * t441 - t563 * t442; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t578 + t563 * t441 + t565 * t442; -mrSges(1,1) * g(2) + mrSges(2,1) * t553 + mrSges(1,2) * g(1) - mrSges(2,2) * t554 + Ifges(2,3) * qJDD(1) + pkin(1) * t568 + pkin(6) * t577 + t564 * t443 + t562 * t444;];
tauB = t1;
