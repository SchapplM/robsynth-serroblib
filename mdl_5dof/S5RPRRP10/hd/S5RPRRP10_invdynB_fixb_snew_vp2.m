% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPRRP10
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPRRP10_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP10_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP10_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP10_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP10_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP10_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP10_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP10_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP10_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:51:11
% EndTime: 2019-12-31 18:51:17
% DurationCPUTime: 3.62s
% Computational Cost: add. (32314->269), mult. (75982->323), div. (0->0), fcn. (53132->8), ass. (0->114)
t598 = Ifges(5,1) + Ifges(6,1);
t592 = Ifges(5,4) + Ifges(6,4);
t591 = Ifges(5,5) + Ifges(6,5);
t597 = Ifges(5,2) + Ifges(6,2);
t596 = Ifges(5,6) + Ifges(6,6);
t595 = Ifges(5,3) + Ifges(6,3);
t559 = qJD(1) ^ 2;
t550 = sin(pkin(8));
t551 = cos(pkin(8));
t553 = sin(qJ(3));
t556 = cos(qJ(3));
t564 = t550 * t553 - t551 * t556;
t535 = t564 * qJD(1);
t565 = t550 * t556 + t551 * t553;
t536 = t565 * qJD(1);
t580 = t536 * qJD(3);
t523 = -t564 * qJDD(1) - t580;
t594 = pkin(2) * t551;
t593 = -mrSges(5,2) - mrSges(6,2);
t589 = mrSges(3,2) * t550;
t549 = t551 ^ 2;
t588 = t549 * t559;
t554 = sin(qJ(1));
t557 = cos(qJ(1));
t541 = -t557 * g(1) - t554 * g(2);
t537 = -t559 * pkin(1) + qJDD(1) * qJ(2) + t541;
t579 = qJD(1) * qJD(2);
t575 = -t551 * g(3) - 0.2e1 * t550 * t579;
t511 = (-pkin(6) * qJDD(1) + t559 * t594 - t537) * t550 + t575;
t526 = -t550 * g(3) + (t537 + 0.2e1 * t579) * t551;
t578 = qJDD(1) * t551;
t512 = -pkin(2) * t588 + pkin(6) * t578 + t526;
t484 = t553 * t511 + t556 * t512;
t518 = t535 * mrSges(4,1) + t536 * mrSges(4,2);
t531 = qJD(3) * mrSges(4,1) - t536 * mrSges(4,3);
t521 = t535 * pkin(3) - t536 * pkin(7);
t558 = qJD(3) ^ 2;
t479 = -t558 * pkin(3) + qJDD(3) * pkin(7) - t535 * t521 + t484;
t548 = t550 ^ 2;
t540 = t554 * g(1) - t557 * g(2);
t569 = qJDD(2) - t540;
t522 = (-pkin(1) - t594) * qJDD(1) + (-qJ(2) + (-t548 - t549) * pkin(6)) * t559 + t569;
t581 = t535 * qJD(3);
t524 = t565 * qJDD(1) - t581;
t482 = (-t524 + t581) * pkin(7) + (-t523 + t580) * pkin(3) + t522;
t552 = sin(qJ(4));
t555 = cos(qJ(4));
t475 = -t552 * t479 + t555 * t482;
t528 = t555 * qJD(3) - t552 * t536;
t498 = t528 * qJD(4) + t552 * qJDD(3) + t555 * t524;
t529 = t552 * qJD(3) + t555 * t536;
t500 = -t528 * mrSges(6,1) + t529 * mrSges(6,2);
t501 = -t528 * mrSges(5,1) + t529 * mrSges(5,2);
t533 = qJD(4) + t535;
t505 = -t533 * mrSges(5,2) + t528 * mrSges(5,3);
t520 = qJDD(4) - t523;
t471 = -0.2e1 * qJD(5) * t529 + (t528 * t533 - t498) * qJ(5) + (t528 * t529 + t520) * pkin(4) + t475;
t504 = -t533 * mrSges(6,2) + t528 * mrSges(6,3);
t577 = m(6) * t471 + t520 * mrSges(6,1) + t533 * t504;
t463 = m(5) * t475 + t520 * mrSges(5,1) + t533 * t505 + (-t500 - t501) * t529 + (-mrSges(5,3) - mrSges(6,3)) * t498 + t577;
t476 = t555 * t479 + t552 * t482;
t497 = -t529 * qJD(4) + t555 * qJDD(3) - t552 * t524;
t506 = t533 * pkin(4) - t529 * qJ(5);
t527 = t528 ^ 2;
t473 = -t527 * pkin(4) + t497 * qJ(5) + 0.2e1 * qJD(5) * t528 - t533 * t506 + t476;
t576 = m(6) * t473 + t497 * mrSges(6,3) + t528 * t500;
t507 = t533 * mrSges(6,1) - t529 * mrSges(6,3);
t583 = -t533 * mrSges(5,1) + t529 * mrSges(5,3) - t507;
t466 = m(5) * t476 + t497 * mrSges(5,3) + t528 * t501 + t593 * t520 + t583 * t533 + t576;
t571 = -t552 * t463 + t555 * t466;
t459 = m(4) * t484 - qJDD(3) * mrSges(4,2) + t523 * mrSges(4,3) - qJD(3) * t531 - t535 * t518 + t571;
t483 = t556 * t511 - t553 * t512;
t530 = -qJD(3) * mrSges(4,2) - t535 * mrSges(4,3);
t478 = -qJDD(3) * pkin(3) - t558 * pkin(7) + t536 * t521 - t483;
t474 = -t497 * pkin(4) - t527 * qJ(5) + t529 * t506 + qJDD(5) + t478;
t570 = m(6) * t474 - t497 * mrSges(6,1) - t528 * t504;
t561 = -m(5) * t478 + t497 * mrSges(5,1) + t593 * t498 + t528 * t505 + t583 * t529 - t570;
t468 = m(4) * t483 + qJDD(3) * mrSges(4,1) - t524 * mrSges(4,3) + qJD(3) * t530 - t536 * t518 + t561;
t453 = t553 * t459 + t556 * t468;
t525 = -t550 * t537 + t575;
t563 = mrSges(3,3) * qJDD(1) + t559 * (-mrSges(3,1) * t551 + t589);
t451 = m(3) * t525 - t563 * t550 + t453;
t572 = t556 * t459 - t553 * t468;
t452 = m(3) * t526 + t563 * t551 + t572;
t573 = -t550 * t451 + t551 * t452;
t445 = m(2) * t541 - t559 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t573;
t534 = -qJDD(1) * pkin(1) - t559 * qJ(2) + t569;
t461 = t555 * t463 + t552 * t466;
t562 = m(4) * t522 - t523 * mrSges(4,1) + t524 * mrSges(4,2) + t535 * t530 + t536 * t531 + t461;
t560 = -m(3) * t534 + mrSges(3,1) * t578 - t562 + (t548 * t559 + t588) * mrSges(3,3);
t456 = t560 + (mrSges(2,1) - t589) * qJDD(1) - t559 * mrSges(2,2) + m(2) * t540;
t587 = t554 * t445 + t557 * t456;
t446 = t551 * t451 + t550 * t452;
t586 = t596 * t528 + t591 * t529 + t595 * t533;
t585 = -t597 * t528 - t592 * t529 - t596 * t533;
t584 = t592 * t528 + t598 * t529 + t591 * t533;
t566 = Ifges(3,5) * t550 + Ifges(3,6) * t551;
t582 = t559 * t566;
t574 = t557 * t445 - t554 * t456;
t568 = Ifges(3,1) * t550 + Ifges(3,4) * t551;
t567 = Ifges(3,4) * t550 + Ifges(3,2) * t551;
t515 = Ifges(4,1) * t536 - Ifges(4,4) * t535 + Ifges(4,5) * qJD(3);
t514 = Ifges(4,4) * t536 - Ifges(4,2) * t535 + Ifges(4,6) * qJD(3);
t513 = Ifges(4,5) * t536 - Ifges(4,6) * t535 + Ifges(4,3) * qJD(3);
t469 = -t498 * mrSges(6,3) - t529 * t500 + t577;
t460 = mrSges(5,2) * t478 + mrSges(6,2) * t474 - mrSges(5,3) * t475 - mrSges(6,3) * t471 - qJ(5) * t469 + t592 * t497 + t598 * t498 + t591 * t520 + t586 * t528 + t585 * t533;
t454 = -mrSges(5,1) * t478 + mrSges(5,3) * t476 - mrSges(6,1) * t474 + mrSges(6,3) * t473 - pkin(4) * t570 + qJ(5) * t576 + (-qJ(5) * t507 + t584) * t533 + (-pkin(4) * t507 - t586) * t529 + (-qJ(5) * mrSges(6,2) + t596) * t520 + (-pkin(4) * mrSges(6,2) + t592) * t498 + t597 * t497;
t447 = -mrSges(4,1) * t522 - mrSges(5,1) * t475 - mrSges(6,1) * t471 + mrSges(5,2) * t476 + mrSges(6,2) * t473 + mrSges(4,3) * t484 + Ifges(4,4) * t524 + Ifges(4,2) * t523 + Ifges(4,6) * qJDD(3) - pkin(3) * t461 - pkin(4) * t469 + qJD(3) * t515 - t536 * t513 + t585 * t529 + t584 * t528 - t595 * t520 - t591 * t498 - t596 * t497;
t442 = mrSges(4,2) * t522 - mrSges(4,3) * t483 + Ifges(4,1) * t524 + Ifges(4,4) * t523 + Ifges(4,5) * qJDD(3) - pkin(7) * t461 - qJD(3) * t514 - t552 * t454 + t555 * t460 - t535 * t513;
t441 = mrSges(3,2) * t534 - mrSges(3,3) * t525 - pkin(6) * t453 + t568 * qJDD(1) + t556 * t442 - t553 * t447 + t551 * t582;
t440 = -mrSges(3,1) * t534 + mrSges(3,3) * t526 - pkin(2) * t562 + pkin(6) * t572 + t567 * qJDD(1) + t553 * t442 + t556 * t447 - t550 * t582;
t439 = -pkin(1) * t446 + mrSges(2,3) * t541 + mrSges(2,1) * g(3) - pkin(2) * t453 - mrSges(3,1) * t525 + mrSges(3,2) * t526 - t552 * t460 - t555 * t454 - pkin(3) * t561 - pkin(7) * t571 - mrSges(4,1) * t483 + mrSges(4,2) * t484 - Ifges(4,5) * t524 - Ifges(4,6) * t523 - Ifges(4,3) * qJDD(3) - t536 * t514 - t535 * t515 + (Ifges(2,6) - t566) * qJDD(1) + (-t550 * t567 + t551 * t568 + Ifges(2,5)) * t559;
t438 = -mrSges(2,2) * g(3) - mrSges(2,3) * t540 + Ifges(2,5) * qJDD(1) - t559 * Ifges(2,6) - qJ(2) * t446 - t550 * t440 + t551 * t441;
t1 = [-m(1) * g(1) + t574; -m(1) * g(2) + t587; (-m(1) - m(2)) * g(3) + t446; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t587 + t557 * t438 - t554 * t439; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t574 + t554 * t438 + t557 * t439; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(2,3) * qJDD(1) + mrSges(2,1) * t540 - mrSges(2,2) * t541 + t550 * t441 + t551 * t440 + pkin(1) * (-qJDD(1) * t589 + t560) + qJ(2) * t573;];
tauB = t1;
