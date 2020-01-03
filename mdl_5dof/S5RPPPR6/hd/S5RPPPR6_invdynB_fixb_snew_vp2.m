% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPPPR6
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
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
% Datum: 2019-12-31 17:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPPPR6_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR6_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR6_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR6_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR6_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR6_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR6_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR6_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR6_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:47:45
% EndTime: 2019-12-31 17:47:47
% DurationCPUTime: 2.00s
% Computational Cost: add. (14039->253), mult. (35981->328), div. (0->0), fcn. (20484->8), ass. (0->122)
t613 = -2 * qJD(4);
t555 = sin(qJ(1));
t557 = cos(qJ(1));
t534 = -t557 * g(1) - t555 * g(2);
t558 = qJD(1) ^ 2;
t612 = -t558 * pkin(1) + qJDD(1) * qJ(2) + (2 * qJD(1) * qJD(2)) + t534;
t611 = Ifges(4,4) - Ifges(3,5);
t610 = Ifges(4,5) - Ifges(3,6);
t551 = sin(pkin(7));
t589 = t551 * qJD(1);
t535 = -0.2e1 * qJD(3) * t589;
t548 = t551 ^ 2;
t553 = cos(pkin(7));
t549 = t553 ^ 2;
t533 = t555 * g(1) - t557 * g(2);
t575 = qJDD(2) - t533;
t601 = qJ(3) * t551;
t491 = t535 + (-qJ(2) + (-t548 - t549) * pkin(3)) * t558 + (-t601 - pkin(1) + (-pkin(2) - qJ(4)) * t553) * qJDD(1) + t575;
t588 = t553 * qJD(1);
t609 = t588 * t613 + t491;
t506 = -t551 * g(3) + t612 * t553;
t572 = -pkin(2) * t553 - t601;
t524 = t572 * qJD(1);
t497 = -t524 * t588 - t506;
t598 = t549 * t558;
t599 = t548 * t558;
t608 = t598 + t599;
t505 = -t553 * g(3) - t612 * t551;
t607 = Ifges(3,1) + Ifges(4,2);
t606 = Ifges(3,4) + Ifges(4,6);
t605 = -Ifges(3,2) - Ifges(4,3);
t604 = mrSges(3,2) * t551;
t603 = mrSges(4,2) * t553;
t550 = sin(pkin(8));
t602 = mrSges(5,2) * t550;
t597 = t550 * t553;
t570 = mrSges(5,1) * t551 + mrSges(5,3) * t597;
t521 = t570 * qJD(1);
t600 = t521 * t550;
t596 = t551 * t558;
t496 = t524 * t589 + qJDD(3) - t505;
t487 = (-qJ(4) * t553 * t558 + pkin(3) * qJDD(1)) * t551 + t496;
t552 = cos(pkin(8));
t595 = t552 * t487;
t594 = t552 * t553;
t525 = (-mrSges(4,3) * t551 + t603) * qJD(1);
t526 = (-mrSges(3,1) * t553 + t604) * qJD(1);
t484 = t550 * t487 + t609 * t552;
t515 = (mrSges(5,1) * t552 - t602) * t588;
t569 = -mrSges(5,2) * t551 - mrSges(5,3) * t594;
t516 = (pkin(4) * t552 + pkin(6) * t550) * t588;
t584 = t552 * t588;
t586 = qJDD(1) * t551;
t482 = -pkin(4) * t599 + pkin(6) * t586 - t516 * t584 + t484;
t585 = qJDD(1) * t553;
t489 = pkin(3) * t585 - qJ(4) * t598 + qJDD(4) - t497;
t485 = ((qJDD(1) * t550 + t552 * t596) * pkin(6) + (qJDD(1) * t552 - t550 * t596) * pkin(4)) * t553 + t489;
t554 = sin(qJ(5));
t556 = cos(qJ(5));
t479 = -t554 * t482 + t556 * t485;
t567 = t551 * t556 + t554 * t597;
t513 = t567 * qJD(1);
t568 = t551 * t554 - t556 * t597;
t514 = t568 * qJD(1);
t498 = -t513 * mrSges(6,1) + t514 * mrSges(6,2);
t501 = t513 * qJD(5) + t568 * qJDD(1);
t530 = qJD(5) + t584;
t503 = -t530 * mrSges(6,2) + t513 * mrSges(6,3);
t580 = t552 * t585;
t529 = qJDD(5) + t580;
t477 = m(6) * t479 + t529 * mrSges(6,1) - t501 * mrSges(6,3) - t514 * t498 + t530 * t503;
t480 = t556 * t482 + t554 * t485;
t500 = -t514 * qJD(5) + t567 * qJDD(1);
t504 = t530 * mrSges(6,1) - t514 * mrSges(6,3);
t478 = m(6) * t480 - t529 * mrSges(6,2) + t500 * mrSges(6,3) + t513 * t498 - t530 * t504;
t576 = -t554 * t477 + t556 * t478;
t468 = m(5) * t484 + t569 * qJDD(1) + (-t515 * t594 - t521 * t551) * qJD(1) + t576;
t483 = -t609 * t550 + t595;
t522 = t569 * qJD(1);
t481 = -pkin(4) * t586 - pkin(6) * t599 - t595 + (t491 + (t613 - t516) * t588) * t550;
t563 = -m(6) * t481 + t500 * mrSges(6,1) - t501 * mrSges(6,2) + t513 * t503 - t514 * t504;
t473 = m(5) * t483 + t570 * qJDD(1) + (t515 * t597 + t522 * t551) * qJD(1) + t563;
t463 = t550 * t468 + t552 * t473;
t564 = -m(4) * t496 - t463;
t459 = m(3) * t505 + ((-mrSges(4,1) - mrSges(3,3)) * qJDD(1) + (-t525 - t526) * qJD(1)) * t551 + t564;
t469 = t556 * t477 + t554 * t478;
t571 = -m(5) * t489 - mrSges(5,1) * t580 - t522 * t584 - t469;
t562 = -m(4) * t497 + mrSges(4,1) * t585 + t525 * t588 - t571;
t466 = m(3) * t506 + ((mrSges(3,3) - t602) * qJDD(1) + (t526 - t600) * qJD(1)) * t553 + t562;
t577 = -t551 * t459 + t553 * t466;
t454 = m(2) * t534 - t558 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t577;
t566 = -t558 * qJ(2) + t575;
t520 = -qJDD(1) * pkin(1) + t566;
t502 = t535 + (-pkin(1) + t572) * qJDD(1) + t566;
t592 = t552 * t468 - t550 * t473;
t565 = m(4) * t502 - t608 * mrSges(4,1) - mrSges(4,3) * t586 + t592;
t559 = -m(3) * t520 + mrSges(3,1) * t585 + t608 * mrSges(3,3) - t565;
t574 = -t603 - t604;
t461 = m(2) * t533 - t558 * mrSges(2,2) + (mrSges(2,1) + t574) * qJDD(1) + t559;
t593 = t555 * t454 + t557 * t461;
t455 = t553 * t459 + t551 * t466;
t590 = (t611 * t551 + t610 * t553) * qJD(1);
t578 = t557 * t454 - t555 * t461;
t573 = -Ifges(5,5) * t550 - Ifges(5,6) * t552;
t561 = Ifges(5,5) * t551 + (-Ifges(5,1) * t550 - Ifges(5,4) * t552) * t553;
t560 = Ifges(5,6) * t551 + (-Ifges(5,4) * t550 - Ifges(5,2) * t552) * t553;
t509 = t561 * qJD(1);
t508 = t560 * qJD(1);
t507 = (Ifges(5,3) * t551 + t573 * t553) * qJD(1);
t494 = Ifges(6,1) * t514 + Ifges(6,4) * t513 + Ifges(6,5) * t530;
t493 = Ifges(6,4) * t514 + Ifges(6,2) * t513 + Ifges(6,6) * t530;
t492 = Ifges(6,5) * t514 + Ifges(6,6) * t513 + Ifges(6,3) * t530;
t471 = mrSges(6,2) * t481 - mrSges(6,3) * t479 + Ifges(6,1) * t501 + Ifges(6,4) * t500 + Ifges(6,5) * t529 + t513 * t492 - t530 * t493;
t470 = -mrSges(6,1) * t481 + mrSges(6,3) * t480 + Ifges(6,4) * t501 + Ifges(6,2) * t500 + Ifges(6,6) * t529 - t514 * t492 + t530 * t494;
t462 = mrSges(4,2) * t585 + t565;
t457 = -mrSges(5,1) * t489 - mrSges(6,1) * t479 + mrSges(6,2) * t480 + mrSges(5,3) * t484 - Ifges(6,5) * t501 - Ifges(6,6) * t500 - Ifges(6,3) * t529 - pkin(4) * t469 - t514 * t493 + t513 * t494 + (t507 * t597 + t551 * t509) * qJD(1) + t560 * qJDD(1);
t456 = mrSges(5,2) * t489 - mrSges(5,3) * t483 - pkin(6) * t469 - t554 * t470 + t556 * t471 + (-t507 * t594 - t508 * t551) * qJD(1) + t561 * qJDD(1);
t451 = -qJ(3) * t462 - mrSges(3,3) * t505 + mrSges(3,2) * t520 + pkin(3) * t463 + mrSges(4,1) * t496 - mrSges(4,3) * t502 + pkin(4) * t563 + pkin(6) * t576 + mrSges(5,1) * t483 - mrSges(5,2) * t484 + t554 * t471 + t556 * t470 + (Ifges(5,3) + t607) * t586 + ((t573 + t606) * qJDD(1) + (-t508 * t550 + t509 * t552 - t590) * qJD(1)) * t553;
t450 = -mrSges(3,1) * t520 + mrSges(3,3) * t506 - mrSges(4,1) * t497 + mrSges(4,2) * t502 - t550 * t456 - t552 * t457 - pkin(3) * t571 - qJ(4) * t592 - pkin(2) * t462 + (-pkin(3) * t521 * t597 + t590 * t551) * qJD(1) + (t606 * t551 + (-pkin(3) * t602 - t605) * t553) * qJDD(1);
t449 = mrSges(2,1) * g(3) - pkin(1) * t455 + mrSges(2,3) * t534 + t558 * Ifges(2,5) - pkin(2) * t564 - qJ(3) * t562 + t550 * t457 + qJ(4) * t463 - mrSges(3,1) * t505 + mrSges(3,2) * t506 - t552 * t456 - mrSges(4,2) * t496 + mrSges(4,3) * t497 + (Ifges(2,6) + (qJ(3) * t602 + t610) * t553 + (mrSges(4,1) * pkin(2) + t611) * t551) * qJDD(1) + ((qJ(3) * t600 + t606 * t588) * t553 + (pkin(2) * t525 - t606 * t589 + (t605 + t607) * t588) * t551) * qJD(1);
t448 = -mrSges(2,2) * g(3) - mrSges(2,3) * t533 + Ifges(2,5) * qJDD(1) - t558 * Ifges(2,6) - qJ(2) * t455 - t551 * t450 + t553 * t451;
t1 = [-m(1) * g(1) + t578; -m(1) * g(2) + t593; (-m(1) - m(2)) * g(3) + t455; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t593 + t557 * t448 - t555 * t449; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t578 + t555 * t448 + t557 * t449; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + mrSges(2,1) * t533 - mrSges(2,2) * t534 + t551 * t451 + t553 * t450 + pkin(1) * t559 + qJ(2) * t577 + (pkin(1) * t574 + Ifges(2,3)) * qJDD(1);];
tauB = t1;
