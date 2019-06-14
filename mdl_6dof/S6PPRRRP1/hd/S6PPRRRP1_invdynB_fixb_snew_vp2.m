% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PPRRRP1
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 20:30
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PPRRRP1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRP1_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRRP1_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRP1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_invdynB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRP1_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRP1_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRRP1_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 20:27:25
% EndTime: 2019-05-04 20:27:37
% DurationCPUTime: 10.12s
% Computational Cost: add. (179406->272), mult. (322302->344), div. (0->0), fcn. (245965->14), ass. (0->125)
t624 = Ifges(6,1) + Ifges(7,1);
t618 = Ifges(6,4) + Ifges(7,4);
t617 = Ifges(6,5) + Ifges(7,5);
t623 = Ifges(6,2) + Ifges(7,2);
t622 = Ifges(6,6) + Ifges(7,6);
t621 = Ifges(6,3) + Ifges(7,3);
t574 = sin(pkin(11));
t578 = cos(pkin(11));
t568 = -t578 * g(1) - t574 * g(2);
t573 = sin(pkin(12));
t577 = cos(pkin(12));
t567 = t574 * g(1) - t578 * g(2);
t572 = -g(3) + qJDD(1);
t576 = sin(pkin(6));
t580 = cos(pkin(6));
t593 = t567 * t580 + t572 * t576;
t527 = -t573 * t568 + t577 * t593;
t528 = t577 * t568 + t573 * t593;
t551 = -t576 * t567 + t580 * t572 + qJDD(2);
t586 = cos(qJ(3));
t579 = cos(pkin(7));
t583 = sin(qJ(3));
t611 = t579 * t583;
t575 = sin(pkin(7));
t612 = t575 * t583;
t520 = t527 * t611 + t586 * t528 + t551 * t612;
t588 = qJD(3) ^ 2;
t518 = -t588 * pkin(3) + qJDD(3) * pkin(9) + t520;
t522 = -t575 * t527 + t579 * t551;
t582 = sin(qJ(4));
t585 = cos(qJ(4));
t511 = t585 * t518 + t582 * t522;
t563 = (-mrSges(5,1) * t585 + mrSges(5,2) * t582) * qJD(3);
t603 = qJD(3) * qJD(4);
t600 = t582 * t603;
t566 = t585 * qJDD(3) - t600;
t605 = qJD(3) * t582;
t569 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t605;
t564 = (-pkin(4) * t585 - pkin(10) * t582) * qJD(3);
t587 = qJD(4) ^ 2;
t604 = t585 * qJD(3);
t509 = -t587 * pkin(4) + qJDD(4) * pkin(10) + t564 * t604 + t511;
t519 = -t583 * t528 + (t527 * t579 + t551 * t575) * t586;
t517 = -qJDD(3) * pkin(3) - t588 * pkin(9) - t519;
t599 = t585 * t603;
t565 = t582 * qJDD(3) + t599;
t514 = (-t565 - t599) * pkin(10) + (-t566 + t600) * pkin(4) + t517;
t581 = sin(qJ(5));
t584 = cos(qJ(5));
t504 = -t581 * t509 + t584 * t514;
t561 = t584 * qJD(4) - t581 * t605;
t541 = t561 * qJD(5) + t581 * qJDD(4) + t584 * t565;
t562 = t581 * qJD(4) + t584 * t605;
t543 = -t561 * mrSges(7,1) + t562 * mrSges(7,2);
t544 = -t561 * mrSges(6,1) + t562 * mrSges(6,2);
t571 = qJD(5) - t604;
t547 = -t571 * mrSges(6,2) + t561 * mrSges(6,3);
t560 = qJDD(5) - t566;
t501 = -0.2e1 * qJD(6) * t562 + (t561 * t571 - t541) * qJ(6) + (t561 * t562 + t560) * pkin(5) + t504;
t546 = -t571 * mrSges(7,2) + t561 * mrSges(7,3);
t602 = m(7) * t501 + t560 * mrSges(7,1) + t571 * t546;
t494 = m(6) * t504 + t560 * mrSges(6,1) + t571 * t547 + (-t543 - t544) * t562 + (-mrSges(6,3) - mrSges(7,3)) * t541 + t602;
t505 = t584 * t509 + t581 * t514;
t540 = -t562 * qJD(5) + t584 * qJDD(4) - t581 * t565;
t548 = t571 * pkin(5) - t562 * qJ(6);
t559 = t561 ^ 2;
t503 = -t559 * pkin(5) + t540 * qJ(6) + 0.2e1 * qJD(6) * t561 - t571 * t548 + t505;
t601 = m(7) * t503 + t540 * mrSges(7,3) + t561 * t543;
t549 = t571 * mrSges(7,1) - t562 * mrSges(7,3);
t606 = -t571 * mrSges(6,1) + t562 * mrSges(6,3) - t549;
t619 = -mrSges(6,2) - mrSges(7,2);
t496 = m(6) * t505 + t540 * mrSges(6,3) + t561 * t544 + t560 * t619 + t606 * t571 + t601;
t596 = -t581 * t494 + t584 * t496;
t491 = m(5) * t511 - qJDD(4) * mrSges(5,2) + t566 * mrSges(5,3) - qJD(4) * t569 + t563 * t604 + t596;
t510 = -t582 * t518 + t585 * t522;
t570 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t604;
t508 = -qJDD(4) * pkin(4) - t587 * pkin(10) + t564 * t605 - t510;
t506 = -t540 * pkin(5) - t559 * qJ(6) + t562 * t548 + qJDD(6) + t508;
t595 = m(7) * t506 - t540 * mrSges(7,1) - t561 * t546;
t589 = -m(6) * t508 + t540 * mrSges(6,1) + t541 * t619 + t561 * t547 + t606 * t562 - t595;
t498 = m(5) * t510 + qJDD(4) * mrSges(5,1) - t565 * mrSges(5,3) + qJD(4) * t570 - t563 * t605 + t589;
t597 = t585 * t491 - t582 * t498;
t481 = m(4) * t520 - t588 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t597;
t484 = t582 * t491 + t585 * t498;
t483 = m(4) * t522 + t484;
t493 = t584 * t494 + t581 * t496;
t590 = -m(5) * t517 + t566 * mrSges(5,1) - t565 * mrSges(5,2) - t569 * t605 + t570 * t604 - t493;
t488 = m(4) * t519 + qJDD(3) * mrSges(4,1) - t588 * mrSges(4,2) + t590;
t613 = t488 * t586;
t470 = t481 * t611 - t575 * t483 + t579 * t613;
t466 = m(3) * t527 + t470;
t476 = t586 * t481 - t583 * t488;
t475 = m(3) * t528 + t476;
t620 = t466 * t577 + t475 * t573;
t469 = t481 * t612 + t579 * t483 + t575 * t613;
t468 = m(3) * t551 + t469;
t456 = -t576 * t468 + t580 * t620;
t454 = m(2) * t567 + t456;
t462 = -t573 * t466 + t577 * t475;
t461 = m(2) * t568 + t462;
t610 = t578 * t454 + t574 * t461;
t609 = t561 * t622 + t562 * t617 + t571 * t621;
t608 = -t561 * t623 - t562 * t618 - t571 * t622;
t607 = t618 * t561 + t562 * t624 + t617 * t571;
t455 = t580 * t468 + t576 * t620;
t598 = -t574 * t454 + t578 * t461;
t485 = -mrSges(6,1) * t508 + mrSges(6,3) * t505 - mrSges(7,1) * t506 + mrSges(7,3) * t503 - pkin(5) * t595 + qJ(6) * t601 + (-qJ(6) * t549 + t607) * t571 + (-pkin(5) * t549 - t609) * t562 + (-qJ(6) * mrSges(7,2) + t622) * t560 + (-pkin(5) * mrSges(7,2) + t618) * t541 + t623 * t540;
t499 = -t541 * mrSges(7,3) - t562 * t543 + t602;
t492 = mrSges(6,2) * t508 + mrSges(7,2) * t506 - mrSges(6,3) * t504 - mrSges(7,3) * t501 - qJ(6) * t499 + t618 * t540 + t541 * t624 + t617 * t560 + t609 * t561 + t608 * t571;
t553 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t582 + Ifges(5,6) * t585) * qJD(3);
t554 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t582 + Ifges(5,2) * t585) * qJD(3);
t471 = mrSges(5,2) * t517 - mrSges(5,3) * t510 + Ifges(5,1) * t565 + Ifges(5,4) * t566 + Ifges(5,5) * qJDD(4) - pkin(10) * t493 - qJD(4) * t554 - t581 * t485 + t584 * t492 + t553 * t604;
t555 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t582 + Ifges(5,4) * t585) * qJD(3);
t477 = -t553 * t605 - mrSges(5,1) * t517 - mrSges(6,1) * t504 - mrSges(7,1) * t501 + mrSges(6,2) * t505 + mrSges(7,2) * t503 + mrSges(5,3) * t511 + Ifges(5,4) * t565 + Ifges(5,2) * t566 + Ifges(5,6) * qJDD(4) - pkin(4) * t493 - pkin(5) * t499 + qJD(4) * t555 + t608 * t562 + t607 * t561 - t621 * t560 - t617 * t541 - t622 * t540;
t458 = mrSges(4,2) * t522 - mrSges(4,3) * t519 + Ifges(4,5) * qJDD(3) - t588 * Ifges(4,6) - pkin(9) * t484 + t585 * t471 - t582 * t477;
t463 = Ifges(4,6) * qJDD(3) + t588 * Ifges(4,5) - mrSges(4,1) * t522 + mrSges(4,3) * t520 - Ifges(5,5) * t565 - Ifges(5,6) * t566 - Ifges(5,3) * qJDD(4) - mrSges(5,1) * t510 + mrSges(5,2) * t511 - t581 * t492 - t584 * t485 - pkin(4) * t589 - pkin(10) * t596 - pkin(3) * t484 + (-t582 * t554 + t585 * t555) * qJD(3);
t592 = pkin(8) * t476 + t458 * t583 + t463 * t586;
t457 = mrSges(4,1) * t519 - mrSges(4,2) * t520 + Ifges(4,3) * qJDD(3) + pkin(3) * t590 + pkin(9) * t597 + t582 * t471 + t585 * t477;
t451 = -mrSges(3,1) * t551 + mrSges(3,3) * t528 - pkin(2) * t469 - t575 * t457 + t579 * t592;
t452 = mrSges(3,2) * t551 - mrSges(3,3) * t527 + t586 * t458 - t583 * t463 + (-t469 * t575 - t470 * t579) * pkin(8);
t591 = qJ(2) * t462 + t451 * t577 + t452 * t573;
t450 = mrSges(3,1) * t527 - mrSges(3,2) * t528 + pkin(2) * t470 + t579 * t457 + t575 * t592;
t449 = mrSges(2,2) * t572 - mrSges(2,3) * t567 - t573 * t451 + t577 * t452 + (-t455 * t576 - t456 * t580) * qJ(2);
t448 = -mrSges(2,1) * t572 + mrSges(2,3) * t568 - pkin(1) * t455 - t576 * t450 + t580 * t591;
t1 = [-m(1) * g(1) + t598; -m(1) * g(2) + t610; -m(1) * g(3) + m(2) * t572 + t455; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t610 - t574 * t448 + t578 * t449; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t598 + t578 * t448 + t574 * t449; -mrSges(1,1) * g(2) + mrSges(2,1) * t567 + mrSges(1,2) * g(1) - mrSges(2,2) * t568 + pkin(1) * t456 + t580 * t450 + t576 * t591;];
tauB  = t1;
