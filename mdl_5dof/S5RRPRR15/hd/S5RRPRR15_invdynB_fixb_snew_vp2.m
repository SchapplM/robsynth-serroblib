% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRPRR15
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRPRR15_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR15_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR15_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR15_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR15_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR15_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR15_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR15_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR15_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:41:30
% EndTime: 2019-12-31 20:41:36
% DurationCPUTime: 3.25s
% Computational Cost: add. (30542->295), mult. (62710->356), div. (0->0), fcn. (36246->8), ass. (0->116)
t624 = -2 * qJD(3);
t623 = Ifges(3,1) + Ifges(4,2);
t619 = Ifges(3,4) + Ifges(4,6);
t618 = Ifges(3,5) - Ifges(4,4);
t622 = Ifges(3,2) + Ifges(4,3);
t617 = Ifges(3,6) - Ifges(4,5);
t621 = (Ifges(3,3) + Ifges(4,1));
t588 = sin(qJ(1));
t592 = cos(qJ(1));
t572 = -t592 * g(1) - t588 * g(2);
t594 = qJD(1) ^ 2;
t548 = -t594 * pkin(1) + qJDD(1) * pkin(6) + t572;
t587 = sin(qJ(2));
t591 = cos(qJ(2));
t534 = -t587 * g(3) + t591 * t548;
t559 = (-pkin(2) * t591 - qJ(3) * t587) * qJD(1);
t593 = qJD(2) ^ 2;
t611 = qJD(1) * t591;
t517 = t593 * pkin(2) - qJDD(2) * qJ(3) + (qJD(2) * t624) - t559 * t611 - t534;
t620 = t594 * pkin(6);
t533 = -t591 * g(3) - t587 * t548;
t560 = (mrSges(4,2) * t591 - mrSges(4,3) * t587) * qJD(1);
t561 = (-mrSges(3,1) * t591 + mrSges(3,2) * t587) * qJD(1);
t610 = qJD(1) * qJD(2);
t607 = t591 * t610;
t562 = t587 * qJDD(1) + t607;
t567 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t611;
t568 = -mrSges(4,1) * t611 - qJD(2) * mrSges(4,3);
t578 = t587 * qJD(1);
t608 = t587 * t610;
t563 = t591 * qJDD(1) - t608;
t570 = pkin(3) * t578 - (qJD(2) * pkin(7));
t584 = t591 ^ 2;
t571 = t588 * g(1) - t592 * g(2);
t603 = -qJDD(1) * pkin(1) - t571;
t598 = pkin(2) * t608 + t578 * t624 + (-t562 - t607) * qJ(3) + t603;
t505 = -t570 * t578 + (-pkin(3) * t584 - pkin(6)) * t594 + (-pkin(2) - pkin(7)) * t563 + t598;
t518 = -qJDD(2) * pkin(2) - t593 * qJ(3) + t559 * t578 + qJDD(3) - t533;
t510 = (-t587 * t591 * t594 - qJDD(2)) * pkin(7) + (t562 - t607) * pkin(3) + t518;
t586 = sin(qJ(4));
t590 = cos(qJ(4));
t497 = -t586 * t505 + t590 * t510;
t557 = -t586 * qJD(2) - t590 * t611;
t527 = t557 * qJD(4) + t590 * qJDD(2) - t586 * t563;
t556 = qJDD(4) + t562;
t558 = t590 * qJD(2) - t586 * t611;
t575 = t578 + qJD(4);
t495 = (t557 * t575 - t527) * pkin(8) + (t557 * t558 + t556) * pkin(4) + t497;
t498 = t590 * t505 + t586 * t510;
t526 = -t558 * qJD(4) - t586 * qJDD(2) - t590 * t563;
t535 = t575 * pkin(4) - t558 * pkin(8);
t555 = t557 ^ 2;
t496 = -t555 * pkin(4) + t526 * pkin(8) - t575 * t535 + t498;
t585 = sin(qJ(5));
t589 = cos(qJ(5));
t493 = t589 * t495 - t585 * t496;
t528 = t589 * t557 - t585 * t558;
t503 = t528 * qJD(5) + t585 * t526 + t589 * t527;
t529 = t585 * t557 + t589 * t558;
t515 = -t528 * mrSges(6,1) + t529 * mrSges(6,2);
t573 = qJD(5) + t575;
t519 = -t573 * mrSges(6,2) + t528 * mrSges(6,3);
t549 = qJDD(5) + t556;
t491 = m(6) * t493 + t549 * mrSges(6,1) - t503 * mrSges(6,3) - t529 * t515 + t573 * t519;
t494 = t585 * t495 + t589 * t496;
t502 = -t529 * qJD(5) + t589 * t526 - t585 * t527;
t520 = t573 * mrSges(6,1) - t529 * mrSges(6,3);
t492 = m(6) * t494 - t549 * mrSges(6,2) + t502 * mrSges(6,3) + t528 * t515 - t573 * t520;
t482 = t589 * t491 + t585 * t492;
t530 = -t557 * mrSges(5,1) + t558 * mrSges(5,2);
t531 = -t575 * mrSges(5,2) + t557 * mrSges(5,3);
t480 = m(5) * t497 + t556 * mrSges(5,1) - t527 * mrSges(5,3) - t558 * t530 + t575 * t531 + t482;
t532 = t575 * mrSges(5,1) - t558 * mrSges(5,3);
t604 = -t585 * t491 + t589 * t492;
t481 = m(5) * t498 - t556 * mrSges(5,2) + t526 * mrSges(5,3) + t557 * t530 - t575 * t532 + t604;
t477 = t590 * t480 + t586 * t481;
t599 = -m(4) * t518 - t562 * mrSges(4,1) - t477;
t475 = m(3) * t533 - t562 * mrSges(3,3) + (mrSges(3,1) - mrSges(4,2)) * qJDD(2) + (t567 - t568) * qJD(2) + (-t560 - t561) * t578 + t599;
t566 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t578;
t569 = mrSges(4,1) * t578 + qJD(2) * mrSges(4,2);
t509 = -t584 * t594 * pkin(7) + t563 * pkin(3) + qJD(2) * t570 - t517;
t500 = -t526 * pkin(4) - t555 * pkin(8) + t558 * t535 + t509;
t600 = m(6) * t500 - t502 * mrSges(6,1) + t503 * mrSges(6,2) - t528 * t519 + t529 * t520;
t597 = -m(5) * t509 + t526 * mrSges(5,1) - t527 * mrSges(5,2) + t557 * t531 - t558 * t532 - t600;
t596 = -m(4) * t517 + qJDD(2) * mrSges(4,3) + qJD(2) * t569 + t560 * t611 - t597;
t487 = t596 + (mrSges(3,3) + mrSges(4,1)) * t563 - qJD(2) * t566 + m(3) * t534 - qJDD(2) * mrSges(3,2) + t561 * t611;
t605 = -t587 * t475 + t591 * t487;
t468 = m(2) * t572 - t594 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t605;
t547 = t603 - t620;
t516 = -t563 * pkin(2) + t598 - t620;
t615 = -t586 * t480 + t590 * t481;
t602 = -m(4) * t516 - t563 * mrSges(4,2) + t569 * t578 - t615;
t595 = -m(3) * t547 + t567 * t611 + t563 * mrSges(3,1) + (-mrSges(3,2) + mrSges(4,3)) * t562 + (-t566 * t587 - t568 * t591) * qJD(1) + t602;
t473 = m(2) * t571 + qJDD(1) * mrSges(2,1) - t594 * mrSges(2,2) + t595;
t616 = t588 * t468 + t592 * t473;
t469 = t591 * t475 + t587 * t487;
t614 = (t621 * qJD(2)) + (t618 * t587 + t617 * t591) * qJD(1);
t613 = -t617 * qJD(2) + (-t619 * t587 - t622 * t591) * qJD(1);
t612 = t618 * qJD(2) + (t623 * t587 + t619 * t591) * qJD(1);
t606 = t592 * t468 - t588 * t473;
t523 = Ifges(5,1) * t558 + Ifges(5,4) * t557 + Ifges(5,5) * t575;
t522 = Ifges(5,4) * t558 + Ifges(5,2) * t557 + Ifges(5,6) * t575;
t521 = Ifges(5,5) * t558 + Ifges(5,6) * t557 + Ifges(5,3) * t575;
t513 = Ifges(6,1) * t529 + Ifges(6,4) * t528 + Ifges(6,5) * t573;
t512 = Ifges(6,4) * t529 + Ifges(6,2) * t528 + Ifges(6,6) * t573;
t511 = Ifges(6,5) * t529 + Ifges(6,6) * t528 + Ifges(6,3) * t573;
t484 = mrSges(6,2) * t500 - mrSges(6,3) * t493 + Ifges(6,1) * t503 + Ifges(6,4) * t502 + Ifges(6,5) * t549 + t528 * t511 - t573 * t512;
t483 = -mrSges(6,1) * t500 + mrSges(6,3) * t494 + Ifges(6,4) * t503 + Ifges(6,2) * t502 + Ifges(6,6) * t549 - t529 * t511 + t573 * t513;
t476 = -t562 * mrSges(4,3) + t568 * t611 - t602;
t471 = mrSges(5,2) * t509 - mrSges(5,3) * t497 + Ifges(5,1) * t527 + Ifges(5,4) * t526 + Ifges(5,5) * t556 - pkin(8) * t482 - t585 * t483 + t589 * t484 + t557 * t521 - t575 * t522;
t470 = -mrSges(5,1) * t509 + mrSges(5,3) * t498 + Ifges(5,4) * t527 + Ifges(5,2) * t526 + Ifges(5,6) * t556 - pkin(4) * t600 + pkin(8) * t604 + t589 * t483 + t585 * t484 - t558 * t521 + t575 * t523;
t465 = t618 * qJDD(2) + t619 * t563 + t613 * qJD(2) + t614 * t611 + Ifges(5,3) * t556 - t557 * t523 + t558 * t522 - mrSges(3,3) * t533 + mrSges(3,2) * t547 + Ifges(6,3) * t549 + Ifges(5,6) * t526 + Ifges(5,5) * t527 - t528 * t513 + t529 * t512 - mrSges(4,3) * t516 + mrSges(4,1) * t518 + Ifges(6,5) * t503 + mrSges(5,1) * t497 - mrSges(5,2) * t498 + Ifges(6,6) * t502 + mrSges(6,1) * t493 - mrSges(6,2) * t494 + pkin(4) * t482 + pkin(3) * t477 - qJ(3) * t476 + t623 * t562;
t464 = -mrSges(3,1) * t547 - mrSges(4,1) * t517 + mrSges(4,2) * t516 + mrSges(3,3) * t534 - pkin(2) * t476 - pkin(3) * t597 - pkin(7) * t615 + t612 * qJD(2) + t617 * qJDD(2) - t590 * t470 - t586 * t471 + t619 * t562 + t622 * t563 - t614 * t578;
t463 = -pkin(1) * t469 + mrSges(2,3) * t572 - pkin(2) * (-qJD(2) * t568 + t599) - qJ(3) * t596 - mrSges(3,1) * t533 + mrSges(3,2) * t534 + mrSges(4,3) * t517 - t590 * t471 + t586 * t470 + pkin(7) * t477 - mrSges(4,2) * t518 + mrSges(2,1) * g(3) + t594 * Ifges(2,5) + Ifges(2,6) * qJDD(1) + (-qJ(3) * mrSges(4,1) - t617) * t563 - t618 * t562 + (pkin(2) * mrSges(4,2) - t621) * qJDD(2) + (t612 * t591 + (pkin(2) * t560 + t613) * t587) * qJD(1);
t462 = -mrSges(2,2) * g(3) - mrSges(2,3) * t571 + Ifges(2,5) * qJDD(1) - t594 * Ifges(2,6) - pkin(6) * t469 - t587 * t464 + t591 * t465;
t1 = [-m(1) * g(1) + t606; -m(1) * g(2) + t616; (-m(1) - m(2)) * g(3) + t469; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t616 + t592 * t462 - t588 * t463; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t606 + t588 * t462 + t592 * t463; -mrSges(1,1) * g(2) + mrSges(2,1) * t571 + mrSges(1,2) * g(1) - mrSges(2,2) * t572 + Ifges(2,3) * qJDD(1) + pkin(1) * t595 + pkin(6) * t605 + t591 * t464 + t587 * t465;];
tauB = t1;
