% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PRRRR12
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d2,d3,d4,d5,theta1]';
% m [6x1]
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
% tauJB [(6+5)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-28 18:09
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PRRRR12_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR12_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR12_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR12_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR12_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR12_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR12_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR12_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR12_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-28 18:07:58
% EndTime: 2024-09-28 18:08:06
% DurationCPUTime: 3.06s
% Computational Cost: add. (135606->195), mult. (184009->268), div. (0->0), fcn. (140574->14), ass. (0->105)
t590 = sin(pkin(6));
t625 = pkin(10) * t590;
t587 = qJD(2) + qJD(3);
t583 = qJD(4) + t587;
t624 = t583 * t590;
t595 = sin(qJ(5));
t623 = t590 * t595;
t599 = cos(qJ(5));
t622 = t590 * t599;
t591 = sin(pkin(5));
t598 = sin(qJ(2));
t621 = t591 * t598;
t602 = cos(qJ(2));
t620 = t591 * t602;
t594 = cos(pkin(5));
t619 = t594 * t598;
t618 = t594 * t602;
t589 = sin(pkin(11));
t592 = cos(pkin(11));
t575 = t589 * g(1) - t592 * g(2);
t576 = -t592 * g(1) - t589 * g(2);
t588 = -g(3) + qJDD(1);
t557 = t575 * t618 - t598 * t576 + t588 * t620;
t555 = qJDD(2) * pkin(2) + t557;
t558 = t575 * t619 + t602 * t576 + t588 * t621;
t603 = qJD(2) ^ 2;
t556 = -t603 * pkin(2) + t558;
t597 = sin(qJ(3));
t601 = cos(qJ(3));
t550 = t601 * t555 - t597 * t556;
t586 = qJDD(2) + qJDD(3);
t547 = t586 * pkin(3) + t550;
t551 = t597 * t555 + t601 * t556;
t585 = t587 ^ 2;
t548 = -t585 * pkin(3) + t551;
t596 = sin(qJ(4));
t600 = cos(qJ(4));
t543 = t596 * t547 + t600 * t548;
t581 = t583 ^ 2;
t582 = qJDD(4) + t586;
t540 = -t581 * pkin(4) + t582 * t625 + t543;
t542 = t600 * t547 - t596 * t548;
t539 = t582 * pkin(4) + t581 * t625 + t542;
t569 = -t591 * t575 + t594 * t588;
t593 = cos(pkin(6));
t607 = t539 * t593 + t569 * t590;
t534 = -t595 * t540 + t607 * t599;
t574 = t593 * t583 + qJD(5);
t614 = t583 * t622;
t563 = -t574 * mrSges(6,2) + mrSges(6,3) * t614;
t564 = (-mrSges(6,1) * t599 + mrSges(6,2) * t595) * t624;
t616 = qJD(5) * t583;
t565 = (t582 * t595 + t599 * t616) * t590;
t573 = t593 * t582 + qJDD(5);
t615 = t583 * t623;
t532 = m(6) * t534 + t573 * mrSges(6,1) - t565 * mrSges(6,3) + t574 * t563 - t564 * t615;
t535 = t599 * t540 + t607 * t595;
t562 = t574 * mrSges(6,1) - mrSges(6,3) * t615;
t566 = (t582 * t599 - t595 * t616) * t590;
t533 = m(6) * t535 - t573 * mrSges(6,2) + t566 * mrSges(6,3) - t574 * t562 + t564 * t614;
t538 = -t590 * t539 + t593 * t569;
t537 = m(6) * t538 - t566 * mrSges(6,1) + t565 * mrSges(6,2) + (t562 * t595 - t563 * t599) * t624;
t515 = -t590 * t537 + (t532 * t599 + t533 * t595) * t593;
t510 = m(5) * t542 + t582 * mrSges(5,1) - t581 * mrSges(5,2) + t515;
t520 = -t595 * t532 + t599 * t533;
t518 = m(5) * t543 - t581 * mrSges(5,1) - t582 * mrSges(5,2) + t520;
t508 = t600 * t510 + t596 * t518;
t505 = m(4) * t550 + t586 * mrSges(4,1) - t585 * mrSges(4,2) + t508;
t611 = -t596 * t510 + t600 * t518;
t506 = m(4) * t551 - t585 * mrSges(4,1) - t586 * mrSges(4,2) + t611;
t499 = t601 * t505 + t597 * t506;
t497 = m(3) * t557 + qJDD(2) * mrSges(3,1) - t603 * mrSges(3,2) + t499;
t612 = -t597 * t505 + t601 * t506;
t498 = m(3) * t558 - t603 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t612;
t514 = t532 * t622 + t533 * t623 + t593 * t537;
t610 = m(5) * t569 + t514;
t608 = m(4) * t569 + t610;
t512 = m(3) * t569 + t608;
t486 = t497 * t618 + t498 * t619 - t591 * t512;
t484 = m(2) * t575 + t486;
t491 = -t598 * t497 + t602 * t498;
t490 = m(2) * t576 + t491;
t617 = t592 * t484 + t589 * t490;
t485 = t497 * t620 + t498 * t621 + t594 * t512;
t613 = -t589 * t484 + t592 * t490;
t609 = m(2) * t588 + t485;
t560 = Ifges(6,6) * t574 + (Ifges(6,4) * t595 + Ifges(6,2) * t599) * t624;
t561 = Ifges(6,5) * t574 + (Ifges(6,1) * t595 + Ifges(6,4) * t599) * t624;
t522 = mrSges(6,1) * t534 - mrSges(6,2) * t535 + Ifges(6,5) * t565 + Ifges(6,6) * t566 + Ifges(6,3) * t573 + (t560 * t595 - t561 * t599) * t624;
t559 = Ifges(6,3) * t574 + (Ifges(6,5) * t595 + Ifges(6,6) * t599) * t624;
t525 = -mrSges(6,1) * t538 + mrSges(6,3) * t535 + Ifges(6,4) * t565 + Ifges(6,2) * t566 + Ifges(6,6) * t573 - t559 * t615 + t574 * t561;
t526 = mrSges(6,2) * t538 - mrSges(6,3) * t534 + Ifges(6,1) * t565 + Ifges(6,4) * t566 + Ifges(6,5) * t573 + t559 * t614 - t574 * t560;
t500 = -mrSges(5,1) * t569 + mrSges(5,3) * t543 + t581 * Ifges(5,5) + Ifges(5,6) * t582 - pkin(4) * t514 - t590 * t522 + (pkin(10) * t520 + t525 * t599 + t526 * t595) * t593;
t501 = mrSges(5,2) * t569 - mrSges(5,3) * t542 + Ifges(5,5) * t582 - t581 * Ifges(5,6) - t595 * t525 + t599 * t526 + (-t514 * t590 - t515 * t593) * pkin(10);
t482 = -mrSges(4,1) * t569 + mrSges(4,3) * t551 + t585 * Ifges(4,5) + Ifges(4,6) * t586 - pkin(3) * t610 + pkin(9) * t611 + t600 * t500 + t596 * t501;
t487 = mrSges(4,2) * t569 - mrSges(4,3) * t550 + Ifges(4,5) * t586 - t585 * Ifges(4,6) - pkin(9) * t508 - t596 * t500 + t600 * t501;
t479 = -mrSges(3,1) * t569 + mrSges(3,3) * t558 + t603 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t608 + pkin(8) * t612 + t601 * t482 + t597 * t487;
t480 = mrSges(3,2) * t569 - mrSges(3,3) * t557 + Ifges(3,5) * qJDD(2) - t603 * Ifges(3,6) - pkin(8) * t499 - t597 * t482 + t601 * t487;
t606 = pkin(7) * t491 + t479 * t602 + t480 * t598;
t605 = mrSges(5,1) * t542 - mrSges(5,2) * t543 + Ifges(5,3) * t582 + pkin(4) * t515 + t520 * t625 + t593 * t522 + t525 * t622 + t526 * t623;
t604 = mrSges(4,1) * t550 - mrSges(4,2) * t551 + Ifges(4,3) * t586 + pkin(3) * t508 + t605;
t481 = mrSges(3,1) * t557 - mrSges(3,2) * t558 + Ifges(3,3) * qJDD(2) + pkin(2) * t499 + t604;
t478 = mrSges(2,2) * t588 - mrSges(2,3) * t575 - t598 * t479 + t602 * t480 + (-t485 * t591 - t486 * t594) * pkin(7);
t477 = -mrSges(2,1) * t588 + mrSges(2,3) * t576 - pkin(1) * t485 - t591 * t481 + t606 * t594;
t1 = [-m(1) * g(1) + t613; -m(1) * g(2) + t617; -m(1) * g(3) + t609; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t617 - t589 * t477 + t592 * t478; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t613 + t592 * t477 + t589 * t478; -mrSges(1,1) * g(2) + mrSges(2,1) * t575 + mrSges(1,2) * g(1) - mrSges(2,2) * t576 + pkin(1) * t486 + t594 * t481 + t606 * t591; t609; t481; t604; t605; t522;];
tauJB = t1;
