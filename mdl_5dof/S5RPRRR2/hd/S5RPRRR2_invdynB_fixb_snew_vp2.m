% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPRRR2
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-05 18:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPRRR2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR2_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR2_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR2_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR2_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR2_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR2_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR2_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:11:54
% EndTime: 2019-12-05 18:12:03
% DurationCPUTime: 8.21s
% Computational Cost: add. (111534->291), mult. (276664->366), div. (0->0), fcn. (209796->10), ass. (0->122)
t587 = qJD(1) ^ 2;
t578 = cos(pkin(9));
t612 = pkin(2) * t578;
t577 = sin(pkin(9));
t611 = mrSges(3,2) * t577;
t575 = t578 ^ 2;
t610 = t575 * t587;
t582 = sin(qJ(1));
t586 = cos(qJ(1));
t564 = -t586 * g(1) - t582 * g(2);
t560 = -t587 * pkin(1) + qJDD(1) * qJ(2) + t564;
t606 = qJD(1) * qJD(2);
t604 = -t578 * g(3) - 0.2e1 * t577 * t606;
t535 = (-pkin(6) * qJDD(1) + t587 * t612 - t560) * t577 + t604;
t551 = -g(3) * t577 + (t560 + 0.2e1 * t606) * t578;
t605 = qJDD(1) * t578;
t536 = -pkin(2) * t610 + pkin(6) * t605 + t551;
t581 = sin(qJ(3));
t585 = cos(qJ(3));
t517 = t585 * t535 - t581 * t536;
t593 = t577 * t585 + t578 * t581;
t592 = -t577 * t581 + t578 * t585;
t558 = t592 * qJD(1);
t607 = t558 * qJD(3);
t549 = t593 * qJDD(1) + t607;
t559 = t593 * qJD(1);
t502 = (-t549 + t607) * pkin(7) + (t558 * t559 + qJDD(3)) * pkin(3) + t517;
t518 = t581 * t535 + t585 * t536;
t548 = -t559 * qJD(3) + t592 * qJDD(1);
t554 = qJD(3) * pkin(3) - pkin(7) * t559;
t557 = t558 ^ 2;
t510 = -pkin(3) * t557 + pkin(7) * t548 - qJD(3) * t554 + t518;
t580 = sin(qJ(4));
t584 = cos(qJ(4));
t494 = t584 * t502 - t580 * t510;
t541 = t558 * t584 - t559 * t580;
t516 = qJD(4) * t541 + t548 * t580 + t549 * t584;
t542 = t558 * t580 + t559 * t584;
t573 = qJDD(3) + qJDD(4);
t576 = qJD(3) + qJD(4);
t490 = (t541 * t576 - t516) * pkin(8) + (t541 * t542 + t573) * pkin(4) + t494;
t495 = t580 * t502 + t584 * t510;
t515 = -qJD(4) * t542 + t548 * t584 - t549 * t580;
t533 = pkin(4) * t576 - pkin(8) * t542;
t537 = t541 ^ 2;
t491 = -pkin(4) * t537 + pkin(8) * t515 - t533 * t576 + t495;
t579 = sin(qJ(5));
t583 = cos(qJ(5));
t488 = t490 * t583 - t491 * t579;
t526 = t541 * t583 - t542 * t579;
t499 = qJD(5) * t526 + t515 * t579 + t516 * t583;
t527 = t541 * t579 + t542 * t583;
t508 = -mrSges(6,1) * t526 + mrSges(6,2) * t527;
t571 = qJD(5) + t576;
t519 = -mrSges(6,2) * t571 + mrSges(6,3) * t526;
t570 = qJDD(5) + t573;
t486 = m(6) * t488 + mrSges(6,1) * t570 - mrSges(6,3) * t499 - t508 * t527 + t519 * t571;
t489 = t490 * t579 + t491 * t583;
t498 = -qJD(5) * t527 + t515 * t583 - t516 * t579;
t520 = mrSges(6,1) * t571 - mrSges(6,3) * t527;
t487 = m(6) * t489 - mrSges(6,2) * t570 + mrSges(6,3) * t498 + t508 * t526 - t520 * t571;
t478 = t583 * t486 + t579 * t487;
t528 = -mrSges(5,1) * t541 + mrSges(5,2) * t542;
t531 = -mrSges(5,2) * t576 + mrSges(5,3) * t541;
t476 = m(5) * t494 + mrSges(5,1) * t573 - mrSges(5,3) * t516 - t528 * t542 + t531 * t576 + t478;
t532 = mrSges(5,1) * t576 - mrSges(5,3) * t542;
t599 = -t486 * t579 + t583 * t487;
t477 = m(5) * t495 - mrSges(5,2) * t573 + mrSges(5,3) * t515 + t528 * t541 - t532 * t576 + t599;
t472 = t584 * t476 + t580 * t477;
t545 = -mrSges(4,1) * t558 + mrSges(4,2) * t559;
t552 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t558;
t470 = m(4) * t517 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t549 + qJD(3) * t552 - t545 * t559 + t472;
t553 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t559;
t600 = -t476 * t580 + t584 * t477;
t471 = m(4) * t518 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t548 - qJD(3) * t553 + t545 * t558 + t600;
t464 = t585 * t470 + t581 * t471;
t550 = -t577 * t560 + t604;
t591 = mrSges(3,3) * qJDD(1) + t587 * (-mrSges(3,1) * t578 + t611);
t462 = m(3) * t550 - t591 * t577 + t464;
t601 = -t581 * t470 + t585 * t471;
t463 = m(3) * t551 + t591 * t578 + t601;
t602 = -t462 * t577 + t578 * t463;
t455 = m(2) * t564 - mrSges(2,1) * t587 - qJDD(1) * mrSges(2,2) + t602;
t563 = t582 * g(1) - t586 * g(2);
t598 = qJDD(2) - t563;
t556 = -qJDD(1) * pkin(1) - t587 * qJ(2) + t598;
t574 = t577 ^ 2;
t547 = (-pkin(1) - t612) * qJDD(1) + (-qJ(2) + (-t574 - t575) * pkin(6)) * t587 + t598;
t512 = -t548 * pkin(3) - t557 * pkin(7) + t559 * t554 + t547;
t493 = -t515 * pkin(4) - t537 * pkin(8) + t542 * t533 + t512;
t594 = m(6) * t493 - t498 * mrSges(6,1) + t499 * mrSges(6,2) - t526 * t519 + t527 * t520;
t590 = m(5) * t512 - t515 * mrSges(5,1) + t516 * mrSges(5,2) - t541 * t531 + t542 * t532 + t594;
t589 = m(4) * t547 - t548 * mrSges(4,1) + t549 * mrSges(4,2) - t558 * t552 + t559 * t553 + t590;
t588 = -m(3) * t556 + mrSges(3,1) * t605 - t589 + (t574 * t587 + t610) * mrSges(3,3);
t482 = t588 - t587 * mrSges(2,2) + m(2) * t563 + (mrSges(2,1) - t611) * qJDD(1);
t609 = t582 * t455 + t586 * t482;
t456 = t578 * t462 + t577 * t463;
t595 = Ifges(3,5) * t577 + Ifges(3,6) * t578;
t608 = t587 * t595;
t603 = t586 * t455 - t482 * t582;
t597 = Ifges(3,1) * t577 + Ifges(3,4) * t578;
t596 = Ifges(3,4) * t577 + Ifges(3,2) * t578;
t540 = Ifges(4,1) * t559 + Ifges(4,4) * t558 + Ifges(4,5) * qJD(3);
t539 = Ifges(4,4) * t559 + Ifges(4,2) * t558 + Ifges(4,6) * qJD(3);
t538 = Ifges(4,5) * t559 + Ifges(4,6) * t558 + Ifges(4,3) * qJD(3);
t523 = Ifges(5,1) * t542 + Ifges(5,4) * t541 + Ifges(5,5) * t576;
t522 = Ifges(5,4) * t542 + Ifges(5,2) * t541 + Ifges(5,6) * t576;
t521 = Ifges(5,5) * t542 + Ifges(5,6) * t541 + Ifges(5,3) * t576;
t505 = Ifges(6,1) * t527 + Ifges(6,4) * t526 + Ifges(6,5) * t571;
t504 = Ifges(6,4) * t527 + Ifges(6,2) * t526 + Ifges(6,6) * t571;
t503 = Ifges(6,5) * t527 + Ifges(6,6) * t526 + Ifges(6,3) * t571;
t480 = mrSges(6,2) * t493 - mrSges(6,3) * t488 + Ifges(6,1) * t499 + Ifges(6,4) * t498 + Ifges(6,5) * t570 + t503 * t526 - t504 * t571;
t479 = -mrSges(6,1) * t493 + mrSges(6,3) * t489 + Ifges(6,4) * t499 + Ifges(6,2) * t498 + Ifges(6,6) * t570 - t503 * t527 + t505 * t571;
t466 = mrSges(5,2) * t512 - mrSges(5,3) * t494 + Ifges(5,1) * t516 + Ifges(5,4) * t515 + Ifges(5,5) * t573 - pkin(8) * t478 - t479 * t579 + t480 * t583 + t521 * t541 - t522 * t576;
t465 = -mrSges(5,1) * t512 + mrSges(5,3) * t495 + Ifges(5,4) * t516 + Ifges(5,2) * t515 + Ifges(5,6) * t573 - pkin(4) * t594 + pkin(8) * t599 + t583 * t479 + t579 * t480 - t542 * t521 + t576 * t523;
t458 = mrSges(4,2) * t547 - mrSges(4,3) * t517 + Ifges(4,1) * t549 + Ifges(4,4) * t548 + Ifges(4,5) * qJDD(3) - pkin(7) * t472 - qJD(3) * t539 - t465 * t580 + t466 * t584 + t538 * t558;
t457 = -mrSges(4,1) * t547 + mrSges(4,3) * t518 + Ifges(4,4) * t549 + Ifges(4,2) * t548 + Ifges(4,6) * qJDD(3) - pkin(3) * t590 + pkin(7) * t600 + qJD(3) * t540 + t584 * t465 + t580 * t466 - t559 * t538;
t452 = -Ifges(4,3) * qJDD(3) + mrSges(2,1) * g(3) - Ifges(5,3) * t573 - t559 * t539 + mrSges(2,3) * t564 - Ifges(6,3) * t570 - Ifges(4,5) * t549 - mrSges(3,1) * t550 + mrSges(3,2) * t551 + t558 * t540 + t541 * t523 - t542 * t522 - Ifges(4,6) * t548 + mrSges(4,2) * t518 + t526 * t505 - t527 * t504 + (-t577 * t596 + t578 * t597 + Ifges(2,5)) * t587 - Ifges(5,6) * t515 - Ifges(5,5) * t516 - mrSges(4,1) * t517 - Ifges(6,6) * t498 - Ifges(6,5) * t499 - mrSges(5,1) * t494 + mrSges(5,2) * t495 + mrSges(6,2) * t489 - mrSges(6,1) * t488 - pkin(4) * t478 - pkin(3) * t472 - pkin(2) * t464 + (Ifges(2,6) - t595) * qJDD(1) - pkin(1) * t456;
t451 = mrSges(3,2) * t556 - mrSges(3,3) * t550 - pkin(6) * t464 + t597 * qJDD(1) - t581 * t457 + t585 * t458 + t578 * t608;
t450 = -mrSges(3,1) * t556 + mrSges(3,3) * t551 - pkin(2) * t589 + pkin(6) * t601 + t596 * qJDD(1) + t585 * t457 + t581 * t458 - t577 * t608;
t449 = -mrSges(2,2) * g(3) - mrSges(2,3) * t563 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t587 - qJ(2) * t456 - t450 * t577 + t451 * t578;
t1 = [-m(1) * g(1) + t603; -m(1) * g(2) + t609; (-m(1) - m(2)) * g(3) + t456; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t609 + t586 * t449 - t582 * t452; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t603 + t582 * t449 + t586 * t452; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(2,3) * qJDD(1) + mrSges(2,1) * t563 - mrSges(2,2) * t564 + t577 * t451 + t578 * t450 + pkin(1) * (-qJDD(1) * t611 + t588) + qJ(2) * t602;];
tauB = t1;
