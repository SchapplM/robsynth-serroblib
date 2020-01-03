% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRRPR12
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-31 21:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRRPR12_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR12_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR12_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR12_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR12_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR12_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR12_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR12_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR12_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:38:07
% EndTime: 2019-12-31 21:38:20
% DurationCPUTime: 13.19s
% Computational Cost: add. (206516->324), mult. (452898->427), div. (0->0), fcn. (351986->12), ass. (0->135)
t594 = sin(pkin(5));
t599 = sin(qJ(2));
t602 = cos(qJ(2));
t616 = qJD(1) * qJD(2);
t581 = (-qJDD(1) * t602 + t599 * t616) * t594;
t627 = cos(qJ(3));
t626 = pkin(7) * t594;
t596 = cos(pkin(5));
t625 = t596 * g(3);
t624 = t594 * t599;
t623 = t594 * t602;
t622 = t596 * t599;
t621 = t596 * t602;
t600 = sin(qJ(1));
t603 = cos(qJ(1));
t586 = t600 * g(1) - t603 * g(2);
t604 = qJD(1) ^ 2;
t576 = qJDD(1) * pkin(1) + t604 * t626 + t586;
t587 = -t603 * g(1) - t600 * g(2);
t577 = -t604 * pkin(1) + qJDD(1) * t626 + t587;
t619 = t576 * t622 + t602 * t577;
t552 = -g(3) * t624 + t619;
t590 = t596 * qJD(1) + qJD(2);
t618 = qJD(1) * t594;
t615 = t599 * t618;
t574 = t590 * mrSges(3,1) - mrSges(3,3) * t615;
t578 = (-mrSges(3,1) * t602 + mrSges(3,2) * t599) * t618;
t589 = t596 * qJDD(1) + qJDD(2);
t579 = (-pkin(2) * t602 - pkin(8) * t599) * t618;
t588 = t590 ^ 2;
t617 = qJD(1) * t602;
t531 = -t588 * pkin(2) + t589 * pkin(8) + (-g(3) * t599 + t579 * t617) * t594 + t619;
t580 = (qJDD(1) * t599 + t602 * t616) * t594;
t532 = t581 * pkin(2) - t580 * pkin(8) - t625 + (-t576 + (pkin(2) * t599 - pkin(8) * t602) * t590 * qJD(1)) * t594;
t598 = sin(qJ(3));
t515 = t627 * t531 + t598 * t532;
t570 = t598 * t590 + t627 * t615;
t549 = t570 * qJD(3) + t598 * t580 - t627 * t589;
t569 = -t627 * t590 + t598 * t615;
t554 = t569 * mrSges(4,1) + t570 * mrSges(4,2);
t614 = t594 * t617;
t585 = qJD(3) - t614;
t561 = t585 * mrSges(4,1) - t570 * mrSges(4,3);
t573 = qJDD(3) + t581;
t553 = t569 * pkin(3) - t570 * qJ(4);
t584 = t585 ^ 2;
t507 = -t584 * pkin(3) + t573 * qJ(4) - t569 * t553 + t515;
t551 = -g(3) * t623 + t576 * t621 - t599 * t577;
t530 = -t589 * pkin(2) - t588 * pkin(8) + t579 * t615 - t551;
t550 = -t569 * qJD(3) + t627 * t580 + t598 * t589;
t510 = (t569 * t585 - t550) * qJ(4) + (t570 * t585 + t549) * pkin(3) + t530;
t593 = sin(pkin(10));
t595 = cos(pkin(10));
t559 = t595 * t570 + t593 * t585;
t502 = -0.2e1 * qJD(4) * t559 - t593 * t507 + t595 * t510;
t537 = t595 * t550 + t593 * t573;
t558 = -t593 * t570 + t595 * t585;
t500 = (t558 * t569 - t537) * pkin(9) + (t558 * t559 + t549) * pkin(4) + t502;
t503 = 0.2e1 * qJD(4) * t558 + t595 * t507 + t593 * t510;
t536 = -t593 * t550 + t595 * t573;
t542 = t569 * pkin(4) - t559 * pkin(9);
t557 = t558 ^ 2;
t501 = -t557 * pkin(4) + t536 * pkin(9) - t569 * t542 + t503;
t597 = sin(qJ(5));
t601 = cos(qJ(5));
t498 = t601 * t500 - t597 * t501;
t533 = t601 * t558 - t597 * t559;
t513 = t533 * qJD(5) + t597 * t536 + t601 * t537;
t534 = t597 * t558 + t601 * t559;
t520 = -t533 * mrSges(6,1) + t534 * mrSges(6,2);
t568 = qJD(5) + t569;
t521 = -t568 * mrSges(6,2) + t533 * mrSges(6,3);
t547 = qJDD(5) + t549;
t496 = m(6) * t498 + t547 * mrSges(6,1) - t513 * mrSges(6,3) - t534 * t520 + t568 * t521;
t499 = t597 * t500 + t601 * t501;
t512 = -t534 * qJD(5) + t601 * t536 - t597 * t537;
t522 = t568 * mrSges(6,1) - t534 * mrSges(6,3);
t497 = m(6) * t499 - t547 * mrSges(6,2) + t512 * mrSges(6,3) + t533 * t520 - t568 * t522;
t488 = t601 * t496 + t597 * t497;
t538 = -t558 * mrSges(5,1) + t559 * mrSges(5,2);
t540 = -t569 * mrSges(5,2) + t558 * mrSges(5,3);
t486 = m(5) * t502 + t549 * mrSges(5,1) - t537 * mrSges(5,3) - t559 * t538 + t569 * t540 + t488;
t541 = t569 * mrSges(5,1) - t559 * mrSges(5,3);
t610 = -t597 * t496 + t601 * t497;
t487 = m(5) * t503 - t549 * mrSges(5,2) + t536 * mrSges(5,3) + t558 * t538 - t569 * t541 + t610;
t611 = -t593 * t486 + t595 * t487;
t483 = m(4) * t515 - t573 * mrSges(4,2) - t549 * mrSges(4,3) - t569 * t554 - t585 * t561 + t611;
t514 = -t598 * t531 + t627 * t532;
t560 = -t585 * mrSges(4,2) - t569 * mrSges(4,3);
t506 = -t573 * pkin(3) - t584 * qJ(4) + t570 * t553 + qJDD(4) - t514;
t504 = -t536 * pkin(4) - t557 * pkin(9) + t559 * t542 + t506;
t607 = m(6) * t504 - t512 * mrSges(6,1) + t513 * mrSges(6,2) - t533 * t521 + t534 * t522;
t605 = -m(5) * t506 + t536 * mrSges(5,1) - t537 * mrSges(5,2) + t558 * t540 - t559 * t541 - t607;
t492 = m(4) * t514 + t573 * mrSges(4,1) - t550 * mrSges(4,3) - t570 * t554 + t585 * t560 + t605;
t612 = t627 * t483 - t598 * t492;
t472 = m(3) * t552 - t589 * mrSges(3,2) - t581 * mrSges(3,3) - t590 * t574 + t578 * t614 + t612;
t475 = t598 * t483 + t627 * t492;
t565 = -t594 * t576 - t625;
t575 = -t590 * mrSges(3,2) + mrSges(3,3) * t614;
t474 = m(3) * t565 + t581 * mrSges(3,1) + t580 * mrSges(3,2) + (t574 * t599 - t575 * t602) * t618 + t475;
t484 = t595 * t486 + t593 * t487;
t606 = -m(4) * t530 - t549 * mrSges(4,1) - t550 * mrSges(4,2) - t569 * t560 - t570 * t561 - t484;
t480 = m(3) * t551 + t589 * mrSges(3,1) - t580 * mrSges(3,3) + t590 * t575 - t578 * t615 + t606;
t462 = t472 * t622 - t594 * t474 + t480 * t621;
t460 = m(2) * t586 + qJDD(1) * mrSges(2,1) - t604 * mrSges(2,2) + t462;
t467 = t602 * t472 - t599 * t480;
t466 = m(2) * t587 - t604 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t467;
t620 = t603 * t460 + t600 * t466;
t461 = t472 * t624 + t596 * t474 + t480 * t623;
t613 = -t600 * t460 + t603 * t466;
t516 = Ifges(6,5) * t534 + Ifges(6,6) * t533 + Ifges(6,3) * t568;
t518 = Ifges(6,1) * t534 + Ifges(6,4) * t533 + Ifges(6,5) * t568;
t489 = -mrSges(6,1) * t504 + mrSges(6,3) * t499 + Ifges(6,4) * t513 + Ifges(6,2) * t512 + Ifges(6,6) * t547 - t534 * t516 + t568 * t518;
t517 = Ifges(6,4) * t534 + Ifges(6,2) * t533 + Ifges(6,6) * t568;
t490 = mrSges(6,2) * t504 - mrSges(6,3) * t498 + Ifges(6,1) * t513 + Ifges(6,4) * t512 + Ifges(6,5) * t547 + t533 * t516 - t568 * t517;
t523 = Ifges(5,5) * t559 + Ifges(5,6) * t558 + Ifges(5,3) * t569;
t525 = Ifges(5,1) * t559 + Ifges(5,4) * t558 + Ifges(5,5) * t569;
t476 = -mrSges(5,1) * t506 + mrSges(5,3) * t503 + Ifges(5,4) * t537 + Ifges(5,2) * t536 + Ifges(5,6) * t549 - pkin(4) * t607 + pkin(9) * t610 + t601 * t489 + t597 * t490 - t559 * t523 + t569 * t525;
t524 = Ifges(5,4) * t559 + Ifges(5,2) * t558 + Ifges(5,6) * t569;
t477 = mrSges(5,2) * t506 - mrSges(5,3) * t502 + Ifges(5,1) * t537 + Ifges(5,4) * t536 + Ifges(5,5) * t549 - pkin(9) * t488 - t597 * t489 + t601 * t490 + t558 * t523 - t569 * t524;
t543 = Ifges(4,5) * t570 - Ifges(4,6) * t569 + Ifges(4,3) * t585;
t544 = Ifges(4,4) * t570 - Ifges(4,2) * t569 + Ifges(4,6) * t585;
t463 = mrSges(4,2) * t530 - mrSges(4,3) * t514 + Ifges(4,1) * t550 - Ifges(4,4) * t549 + Ifges(4,5) * t573 - qJ(4) * t484 - t593 * t476 + t595 * t477 - t569 * t543 - t585 * t544;
t545 = Ifges(4,1) * t570 - Ifges(4,4) * t569 + Ifges(4,5) * t585;
t468 = Ifges(4,4) * t550 + Ifges(4,6) * t573 - t570 * t543 + t585 * t545 - mrSges(4,1) * t530 + mrSges(4,3) * t515 - Ifges(5,5) * t537 - Ifges(5,6) * t536 - t559 * t524 + t558 * t525 - mrSges(5,1) * t502 + mrSges(5,2) * t503 - Ifges(6,5) * t513 - Ifges(6,6) * t512 - Ifges(6,3) * t547 - t534 * t517 + t533 * t518 - mrSges(6,1) * t498 + mrSges(6,2) * t499 - pkin(4) * t488 - pkin(3) * t484 + (-Ifges(4,2) - Ifges(5,3)) * t549;
t562 = Ifges(3,3) * t590 + (Ifges(3,5) * t599 + Ifges(3,6) * t602) * t618;
t563 = Ifges(3,6) * t590 + (Ifges(3,4) * t599 + Ifges(3,2) * t602) * t618;
t457 = mrSges(3,2) * t565 - mrSges(3,3) * t551 + Ifges(3,1) * t580 - Ifges(3,4) * t581 + Ifges(3,5) * t589 - pkin(8) * t475 + t627 * t463 - t598 * t468 + t562 * t614 - t590 * t563;
t564 = Ifges(3,5) * t590 + (Ifges(3,1) * t599 + Ifges(3,4) * t602) * t618;
t458 = Ifges(3,4) * t580 - Ifges(3,2) * t581 + Ifges(3,6) * t589 - t562 * t615 + t590 * t564 - mrSges(3,1) * t565 + mrSges(3,3) * t552 - Ifges(4,5) * t550 + Ifges(4,6) * t549 - Ifges(4,3) * t573 - t570 * t544 - t569 * t545 - mrSges(4,1) * t514 + mrSges(4,2) * t515 - t593 * t477 - t595 * t476 - pkin(3) * t605 - qJ(4) * t611 - pkin(2) * t475;
t608 = pkin(7) * t467 + t457 * t599 + t458 * t602;
t456 = Ifges(3,5) * t580 - Ifges(3,6) * t581 + Ifges(3,3) * t589 + mrSges(3,1) * t551 - mrSges(3,2) * t552 + t598 * t463 + t627 * t468 + pkin(2) * t606 + pkin(8) * t612 + (t563 * t599 - t564 * t602) * t618;
t455 = -mrSges(2,2) * g(3) - mrSges(2,3) * t586 + Ifges(2,5) * qJDD(1) - t604 * Ifges(2,6) + t602 * t457 - t599 * t458 + (-t461 * t594 - t462 * t596) * pkin(7);
t454 = mrSges(2,1) * g(3) + mrSges(2,3) * t587 + t604 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t461 - t594 * t456 + t608 * t596;
t1 = [-m(1) * g(1) + t613; -m(1) * g(2) + t620; (-m(1) - m(2)) * g(3) + t461; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t620 - t600 * t454 + t603 * t455; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t613 + t603 * t454 + t600 * t455; -mrSges(1,1) * g(2) + mrSges(2,1) * t586 + mrSges(1,2) * g(1) - mrSges(2,2) * t587 + Ifges(2,3) * qJDD(1) + pkin(1) * t462 + t596 * t456 + t608 * t594;];
tauB = t1;
