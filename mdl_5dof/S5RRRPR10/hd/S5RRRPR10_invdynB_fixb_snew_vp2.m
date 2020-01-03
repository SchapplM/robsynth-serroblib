% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRRPR10
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
% Datum: 2019-12-31 21:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRRPR10_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR10_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR10_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR10_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR10_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR10_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR10_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR10_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR10_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:27:38
% EndTime: 2019-12-31 21:27:51
% DurationCPUTime: 12.82s
% Computational Cost: add. (197985->325), mult. (438560->427), div. (0->0), fcn. (342833->12), ass. (0->137)
t600 = sin(pkin(5));
t605 = sin(qJ(2));
t609 = cos(qJ(2));
t624 = qJD(1) * qJD(2);
t589 = (-qJDD(1) * t609 + t605 * t624) * t600;
t635 = 2 * qJD(4);
t634 = pkin(7) * t600;
t602 = cos(pkin(5));
t633 = t602 * g(3);
t632 = t600 * t605;
t631 = t600 * t609;
t630 = t602 * t605;
t629 = t602 * t609;
t606 = sin(qJ(1));
t610 = cos(qJ(1));
t592 = t606 * g(1) - g(2) * t610;
t611 = qJD(1) ^ 2;
t584 = qJDD(1) * pkin(1) + t611 * t634 + t592;
t593 = -g(1) * t610 - g(2) * t606;
t585 = -pkin(1) * t611 + qJDD(1) * t634 + t593;
t627 = t584 * t630 + t609 * t585;
t561 = -g(3) * t632 + t627;
t596 = qJD(1) * t602 + qJD(2);
t626 = qJD(1) * t600;
t623 = t605 * t626;
t582 = mrSges(3,1) * t596 - mrSges(3,3) * t623;
t586 = (-mrSges(3,1) * t609 + mrSges(3,2) * t605) * t626;
t595 = qJDD(1) * t602 + qJDD(2);
t587 = (-pkin(2) * t609 - pkin(8) * t605) * t626;
t594 = t596 ^ 2;
t625 = qJD(1) * t609;
t546 = -t594 * pkin(2) + t595 * pkin(8) + (-g(3) * t605 + t587 * t625) * t600 + t627;
t588 = (qJDD(1) * t605 + t609 * t624) * t600;
t547 = t589 * pkin(2) - t588 * pkin(8) - t633 + (-t584 + (pkin(2) * t605 - pkin(8) * t609) * t596 * qJD(1)) * t600;
t604 = sin(qJ(3));
t608 = cos(qJ(3));
t523 = -t604 * t546 + t608 * t547;
t577 = t596 * t608 - t604 * t623;
t559 = qJD(3) * t577 + t588 * t608 + t595 * t604;
t578 = t596 * t604 + t608 * t623;
t581 = qJDD(3) + t589;
t622 = t600 * t625;
t591 = qJD(3) - t622;
t516 = (t577 * t591 - t559) * qJ(4) + (t577 * t578 + t581) * pkin(3) + t523;
t524 = t608 * t546 + t604 * t547;
t558 = -qJD(3) * t578 - t588 * t604 + t595 * t608;
t568 = pkin(3) * t591 - qJ(4) * t578;
t576 = t577 ^ 2;
t518 = -pkin(3) * t576 + qJ(4) * t558 - t568 * t591 + t524;
t599 = sin(pkin(10));
t601 = cos(pkin(10));
t564 = t577 * t601 - t578 * t599;
t513 = t599 * t516 + t601 * t518 + t564 * t635;
t534 = t558 * t601 - t559 * t599;
t565 = t577 * t599 + t578 * t601;
t540 = -mrSges(5,1) * t564 + mrSges(5,2) * t565;
t551 = mrSges(5,1) * t591 - mrSges(5,3) * t565;
t541 = -pkin(4) * t564 - pkin(9) * t565;
t590 = t591 ^ 2;
t511 = -pkin(4) * t590 + pkin(9) * t581 + t541 * t564 + t513;
t560 = -g(3) * t631 + t584 * t629 - t605 * t585;
t545 = -t595 * pkin(2) - t594 * pkin(8) + t587 * t623 - t560;
t519 = -t558 * pkin(3) - t576 * qJ(4) + t578 * t568 + qJDD(4) + t545;
t535 = t558 * t599 + t559 * t601;
t514 = (-t564 * t591 - t535) * pkin(9) + (t565 * t591 - t534) * pkin(4) + t519;
t603 = sin(qJ(5));
t607 = cos(qJ(5));
t508 = -t511 * t603 + t514 * t607;
t548 = -t565 * t603 + t591 * t607;
t522 = qJD(5) * t548 + t535 * t607 + t581 * t603;
t549 = t565 * t607 + t591 * t603;
t529 = -mrSges(6,1) * t548 + mrSges(6,2) * t549;
t563 = qJD(5) - t564;
t530 = -mrSges(6,2) * t563 + mrSges(6,3) * t548;
t533 = qJDD(5) - t534;
t506 = m(6) * t508 + mrSges(6,1) * t533 - mrSges(6,3) * t522 - t529 * t549 + t530 * t563;
t509 = t511 * t607 + t514 * t603;
t521 = -qJD(5) * t549 - t535 * t603 + t581 * t607;
t531 = mrSges(6,1) * t563 - mrSges(6,3) * t549;
t507 = m(6) * t509 - mrSges(6,2) * t533 + mrSges(6,3) * t521 + t529 * t548 - t531 * t563;
t618 = -t506 * t603 + t507 * t607;
t497 = m(5) * t513 - mrSges(5,2) * t581 + mrSges(5,3) * t534 + t540 * t564 - t551 * t591 + t618;
t617 = -t601 * t516 + t599 * t518;
t512 = -0.2e1 * qJD(4) * t565 - t617;
t550 = -mrSges(5,2) * t591 + mrSges(5,3) * t564;
t510 = -t581 * pkin(4) - t590 * pkin(9) + (t635 + t541) * t565 + t617;
t614 = -m(6) * t510 + t521 * mrSges(6,1) - mrSges(6,2) * t522 + t548 * t530 - t531 * t549;
t502 = m(5) * t512 + mrSges(5,1) * t581 - mrSges(5,3) * t535 - t540 * t565 + t550 * t591 + t614;
t491 = t497 * t599 + t502 * t601;
t566 = -mrSges(4,1) * t577 + mrSges(4,2) * t578;
t567 = -mrSges(4,2) * t591 + mrSges(4,3) * t577;
t489 = m(4) * t523 + mrSges(4,1) * t581 - mrSges(4,3) * t559 - t566 * t578 + t567 * t591 + t491;
t569 = mrSges(4,1) * t591 - mrSges(4,3) * t578;
t619 = t497 * t601 - t502 * t599;
t490 = m(4) * t524 - mrSges(4,2) * t581 + mrSges(4,3) * t558 + t566 * t577 - t569 * t591 + t619;
t620 = -t489 * t604 + t490 * t608;
t480 = m(3) * t561 - mrSges(3,2) * t595 - mrSges(3,3) * t589 - t582 * t596 + t586 * t622 + t620;
t483 = t489 * t608 + t490 * t604;
t573 = -t600 * t584 - t633;
t583 = -mrSges(3,2) * t596 + mrSges(3,3) * t622;
t482 = m(3) * t573 + t589 * mrSges(3,1) + t588 * mrSges(3,2) + (t582 * t605 - t583 * t609) * t626 + t483;
t498 = t506 * t607 + t507 * t603;
t613 = m(5) * t519 - t534 * mrSges(5,1) + mrSges(5,2) * t535 - t564 * t550 + t551 * t565 + t498;
t612 = -m(4) * t545 + t558 * mrSges(4,1) - mrSges(4,2) * t559 + t577 * t567 - t569 * t578 - t613;
t494 = m(3) * t560 + mrSges(3,1) * t595 - mrSges(3,3) * t588 + t583 * t596 - t586 * t623 + t612;
t470 = t480 * t630 - t482 * t600 + t494 * t629;
t468 = m(2) * t592 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t611 + t470;
t476 = t480 * t609 - t494 * t605;
t475 = m(2) * t593 - mrSges(2,1) * t611 - qJDD(1) * mrSges(2,2) + t476;
t628 = t468 * t610 + t475 * t606;
t469 = t480 * t632 + t482 * t602 + t494 * t631;
t621 = -t468 * t606 + t475 * t610;
t525 = Ifges(6,5) * t549 + Ifges(6,6) * t548 + Ifges(6,3) * t563;
t527 = Ifges(6,1) * t549 + Ifges(6,4) * t548 + Ifges(6,5) * t563;
t499 = -mrSges(6,1) * t510 + mrSges(6,3) * t509 + Ifges(6,4) * t522 + Ifges(6,2) * t521 + Ifges(6,6) * t533 - t525 * t549 + t527 * t563;
t526 = Ifges(6,4) * t549 + Ifges(6,2) * t548 + Ifges(6,6) * t563;
t500 = mrSges(6,2) * t510 - mrSges(6,3) * t508 + Ifges(6,1) * t522 + Ifges(6,4) * t521 + Ifges(6,5) * t533 + t525 * t548 - t526 * t563;
t536 = Ifges(5,5) * t565 + Ifges(5,6) * t564 + Ifges(5,3) * t591;
t537 = Ifges(5,4) * t565 + Ifges(5,2) * t564 + Ifges(5,6) * t591;
t484 = mrSges(5,2) * t519 - mrSges(5,3) * t512 + Ifges(5,1) * t535 + Ifges(5,4) * t534 + Ifges(5,5) * t581 - pkin(9) * t498 - t499 * t603 + t500 * t607 + t536 * t564 - t537 * t591;
t538 = Ifges(5,1) * t565 + Ifges(5,4) * t564 + Ifges(5,5) * t591;
t485 = -mrSges(5,1) * t519 - mrSges(6,1) * t508 + mrSges(6,2) * t509 + mrSges(5,3) * t513 + Ifges(5,4) * t535 - Ifges(6,5) * t522 + Ifges(5,2) * t534 + Ifges(5,6) * t581 - Ifges(6,6) * t521 - Ifges(6,3) * t533 - pkin(4) * t498 - t526 * t549 + t527 * t548 - t536 * t565 + t538 * t591;
t552 = Ifges(4,5) * t578 + Ifges(4,6) * t577 + Ifges(4,3) * t591;
t554 = Ifges(4,1) * t578 + Ifges(4,4) * t577 + Ifges(4,5) * t591;
t471 = -mrSges(4,1) * t545 + mrSges(4,3) * t524 + Ifges(4,4) * t559 + Ifges(4,2) * t558 + Ifges(4,6) * t581 - pkin(3) * t613 + qJ(4) * t619 + t599 * t484 + t601 * t485 - t578 * t552 + t591 * t554;
t553 = Ifges(4,4) * t578 + Ifges(4,2) * t577 + Ifges(4,6) * t591;
t472 = mrSges(4,2) * t545 - mrSges(4,3) * t523 + Ifges(4,1) * t559 + Ifges(4,4) * t558 + Ifges(4,5) * t581 - qJ(4) * t491 + t484 * t601 - t485 * t599 + t552 * t577 - t553 * t591;
t570 = Ifges(3,3) * t596 + (Ifges(3,5) * t605 + Ifges(3,6) * t609) * t626;
t571 = Ifges(3,6) * t596 + (Ifges(3,4) * t605 + Ifges(3,2) * t609) * t626;
t465 = mrSges(3,2) * t573 - mrSges(3,3) * t560 + Ifges(3,1) * t588 - Ifges(3,4) * t589 + Ifges(3,5) * t595 - pkin(8) * t483 - t471 * t604 + t472 * t608 + t570 * t622 - t571 * t596;
t572 = Ifges(3,5) * t596 + (Ifges(3,1) * t605 + Ifges(3,4) * t609) * t626;
t466 = Ifges(3,6) * t595 + t596 * t572 - Ifges(4,6) * t558 - Ifges(4,5) * t559 + mrSges(3,3) * t561 + Ifges(3,4) * t588 - Ifges(3,2) * t589 - pkin(3) * t491 - pkin(2) * t483 - t607 * t499 - mrSges(4,1) * t523 + mrSges(4,2) * t524 - mrSges(5,1) * t512 + mrSges(5,2) * t513 + (-Ifges(4,3) - Ifges(5,3)) * t581 - Ifges(5,5) * t535 - t603 * t500 - pkin(9) * t618 + t564 * t538 - t565 * t537 - mrSges(3,1) * t573 - t570 * t623 - pkin(4) * t614 + t577 * t554 - t578 * t553 - Ifges(5,6) * t534;
t615 = pkin(7) * t476 + t465 * t605 + t466 * t609;
t464 = Ifges(3,5) * t588 - Ifges(3,6) * t589 + Ifges(3,3) * t595 + mrSges(3,1) * t560 - mrSges(3,2) * t561 + t604 * t472 + t608 * t471 + pkin(2) * t612 + pkin(8) * t620 + (t571 * t605 - t572 * t609) * t626;
t463 = -mrSges(2,2) * g(3) - mrSges(2,3) * t592 + Ifges(2,5) * qJDD(1) - t611 * Ifges(2,6) + t609 * t465 - t605 * t466 + (-t469 * t600 - t470 * t602) * pkin(7);
t462 = mrSges(2,1) * g(3) + mrSges(2,3) * t593 + t611 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t469 - t600 * t464 + t602 * t615;
t1 = [-m(1) * g(1) + t621; -m(1) * g(2) + t628; (-m(1) - m(2)) * g(3) + t469; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t628 - t462 * t606 + t463 * t610; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t621 + t610 * t462 + t606 * t463; -mrSges(1,1) * g(2) + mrSges(2,1) * t592 + mrSges(1,2) * g(1) - mrSges(2,2) * t593 + Ifges(2,3) * qJDD(1) + pkin(1) * t470 + t602 * t464 + t600 * t615;];
tauB = t1;
