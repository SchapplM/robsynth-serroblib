% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRRRR11
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
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-31 22:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRRRR11_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR11_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR11_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR11_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR11_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR11_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR11_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR11_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR11_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:39:54
% EndTime: 2019-12-31 22:40:08
% DurationCPUTime: 13.55s
% Computational Cost: add. (218519->325), mult. (467294->426), div. (0->0), fcn. (366008->12), ass. (0->137)
t603 = sin(pkin(5));
t608 = sin(qJ(2));
t613 = cos(qJ(2));
t627 = qJD(1) * qJD(2);
t592 = (-qJDD(1) * t613 + t608 * t627) * t603;
t637 = pkin(7) * t603;
t604 = cos(pkin(5));
t636 = t604 * g(3);
t635 = t603 * t608;
t634 = t603 * t613;
t633 = t604 * t608;
t632 = t604 * t613;
t609 = sin(qJ(1));
t614 = cos(qJ(1));
t596 = t609 * g(1) - t614 * g(2);
t615 = qJD(1) ^ 2;
t587 = qJDD(1) * pkin(1) + t615 * t637 + t596;
t597 = -t614 * g(1) - t609 * g(2);
t588 = -t615 * pkin(1) + qJDD(1) * t637 + t597;
t630 = t587 * t633 + t613 * t588;
t562 = -g(3) * t635 + t630;
t600 = t604 * qJD(1) + qJD(2);
t629 = qJD(1) * t603;
t626 = t608 * t629;
t585 = t600 * mrSges(3,1) - mrSges(3,3) * t626;
t589 = (-mrSges(3,1) * t613 + mrSges(3,2) * t608) * t629;
t599 = t604 * qJDD(1) + qJDD(2);
t590 = (-pkin(2) * t613 - pkin(8) * t608) * t629;
t598 = t600 ^ 2;
t628 = qJD(1) * t613;
t543 = -t598 * pkin(2) + t599 * pkin(8) + (-g(3) * t608 + t590 * t628) * t603 + t630;
t591 = (qJDD(1) * t608 + t613 * t627) * t603;
t544 = t592 * pkin(2) - t591 * pkin(8) - t636 + (-t587 + (pkin(2) * t608 - pkin(8) * t613) * t600 * qJD(1)) * t603;
t607 = sin(qJ(3));
t612 = cos(qJ(3));
t524 = t612 * t543 + t607 * t544;
t580 = t607 * t600 + t612 * t626;
t559 = -t580 * qJD(3) - t607 * t591 + t612 * t599;
t579 = t612 * t600 - t607 * t626;
t563 = -t579 * mrSges(4,1) + t580 * mrSges(4,2);
t625 = t603 * t628;
t595 = qJD(3) - t625;
t569 = t595 * mrSges(4,1) - t580 * mrSges(4,3);
t584 = qJDD(3) + t592;
t564 = -t579 * pkin(3) - t580 * pkin(9);
t593 = t595 ^ 2;
t519 = -t593 * pkin(3) + t584 * pkin(9) + t579 * t564 + t524;
t561 = -g(3) * t634 + t587 * t632 - t608 * t588;
t542 = -t599 * pkin(2) - t598 * pkin(8) + t590 * t626 - t561;
t560 = t579 * qJD(3) + t612 * t591 + t607 * t599;
t522 = (-t579 * t595 - t560) * pkin(9) + (t580 * t595 - t559) * pkin(3) + t542;
t606 = sin(qJ(4));
t611 = cos(qJ(4));
t511 = -t606 * t519 + t611 * t522;
t566 = -t606 * t580 + t611 * t595;
t532 = t566 * qJD(4) + t611 * t560 + t606 * t584;
t557 = qJDD(4) - t559;
t567 = t611 * t580 + t606 * t595;
t578 = qJD(4) - t579;
t509 = (t566 * t578 - t532) * pkin(10) + (t566 * t567 + t557) * pkin(4) + t511;
t512 = t611 * t519 + t606 * t522;
t531 = -t567 * qJD(4) - t606 * t560 + t611 * t584;
t551 = t578 * pkin(4) - t567 * pkin(10);
t565 = t566 ^ 2;
t510 = -t565 * pkin(4) + t531 * pkin(10) - t578 * t551 + t512;
t605 = sin(qJ(5));
t610 = cos(qJ(5));
t507 = t610 * t509 - t605 * t510;
t545 = t610 * t566 - t605 * t567;
t516 = t545 * qJD(5) + t605 * t531 + t610 * t532;
t546 = t605 * t566 + t610 * t567;
t529 = -t545 * mrSges(6,1) + t546 * mrSges(6,2);
t576 = qJD(5) + t578;
t533 = -t576 * mrSges(6,2) + t545 * mrSges(6,3);
t552 = qJDD(5) + t557;
t505 = m(6) * t507 + t552 * mrSges(6,1) - t516 * mrSges(6,3) - t546 * t529 + t576 * t533;
t508 = t605 * t509 + t610 * t510;
t515 = -t546 * qJD(5) + t610 * t531 - t605 * t532;
t534 = t576 * mrSges(6,1) - t546 * mrSges(6,3);
t506 = m(6) * t508 - t552 * mrSges(6,2) + t515 * mrSges(6,3) + t545 * t529 - t576 * t534;
t497 = t610 * t505 + t605 * t506;
t547 = -t566 * mrSges(5,1) + t567 * mrSges(5,2);
t549 = -t578 * mrSges(5,2) + t566 * mrSges(5,3);
t495 = m(5) * t511 + t557 * mrSges(5,1) - t532 * mrSges(5,3) - t567 * t547 + t578 * t549 + t497;
t550 = t578 * mrSges(5,1) - t567 * mrSges(5,3);
t621 = -t605 * t505 + t610 * t506;
t496 = m(5) * t512 - t557 * mrSges(5,2) + t531 * mrSges(5,3) + t566 * t547 - t578 * t550 + t621;
t622 = -t606 * t495 + t611 * t496;
t492 = m(4) * t524 - t584 * mrSges(4,2) + t559 * mrSges(4,3) + t579 * t563 - t595 * t569 + t622;
t523 = -t607 * t543 + t612 * t544;
t568 = -t595 * mrSges(4,2) + t579 * mrSges(4,3);
t518 = -t584 * pkin(3) - t593 * pkin(9) + t580 * t564 - t523;
t513 = -t531 * pkin(4) - t565 * pkin(10) + t567 * t551 + t518;
t618 = m(6) * t513 - t515 * mrSges(6,1) + t516 * mrSges(6,2) - t545 * t533 + t546 * t534;
t616 = -m(5) * t518 + t531 * mrSges(5,1) - t532 * mrSges(5,2) + t566 * t549 - t567 * t550 - t618;
t501 = m(4) * t523 + t584 * mrSges(4,1) - t560 * mrSges(4,3) - t580 * t563 + t595 * t568 + t616;
t623 = t612 * t492 - t607 * t501;
t481 = m(3) * t562 - t599 * mrSges(3,2) - t592 * mrSges(3,3) - t600 * t585 + t589 * t625 + t623;
t484 = t607 * t492 + t612 * t501;
t573 = -t603 * t587 - t636;
t586 = -t600 * mrSges(3,2) + mrSges(3,3) * t625;
t483 = m(3) * t573 + t592 * mrSges(3,1) + t591 * mrSges(3,2) + (t585 * t608 - t586 * t613) * t629 + t484;
t493 = t611 * t495 + t606 * t496;
t617 = -m(4) * t542 + t559 * mrSges(4,1) - t560 * mrSges(4,2) + t579 * t568 - t580 * t569 - t493;
t489 = m(3) * t561 + t599 * mrSges(3,1) - t591 * mrSges(3,3) + t600 * t586 - t589 * t626 + t617;
t471 = t481 * t633 - t603 * t483 + t489 * t632;
t469 = m(2) * t596 + qJDD(1) * mrSges(2,1) - t615 * mrSges(2,2) + t471;
t476 = t613 * t481 - t608 * t489;
t475 = m(2) * t597 - t615 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t476;
t631 = t614 * t469 + t609 * t475;
t470 = t481 * t635 + t604 * t483 + t489 * t634;
t624 = -t609 * t469 + t614 * t475;
t525 = Ifges(6,5) * t546 + Ifges(6,6) * t545 + Ifges(6,3) * t576;
t527 = Ifges(6,1) * t546 + Ifges(6,4) * t545 + Ifges(6,5) * t576;
t498 = -mrSges(6,1) * t513 + mrSges(6,3) * t508 + Ifges(6,4) * t516 + Ifges(6,2) * t515 + Ifges(6,6) * t552 - t546 * t525 + t576 * t527;
t526 = Ifges(6,4) * t546 + Ifges(6,2) * t545 + Ifges(6,6) * t576;
t499 = mrSges(6,2) * t513 - mrSges(6,3) * t507 + Ifges(6,1) * t516 + Ifges(6,4) * t515 + Ifges(6,5) * t552 + t545 * t525 - t576 * t526;
t535 = Ifges(5,5) * t567 + Ifges(5,6) * t566 + Ifges(5,3) * t578;
t537 = Ifges(5,1) * t567 + Ifges(5,4) * t566 + Ifges(5,5) * t578;
t485 = -mrSges(5,1) * t518 + mrSges(5,3) * t512 + Ifges(5,4) * t532 + Ifges(5,2) * t531 + Ifges(5,6) * t557 - pkin(4) * t618 + pkin(10) * t621 + t610 * t498 + t605 * t499 - t567 * t535 + t578 * t537;
t536 = Ifges(5,4) * t567 + Ifges(5,2) * t566 + Ifges(5,6) * t578;
t486 = mrSges(5,2) * t518 - mrSges(5,3) * t511 + Ifges(5,1) * t532 + Ifges(5,4) * t531 + Ifges(5,5) * t557 - pkin(10) * t497 - t605 * t498 + t610 * t499 + t566 * t535 - t578 * t536;
t553 = Ifges(4,5) * t580 + Ifges(4,6) * t579 + Ifges(4,3) * t595;
t554 = Ifges(4,4) * t580 + Ifges(4,2) * t579 + Ifges(4,6) * t595;
t472 = mrSges(4,2) * t542 - mrSges(4,3) * t523 + Ifges(4,1) * t560 + Ifges(4,4) * t559 + Ifges(4,5) * t584 - pkin(9) * t493 - t606 * t485 + t611 * t486 + t579 * t553 - t595 * t554;
t555 = Ifges(4,1) * t580 + Ifges(4,4) * t579 + Ifges(4,5) * t595;
t477 = Ifges(4,4) * t560 + Ifges(4,2) * t559 + Ifges(4,6) * t584 - t580 * t553 + t595 * t555 - mrSges(4,1) * t542 + mrSges(4,3) * t524 - Ifges(5,5) * t532 - Ifges(5,6) * t531 - Ifges(5,3) * t557 - t567 * t536 + t566 * t537 - mrSges(5,1) * t511 + mrSges(5,2) * t512 - Ifges(6,5) * t516 - Ifges(6,6) * t515 - Ifges(6,3) * t552 - t546 * t526 + t545 * t527 - mrSges(6,1) * t507 + mrSges(6,2) * t508 - pkin(4) * t497 - pkin(3) * t493;
t570 = Ifges(3,3) * t600 + (Ifges(3,5) * t608 + Ifges(3,6) * t613) * t629;
t571 = Ifges(3,6) * t600 + (Ifges(3,4) * t608 + Ifges(3,2) * t613) * t629;
t466 = mrSges(3,2) * t573 - mrSges(3,3) * t561 + Ifges(3,1) * t591 - Ifges(3,4) * t592 + Ifges(3,5) * t599 - pkin(8) * t484 + t612 * t472 - t607 * t477 + t570 * t625 - t600 * t571;
t572 = Ifges(3,5) * t600 + (Ifges(3,1) * t608 + Ifges(3,4) * t613) * t629;
t467 = Ifges(3,4) * t591 - Ifges(3,2) * t592 + Ifges(3,6) * t599 - t570 * t626 + t600 * t572 - mrSges(3,1) * t573 + mrSges(3,3) * t562 - Ifges(4,5) * t560 - Ifges(4,6) * t559 - Ifges(4,3) * t584 - t580 * t554 + t579 * t555 - mrSges(4,1) * t523 + mrSges(4,2) * t524 - t606 * t486 - t611 * t485 - pkin(3) * t616 - pkin(9) * t622 - pkin(2) * t484;
t619 = pkin(7) * t476 + t466 * t608 + t467 * t613;
t465 = Ifges(3,5) * t591 - Ifges(3,6) * t592 + Ifges(3,3) * t599 + mrSges(3,1) * t561 - mrSges(3,2) * t562 + t607 * t472 + t612 * t477 + pkin(2) * t617 + pkin(8) * t623 + (t571 * t608 - t572 * t613) * t629;
t464 = -mrSges(2,2) * g(3) - mrSges(2,3) * t596 + Ifges(2,5) * qJDD(1) - t615 * Ifges(2,6) + t613 * t466 - t608 * t467 + (-t470 * t603 - t471 * t604) * pkin(7);
t463 = mrSges(2,1) * g(3) + mrSges(2,3) * t597 + t615 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t470 - t603 * t465 + t619 * t604;
t1 = [-m(1) * g(1) + t624; -m(1) * g(2) + t631; (-m(1) - m(2)) * g(3) + t470; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t631 - t609 * t463 + t614 * t464; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t624 + t614 * t463 + t609 * t464; -mrSges(1,1) * g(2) + mrSges(2,1) * t596 + mrSges(1,2) * g(1) - mrSges(2,2) * t597 + Ifges(2,3) * qJDD(1) + pkin(1) * t471 + t604 * t465 + t619 * t603;];
tauB = t1;
