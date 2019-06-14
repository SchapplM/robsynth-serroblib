% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRRRR2
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-05-06 02:47
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRRRR2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR2_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR2_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR2_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR2_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR2_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR2_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR2_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 02:44:32
% EndTime: 2019-05-06 02:44:45
% DurationCPUTime: 12.90s
% Computational Cost: add. (225071->344), mult. (439675->432), div. (0->0), fcn. (303909->12), ass. (0->136)
t664 = sin(qJ(1));
t669 = cos(qJ(1));
t646 = t664 * g(1) - g(2) * t669;
t638 = qJDD(1) * pkin(1) + t646;
t647 = -g(1) * t669 - g(2) * t664;
t670 = qJD(1) ^ 2;
t640 = -pkin(1) * t670 + t647;
t658 = sin(pkin(11));
t659 = cos(pkin(11));
t617 = t658 * t638 + t659 * t640;
t613 = -pkin(2) * t670 + qJDD(1) * pkin(7) + t617;
t657 = -g(3) + qJDD(2);
t663 = sin(qJ(3));
t668 = cos(qJ(3));
t600 = -t663 * t613 + t668 * t657;
t684 = qJD(1) * qJD(3);
t682 = t668 * t684;
t641 = qJDD(1) * t663 + t682;
t586 = (-t641 + t682) * pkin(8) + (t663 * t668 * t670 + qJDD(3)) * pkin(3) + t600;
t601 = t668 * t613 + t663 * t657;
t642 = qJDD(1) * t668 - t663 * t684;
t686 = qJD(1) * t663;
t645 = qJD(3) * pkin(3) - pkin(8) * t686;
t656 = t668 ^ 2;
t587 = -pkin(3) * t656 * t670 + pkin(8) * t642 - qJD(3) * t645 + t601;
t662 = sin(qJ(4));
t667 = cos(qJ(4));
t572 = t662 * t586 + t667 * t587;
t633 = (t662 * t668 + t663 * t667) * qJD(1);
t602 = -t633 * qJD(4) - t662 * t641 + t642 * t667;
t685 = qJD(1) * t668;
t632 = -t662 * t686 + t667 * t685;
t614 = -mrSges(5,1) * t632 + mrSges(5,2) * t633;
t655 = qJD(3) + qJD(4);
t622 = mrSges(5,1) * t655 - mrSges(5,3) * t633;
t654 = qJDD(3) + qJDD(4);
t615 = -pkin(4) * t632 - pkin(9) * t633;
t653 = t655 ^ 2;
t567 = -pkin(4) * t653 + pkin(9) * t654 + t615 * t632 + t572;
t616 = t659 * t638 - t658 * t640;
t675 = -qJDD(1) * pkin(2) - t616;
t592 = -t642 * pkin(3) + t645 * t686 + (-pkin(8) * t656 - pkin(7)) * t670 + t675;
t603 = qJD(4) * t632 + t641 * t667 + t642 * t662;
t570 = (-t632 * t655 - t603) * pkin(9) + (t633 * t655 - t602) * pkin(4) + t592;
t661 = sin(qJ(5));
t666 = cos(qJ(5));
t559 = -t661 * t567 + t666 * t570;
t619 = -t633 * t661 + t655 * t666;
t580 = qJD(5) * t619 + t603 * t666 + t654 * t661;
t599 = qJDD(5) - t602;
t620 = t633 * t666 + t655 * t661;
t625 = qJD(5) - t632;
t557 = (t619 * t625 - t580) * pkin(10) + (t619 * t620 + t599) * pkin(5) + t559;
t560 = t666 * t567 + t661 * t570;
t579 = -qJD(5) * t620 - t603 * t661 + t654 * t666;
t606 = pkin(5) * t625 - pkin(10) * t620;
t618 = t619 ^ 2;
t558 = -pkin(5) * t618 + pkin(10) * t579 - t606 * t625 + t560;
t660 = sin(qJ(6));
t665 = cos(qJ(6));
t555 = t557 * t665 - t558 * t660;
t593 = t619 * t665 - t620 * t660;
t564 = qJD(6) * t593 + t579 * t660 + t580 * t665;
t594 = t619 * t660 + t620 * t665;
t577 = -mrSges(7,1) * t593 + mrSges(7,2) * t594;
t623 = qJD(6) + t625;
t584 = -mrSges(7,2) * t623 + mrSges(7,3) * t593;
t596 = qJDD(6) + t599;
t553 = m(7) * t555 + mrSges(7,1) * t596 - mrSges(7,3) * t564 - t577 * t594 + t584 * t623;
t556 = t557 * t660 + t558 * t665;
t563 = -qJD(6) * t594 + t579 * t665 - t580 * t660;
t585 = mrSges(7,1) * t623 - mrSges(7,3) * t594;
t554 = m(7) * t556 - mrSges(7,2) * t596 + mrSges(7,3) * t563 + t577 * t593 - t585 * t623;
t545 = t665 * t553 + t660 * t554;
t595 = -mrSges(6,1) * t619 + mrSges(6,2) * t620;
t604 = -mrSges(6,2) * t625 + mrSges(6,3) * t619;
t543 = m(6) * t559 + mrSges(6,1) * t599 - mrSges(6,3) * t580 - t595 * t620 + t604 * t625 + t545;
t605 = mrSges(6,1) * t625 - mrSges(6,3) * t620;
t676 = -t553 * t660 + t665 * t554;
t544 = m(6) * t560 - mrSges(6,2) * t599 + mrSges(6,3) * t579 + t595 * t619 - t605 * t625 + t676;
t677 = -t543 * t661 + t666 * t544;
t538 = m(5) * t572 - mrSges(5,2) * t654 + mrSges(5,3) * t602 + t614 * t632 - t622 * t655 + t677;
t571 = t586 * t667 - t662 * t587;
t621 = -mrSges(5,2) * t655 + mrSges(5,3) * t632;
t566 = -pkin(4) * t654 - pkin(9) * t653 + t633 * t615 - t571;
t561 = -pkin(5) * t579 - pkin(10) * t618 + t606 * t620 + t566;
t674 = m(7) * t561 - t563 * mrSges(7,1) + mrSges(7,2) * t564 - t593 * t584 + t585 * t594;
t672 = -m(6) * t566 + t579 * mrSges(6,1) - mrSges(6,2) * t580 + t619 * t604 - t605 * t620 - t674;
t549 = m(5) * t571 + mrSges(5,1) * t654 - mrSges(5,3) * t603 - t614 * t633 + t621 * t655 + t672;
t531 = t662 * t538 + t667 * t549;
t639 = (-mrSges(4,1) * t668 + mrSges(4,2) * t663) * qJD(1);
t644 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t685;
t529 = m(4) * t600 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t641 + qJD(3) * t644 - t639 * t686 + t531;
t643 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t686;
t678 = t667 * t538 - t549 * t662;
t530 = m(4) * t601 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t642 - qJD(3) * t643 + t639 * t685 + t678;
t679 = -t529 * t663 + t668 * t530;
t523 = m(3) * t617 - mrSges(3,1) * t670 - qJDD(1) * mrSges(3,2) + t679;
t612 = -t670 * pkin(7) + t675;
t539 = t666 * t543 + t661 * t544;
t673 = m(5) * t592 - t602 * mrSges(5,1) + mrSges(5,2) * t603 - t632 * t621 + t622 * t633 + t539;
t671 = -m(4) * t612 + t642 * mrSges(4,1) - mrSges(4,2) * t641 - t643 * t686 + t644 * t685 - t673;
t535 = m(3) * t616 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t670 + t671;
t519 = t658 * t523 + t659 * t535;
t517 = m(2) * t646 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t670 + t519;
t680 = t659 * t523 - t535 * t658;
t518 = m(2) * t647 - mrSges(2,1) * t670 - qJDD(1) * mrSges(2,2) + t680;
t687 = t669 * t517 + t664 * t518;
t524 = t668 * t529 + t663 * t530;
t683 = m(3) * t657 + t524;
t681 = -t517 * t664 + t669 * t518;
t631 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t663 + Ifges(4,4) * t668) * qJD(1);
t630 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t663 + Ifges(4,2) * t668) * qJD(1);
t629 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t663 + Ifges(4,6) * t668) * qJD(1);
t610 = Ifges(5,1) * t633 + Ifges(5,4) * t632 + Ifges(5,5) * t655;
t609 = Ifges(5,4) * t633 + Ifges(5,2) * t632 + Ifges(5,6) * t655;
t608 = Ifges(5,5) * t633 + Ifges(5,6) * t632 + Ifges(5,3) * t655;
t590 = Ifges(6,1) * t620 + Ifges(6,4) * t619 + Ifges(6,5) * t625;
t589 = Ifges(6,4) * t620 + Ifges(6,2) * t619 + Ifges(6,6) * t625;
t588 = Ifges(6,5) * t620 + Ifges(6,6) * t619 + Ifges(6,3) * t625;
t575 = Ifges(7,1) * t594 + Ifges(7,4) * t593 + Ifges(7,5) * t623;
t574 = Ifges(7,4) * t594 + Ifges(7,2) * t593 + Ifges(7,6) * t623;
t573 = Ifges(7,5) * t594 + Ifges(7,6) * t593 + Ifges(7,3) * t623;
t547 = mrSges(7,2) * t561 - mrSges(7,3) * t555 + Ifges(7,1) * t564 + Ifges(7,4) * t563 + Ifges(7,5) * t596 + t573 * t593 - t574 * t623;
t546 = -mrSges(7,1) * t561 + mrSges(7,3) * t556 + Ifges(7,4) * t564 + Ifges(7,2) * t563 + Ifges(7,6) * t596 - t573 * t594 + t575 * t623;
t533 = mrSges(6,2) * t566 - mrSges(6,3) * t559 + Ifges(6,1) * t580 + Ifges(6,4) * t579 + Ifges(6,5) * t599 - pkin(10) * t545 - t546 * t660 + t547 * t665 + t588 * t619 - t589 * t625;
t532 = -mrSges(6,1) * t566 + mrSges(6,3) * t560 + Ifges(6,4) * t580 + Ifges(6,2) * t579 + Ifges(6,6) * t599 - pkin(5) * t674 + pkin(10) * t676 + t665 * t546 + t660 * t547 - t620 * t588 + t625 * t590;
t525 = Ifges(5,4) * t603 + Ifges(5,2) * t602 + Ifges(5,6) * t654 - t633 * t608 + t655 * t610 - mrSges(5,1) * t592 + mrSges(5,3) * t572 - Ifges(6,5) * t580 - Ifges(6,6) * t579 - Ifges(6,3) * t599 - t620 * t589 + t619 * t590 - mrSges(6,1) * t559 + mrSges(6,2) * t560 - Ifges(7,5) * t564 - Ifges(7,6) * t563 - Ifges(7,3) * t596 - t594 * t574 + t593 * t575 - mrSges(7,1) * t555 + mrSges(7,2) * t556 - pkin(5) * t545 - pkin(4) * t539;
t520 = mrSges(5,2) * t592 - mrSges(5,3) * t571 + Ifges(5,1) * t603 + Ifges(5,4) * t602 + Ifges(5,5) * t654 - pkin(9) * t539 - t532 * t661 + t533 * t666 + t608 * t632 - t609 * t655;
t513 = mrSges(4,2) * t612 - mrSges(4,3) * t600 + Ifges(4,1) * t641 + Ifges(4,4) * t642 + Ifges(4,5) * qJDD(3) - pkin(8) * t531 - qJD(3) * t630 + t520 * t667 - t525 * t662 + t629 * t685;
t512 = -pkin(2) * t524 - mrSges(3,1) * t657 + mrSges(3,3) * t617 - pkin(3) * t531 + mrSges(4,2) * t601 - Ifges(4,5) * t641 - Ifges(4,6) * t642 - Ifges(4,3) * qJDD(3) - mrSges(4,1) * t600 - t661 * t533 - t666 * t532 - pkin(4) * t672 - pkin(9) * t677 - Ifges(5,5) * t603 - Ifges(5,6) * t602 - Ifges(5,3) * t654 - mrSges(5,1) * t571 + mrSges(5,2) * t572 + t670 * Ifges(3,5) - t633 * t609 + t632 * t610 + Ifges(3,6) * qJDD(1) + (-t630 * t663 + t631 * t668) * qJD(1);
t511 = -mrSges(4,1) * t612 + mrSges(4,3) * t601 + Ifges(4,4) * t641 + Ifges(4,2) * t642 + Ifges(4,6) * qJDD(3) - pkin(3) * t673 + pkin(8) * t678 + qJD(3) * t631 + t662 * t520 + t667 * t525 - t629 * t686;
t510 = mrSges(3,2) * t657 - mrSges(3,3) * t616 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t670 - pkin(7) * t524 - t511 * t663 + t513 * t668;
t509 = -mrSges(2,2) * g(3) - mrSges(2,3) * t646 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t670 - qJ(2) * t519 + t510 * t659 - t512 * t658;
t508 = mrSges(2,1) * g(3) + mrSges(2,3) * t647 + t670 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t683 + qJ(2) * t680 + t658 * t510 + t659 * t512;
t1 = [-m(1) * g(1) + t681; -m(1) * g(2) + t687; (-m(1) - m(2)) * g(3) + t683; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t687 - t664 * t508 + t669 * t509; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t681 + t669 * t508 + t664 * t509; pkin(1) * t519 + mrSges(2,1) * t646 - mrSges(2,2) * t647 + t668 * t511 + pkin(2) * t671 + pkin(7) * t679 + t663 * t513 + mrSges(3,1) * t616 - mrSges(3,2) * t617 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);];
tauB  = t1;
