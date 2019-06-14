% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRRRPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
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
% Datum: 2019-05-05 08:30
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRRRPR6_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR6_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR6_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR6_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR6_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR6_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR6_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR6_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR6_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 08:22:43
% EndTime: 2019-05-05 08:22:56
% DurationCPUTime: 7.88s
% Computational Cost: add. (122337->320), mult. (233320->395), div. (0->0), fcn. (158657->12), ass. (0->136)
t727 = Ifges(5,1) + Ifges(6,1);
t720 = Ifges(5,4) - Ifges(6,5);
t719 = Ifges(5,5) + Ifges(6,4);
t726 = Ifges(5,2) + Ifges(6,3);
t718 = Ifges(5,6) - Ifges(6,6);
t725 = -Ifges(5,3) - Ifges(6,2);
t683 = sin(qJ(4));
t684 = sin(qJ(3));
t708 = qJD(2) * t684;
t722 = cos(qJ(4));
t659 = -t722 * qJD(3) + t683 * t708;
t687 = cos(qJ(3));
t706 = qJD(2) * qJD(3);
t704 = t687 * t706;
t663 = t684 * qJDD(2) + t704;
t628 = -t659 * qJD(4) + t683 * qJDD(3) + t722 * t663;
t678 = sin(pkin(11));
t680 = cos(pkin(11));
t665 = t678 * g(1) - t680 * g(2);
t666 = -t680 * g(1) - t678 * g(2);
t677 = -g(3) + qJDD(1);
t688 = cos(qJ(2));
t681 = cos(pkin(6));
t685 = sin(qJ(2));
t714 = t681 * t685;
t679 = sin(pkin(6));
t715 = t679 * t685;
t616 = t665 * t714 + t688 * t666 + t677 * t715;
t690 = qJD(2) ^ 2;
t610 = -t690 * pkin(2) + qJDD(2) * pkin(8) + t616;
t643 = -t679 * t665 + t681 * t677;
t604 = -t684 * t610 + t687 * t643;
t662 = (-pkin(3) * t687 - pkin(9) * t684) * qJD(2);
t689 = qJD(3) ^ 2;
t694 = qJDD(3) * pkin(3) + t689 * pkin(9) - t662 * t708 + t604;
t707 = t687 * qJD(2);
t673 = qJD(4) - t707;
t716 = t659 * t673;
t724 = (-t628 + t716) * qJ(5) - t694;
t615 = -t685 * t666 + (t665 * t681 + t677 * t679) * t688;
t723 = 2 * qJD(5);
t721 = -mrSges(5,3) - mrSges(6,2);
t605 = t687 * t610 + t684 * t643;
t596 = -t689 * pkin(3) + qJDD(3) * pkin(9) + t662 * t707 + t605;
t609 = -qJDD(2) * pkin(2) - t690 * pkin(8) - t615;
t705 = t684 * t706;
t664 = t687 * qJDD(2) - t705;
t598 = (-t663 - t704) * pkin(9) + (-t664 + t705) * pkin(3) + t609;
t588 = t722 * t596 + t683 * t598;
t660 = t683 * qJD(3) + t722 * t708;
t627 = t660 * qJD(4) - t722 * qJDD(3) + t683 * t663;
t639 = t673 * mrSges(5,1) - t660 * mrSges(5,3);
t656 = qJDD(4) - t664;
t633 = t659 * pkin(4) - t660 * qJ(5);
t672 = t673 ^ 2;
t584 = -t672 * pkin(4) + t656 * qJ(5) - t659 * t633 + t673 * t723 + t588;
t640 = -t673 * mrSges(6,1) + t660 * mrSges(6,2);
t587 = -t683 * t596 + t722 * t598;
t585 = -t656 * pkin(4) - t672 * qJ(5) + t660 * t633 + qJDD(5) - t587;
t579 = (-t628 - t716) * pkin(10) + (t659 * t660 - t656) * pkin(5) + t585;
t642 = -t673 * pkin(5) - t660 * pkin(10);
t655 = t659 ^ 2;
t580 = -t655 * pkin(5) + t627 * pkin(10) + t673 * t642 + t584;
t682 = sin(qJ(6));
t686 = cos(qJ(6));
t577 = t686 * t579 - t682 * t580;
t629 = t686 * t659 - t682 * t660;
t592 = t629 * qJD(6) + t682 * t627 + t686 * t628;
t630 = t682 * t659 + t686 * t660;
t606 = -t629 * mrSges(7,1) + t630 * mrSges(7,2);
t671 = qJD(6) - t673;
t611 = -t671 * mrSges(7,2) + t629 * mrSges(7,3);
t652 = qJDD(6) - t656;
t575 = m(7) * t577 + t652 * mrSges(7,1) - t592 * mrSges(7,3) - t630 * t606 + t671 * t611;
t578 = t682 * t579 + t686 * t580;
t591 = -t630 * qJD(6) + t686 * t627 - t682 * t628;
t612 = t671 * mrSges(7,1) - t630 * mrSges(7,3);
t576 = m(7) * t578 - t652 * mrSges(7,2) + t591 * mrSges(7,3) + t629 * t606 - t671 * t612;
t700 = -t682 * t575 + t686 * t576;
t696 = m(6) * t584 + t656 * mrSges(6,3) + t673 * t640 + t700;
t634 = t659 * mrSges(6,1) - t660 * mrSges(6,3);
t709 = -t659 * mrSges(5,1) - t660 * mrSges(5,2) - t634;
t566 = m(5) * t588 - t656 * mrSges(5,2) + t721 * t627 - t673 * t639 + t709 * t659 + t696;
t638 = -t673 * mrSges(5,2) - t659 * mrSges(5,3);
t568 = t686 * t575 + t682 * t576;
t641 = -t659 * mrSges(6,2) + t673 * mrSges(6,3);
t693 = -m(6) * t585 + t656 * mrSges(6,1) + t673 * t641 - t568;
t567 = m(5) * t587 + t656 * mrSges(5,1) + t721 * t628 + t673 * t638 + t709 * t660 + t693;
t564 = t683 * t566 + t722 * t567;
t667 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t708;
t668 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t707;
t692 = -m(4) * t609 + t664 * mrSges(4,1) - t663 * mrSges(4,2) - t667 * t708 + t668 * t707 - t564;
t560 = m(3) * t615 + qJDD(2) * mrSges(3,1) - t690 * mrSges(3,2) + t692;
t717 = t560 * t688;
t661 = (-mrSges(4,1) * t687 + mrSges(4,2) * t684) * qJD(2);
t701 = t722 * t566 - t683 * t567;
t563 = m(4) * t605 - qJDD(3) * mrSges(4,2) + t664 * mrSges(4,3) - qJD(3) * t667 + t661 * t707 + t701;
t586 = -0.2e1 * qJD(5) * t660 + (t660 * t673 + t627) * pkin(4) + t724;
t582 = -t655 * pkin(10) + (-pkin(4) - pkin(5)) * t627 + (-pkin(4) * t673 + t642 + t723) * t660 - t724;
t698 = -m(7) * t582 + t591 * mrSges(7,1) - t592 * mrSges(7,2) + t629 * t611 - t630 * t612;
t573 = m(6) * t586 + t627 * mrSges(6,1) - t628 * mrSges(6,3) - t660 * t640 + t659 * t641 + t698;
t691 = m(5) * t694 - t627 * mrSges(5,1) - t628 * mrSges(5,2) - t659 * t638 - t660 * t639 - t573;
t572 = m(4) * t604 + qJDD(3) * mrSges(4,1) - t663 * mrSges(4,3) + qJD(3) * t668 - t661 * t708 + t691;
t702 = t687 * t563 - t684 * t572;
t554 = m(3) * t616 - t690 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t702;
t557 = t684 * t563 + t687 * t572;
t556 = m(3) * t643 + t557;
t543 = t554 * t714 - t679 * t556 + t681 * t717;
t541 = m(2) * t665 + t543;
t548 = t688 * t554 - t685 * t560;
t547 = m(2) * t666 + t548;
t713 = t680 * t541 + t678 * t547;
t712 = t659 * t726 - t660 * t720 - t673 * t718;
t711 = t659 * t718 - t660 * t719 + t673 * t725;
t710 = -t720 * t659 + t660 * t727 + t719 * t673;
t542 = t554 * t715 + t681 * t556 + t679 * t717;
t703 = -t678 * t541 + t680 * t547;
t599 = Ifges(7,5) * t630 + Ifges(7,6) * t629 + Ifges(7,3) * t671;
t601 = Ifges(7,1) * t630 + Ifges(7,4) * t629 + Ifges(7,5) * t671;
t569 = -mrSges(7,1) * t582 + mrSges(7,3) * t578 + Ifges(7,4) * t592 + Ifges(7,2) * t591 + Ifges(7,6) * t652 - t630 * t599 + t671 * t601;
t600 = Ifges(7,4) * t630 + Ifges(7,2) * t629 + Ifges(7,6) * t671;
t570 = mrSges(7,2) * t582 - mrSges(7,3) * t577 + Ifges(7,1) * t592 + Ifges(7,4) * t591 + Ifges(7,5) * t652 + t629 * t599 - t671 * t600;
t549 = mrSges(5,1) * t694 - mrSges(6,1) * t586 + mrSges(6,2) * t584 + mrSges(5,3) * t588 - pkin(4) * t573 - pkin(5) * t698 - pkin(10) * t700 - t686 * t569 - t682 * t570 - t627 * t726 + t720 * t628 + t718 * t656 + t711 * t660 + t710 * t673;
t550 = -mrSges(5,2) * t694 + mrSges(6,2) * t585 - mrSges(5,3) * t587 - mrSges(6,3) * t586 - pkin(10) * t568 - qJ(5) * t573 - t682 * t569 + t686 * t570 - t720 * t627 + t628 * t727 + t719 * t656 + t711 * t659 + t712 * t673;
t646 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t684 + Ifges(4,6) * t687) * qJD(2);
t647 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t684 + Ifges(4,2) * t687) * qJD(2);
t539 = mrSges(4,2) * t609 - mrSges(4,3) * t604 + Ifges(4,1) * t663 + Ifges(4,4) * t664 + Ifges(4,5) * qJDD(3) - pkin(9) * t564 - qJD(3) * t647 - t683 * t549 + t722 * t550 + t646 * t707;
t648 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t684 + Ifges(4,4) * t687) * qJD(2);
t544 = -pkin(3) * t564 + t725 * t656 + mrSges(7,1) * t577 - mrSges(7,2) * t578 + pkin(5) * t568 - qJ(5) * t696 + Ifges(4,4) * t663 + Ifges(4,2) * t664 + qJD(3) * t648 + Ifges(7,3) * t652 - t629 * t601 + t630 * t600 - mrSges(4,1) * t609 + Ifges(7,6) * t591 + Ifges(7,5) * t592 + mrSges(4,3) * t605 - mrSges(6,3) * t584 + mrSges(6,1) * t585 - mrSges(5,1) * t587 + mrSges(5,2) * t588 + Ifges(4,6) * qJDD(3) - pkin(4) * t693 - t646 * t708 + (qJ(5) * t634 - t710) * t659 + (pkin(4) * t634 + t712) * t660 + (qJ(5) * mrSges(6,2) + t718) * t627 + (pkin(4) * mrSges(6,2) - t719) * t628;
t537 = mrSges(3,2) * t643 - mrSges(3,3) * t615 + Ifges(3,5) * qJDD(2) - t690 * Ifges(3,6) - pkin(8) * t557 + t687 * t539 - t684 * t544;
t538 = Ifges(3,6) * qJDD(2) + t690 * Ifges(3,5) - mrSges(3,1) * t643 + mrSges(3,3) * t616 - Ifges(4,5) * t663 - Ifges(4,6) * t664 - Ifges(4,3) * qJDD(3) - mrSges(4,1) * t604 + mrSges(4,2) * t605 - t683 * t550 - t722 * t549 - pkin(3) * t691 - pkin(9) * t701 - pkin(2) * t557 + (-t684 * t647 + t687 * t648) * qJD(2);
t695 = pkin(7) * t548 + t537 * t685 + t538 * t688;
t536 = mrSges(3,1) * t615 - mrSges(3,2) * t616 + Ifges(3,3) * qJDD(2) + pkin(2) * t692 + pkin(8) * t702 + t684 * t539 + t687 * t544;
t535 = mrSges(2,2) * t677 - mrSges(2,3) * t665 + t688 * t537 - t685 * t538 + (-t542 * t679 - t543 * t681) * pkin(7);
t534 = -mrSges(2,1) * t677 + mrSges(2,3) * t666 - pkin(1) * t542 - t679 * t536 + t695 * t681;
t1 = [-m(1) * g(1) + t703; -m(1) * g(2) + t713; -m(1) * g(3) + m(2) * t677 + t542; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t713 - t678 * t534 + t680 * t535; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t703 + t680 * t534 + t678 * t535; -mrSges(1,1) * g(2) + mrSges(2,1) * t665 + mrSges(1,2) * g(1) - mrSges(2,2) * t666 + pkin(1) * t543 + t681 * t536 + t695 * t679;];
tauB  = t1;
