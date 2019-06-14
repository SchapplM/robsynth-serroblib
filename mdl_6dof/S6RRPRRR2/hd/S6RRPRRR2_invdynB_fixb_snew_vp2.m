% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPRRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-05-06 19:59
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPRRR2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR2_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR2_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR2_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR2_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR2_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR2_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR2_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 19:53:51
% EndTime: 2019-05-06 19:54:35
% DurationCPUTime: 27.76s
% Computational Cost: add. (412317->387), mult. (952706->490), div. (0->0), fcn. (720896->12), ass. (0->149)
t740 = qJD(1) ^ 2;
t758 = pkin(2) * t740;
t734 = sin(qJ(1));
t739 = cos(qJ(1));
t722 = -g(1) * t739 - g(2) * t734;
t711 = -pkin(1) * t740 + qJDD(1) * pkin(7) + t722;
t733 = sin(qJ(2));
t757 = t733 * t711;
t738 = cos(qJ(2));
t753 = qJD(1) * qJD(2);
t716 = qJDD(1) * t733 + t738 * t753;
t672 = qJDD(2) * pkin(2) - t716 * qJ(3) - t757 + (qJ(3) * t753 + t733 * t758 - g(3)) * t738;
t696 = -g(3) * t733 + t738 * t711;
t717 = qJDD(1) * t738 - t733 * t753;
t755 = qJD(1) * t733;
t718 = qJD(2) * pkin(2) - qJ(3) * t755;
t727 = t738 ^ 2;
t673 = qJ(3) * t717 - qJD(2) * t718 - t727 * t758 + t696;
t728 = sin(pkin(11));
t729 = cos(pkin(11));
t706 = (t728 * t738 + t729 * t733) * qJD(1);
t647 = -0.2e1 * qJD(3) * t706 + t729 * t672 - t728 * t673;
t694 = t716 * t729 + t717 * t728;
t705 = (-t728 * t733 + t729 * t738) * qJD(1);
t631 = (qJD(2) * t705 - t694) * pkin(8) + (t705 * t706 + qJDD(2)) * pkin(3) + t647;
t648 = 0.2e1 * qJD(3) * t705 + t728 * t672 + t729 * t673;
t693 = -t716 * t728 + t717 * t729;
t699 = qJD(2) * pkin(3) - pkin(8) * t706;
t704 = t705 ^ 2;
t635 = -pkin(3) * t704 + pkin(8) * t693 - qJD(2) * t699 + t648;
t732 = sin(qJ(4));
t737 = cos(qJ(4));
t625 = t732 * t631 + t737 * t635;
t687 = t705 * t732 + t706 * t737;
t654 = -t687 * qJD(4) + t693 * t737 - t732 * t694;
t686 = t705 * t737 - t732 * t706;
t667 = -mrSges(5,1) * t686 + mrSges(5,2) * t687;
t726 = qJD(2) + qJD(4);
t679 = mrSges(5,1) * t726 - mrSges(5,3) * t687;
t725 = qJDD(2) + qJDD(4);
t668 = -pkin(4) * t686 - pkin(9) * t687;
t724 = t726 ^ 2;
t617 = -pkin(4) * t724 + pkin(9) * t725 + t668 * t686 + t625;
t721 = t734 * g(1) - t739 * g(2);
t746 = -qJDD(1) * pkin(1) - t721;
t675 = -t717 * pkin(2) + qJDD(3) + t718 * t755 + (-qJ(3) * t727 - pkin(7)) * t740 + t746;
t646 = -t693 * pkin(3) - t704 * pkin(8) + t706 * t699 + t675;
t655 = qJD(4) * t686 + t693 * t732 + t694 * t737;
t623 = (-t686 * t726 - t655) * pkin(9) + (t687 * t726 - t654) * pkin(4) + t646;
t731 = sin(qJ(5));
t736 = cos(qJ(5));
t612 = -t731 * t617 + t736 * t623;
t676 = -t687 * t731 + t726 * t736;
t638 = qJD(5) * t676 + t655 * t736 + t725 * t731;
t653 = qJDD(5) - t654;
t677 = t687 * t736 + t726 * t731;
t682 = qJD(5) - t686;
t610 = (t676 * t682 - t638) * pkin(10) + (t676 * t677 + t653) * pkin(5) + t612;
t613 = t736 * t617 + t731 * t623;
t637 = -qJD(5) * t677 - t655 * t731 + t725 * t736;
t661 = pkin(5) * t682 - pkin(10) * t677;
t674 = t676 ^ 2;
t611 = -pkin(5) * t674 + pkin(10) * t637 - t661 * t682 + t613;
t730 = sin(qJ(6));
t735 = cos(qJ(6));
t608 = t610 * t735 - t611 * t730;
t656 = t676 * t735 - t677 * t730;
t620 = qJD(6) * t656 + t637 * t730 + t638 * t735;
t657 = t676 * t730 + t677 * t735;
t632 = -mrSges(7,1) * t656 + mrSges(7,2) * t657;
t680 = qJD(6) + t682;
t639 = -mrSges(7,2) * t680 + mrSges(7,3) * t656;
t650 = qJDD(6) + t653;
t606 = m(7) * t608 + mrSges(7,1) * t650 - mrSges(7,3) * t620 - t632 * t657 + t639 * t680;
t609 = t610 * t730 + t611 * t735;
t619 = -qJD(6) * t657 + t637 * t735 - t638 * t730;
t640 = mrSges(7,1) * t680 - mrSges(7,3) * t657;
t607 = m(7) * t609 - mrSges(7,2) * t650 + mrSges(7,3) * t619 + t632 * t656 - t640 * t680;
t598 = t735 * t606 + t730 * t607;
t658 = -mrSges(6,1) * t676 + mrSges(6,2) * t677;
t659 = -mrSges(6,2) * t682 + mrSges(6,3) * t676;
t596 = m(6) * t612 + mrSges(6,1) * t653 - mrSges(6,3) * t638 - t658 * t677 + t659 * t682 + t598;
t660 = mrSges(6,1) * t682 - mrSges(6,3) * t677;
t747 = -t606 * t730 + t735 * t607;
t597 = m(6) * t613 - mrSges(6,2) * t653 + mrSges(6,3) * t637 + t658 * t676 - t660 * t682 + t747;
t748 = -t596 * t731 + t736 * t597;
t591 = m(5) * t625 - mrSges(5,2) * t725 + mrSges(5,3) * t654 + t667 * t686 - t679 * t726 + t748;
t624 = t631 * t737 - t732 * t635;
t678 = -mrSges(5,2) * t726 + mrSges(5,3) * t686;
t616 = -pkin(4) * t725 - pkin(9) * t724 + t687 * t668 - t624;
t614 = -pkin(5) * t637 - pkin(10) * t674 + t661 * t677 + t616;
t744 = m(7) * t614 - t619 * mrSges(7,1) + mrSges(7,2) * t620 - t656 * t639 + t640 * t657;
t742 = -m(6) * t616 + t637 * mrSges(6,1) - mrSges(6,2) * t638 + t676 * t659 - t660 * t677 - t744;
t602 = m(5) * t624 + mrSges(5,1) * t725 - mrSges(5,3) * t655 - t667 * t687 + t678 * t726 + t742;
t584 = t732 * t591 + t737 * t602;
t690 = -mrSges(4,1) * t705 + mrSges(4,2) * t706;
t697 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t705;
t582 = m(4) * t647 + qJDD(2) * mrSges(4,1) - mrSges(4,3) * t694 + qJD(2) * t697 - t690 * t706 + t584;
t698 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t706;
t749 = t737 * t591 - t602 * t732;
t583 = m(4) * t648 - qJDD(2) * mrSges(4,2) + mrSges(4,3) * t693 - qJD(2) * t698 + t690 * t705 + t749;
t577 = t729 * t582 + t728 * t583;
t695 = -t738 * g(3) - t757;
t715 = (-mrSges(3,1) * t738 + mrSges(3,2) * t733) * qJD(1);
t754 = qJD(1) * t738;
t720 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t754;
t575 = m(3) * t695 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t716 + qJD(2) * t720 - t715 * t755 + t577;
t719 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t755;
t750 = -t582 * t728 + t729 * t583;
t576 = m(3) * t696 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t717 - qJD(2) * t719 + t715 * t754 + t750;
t751 = -t575 * t733 + t738 * t576;
t569 = m(2) * t722 - mrSges(2,1) * t740 - qJDD(1) * mrSges(2,2) + t751;
t710 = -t740 * pkin(7) + t746;
t592 = t736 * t596 + t731 * t597;
t745 = m(5) * t646 - t654 * mrSges(5,1) + t655 * mrSges(5,2) - t686 * t678 + t687 * t679 + t592;
t743 = m(4) * t675 - t693 * mrSges(4,1) + mrSges(4,2) * t694 - t705 * t697 + t698 * t706 + t745;
t741 = -m(3) * t710 + t717 * mrSges(3,1) - mrSges(3,2) * t716 - t719 * t755 + t720 * t754 - t743;
t588 = m(2) * t721 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t740 + t741;
t756 = t734 * t569 + t739 * t588;
t570 = t738 * t575 + t733 * t576;
t752 = t739 * t569 - t588 * t734;
t709 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t733 + Ifges(3,4) * t738) * qJD(1);
t708 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t733 + Ifges(3,2) * t738) * qJD(1);
t707 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t733 + Ifges(3,6) * t738) * qJD(1);
t685 = Ifges(4,1) * t706 + Ifges(4,4) * t705 + Ifges(4,5) * qJD(2);
t684 = Ifges(4,4) * t706 + Ifges(4,2) * t705 + Ifges(4,6) * qJD(2);
t683 = Ifges(4,5) * t706 + Ifges(4,6) * t705 + Ifges(4,3) * qJD(2);
t664 = Ifges(5,1) * t687 + Ifges(5,4) * t686 + Ifges(5,5) * t726;
t663 = Ifges(5,4) * t687 + Ifges(5,2) * t686 + Ifges(5,6) * t726;
t662 = Ifges(5,5) * t687 + Ifges(5,6) * t686 + Ifges(5,3) * t726;
t643 = Ifges(6,1) * t677 + Ifges(6,4) * t676 + Ifges(6,5) * t682;
t642 = Ifges(6,4) * t677 + Ifges(6,2) * t676 + Ifges(6,6) * t682;
t641 = Ifges(6,5) * t677 + Ifges(6,6) * t676 + Ifges(6,3) * t682;
t628 = Ifges(7,1) * t657 + Ifges(7,4) * t656 + Ifges(7,5) * t680;
t627 = Ifges(7,4) * t657 + Ifges(7,2) * t656 + Ifges(7,6) * t680;
t626 = Ifges(7,5) * t657 + Ifges(7,6) * t656 + Ifges(7,3) * t680;
t600 = mrSges(7,2) * t614 - mrSges(7,3) * t608 + Ifges(7,1) * t620 + Ifges(7,4) * t619 + Ifges(7,5) * t650 + t626 * t656 - t627 * t680;
t599 = -mrSges(7,1) * t614 + mrSges(7,3) * t609 + Ifges(7,4) * t620 + Ifges(7,2) * t619 + Ifges(7,6) * t650 - t626 * t657 + t628 * t680;
t586 = mrSges(6,2) * t616 - mrSges(6,3) * t612 + Ifges(6,1) * t638 + Ifges(6,4) * t637 + Ifges(6,5) * t653 - pkin(10) * t598 - t599 * t730 + t600 * t735 + t641 * t676 - t642 * t682;
t585 = -mrSges(6,1) * t616 + mrSges(6,3) * t613 + Ifges(6,4) * t638 + Ifges(6,2) * t637 + Ifges(6,6) * t653 - pkin(5) * t744 + pkin(10) * t747 + t735 * t599 + t730 * t600 - t677 * t641 + t682 * t643;
t578 = Ifges(5,4) * t655 + Ifges(5,2) * t654 + Ifges(5,6) * t725 - t687 * t662 + t726 * t664 - mrSges(5,1) * t646 + mrSges(5,3) * t625 - Ifges(6,5) * t638 - Ifges(6,6) * t637 - Ifges(6,3) * t653 - t677 * t642 + t676 * t643 - mrSges(6,1) * t612 + mrSges(6,2) * t613 - Ifges(7,5) * t620 - Ifges(7,6) * t619 - Ifges(7,3) * t650 - t657 * t627 + t656 * t628 - mrSges(7,1) * t608 + mrSges(7,2) * t609 - pkin(5) * t598 - pkin(4) * t592;
t571 = mrSges(5,2) * t646 - mrSges(5,3) * t624 + Ifges(5,1) * t655 + Ifges(5,4) * t654 + Ifges(5,5) * t725 - pkin(9) * t592 - t585 * t731 + t586 * t736 + t662 * t686 - t663 * t726;
t566 = mrSges(4,2) * t675 - mrSges(4,3) * t647 + Ifges(4,1) * t694 + Ifges(4,4) * t693 + Ifges(4,5) * qJDD(2) - pkin(8) * t584 - qJD(2) * t684 + t571 * t737 - t578 * t732 + t683 * t705;
t565 = -mrSges(4,1) * t675 + mrSges(4,3) * t648 + Ifges(4,4) * t694 + Ifges(4,2) * t693 + Ifges(4,6) * qJDD(2) - pkin(3) * t745 + pkin(8) * t749 + qJD(2) * t685 + t732 * t571 + t737 * t578 - t706 * t683;
t564 = mrSges(2,1) * g(3) + Ifges(2,6) * qJDD(1) - mrSges(5,1) * t624 + mrSges(5,2) * t625 + t705 * t685 - t706 * t684 - Ifges(3,5) * t716 - Ifges(3,6) * t717 - pkin(1) * t570 - t736 * t585 - pkin(4) * t742 + mrSges(2,3) * t722 - Ifges(5,3) * t725 - Ifges(5,6) * t654 - Ifges(5,5) * t655 + t686 * t664 - t687 * t663 - Ifges(4,6) * t693 - pkin(9) * t748 - pkin(3) * t584 - mrSges(4,1) * t647 + mrSges(4,2) * t648 + (-t708 * t733 + t709 * t738) * qJD(1) - t731 * t586 + t740 * Ifges(2,5) - pkin(2) * t577 + (-Ifges(3,3) - Ifges(4,3)) * qJDD(2) - Ifges(4,5) * t694 - mrSges(3,1) * t695 + mrSges(3,2) * t696;
t563 = mrSges(3,2) * t710 - mrSges(3,3) * t695 + Ifges(3,1) * t716 + Ifges(3,4) * t717 + Ifges(3,5) * qJDD(2) - qJ(3) * t577 - qJD(2) * t708 - t565 * t728 + t566 * t729 + t707 * t754;
t562 = -mrSges(3,1) * t710 + mrSges(3,3) * t696 + Ifges(3,4) * t716 + Ifges(3,2) * t717 + Ifges(3,6) * qJDD(2) - pkin(2) * t743 + qJ(3) * t750 + qJD(2) * t709 + t729 * t565 + t728 * t566 - t707 * t755;
t561 = -mrSges(2,2) * g(3) - mrSges(2,3) * t721 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t740 - pkin(7) * t570 - t562 * t733 + t563 * t738;
t1 = [-m(1) * g(1) + t752; -m(1) * g(2) + t756; (-m(1) - m(2)) * g(3) + t570; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t756 + t739 * t561 - t734 * t564; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t752 + t734 * t561 + t739 * t564; -mrSges(1,1) * g(2) + mrSges(2,1) * t721 + mrSges(1,2) * g(1) - mrSges(2,2) * t722 + Ifges(2,3) * qJDD(1) + pkin(1) * t741 + pkin(7) * t751 + t738 * t562 + t733 * t563;];
tauB  = t1;
