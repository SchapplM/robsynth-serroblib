% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPRPR7
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
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
% Datum: 2019-05-06 14:47
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPRPR7_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR7_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR7_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR7_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR7_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR7_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR7_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR7_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR7_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 14:43:45
% EndTime: 2019-05-06 14:43:56
% DurationCPUTime: 8.56s
% Computational Cost: add. (119442->364), mult. (264973->447), div. (0->0), fcn. (175227->10), ass. (0->142)
t794 = Ifges(3,1) + Ifges(4,1);
t788 = Ifges(3,4) - Ifges(4,5);
t787 = Ifges(3,5) + Ifges(4,4);
t793 = Ifges(3,2) + Ifges(4,3);
t786 = Ifges(3,6) - Ifges(4,6);
t792 = Ifges(3,3) + Ifges(4,2);
t791 = 2 * qJD(3);
t790 = 2 * qJD(5);
t789 = mrSges(3,3) + mrSges(4,2);
t759 = cos(qJ(2));
t762 = qJD(1) ^ 2;
t785 = t759 ^ 2 * t762;
t756 = sin(qJ(1));
t760 = cos(qJ(1));
t734 = -g(1) * t760 - g(2) * t756;
t712 = -pkin(1) * t762 + qJDD(1) * pkin(7) + t734;
t755 = sin(qJ(2));
t693 = -g(3) * t755 + t759 * t712;
t723 = (-mrSges(3,1) * t759 + mrSges(3,2) * t755) * qJD(1);
t778 = qJD(1) * qJD(2);
t776 = t755 * t778;
t725 = qJDD(1) * t759 - t776;
t780 = qJD(1) * t755;
t728 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t780;
t721 = (-pkin(2) * t759 - qJ(3) * t755) * qJD(1);
t761 = qJD(2) ^ 2;
t779 = qJD(1) * t759;
t671 = -pkin(2) * t761 + qJDD(2) * qJ(3) + qJD(2) * t791 + t721 * t779 + t693;
t722 = (-mrSges(4,1) * t759 - mrSges(4,3) * t755) * qJD(1);
t729 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t780;
t732 = -qJD(2) * pkin(3) - pkin(8) * t780;
t667 = -pkin(3) * t785 - pkin(8) * t725 + qJD(2) * t732 + t671;
t692 = -t759 * g(3) - t755 * t712;
t676 = -qJDD(2) * pkin(2) - t761 * qJ(3) + t721 * t780 + qJDD(3) - t692;
t777 = t759 * t778;
t724 = qJDD(1) * t755 + t777;
t668 = (-t724 + t777) * pkin(8) + (-t755 * t759 * t762 - qJDD(2)) * pkin(3) + t676;
t754 = sin(qJ(4));
t758 = cos(qJ(4));
t639 = -t754 * t667 + t758 * t668;
t709 = (-t754 * t755 - t758 * t759) * qJD(1);
t678 = qJD(4) * t709 + t724 * t758 - t725 * t754;
t710 = (-t754 * t759 + t755 * t758) * qJD(1);
t743 = -qJDD(2) + qJDD(4);
t744 = -qJD(2) + qJD(4);
t633 = (t709 * t744 - t678) * qJ(5) + (t709 * t710 + t743) * pkin(4) + t639;
t640 = t758 * t667 + t754 * t668;
t677 = -qJD(4) * t710 - t724 * t754 - t725 * t758;
t695 = pkin(4) * t744 - qJ(5) * t710;
t702 = t709 ^ 2;
t635 = -pkin(4) * t702 + qJ(5) * t677 - t695 * t744 + t640;
t750 = sin(pkin(10));
t751 = cos(pkin(10));
t689 = t709 * t751 - t710 * t750;
t630 = t750 * t633 + t751 * t635 + t689 * t790;
t650 = t677 * t751 - t678 * t750;
t690 = t709 * t750 + t710 * t751;
t665 = -mrSges(6,1) * t689 + mrSges(6,2) * t690;
t680 = mrSges(6,1) * t744 - mrSges(6,3) * t690;
t666 = -pkin(5) * t689 - pkin(9) * t690;
t742 = t744 ^ 2;
t628 = -pkin(5) * t742 + pkin(9) * t743 + t666 * t689 + t630;
t733 = t756 * g(1) - t760 * g(2);
t711 = -qJDD(1) * pkin(1) - t762 * pkin(7) - t733;
t768 = -t725 * pkin(2) + t711 + (-t724 - t777) * qJ(3);
t654 = -pkin(2) * t776 + t725 * pkin(3) - pkin(8) * t785 - t768 + (t732 + t791) * t780;
t637 = -t677 * pkin(4) - t702 * qJ(5) + t710 * t695 + qJDD(5) + t654;
t651 = t677 * t750 + t678 * t751;
t631 = t637 + (t690 * t744 - t650) * pkin(5) + (-t689 * t744 - t651) * pkin(9);
t753 = sin(qJ(6));
t757 = cos(qJ(6));
t625 = -t628 * t753 + t631 * t757;
t674 = -t690 * t753 + t744 * t757;
t642 = qJD(6) * t674 + t651 * t757 + t743 * t753;
t649 = qJDD(6) - t650;
t675 = t690 * t757 + t744 * t753;
t652 = -mrSges(7,1) * t674 + mrSges(7,2) * t675;
t682 = qJD(6) - t689;
t655 = -mrSges(7,2) * t682 + mrSges(7,3) * t674;
t623 = m(7) * t625 + mrSges(7,1) * t649 - mrSges(7,3) * t642 - t652 * t675 + t655 * t682;
t626 = t628 * t757 + t631 * t753;
t641 = -qJD(6) * t675 - t651 * t753 + t743 * t757;
t656 = mrSges(7,1) * t682 - mrSges(7,3) * t675;
t624 = m(7) * t626 - mrSges(7,2) * t649 + mrSges(7,3) * t641 + t652 * t674 - t656 * t682;
t771 = -t623 * t753 + t757 * t624;
t614 = m(6) * t630 - mrSges(6,2) * t743 + mrSges(6,3) * t650 + t665 * t689 - t680 * t744 + t771;
t770 = -t751 * t633 + t750 * t635;
t629 = -0.2e1 * qJD(5) * t690 - t770;
t679 = -mrSges(6,2) * t744 + mrSges(6,3) * t689;
t627 = -t743 * pkin(5) - t742 * pkin(9) + (t790 + t666) * t690 + t770;
t766 = -m(7) * t627 + t641 * mrSges(7,1) - mrSges(7,2) * t642 + t674 * t655 - t656 * t675;
t619 = m(6) * t629 + mrSges(6,1) * t743 - mrSges(6,3) * t651 - t665 * t690 + t679 * t744 + t766;
t608 = t750 * t614 + t751 * t619;
t691 = -mrSges(5,1) * t709 + mrSges(5,2) * t710;
t694 = -mrSges(5,2) * t744 + mrSges(5,3) * t709;
t606 = m(5) * t639 + mrSges(5,1) * t743 - mrSges(5,3) * t678 - t691 * t710 + t694 * t744 + t608;
t696 = mrSges(5,1) * t744 - mrSges(5,3) * t710;
t772 = t751 * t614 - t619 * t750;
t607 = m(5) * t640 - mrSges(5,2) * t743 + mrSges(5,3) * t677 + t691 * t709 - t696 * t744 + t772;
t773 = -t754 * t606 + t758 * t607;
t767 = m(4) * t671 + qJDD(2) * mrSges(4,3) + qJD(2) * t729 + t722 * t779 + t773;
t600 = m(3) * t693 - qJDD(2) * mrSges(3,2) - qJD(2) * t728 + t723 * t779 + t725 * t789 + t767;
t730 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t779;
t602 = t758 * t606 + t754 * t607;
t731 = mrSges(4,2) * t779 + qJD(2) * mrSges(4,3);
t765 = -m(4) * t676 + qJDD(2) * mrSges(4,1) + qJD(2) * t731 - t602;
t601 = m(3) * t692 + qJDD(2) * mrSges(3,1) + qJD(2) * t730 - t789 * t724 + (-t722 - t723) * t780 + t765;
t774 = t759 * t600 - t601 * t755;
t594 = m(2) * t734 - mrSges(2,1) * t762 - qJDD(1) * mrSges(2,2) + t774;
t669 = (pkin(2) * qJD(2) - (2 * qJD(3))) * t780 + t768;
t615 = t757 * t623 + t753 * t624;
t769 = m(6) * t637 - t650 * mrSges(6,1) + t651 * mrSges(6,2) - t689 * t679 + t690 * t680 + t615;
t764 = -m(5) * t654 + t677 * mrSges(5,1) - t678 * mrSges(5,2) + t709 * t694 - t710 * t696 - t769;
t611 = m(4) * t669 - t725 * mrSges(4,1) - t724 * mrSges(4,3) - t729 * t780 - t731 * t779 + t764;
t763 = -m(3) * t711 + t725 * mrSges(3,1) - t724 * mrSges(3,2) - t728 * t780 + t730 * t779 - t611;
t610 = m(2) * t733 + qJDD(1) * mrSges(2,1) - t762 * mrSges(2,2) + t763;
t784 = t756 * t594 + t760 * t610;
t595 = t755 * t600 + t759 * t601;
t783 = t792 * qJD(2) + (t787 * t755 + t786 * t759) * qJD(1);
t782 = -t786 * qJD(2) + (-t788 * t755 - t793 * t759) * qJD(1);
t781 = t787 * qJD(2) + (t794 * t755 + t788 * t759) * qJD(1);
t775 = t760 * t594 - t610 * t756;
t685 = Ifges(5,1) * t710 + Ifges(5,4) * t709 + Ifges(5,5) * t744;
t684 = Ifges(5,4) * t710 + Ifges(5,2) * t709 + Ifges(5,6) * t744;
t683 = Ifges(5,5) * t710 + Ifges(5,6) * t709 + Ifges(5,3) * t744;
t659 = Ifges(6,1) * t690 + Ifges(6,4) * t689 + Ifges(6,5) * t744;
t658 = Ifges(6,4) * t690 + Ifges(6,2) * t689 + Ifges(6,6) * t744;
t657 = Ifges(6,5) * t690 + Ifges(6,6) * t689 + Ifges(6,3) * t744;
t645 = Ifges(7,1) * t675 + Ifges(7,4) * t674 + Ifges(7,5) * t682;
t644 = Ifges(7,4) * t675 + Ifges(7,2) * t674 + Ifges(7,6) * t682;
t643 = Ifges(7,5) * t675 + Ifges(7,6) * t674 + Ifges(7,3) * t682;
t617 = mrSges(7,2) * t627 - mrSges(7,3) * t625 + Ifges(7,1) * t642 + Ifges(7,4) * t641 + Ifges(7,5) * t649 + t643 * t674 - t644 * t682;
t616 = -mrSges(7,1) * t627 + mrSges(7,3) * t626 + Ifges(7,4) * t642 + Ifges(7,2) * t641 + Ifges(7,6) * t649 - t643 * t675 + t645 * t682;
t604 = -mrSges(6,1) * t637 - mrSges(7,1) * t625 + mrSges(7,2) * t626 + mrSges(6,3) * t630 + Ifges(6,4) * t651 - Ifges(7,5) * t642 + Ifges(6,2) * t650 + Ifges(6,6) * t743 - Ifges(7,6) * t641 - Ifges(7,3) * t649 - pkin(5) * t615 - t644 * t675 + t645 * t674 - t657 * t690 + t659 * t744;
t603 = mrSges(6,2) * t637 - mrSges(6,3) * t629 + Ifges(6,1) * t651 + Ifges(6,4) * t650 + Ifges(6,5) * t743 - pkin(9) * t615 - t616 * t753 + t617 * t757 + t657 * t689 - t658 * t744;
t596 = mrSges(5,2) * t654 - mrSges(5,3) * t639 + Ifges(5,1) * t678 + Ifges(5,4) * t677 + Ifges(5,5) * t743 - qJ(5) * t608 + t603 * t751 - t604 * t750 + t683 * t709 - t684 * t744;
t591 = -mrSges(5,1) * t654 + mrSges(5,3) * t640 + Ifges(5,4) * t678 + Ifges(5,2) * t677 + Ifges(5,6) * t743 - pkin(4) * t769 + qJ(5) * t772 + t750 * t603 + t751 * t604 - t710 * t683 + t744 * t685;
t590 = mrSges(3,2) * t711 + mrSges(4,2) * t676 - mrSges(3,3) * t692 - mrSges(4,3) * t669 - pkin(8) * t602 - qJ(3) * t611 + t782 * qJD(2) + t787 * qJDD(2) - t754 * t591 + t758 * t596 + t794 * t724 + t788 * t725 + t783 * t779;
t589 = -mrSges(3,1) * t711 - mrSges(4,1) * t669 + mrSges(4,2) * t671 + mrSges(3,3) * t693 - pkin(2) * t611 - pkin(3) * t764 - pkin(8) * t773 + t781 * qJD(2) + t786 * qJDD(2) - t758 * t591 - t754 * t596 + t788 * t724 + t793 * t725 - t783 * t780;
t588 = -pkin(2) * t765 + pkin(5) * t766 - qJ(3) * t767 + (Ifges(5,3) + Ifges(6,3)) * t743 + (-mrSges(4,2) * qJ(3) - t786) * t725 + (mrSges(4,2) * pkin(2) - t787) * t724 + mrSges(2,1) * g(3) + pkin(9) * t771 - t792 * qJDD(2) + t757 * t616 + t762 * Ifges(2,5) + t753 * t617 + mrSges(2,3) * t734 - t709 * t685 + t710 * t684 - t689 * t659 + t690 * t658 - mrSges(3,1) * t692 + mrSges(3,2) * t693 + mrSges(4,1) * t676 + Ifges(5,6) * t677 + Ifges(5,5) * t678 - mrSges(4,3) * t671 + Ifges(6,6) * t650 + Ifges(6,5) * t651 - pkin(1) * t595 - mrSges(5,2) * t640 + mrSges(5,1) * t639 + mrSges(6,1) * t629 - mrSges(6,2) * t630 + (t781 * t759 + (pkin(2) * t722 + t782) * t755) * qJD(1) + pkin(4) * t608 + pkin(3) * t602 + Ifges(2,6) * qJDD(1);
t587 = -mrSges(2,2) * g(3) - mrSges(2,3) * t733 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t762 - pkin(7) * t595 - t589 * t755 + t590 * t759;
t1 = [-m(1) * g(1) + t775; -m(1) * g(2) + t784; (-m(1) - m(2)) * g(3) + t595; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t784 + t760 * t587 - t756 * t588; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t775 + t756 * t587 + t760 * t588; -mrSges(1,1) * g(2) + mrSges(2,1) * t733 + mrSges(1,2) * g(1) - mrSges(2,2) * t734 + Ifges(2,3) * qJDD(1) + pkin(1) * t763 + pkin(7) * t774 + t759 * t589 + t755 * t590;];
tauB  = t1;
