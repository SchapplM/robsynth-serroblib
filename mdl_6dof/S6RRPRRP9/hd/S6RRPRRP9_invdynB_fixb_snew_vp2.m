% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPRRP9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-05-06 18:37
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPRRP9_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP9_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP9_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP9_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP9_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP9_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP9_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP9_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP9_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 18:29:01
% EndTime: 2019-05-06 18:29:27
% DurationCPUTime: 25.20s
% Computational Cost: add. (391505->376), mult. (887872->474), div. (0->0), fcn. (712182->12), ass. (0->153)
t788 = Ifges(6,1) + Ifges(7,1);
t782 = Ifges(6,4) + Ifges(7,4);
t781 = Ifges(6,5) + Ifges(7,5);
t787 = Ifges(6,2) + Ifges(7,2);
t786 = Ifges(6,6) + Ifges(7,6);
t785 = Ifges(6,3) + Ifges(7,3);
t744 = cos(pkin(6));
t784 = t744 * g(3);
t783 = -mrSges(6,2) - mrSges(7,2);
t742 = sin(pkin(6));
t747 = sin(qJ(2));
t779 = t742 * t747;
t751 = cos(qJ(2));
t778 = t742 * t751;
t777 = t744 * t747;
t776 = t744 * t751;
t748 = sin(qJ(1));
t752 = cos(qJ(1));
t733 = t748 * g(1) - g(2) * t752;
t753 = qJD(1) ^ 2;
t725 = pkin(8) * t742 * t753 + qJDD(1) * pkin(1) + t733;
t734 = -g(1) * t752 - g(2) * t748;
t767 = qJDD(1) * t742;
t726 = -pkin(1) * t753 + pkin(8) * t767 + t734;
t770 = t725 * t777 + t751 * t726;
t698 = -g(3) * t779 + t770;
t738 = qJD(1) * t744 + qJD(2);
t769 = qJD(1) * t742;
t764 = t747 * t769;
t723 = mrSges(3,1) * t738 - mrSges(3,3) * t764;
t728 = (-mrSges(3,1) * t751 + mrSges(3,2) * t747) * t769;
t730 = -qJD(2) * t764 + t751 * t767;
t737 = qJDD(1) * t744 + qJDD(2);
t727 = (-pkin(2) * t751 - qJ(3) * t747) * t769;
t736 = t738 ^ 2;
t768 = qJD(1) * t751;
t685 = -t736 * pkin(2) + t737 * qJ(3) + (-g(3) * t747 + t727 * t768) * t742 + t770;
t729 = (qJD(2) * t768 + qJDD(1) * t747) * t742;
t686 = -t730 * pkin(2) - t784 - t729 * qJ(3) + (-t725 + (pkin(2) * t747 - qJ(3) * t751) * t738 * qJD(1)) * t742;
t741 = sin(pkin(11));
t743 = cos(pkin(11));
t719 = t738 * t741 + t743 * t764;
t644 = -0.2e1 * qJD(3) * t719 - t741 * t685 + t743 * t686;
t707 = t729 * t743 + t737 * t741;
t718 = t738 * t743 - t741 * t764;
t763 = t742 * t768;
t640 = (-t718 * t763 - t707) * pkin(9) + (t718 * t719 - t730) * pkin(3) + t644;
t645 = 0.2e1 * qJD(3) * t718 + t743 * t685 + t741 * t686;
t706 = -t729 * t741 + t737 * t743;
t708 = -pkin(3) * t763 - pkin(9) * t719;
t717 = t718 ^ 2;
t643 = -pkin(3) * t717 + pkin(9) * t706 + t708 * t763 + t645;
t746 = sin(qJ(4));
t750 = cos(qJ(4));
t635 = t746 * t640 + t750 * t643;
t701 = t718 * t746 + t719 * t750;
t668 = -qJD(4) * t701 + t706 * t750 - t707 * t746;
t700 = t718 * t750 - t719 * t746;
t679 = -mrSges(5,1) * t700 + mrSges(5,2) * t701;
t732 = qJD(4) - t763;
t691 = mrSges(5,1) * t732 - mrSges(5,3) * t701;
t722 = qJDD(4) - t730;
t680 = -pkin(4) * t700 - pkin(10) * t701;
t731 = t732 ^ 2;
t633 = -pkin(4) * t731 + pkin(10) * t722 + t680 * t700 + t635;
t697 = -g(3) * t778 + t725 * t776 - t747 * t726;
t684 = -t737 * pkin(2) - t736 * qJ(3) + t727 * t764 + qJDD(3) - t697;
t651 = -t706 * pkin(3) - t717 * pkin(9) + t719 * t708 + t684;
t669 = qJD(4) * t700 + t706 * t746 + t707 * t750;
t638 = (-t700 * t732 - t669) * pkin(10) + (t701 * t732 - t668) * pkin(4) + t651;
t745 = sin(qJ(5));
t749 = cos(qJ(5));
t628 = -t745 * t633 + t749 * t638;
t688 = -t701 * t745 + t732 * t749;
t650 = qJD(5) * t688 + t669 * t749 + t722 * t745;
t689 = t701 * t749 + t732 * t745;
t663 = -mrSges(7,1) * t688 + mrSges(7,2) * t689;
t664 = -mrSges(6,1) * t688 + mrSges(6,2) * t689;
t667 = qJDD(5) - t668;
t699 = qJD(5) - t700;
t671 = -mrSges(6,2) * t699 + mrSges(6,3) * t688;
t625 = -0.2e1 * qJD(6) * t689 + (t688 * t699 - t650) * qJ(6) + (t688 * t689 + t667) * pkin(5) + t628;
t670 = -mrSges(7,2) * t699 + mrSges(7,3) * t688;
t766 = m(7) * t625 + t667 * mrSges(7,1) + t699 * t670;
t617 = m(6) * t628 + t667 * mrSges(6,1) + t699 * t671 + (-t663 - t664) * t689 + (-mrSges(6,3) - mrSges(7,3)) * t650 + t766;
t629 = t749 * t633 + t745 * t638;
t649 = -qJD(5) * t689 - t669 * t745 + t722 * t749;
t672 = pkin(5) * t699 - qJ(6) * t689;
t687 = t688 ^ 2;
t627 = -pkin(5) * t687 + qJ(6) * t649 + 0.2e1 * qJD(6) * t688 - t672 * t699 + t629;
t765 = m(7) * t627 + t649 * mrSges(7,3) + t688 * t663;
t673 = mrSges(7,1) * t699 - mrSges(7,3) * t689;
t771 = -mrSges(6,1) * t699 + mrSges(6,3) * t689 - t673;
t622 = m(6) * t629 + t649 * mrSges(6,3) + t688 * t664 + t667 * t783 + t699 * t771 + t765;
t759 = -t617 * t745 + t749 * t622;
t613 = m(5) * t635 - mrSges(5,2) * t722 + mrSges(5,3) * t668 + t679 * t700 - t691 * t732 + t759;
t634 = t640 * t750 - t746 * t643;
t690 = -mrSges(5,2) * t732 + mrSges(5,3) * t700;
t632 = -pkin(4) * t722 - pkin(10) * t731 + t701 * t680 - t634;
t630 = -pkin(5) * t649 - qJ(6) * t687 + t672 * t689 + qJDD(6) + t632;
t758 = m(7) * t630 - t649 * mrSges(7,1) - t688 * t670;
t755 = -m(6) * t632 + t649 * mrSges(6,1) + t650 * t783 + t688 * t671 + t689 * t771 - t758;
t619 = m(5) * t634 + t722 * mrSges(5,1) - t669 * mrSges(5,3) - t701 * t679 + t732 * t690 + t755;
t606 = t746 * t613 + t750 * t619;
t702 = -mrSges(4,1) * t718 + mrSges(4,2) * t719;
t704 = mrSges(4,2) * t763 + mrSges(4,3) * t718;
t604 = m(4) * t644 - mrSges(4,1) * t730 - mrSges(4,3) * t707 - t702 * t719 - t704 * t763 + t606;
t705 = -mrSges(4,1) * t763 - mrSges(4,3) * t719;
t760 = t750 * t613 - t619 * t746;
t605 = m(4) * t645 + mrSges(4,2) * t730 + mrSges(4,3) * t706 + t702 * t718 + t705 * t763 + t760;
t761 = -t604 * t741 + t743 * t605;
t596 = m(3) * t698 - mrSges(3,2) * t737 + mrSges(3,3) * t730 - t723 * t738 + t728 * t763 + t761;
t599 = t743 * t604 + t741 * t605;
t712 = -t742 * t725 - t784;
t724 = -mrSges(3,2) * t738 + mrSges(3,3) * t763;
t598 = m(3) * t712 - t730 * mrSges(3,1) + t729 * mrSges(3,2) + (t723 * t747 - t724 * t751) * t769 + t599;
t615 = t749 * t617 + t745 * t622;
t756 = m(5) * t651 - t668 * mrSges(5,1) + mrSges(5,2) * t669 - t700 * t690 + t691 * t701 + t615;
t754 = -m(4) * t684 + t706 * mrSges(4,1) - mrSges(4,2) * t707 + t718 * t704 - t705 * t719 - t756;
t610 = m(3) * t697 + mrSges(3,1) * t737 - mrSges(3,3) * t729 + t724 * t738 - t728 * t764 + t754;
t586 = t596 * t777 - t598 * t742 + t610 * t776;
t584 = m(2) * t733 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t753 + t586;
t591 = t751 * t596 - t610 * t747;
t590 = m(2) * t734 - mrSges(2,1) * t753 - qJDD(1) * mrSges(2,2) + t591;
t775 = t752 * t584 + t748 * t590;
t774 = t786 * t688 + t781 * t689 + t785 * t699;
t773 = -t787 * t688 - t782 * t689 - t786 * t699;
t772 = t782 * t688 + t788 * t689 + t781 * t699;
t585 = t596 * t779 + t744 * t598 + t610 * t778;
t762 = -t584 * t748 + t752 * t590;
t607 = -mrSges(6,1) * t632 + mrSges(6,3) * t629 - mrSges(7,1) * t630 + mrSges(7,3) * t627 - pkin(5) * t758 + qJ(6) * t765 + (-qJ(6) * t673 + t772) * t699 + (-pkin(5) * t673 - t774) * t689 + (-mrSges(7,2) * qJ(6) + t786) * t667 + (-mrSges(7,2) * pkin(5) + t782) * t650 + t787 * t649;
t623 = -t650 * mrSges(7,3) - t689 * t663 + t766;
t614 = mrSges(6,2) * t632 + mrSges(7,2) * t630 - mrSges(6,3) * t628 - mrSges(7,3) * t625 - qJ(6) * t623 + t782 * t649 + t788 * t650 + t781 * t667 + t774 * t688 + t773 * t699;
t675 = Ifges(5,5) * t701 + Ifges(5,6) * t700 + Ifges(5,3) * t732;
t676 = Ifges(5,4) * t701 + Ifges(5,2) * t700 + Ifges(5,6) * t732;
t592 = mrSges(5,2) * t651 - mrSges(5,3) * t634 + Ifges(5,1) * t669 + Ifges(5,4) * t668 + Ifges(5,5) * t722 - pkin(10) * t615 - t607 * t745 + t614 * t749 + t675 * t700 - t676 * t732;
t677 = Ifges(5,1) * t701 + Ifges(5,4) * t700 + Ifges(5,5) * t732;
t600 = -mrSges(5,1) * t651 - mrSges(6,1) * t628 - mrSges(7,1) * t625 + mrSges(6,2) * t629 + mrSges(7,2) * t627 + mrSges(5,3) * t635 + Ifges(5,4) * t669 + Ifges(5,2) * t668 + Ifges(5,6) * t722 - pkin(4) * t615 - pkin(5) * t623 - t701 * t675 + t732 * t677 + t773 * t689 + t772 * t688 - t785 * t667 - t781 * t650 - t786 * t649;
t692 = Ifges(4,5) * t719 + Ifges(4,6) * t718 - Ifges(4,3) * t763;
t694 = Ifges(4,1) * t719 + Ifges(4,4) * t718 - Ifges(4,5) * t763;
t582 = -mrSges(4,1) * t684 + mrSges(4,3) * t645 + Ifges(4,4) * t707 + Ifges(4,2) * t706 - Ifges(4,6) * t730 - pkin(3) * t756 + pkin(9) * t760 + t746 * t592 + t750 * t600 - t719 * t692 - t694 * t763;
t693 = Ifges(4,4) * t719 + Ifges(4,2) * t718 - Ifges(4,6) * t763;
t587 = mrSges(4,2) * t684 - mrSges(4,3) * t644 + Ifges(4,1) * t707 + Ifges(4,4) * t706 - Ifges(4,5) * t730 - pkin(9) * t606 + t592 * t750 - t600 * t746 + t692 * t718 + t693 * t763;
t709 = Ifges(3,3) * t738 + (Ifges(3,5) * t747 + Ifges(3,6) * t751) * t769;
t710 = Ifges(3,6) * t738 + (Ifges(3,4) * t747 + Ifges(3,2) * t751) * t769;
t580 = mrSges(3,2) * t712 - mrSges(3,3) * t697 + Ifges(3,1) * t729 + Ifges(3,4) * t730 + Ifges(3,5) * t737 - qJ(3) * t599 - t582 * t741 + t587 * t743 + t709 * t763 - t710 * t738;
t711 = Ifges(3,5) * t738 + (Ifges(3,1) * t747 + Ifges(3,4) * t751) * t769;
t581 = -t709 * t764 - pkin(10) * t759 - pkin(4) * t755 - pkin(2) * t599 + (Ifges(3,2) + Ifges(4,3)) * t730 - pkin(3) * t606 - t749 * t607 - t745 * t614 + Ifges(3,6) * t737 + t738 * t711 - Ifges(5,3) * t722 + Ifges(3,4) * t729 - mrSges(5,1) * t634 + mrSges(5,2) * t635 - mrSges(4,1) * t644 + mrSges(4,2) * t645 - Ifges(5,6) * t668 - Ifges(5,5) * t669 + mrSges(3,3) * t698 + t700 * t677 - t701 * t676 - Ifges(4,6) * t706 - Ifges(4,5) * t707 - mrSges(3,1) * t712 + t718 * t694 - t719 * t693;
t757 = pkin(8) * t591 + t580 * t747 + t581 * t751;
t579 = Ifges(3,5) * t729 + Ifges(3,6) * t730 + Ifges(3,3) * t737 + mrSges(3,1) * t697 - mrSges(3,2) * t698 + t741 * t587 + t743 * t582 + pkin(2) * t754 + qJ(3) * t761 + (t710 * t747 - t711 * t751) * t769;
t578 = -mrSges(2,2) * g(3) - mrSges(2,3) * t733 + Ifges(2,5) * qJDD(1) - t753 * Ifges(2,6) + t751 * t580 - t747 * t581 + (-t585 * t742 - t586 * t744) * pkin(8);
t577 = mrSges(2,1) * g(3) + mrSges(2,3) * t734 + t753 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t585 - t742 * t579 + t744 * t757;
t1 = [-m(1) * g(1) + t762; -m(1) * g(2) + t775; (-m(1) - m(2)) * g(3) + t585; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t775 - t748 * t577 + t752 * t578; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t762 + t752 * t577 + t748 * t578; -mrSges(1,1) * g(2) + mrSges(2,1) * t733 + mrSges(1,2) * g(1) - mrSges(2,2) * t734 + Ifges(2,3) * qJDD(1) + pkin(1) * t586 + t744 * t579 + t742 * t757;];
tauB  = t1;
