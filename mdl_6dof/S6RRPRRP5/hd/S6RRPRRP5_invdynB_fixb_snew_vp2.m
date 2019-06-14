% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPRRP5
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
% Datum: 2019-05-06 17:56
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPRRP5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP5_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP5_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP5_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP5_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP5_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP5_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP5_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 17:46:59
% EndTime: 2019-05-06 17:47:24
% DurationCPUTime: 23.09s
% Computational Cost: add. (338281->375), mult. (884943->474), div. (0->0), fcn. (696330->12), ass. (0->154)
t798 = -2 * qJD(3);
t797 = Ifges(6,1) + Ifges(7,1);
t791 = Ifges(6,4) + Ifges(7,4);
t790 = Ifges(6,5) + Ifges(7,5);
t796 = Ifges(6,2) + Ifges(7,2);
t795 = Ifges(6,6) + Ifges(7,6);
t794 = Ifges(6,3) + Ifges(7,3);
t746 = sin(pkin(11));
t748 = cos(pkin(11));
t752 = sin(qJ(2));
t756 = cos(qJ(2));
t747 = sin(pkin(6));
t778 = qJD(1) * t747;
t727 = (t746 * t752 - t748 * t756) * t778;
t776 = qJD(1) * qJD(2);
t736 = (qJDD(1) * t752 + t756 * t776) * t747;
t749 = cos(pkin(6));
t741 = qJDD(1) * t749 + qJDD(2);
t742 = qJD(1) * t749 + qJD(2);
t753 = sin(qJ(1));
t757 = cos(qJ(1));
t738 = t753 * g(1) - g(2) * t757;
t758 = qJD(1) ^ 2;
t793 = pkin(8) * t747;
t733 = qJDD(1) * pkin(1) + t758 * t793 + t738;
t739 = -g(1) * t757 - g(2) * t753;
t734 = -pkin(1) * t758 + qJDD(1) * t793 + t739;
t784 = t749 * t756;
t765 = t733 * t784 - t734 * t752;
t788 = t747 ^ 2 * t758;
t673 = pkin(2) * t741 - qJ(3) * t736 + (pkin(2) * t752 * t788 + (qJ(3) * qJD(1) * t742 - g(3)) * t747) * t756 + t765;
t785 = t749 * t752;
t787 = t747 * t752;
t704 = -g(3) * t787 + t733 * t785 + t756 * t734;
t772 = t752 * t778;
t730 = pkin(2) * t742 - qJ(3) * t772;
t737 = (qJDD(1) * t756 - t752 * t776) * t747;
t775 = t756 ^ 2 * t788;
t679 = -pkin(2) * t775 + qJ(3) * t737 - t730 * t742 + t704;
t728 = (t746 * t756 + t748 * t752) * t778;
t648 = t673 * t748 - t746 * t679 + t728 * t798;
t792 = -mrSges(6,2) - mrSges(7,2);
t786 = t747 * t756;
t649 = t746 * t673 + t748 * t679 + t727 * t798;
t705 = mrSges(4,1) * t727 + mrSges(4,2) * t728;
t709 = -t736 * t746 + t737 * t748;
t715 = mrSges(4,1) * t742 - mrSges(4,3) * t728;
t706 = pkin(3) * t727 - pkin(9) * t728;
t740 = t742 ^ 2;
t647 = -pkin(3) * t740 + pkin(9) * t741 - t706 * t727 + t649;
t719 = -g(3) * t749 - t733 * t747;
t690 = -pkin(2) * t737 - qJ(3) * t775 + t730 * t772 + qJDD(3) + t719;
t710 = t736 * t748 + t737 * t746;
t651 = (t727 * t742 - t710) * pkin(9) + (t728 * t742 - t709) * pkin(3) + t690;
t751 = sin(qJ(4));
t755 = cos(qJ(4));
t643 = t755 * t647 + t751 * t651;
t713 = t728 * t755 + t742 * t751;
t687 = -qJD(4) * t713 - t710 * t751 + t741 * t755;
t712 = -t728 * t751 + t742 * t755;
t691 = -mrSges(5,1) * t712 + mrSges(5,2) * t713;
t726 = qJD(4) + t727;
t698 = mrSges(5,1) * t726 - mrSges(5,3) * t713;
t708 = qJDD(4) - t709;
t692 = -pkin(4) * t712 - pkin(10) * t713;
t725 = t726 ^ 2;
t638 = -pkin(4) * t725 + pkin(10) * t708 + t692 * t712 + t643;
t646 = -pkin(3) * t741 - pkin(9) * t740 + t728 * t706 - t648;
t688 = qJD(4) * t712 + t710 * t755 + t741 * t751;
t641 = (-t712 * t726 - t688) * pkin(10) + (t713 * t726 - t687) * pkin(4) + t646;
t750 = sin(qJ(5));
t754 = cos(qJ(5));
t633 = -t638 * t750 + t754 * t641;
t695 = -t713 * t750 + t726 * t754;
t656 = qJD(5) * t695 + t688 * t754 + t708 * t750;
t696 = t713 * t754 + t726 * t750;
t667 = -mrSges(7,1) * t695 + mrSges(7,2) * t696;
t668 = -mrSges(6,1) * t695 + mrSges(6,2) * t696;
t711 = qJD(5) - t712;
t675 = -mrSges(6,2) * t711 + mrSges(6,3) * t695;
t686 = qJDD(5) - t687;
t630 = -0.2e1 * qJD(6) * t696 + (t695 * t711 - t656) * qJ(6) + (t695 * t696 + t686) * pkin(5) + t633;
t674 = -mrSges(7,2) * t711 + mrSges(7,3) * t695;
t774 = m(7) * t630 + t686 * mrSges(7,1) + t711 * t674;
t623 = m(6) * t633 + mrSges(6,1) * t686 + t675 * t711 + (-t667 - t668) * t696 + (-mrSges(6,3) - mrSges(7,3)) * t656 + t774;
t634 = t754 * t638 + t750 * t641;
t655 = -qJD(5) * t696 - t688 * t750 + t708 * t754;
t676 = pkin(5) * t711 - qJ(6) * t696;
t694 = t695 ^ 2;
t632 = -pkin(5) * t694 + qJ(6) * t655 + 0.2e1 * qJD(6) * t695 - t676 * t711 + t634;
t773 = m(7) * t632 + t655 * mrSges(7,3) + t695 * t667;
t677 = mrSges(7,1) * t711 - mrSges(7,3) * t696;
t779 = -mrSges(6,1) * t711 + mrSges(6,3) * t696 - t677;
t625 = m(6) * t634 + mrSges(6,3) * t655 + t668 * t695 + t792 * t686 + t779 * t711 + t773;
t767 = -t623 * t750 + t754 * t625;
t620 = m(5) * t643 - mrSges(5,2) * t708 + mrSges(5,3) * t687 + t691 * t712 - t698 * t726 + t767;
t642 = -t751 * t647 + t651 * t755;
t697 = -mrSges(5,2) * t726 + mrSges(5,3) * t712;
t637 = -pkin(4) * t708 - pkin(10) * t725 + t713 * t692 - t642;
t635 = -pkin(5) * t655 - qJ(6) * t694 + t676 * t696 + qJDD(6) + t637;
t764 = m(7) * t635 - t655 * mrSges(7,1) - t695 * t674;
t759 = -m(6) * t637 + t655 * mrSges(6,1) + t792 * t656 + t695 * t675 + t779 * t696 - t764;
t627 = m(5) * t642 + mrSges(5,1) * t708 - mrSges(5,3) * t688 - t691 * t713 + t697 * t726 + t759;
t768 = t755 * t620 - t627 * t751;
t611 = m(4) * t649 - mrSges(4,2) * t741 + mrSges(4,3) * t709 - t705 * t727 - t715 * t742 + t768;
t714 = -mrSges(4,2) * t742 - mrSges(4,3) * t727;
t622 = t623 * t754 + t625 * t750;
t760 = -m(5) * t646 + t687 * mrSges(5,1) - mrSges(5,2) * t688 + t712 * t697 - t698 * t713 - t622;
t617 = m(4) * t648 + mrSges(4,1) * t741 - mrSges(4,3) * t710 - t705 * t728 + t714 * t742 + t760;
t607 = t746 * t611 + t748 * t617;
t703 = -g(3) * t786 + t765;
t771 = t756 * t778;
t732 = -mrSges(3,2) * t742 + mrSges(3,3) * t771;
t735 = (-mrSges(3,1) * t756 + mrSges(3,2) * t752) * t778;
t605 = m(3) * t703 + mrSges(3,1) * t741 - mrSges(3,3) * t736 + t732 * t742 - t735 * t772 + t607;
t731 = mrSges(3,1) * t742 - mrSges(3,3) * t772;
t769 = t748 * t611 - t617 * t746;
t606 = m(3) * t704 - mrSges(3,2) * t741 + mrSges(3,3) * t737 - t731 * t742 + t735 * t771 + t769;
t614 = t751 * t620 + t755 * t627;
t761 = m(4) * t690 - mrSges(4,1) * t709 + t710 * mrSges(4,2) + t714 * t727 + t728 * t715 + t614;
t613 = m(3) * t719 - mrSges(3,1) * t737 + mrSges(3,2) * t736 + (t731 * t752 - t732 * t756) * t778 + t761;
t593 = t605 * t784 + t606 * t785 - t613 * t747;
t591 = m(2) * t738 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t758 + t593;
t598 = -t605 * t752 + t756 * t606;
t597 = m(2) * t739 - mrSges(2,1) * t758 - qJDD(1) * mrSges(2,2) + t598;
t783 = t757 * t591 + t753 * t597;
t782 = t695 * t795 + t696 * t790 + t711 * t794;
t781 = -t695 * t796 - t696 * t791 - t711 * t795;
t780 = t695 * t791 + t696 * t797 + t711 * t790;
t592 = t605 * t786 + t606 * t787 + t749 * t613;
t770 = -t591 * t753 + t757 * t597;
t615 = -mrSges(6,1) * t637 + mrSges(6,3) * t634 - mrSges(7,1) * t635 + mrSges(7,3) * t632 - pkin(5) * t764 + qJ(6) * t773 + (-qJ(6) * t677 + t780) * t711 + (-pkin(5) * t677 - t782) * t696 + (-mrSges(7,2) * qJ(6) + t795) * t686 + (-mrSges(7,2) * pkin(5) + t791) * t656 + t796 * t655;
t628 = -mrSges(7,3) * t656 - t667 * t696 + t774;
t621 = mrSges(6,2) * t637 + mrSges(7,2) * t635 - mrSges(6,3) * t633 - mrSges(7,3) * t630 - qJ(6) * t628 + t791 * t655 + t656 * t797 + t790 * t686 + t782 * t695 + t781 * t711;
t680 = Ifges(5,5) * t713 + Ifges(5,6) * t712 + Ifges(5,3) * t726;
t681 = Ifges(5,4) * t713 + Ifges(5,2) * t712 + Ifges(5,6) * t726;
t599 = mrSges(5,2) * t646 - mrSges(5,3) * t642 + Ifges(5,1) * t688 + Ifges(5,4) * t687 + Ifges(5,5) * t708 - pkin(10) * t622 - t615 * t750 + t621 * t754 + t680 * t712 - t681 * t726;
t682 = Ifges(5,1) * t713 + Ifges(5,4) * t712 + Ifges(5,5) * t726;
t608 = -mrSges(5,1) * t646 - mrSges(6,1) * t633 - mrSges(7,1) * t630 + mrSges(6,2) * t634 + mrSges(7,2) * t632 + mrSges(5,3) * t643 + Ifges(5,4) * t688 + Ifges(5,2) * t687 + Ifges(5,6) * t708 - pkin(4) * t622 - pkin(5) * t628 - t713 * t680 + t726 * t682 + t781 * t696 + t780 * t695 - t794 * t686 - t790 * t656 - t795 * t655;
t699 = Ifges(4,5) * t728 - Ifges(4,6) * t727 + Ifges(4,3) * t742;
t700 = Ifges(4,4) * t728 - Ifges(4,2) * t727 + Ifges(4,6) * t742;
t589 = mrSges(4,2) * t690 - mrSges(4,3) * t648 + Ifges(4,1) * t710 + Ifges(4,4) * t709 + Ifges(4,5) * t741 - pkin(9) * t614 + t599 * t755 - t608 * t751 - t699 * t727 - t700 * t742;
t701 = Ifges(4,1) * t728 - Ifges(4,4) * t727 + Ifges(4,5) * t742;
t594 = Ifges(4,4) * t710 + Ifges(4,2) * t709 + Ifges(4,6) * t741 - t728 * t699 + t742 * t701 - mrSges(4,1) * t690 + mrSges(4,3) * t649 - Ifges(5,5) * t688 - Ifges(5,6) * t687 - Ifges(5,3) * t708 - t713 * t681 + t712 * t682 - mrSges(5,1) * t642 + mrSges(5,2) * t643 - t750 * t621 - t754 * t615 - pkin(4) * t759 - pkin(10) * t767 - pkin(3) * t614;
t716 = Ifges(3,3) * t742 + (Ifges(3,5) * t752 + Ifges(3,6) * t756) * t778;
t718 = Ifges(3,5) * t742 + (Ifges(3,1) * t752 + Ifges(3,4) * t756) * t778;
t586 = -mrSges(3,1) * t719 + mrSges(3,3) * t704 + Ifges(3,4) * t736 + Ifges(3,2) * t737 + Ifges(3,6) * t741 - pkin(2) * t761 + qJ(3) * t769 + t746 * t589 + t748 * t594 - t716 * t772 + t742 * t718;
t717 = Ifges(3,6) * t742 + (Ifges(3,4) * t752 + Ifges(3,2) * t756) * t778;
t587 = mrSges(3,2) * t719 - mrSges(3,3) * t703 + Ifges(3,1) * t736 + Ifges(3,4) * t737 + Ifges(3,5) * t741 - qJ(3) * t607 + t589 * t748 - t594 * t746 + t716 * t771 - t717 * t742;
t762 = pkin(8) * t598 + t586 * t756 + t587 * t752;
t588 = Ifges(3,5) * t736 + Ifges(3,6) * t737 + mrSges(3,1) * t703 - mrSges(3,2) * t704 + Ifges(4,5) * t710 + Ifges(4,6) * t709 + t728 * t700 + t727 * t701 + mrSges(4,1) * t648 - mrSges(4,2) * t649 + t751 * t599 + t755 * t608 + pkin(3) * t760 + pkin(9) * t768 + pkin(2) * t607 + (Ifges(3,3) + Ifges(4,3)) * t741 + (t717 * t752 - t718 * t756) * t778;
t585 = -mrSges(2,2) * g(3) - mrSges(2,3) * t738 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t758 - t586 * t752 + t587 * t756 + (-t592 * t747 - t593 * t749) * pkin(8);
t584 = mrSges(2,1) * g(3) + mrSges(2,3) * t739 + Ifges(2,5) * t758 + Ifges(2,6) * qJDD(1) - pkin(1) * t592 - t588 * t747 + t762 * t749;
t1 = [-m(1) * g(1) + t770; -m(1) * g(2) + t783; (-m(1) - m(2)) * g(3) + t592; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t783 - t753 * t584 + t757 * t585; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t770 + t757 * t584 + t753 * t585; -mrSges(1,1) * g(2) + mrSges(2,1) * t738 + mrSges(1,2) * g(1) - mrSges(2,2) * t739 + Ifges(2,3) * qJDD(1) + pkin(1) * t593 + t588 * t749 + t762 * t747;];
tauB  = t1;
