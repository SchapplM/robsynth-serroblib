% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRPRR12
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-05-07 15:13
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRPRR12_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR12_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR12_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR12_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR12_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR12_invdynB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR12_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR12_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR12_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 15:00:06
% EndTime: 2019-05-07 15:01:22
% DurationCPUTime: 64.37s
% Computational Cost: add. (1070292->398), mult. (2337960->517), div. (0->0), fcn. (1889306->14), ass. (0->163)
t749 = sin(pkin(6));
t755 = sin(qJ(2));
t759 = cos(qJ(2));
t775 = qJD(1) * qJD(2);
t736 = (-qJDD(1) * t759 + t755 * t775) * t749;
t786 = cos(qJ(3));
t785 = pkin(8) * t749;
t751 = cos(pkin(6));
t784 = t751 * g(3);
t783 = t749 * t755;
t782 = t749 * t759;
t781 = t751 * t755;
t780 = t751 * t759;
t756 = sin(qJ(1));
t760 = cos(qJ(1));
t741 = t756 * g(1) - g(2) * t760;
t761 = qJD(1) ^ 2;
t731 = qJDD(1) * pkin(1) + t761 * t785 + t741;
t742 = -g(1) * t760 - g(2) * t756;
t732 = -pkin(1) * t761 + qJDD(1) * t785 + t742;
t778 = t731 * t781 + t759 * t732;
t706 = -g(3) * t783 + t778;
t745 = qJD(1) * t751 + qJD(2);
t777 = qJD(1) * t749;
t774 = t755 * t777;
t729 = mrSges(3,1) * t745 - mrSges(3,3) * t774;
t733 = (-mrSges(3,1) * t759 + mrSges(3,2) * t755) * t777;
t744 = qJDD(1) * t751 + qJDD(2);
t734 = (-pkin(2) * t759 - pkin(9) * t755) * t777;
t743 = t745 ^ 2;
t776 = qJD(1) * t759;
t683 = -t743 * pkin(2) + t744 * pkin(9) + (-g(3) * t755 + t734 * t776) * t749 + t778;
t735 = (qJDD(1) * t755 + t759 * t775) * t749;
t684 = t736 * pkin(2) - t735 * pkin(9) - t784 + (-t731 + (pkin(2) * t755 - pkin(9) * t759) * t745 * qJD(1)) * t749;
t754 = sin(qJ(3));
t663 = t683 * t786 + t754 * t684;
t725 = t754 * t745 + t774 * t786;
t703 = qJD(3) * t725 + t735 * t754 - t744 * t786;
t724 = -t745 * t786 + t754 * t774;
t708 = mrSges(4,1) * t724 + mrSges(4,2) * t725;
t773 = t749 * t776;
t740 = qJD(3) - t773;
t715 = mrSges(4,1) * t740 - mrSges(4,3) * t725;
t728 = qJDD(3) + t736;
t707 = pkin(3) * t724 - qJ(4) * t725;
t739 = t740 ^ 2;
t653 = -pkin(3) * t739 + qJ(4) * t728 - t707 * t724 + t663;
t705 = -g(3) * t782 + t731 * t780 - t755 * t732;
t682 = -t744 * pkin(2) - t743 * pkin(9) + t734 * t774 - t705;
t704 = -t724 * qJD(3) + t735 * t786 + t754 * t744;
t656 = (t724 * t740 - t704) * qJ(4) + (t725 * t740 + t703) * pkin(3) + t682;
t748 = sin(pkin(12));
t750 = cos(pkin(12));
t713 = t725 * t750 + t740 * t748;
t642 = -0.2e1 * qJD(4) * t713 - t748 * t653 + t750 * t656;
t690 = t704 * t750 + t728 * t748;
t712 = -t725 * t748 + t740 * t750;
t635 = (t712 * t724 - t690) * pkin(10) + (t712 * t713 + t703) * pkin(4) + t642;
t643 = 0.2e1 * qJD(4) * t712 + t750 * t653 + t748 * t656;
t689 = -t704 * t748 + t728 * t750;
t695 = pkin(4) * t724 - pkin(10) * t713;
t711 = t712 ^ 2;
t637 = -pkin(4) * t711 + pkin(10) * t689 - t695 * t724 + t643;
t753 = sin(qJ(5));
t758 = cos(qJ(5));
t629 = t758 * t635 - t753 * t637;
t686 = t712 * t758 - t713 * t753;
t659 = qJD(5) * t686 + t689 * t753 + t690 * t758;
t687 = t712 * t753 + t713 * t758;
t701 = qJDD(5) + t703;
t723 = qJD(5) + t724;
t627 = (t686 * t723 - t659) * pkin(11) + (t686 * t687 + t701) * pkin(5) + t629;
t630 = t753 * t635 + t758 * t637;
t658 = -qJD(5) * t687 + t689 * t758 - t690 * t753;
t673 = pkin(5) * t723 - pkin(11) * t687;
t685 = t686 ^ 2;
t628 = -pkin(5) * t685 + pkin(11) * t658 - t673 * t723 + t630;
t752 = sin(qJ(6));
t757 = cos(qJ(6));
t625 = t627 * t757 - t628 * t752;
t668 = t686 * t757 - t687 * t752;
t641 = qJD(6) * t668 + t658 * t752 + t659 * t757;
t669 = t686 * t752 + t687 * t757;
t650 = -mrSges(7,1) * t668 + mrSges(7,2) * t669;
t722 = qJD(6) + t723;
t660 = -mrSges(7,2) * t722 + mrSges(7,3) * t668;
t696 = qJDD(6) + t701;
t621 = m(7) * t625 + mrSges(7,1) * t696 - mrSges(7,3) * t641 - t650 * t669 + t660 * t722;
t626 = t627 * t752 + t628 * t757;
t640 = -qJD(6) * t669 + t658 * t757 - t659 * t752;
t661 = mrSges(7,1) * t722 - mrSges(7,3) * t669;
t622 = m(7) * t626 - mrSges(7,2) * t696 + mrSges(7,3) * t640 + t650 * t668 - t661 * t722;
t615 = t757 * t621 + t752 * t622;
t670 = -mrSges(6,1) * t686 + mrSges(6,2) * t687;
t671 = -mrSges(6,2) * t723 + mrSges(6,3) * t686;
t613 = m(6) * t629 + mrSges(6,1) * t701 - mrSges(6,3) * t659 - t670 * t687 + t671 * t723 + t615;
t672 = mrSges(6,1) * t723 - mrSges(6,3) * t687;
t768 = -t621 * t752 + t757 * t622;
t614 = m(6) * t630 - mrSges(6,2) * t701 + mrSges(6,3) * t658 + t670 * t686 - t672 * t723 + t768;
t609 = t758 * t613 + t753 * t614;
t691 = -mrSges(5,1) * t712 + mrSges(5,2) * t713;
t693 = -mrSges(5,2) * t724 + mrSges(5,3) * t712;
t607 = m(5) * t642 + mrSges(5,1) * t703 - mrSges(5,3) * t690 - t691 * t713 + t693 * t724 + t609;
t694 = mrSges(5,1) * t724 - mrSges(5,3) * t713;
t769 = -t613 * t753 + t758 * t614;
t608 = m(5) * t643 - mrSges(5,2) * t703 + mrSges(5,3) * t689 + t691 * t712 - t694 * t724 + t769;
t770 = -t607 * t748 + t750 * t608;
t602 = m(4) * t663 - mrSges(4,2) * t728 - mrSges(4,3) * t703 - t708 * t724 - t715 * t740 + t770;
t662 = -t754 * t683 + t684 * t786;
t714 = -mrSges(4,2) * t740 - mrSges(4,3) * t724;
t652 = -t728 * pkin(3) - t739 * qJ(4) + t725 * t707 + qJDD(4) - t662;
t644 = -t689 * pkin(4) - t711 * pkin(10) + t713 * t695 + t652;
t632 = -t658 * pkin(5) - t685 * pkin(11) + t687 * t673 + t644;
t767 = m(7) * t632 - t640 * mrSges(7,1) + t641 * mrSges(7,2) - t668 * t660 + t669 * t661;
t764 = m(6) * t644 - t658 * mrSges(6,1) + mrSges(6,2) * t659 - t686 * t671 + t672 * t687 + t767;
t762 = -m(5) * t652 + t689 * mrSges(5,1) - mrSges(5,2) * t690 + t712 * t693 - t694 * t713 - t764;
t624 = m(4) * t662 + mrSges(4,1) * t728 - mrSges(4,3) * t704 - t708 * t725 + t714 * t740 + t762;
t771 = t602 * t786 - t624 * t754;
t593 = m(3) * t706 - mrSges(3,2) * t744 - mrSges(3,3) * t736 - t729 * t745 + t733 * t773 + t771;
t596 = t754 * t602 + t624 * t786;
t719 = -t749 * t731 - t784;
t730 = -mrSges(3,2) * t745 + mrSges(3,3) * t773;
t595 = m(3) * t719 + t736 * mrSges(3,1) + t735 * mrSges(3,2) + (t729 * t755 - t730 * t759) * t777 + t596;
t603 = t607 * t750 + t608 * t748;
t763 = -m(4) * t682 - t703 * mrSges(4,1) - mrSges(4,2) * t704 - t724 * t714 - t715 * t725 - t603;
t599 = m(3) * t705 + mrSges(3,1) * t744 - mrSges(3,3) * t735 + t730 * t745 - t733 * t774 + t763;
t582 = t593 * t781 - t595 * t749 + t599 * t780;
t580 = m(2) * t741 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t761 + t582;
t586 = t759 * t593 - t599 * t755;
t585 = m(2) * t742 - mrSges(2,1) * t761 - qJDD(1) * mrSges(2,2) + t586;
t779 = t760 * t580 + t756 * t585;
t581 = t593 * t783 + t751 * t595 + t599 * t782;
t772 = -t580 * t756 + t760 * t585;
t645 = Ifges(7,5) * t669 + Ifges(7,6) * t668 + Ifges(7,3) * t722;
t647 = Ifges(7,1) * t669 + Ifges(7,4) * t668 + Ifges(7,5) * t722;
t616 = -mrSges(7,1) * t632 + mrSges(7,3) * t626 + Ifges(7,4) * t641 + Ifges(7,2) * t640 + Ifges(7,6) * t696 - t645 * t669 + t647 * t722;
t646 = Ifges(7,4) * t669 + Ifges(7,2) * t668 + Ifges(7,6) * t722;
t617 = mrSges(7,2) * t632 - mrSges(7,3) * t625 + Ifges(7,1) * t641 + Ifges(7,4) * t640 + Ifges(7,5) * t696 + t645 * t668 - t646 * t722;
t664 = Ifges(6,5) * t687 + Ifges(6,6) * t686 + Ifges(6,3) * t723;
t666 = Ifges(6,1) * t687 + Ifges(6,4) * t686 + Ifges(6,5) * t723;
t604 = -mrSges(6,1) * t644 + mrSges(6,3) * t630 + Ifges(6,4) * t659 + Ifges(6,2) * t658 + Ifges(6,6) * t701 - pkin(5) * t767 + pkin(11) * t768 + t757 * t616 + t752 * t617 - t687 * t664 + t723 * t666;
t665 = Ifges(6,4) * t687 + Ifges(6,2) * t686 + Ifges(6,6) * t723;
t605 = mrSges(6,2) * t644 - mrSges(6,3) * t629 + Ifges(6,1) * t659 + Ifges(6,4) * t658 + Ifges(6,5) * t701 - pkin(11) * t615 - t616 * t752 + t617 * t757 + t664 * t686 - t665 * t723;
t674 = Ifges(5,5) * t713 + Ifges(5,6) * t712 + Ifges(5,3) * t724;
t676 = Ifges(5,1) * t713 + Ifges(5,4) * t712 + Ifges(5,5) * t724;
t588 = -mrSges(5,1) * t652 + mrSges(5,3) * t643 + Ifges(5,4) * t690 + Ifges(5,2) * t689 + Ifges(5,6) * t703 - pkin(4) * t764 + pkin(10) * t769 + t758 * t604 + t753 * t605 - t713 * t674 + t724 * t676;
t675 = Ifges(5,4) * t713 + Ifges(5,2) * t712 + Ifges(5,6) * t724;
t589 = mrSges(5,2) * t652 - mrSges(5,3) * t642 + Ifges(5,1) * t690 + Ifges(5,4) * t689 + Ifges(5,5) * t703 - pkin(10) * t609 - t604 * t753 + t605 * t758 + t674 * t712 - t675 * t724;
t697 = Ifges(4,5) * t725 - Ifges(4,6) * t724 + Ifges(4,3) * t740;
t698 = Ifges(4,4) * t725 - Ifges(4,2) * t724 + Ifges(4,6) * t740;
t578 = mrSges(4,2) * t682 - mrSges(4,3) * t662 + Ifges(4,1) * t704 - Ifges(4,4) * t703 + Ifges(4,5) * t728 - qJ(4) * t603 - t588 * t748 + t589 * t750 - t697 * t724 - t698 * t740;
t699 = Ifges(4,1) * t725 - Ifges(4,4) * t724 + Ifges(4,5) * t740;
t587 = -pkin(3) * t603 + t740 * t699 + Ifges(4,6) * t728 - t725 * t697 - t713 * t675 + t712 * t676 - Ifges(6,3) * t701 + Ifges(4,4) * t704 - Ifges(7,3) * t696 - t687 * t665 - Ifges(5,6) * t689 - Ifges(5,5) * t690 - mrSges(4,1) * t682 + t686 * t666 - t669 * t646 + mrSges(4,3) * t663 + t668 * t647 - Ifges(6,5) * t659 - Ifges(6,6) * t658 - Ifges(7,5) * t641 - mrSges(5,1) * t642 + mrSges(5,2) * t643 - Ifges(7,6) * t640 + mrSges(6,2) * t630 - mrSges(6,1) * t629 + mrSges(7,2) * t626 - mrSges(7,1) * t625 - pkin(5) * t615 + (-Ifges(5,3) - Ifges(4,2)) * t703 - pkin(4) * t609;
t716 = Ifges(3,3) * t745 + (Ifges(3,5) * t755 + Ifges(3,6) * t759) * t777;
t717 = Ifges(3,6) * t745 + (Ifges(3,4) * t755 + Ifges(3,2) * t759) * t777;
t576 = mrSges(3,2) * t719 - mrSges(3,3) * t705 + Ifges(3,1) * t735 - Ifges(3,4) * t736 + Ifges(3,5) * t744 - pkin(9) * t596 + t578 * t786 - t754 * t587 + t716 * t773 - t745 * t717;
t718 = Ifges(3,5) * t745 + (Ifges(3,1) * t755 + Ifges(3,4) * t759) * t777;
t577 = Ifges(3,4) * t735 - Ifges(3,2) * t736 + Ifges(3,6) * t744 - t716 * t774 + t745 * t718 - mrSges(3,1) * t719 + mrSges(3,3) * t706 - Ifges(4,5) * t704 + Ifges(4,6) * t703 - Ifges(4,3) * t728 - t725 * t698 - t724 * t699 - mrSges(4,1) * t662 + mrSges(4,2) * t663 - t748 * t589 - t750 * t588 - pkin(3) * t762 - qJ(4) * t770 - pkin(2) * t596;
t765 = pkin(8) * t586 + t576 * t755 + t577 * t759;
t575 = Ifges(3,5) * t735 - Ifges(3,6) * t736 + Ifges(3,3) * t744 + mrSges(3,1) * t705 - mrSges(3,2) * t706 + t754 * t578 + t786 * t587 + pkin(2) * t763 + pkin(9) * t771 + (t717 * t755 - t718 * t759) * t777;
t574 = -mrSges(2,2) * g(3) - mrSges(2,3) * t741 + Ifges(2,5) * qJDD(1) - t761 * Ifges(2,6) + t759 * t576 - t755 * t577 + (-t581 * t749 - t582 * t751) * pkin(8);
t573 = mrSges(2,1) * g(3) + mrSges(2,3) * t742 + t761 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t581 - t749 * t575 + t751 * t765;
t1 = [-m(1) * g(1) + t772; -m(1) * g(2) + t779; (-m(1) - m(2)) * g(3) + t581; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t779 - t756 * t573 + t760 * t574; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t772 + t760 * t573 + t756 * t574; -mrSges(1,1) * g(2) + mrSges(2,1) * t741 + mrSges(1,2) * g(1) - mrSges(2,2) * t742 + Ifges(2,3) * qJDD(1) + pkin(1) * t582 + t751 * t575 + t749 * t765;];
tauB  = t1;
