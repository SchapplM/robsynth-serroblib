% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRRRR14
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
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,d3,d4,d5]';
% m [6x1]
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
% tauJB [(6+5)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 18:44
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRRRR14_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR14_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR14_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR14_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR14_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR14_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR14_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR14_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR14_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 18:43:14
% EndTime: 2024-09-27 18:43:29
% DurationCPUTime: 6.69s
% Computational Cost: add. (270419->285), mult. (414682->375), div. (0->0), fcn. (296334->12), ass. (0->126)
t758 = sin(pkin(5));
t793 = pkin(8) * t758;
t755 = qJD(1) + qJD(2);
t752 = t755 ^ 2;
t792 = t752 * t758 ^ 2;
t791 = t755 * t758;
t762 = sin(qJ(3));
t790 = t758 * t762;
t767 = cos(qJ(3));
t789 = t758 * t767;
t759 = cos(pkin(5));
t788 = t759 * t762;
t787 = t759 * t767;
t753 = qJDD(1) + qJDD(2);
t785 = qJD(3) * t755;
t730 = (t753 * t762 + t767 * t785) * t758;
t744 = t759 * t753 + qJDD(3);
t745 = t759 * t755 + qJD(3);
t764 = sin(qJ(1));
t769 = cos(qJ(1));
t746 = t764 * g(1) - t769 * g(2);
t738 = qJDD(1) * pkin(1) + t746;
t747 = -t769 * g(1) - t764 * g(2);
t770 = qJD(1) ^ 2;
t740 = -t770 * pkin(1) + t747;
t763 = sin(qJ(2));
t768 = cos(qJ(2));
t721 = t768 * t738 - t763 * t740;
t711 = t753 * pkin(2) + t752 * t793 + t721;
t722 = t763 * t738 + t768 * t740;
t712 = -t752 * pkin(2) + t753 * t793 + t722;
t777 = t711 * t787 - t762 * t712;
t681 = t744 * pkin(3) - t730 * pkin(9) + (pkin(3) * t762 * t792 + (pkin(9) * t745 * t755 - g(3)) * t758) * t767 + t777;
t691 = -g(3) * t790 + t711 * t788 + t767 * t712;
t783 = t755 * t790;
t729 = t745 * pkin(3) - pkin(9) * t783;
t731 = (t753 * t767 - t762 * t785) * t758;
t784 = t767 ^ 2 * t792;
t682 = -pkin(3) * t784 + t731 * pkin(9) - t745 * t729 + t691;
t761 = sin(qJ(4));
t766 = cos(qJ(4));
t668 = t766 * t681 - t761 * t682;
t724 = (-t761 * t762 + t766 * t767) * t791;
t697 = t724 * qJD(4) + t766 * t730 + t761 * t731;
t725 = (t761 * t767 + t762 * t766) * t791;
t741 = qJDD(4) + t744;
t743 = qJD(4) + t745;
t665 = (t724 * t743 - t697) * pkin(10) + (t724 * t725 + t741) * pkin(4) + t668;
t669 = t761 * t681 + t766 * t682;
t696 = -t725 * qJD(4) - t761 * t730 + t766 * t731;
t715 = t743 * pkin(4) - t725 * pkin(10);
t723 = t724 ^ 2;
t666 = -t723 * pkin(4) + t696 * pkin(10) - t743 * t715 + t669;
t760 = sin(qJ(5));
t765 = cos(qJ(5));
t663 = t765 * t665 - t760 * t666;
t704 = t765 * t724 - t760 * t725;
t676 = t704 * qJD(5) + t760 * t696 + t765 * t697;
t705 = t760 * t724 + t765 * t725;
t688 = -t704 * mrSges(6,1) + t705 * mrSges(6,2);
t739 = qJD(5) + t743;
t698 = -t739 * mrSges(6,2) + t704 * mrSges(6,3);
t737 = qJDD(5) + t741;
t660 = m(6) * t663 + t737 * mrSges(6,1) - t676 * mrSges(6,3) - t705 * t688 + t739 * t698;
t664 = t760 * t665 + t765 * t666;
t675 = -t705 * qJD(5) + t765 * t696 - t760 * t697;
t699 = t739 * mrSges(6,1) - t705 * mrSges(6,3);
t661 = m(6) * t664 - t737 * mrSges(6,2) + t675 * mrSges(6,3) + t704 * t688 - t739 * t699;
t652 = t765 * t660 + t760 * t661;
t706 = -t724 * mrSges(5,1) + t725 * mrSges(5,2);
t713 = -t743 * mrSges(5,2) + t724 * mrSges(5,3);
t649 = m(5) * t668 + t741 * mrSges(5,1) - t697 * mrSges(5,3) - t725 * t706 + t743 * t713 + t652;
t714 = t743 * mrSges(5,1) - t725 * mrSges(5,3);
t778 = -t760 * t660 + t765 * t661;
t650 = m(5) * t669 - t741 * mrSges(5,2) + t696 * mrSges(5,3) + t724 * t706 - t743 * t714 + t778;
t645 = t766 * t649 + t761 * t650;
t690 = -g(3) * t789 + t777;
t782 = t755 * t789;
t727 = -t745 * mrSges(4,2) + mrSges(4,3) * t782;
t728 = (-mrSges(4,1) * t767 + mrSges(4,2) * t762) * t791;
t643 = m(4) * t690 + t744 * mrSges(4,1) - t730 * mrSges(4,3) + t745 * t727 - t728 * t783 + t645;
t726 = t745 * mrSges(4,1) - mrSges(4,3) * t783;
t779 = -t761 * t649 + t766 * t650;
t644 = m(4) * t691 - t744 * mrSges(4,2) + t731 * mrSges(4,3) - t745 * t726 + t728 * t782 + t779;
t707 = -t759 * g(3) - t758 * t711;
t689 = -t731 * pkin(3) - pkin(9) * t784 + t729 * t783 + t707;
t671 = -t696 * pkin(4) - t723 * pkin(10) + t725 * t715 + t689;
t776 = m(6) * t671 - t675 * mrSges(6,1) + t676 * mrSges(6,2) - t704 * t698 + t705 * t699;
t772 = m(5) * t689 - t696 * mrSges(5,1) + t697 * mrSges(5,2) - t724 * t713 + t725 * t714 + t776;
t656 = -t731 * mrSges(4,1) + t772 + t730 * mrSges(4,2) + m(4) * t707 + (t726 * t762 - t727 * t767) * t791;
t628 = t643 * t787 + t644 * t788 - t758 * t656;
t625 = m(3) * t721 + t753 * mrSges(3,1) - t752 * mrSges(3,2) + t628;
t633 = -t762 * t643 + t767 * t644;
t631 = m(3) * t722 - t752 * mrSges(3,1) - t753 * mrSges(3,2) + t633;
t619 = t768 * t625 + t763 * t631;
t616 = m(2) * t746 + qJDD(1) * mrSges(2,1) - t770 * mrSges(2,2) + t619;
t780 = -t763 * t625 + t768 * t631;
t617 = m(2) * t747 - t770 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t780;
t786 = t769 * t616 + t764 * t617;
t627 = t643 * t789 + t644 * t790 + t759 * t656;
t781 = -t764 * t616 + t769 * t617;
t683 = Ifges(6,5) * t705 + Ifges(6,6) * t704 + Ifges(6,3) * t739;
t685 = Ifges(6,1) * t705 + Ifges(6,4) * t704 + Ifges(6,5) * t739;
t653 = -mrSges(6,1) * t671 + mrSges(6,3) * t664 + Ifges(6,4) * t676 + Ifges(6,2) * t675 + Ifges(6,6) * t737 - t705 * t683 + t739 * t685;
t684 = Ifges(6,4) * t705 + Ifges(6,2) * t704 + Ifges(6,6) * t739;
t654 = mrSges(6,2) * t671 - mrSges(6,3) * t663 + Ifges(6,1) * t676 + Ifges(6,4) * t675 + Ifges(6,5) * t737 + t704 * t683 - t739 * t684;
t700 = Ifges(5,5) * t725 + Ifges(5,6) * t724 + Ifges(5,3) * t743;
t702 = Ifges(5,1) * t725 + Ifges(5,4) * t724 + Ifges(5,5) * t743;
t636 = -mrSges(5,1) * t689 + mrSges(5,3) * t669 + Ifges(5,4) * t697 + Ifges(5,2) * t696 + Ifges(5,6) * t741 - pkin(4) * t776 + pkin(10) * t778 + t765 * t653 + t760 * t654 - t725 * t700 + t743 * t702;
t701 = Ifges(5,4) * t725 + Ifges(5,2) * t724 + Ifges(5,6) * t743;
t637 = mrSges(5,2) * t689 - mrSges(5,3) * t668 + Ifges(5,1) * t697 + Ifges(5,4) * t696 + Ifges(5,5) * t741 - pkin(10) * t652 - t760 * t653 + t765 * t654 + t724 * t700 - t743 * t701;
t716 = Ifges(4,3) * t745 + (Ifges(4,5) * t762 + Ifges(4,6) * t767) * t791;
t718 = Ifges(4,5) * t745 + (Ifges(4,1) * t762 + Ifges(4,4) * t767) * t791;
t621 = -mrSges(4,1) * t707 + mrSges(4,3) * t691 + Ifges(4,4) * t730 + Ifges(4,2) * t731 + Ifges(4,6) * t744 - pkin(3) * t772 + pkin(9) * t779 + t766 * t636 + t761 * t637 - t716 * t783 + t745 * t718;
t717 = Ifges(4,6) * t745 + (Ifges(4,4) * t762 + Ifges(4,2) * t767) * t791;
t623 = mrSges(4,2) * t707 - mrSges(4,3) * t690 + Ifges(4,1) * t730 + Ifges(4,4) * t731 + Ifges(4,5) * t744 - pkin(9) * t645 - t761 * t636 + t766 * t637 + t716 * t782 - t745 * t717;
t774 = mrSges(6,1) * t663 - mrSges(6,2) * t664 + Ifges(6,5) * t676 + Ifges(6,6) * t675 + Ifges(6,3) * t737 + t705 * t684 - t704 * t685;
t771 = mrSges(5,1) * t668 - mrSges(5,2) * t669 + Ifges(5,5) * t697 + Ifges(5,6) * t696 + Ifges(5,3) * t741 + pkin(4) * t652 + t725 * t701 - t724 * t702 + t774;
t635 = t771 + (t717 * t762 - t718 * t767) * t791 + Ifges(4,3) * t744 + Ifges(4,5) * t730 + Ifges(4,6) * t731 + mrSges(4,1) * t690 - mrSges(4,2) * t691 + pkin(3) * t645;
t775 = mrSges(3,1) * t721 - mrSges(3,2) * t722 + Ifges(3,3) * t753 + pkin(2) * t628 + t621 * t789 + t623 * t790 + t633 * t793 + t759 * t635;
t773 = mrSges(2,1) * t746 - mrSges(2,2) * t747 + Ifges(2,3) * qJDD(1) + pkin(1) * t619 + t775;
t612 = -mrSges(3,2) * g(3) - mrSges(3,3) * t721 + Ifges(3,5) * t753 - t752 * Ifges(3,6) - t762 * t621 + t767 * t623 + (-t627 * t758 - t628 * t759) * pkin(8);
t611 = mrSges(3,1) * g(3) + mrSges(3,3) * t722 + t752 * Ifges(3,5) + Ifges(3,6) * t753 - pkin(2) * t627 - t758 * t635 + (pkin(8) * t633 + t621 * t767 + t623 * t762) * t759;
t610 = -mrSges(2,2) * g(3) - mrSges(2,3) * t746 + Ifges(2,5) * qJDD(1) - t770 * Ifges(2,6) - pkin(7) * t619 - t763 * t611 + t768 * t612;
t609 = Ifges(2,6) * qJDD(1) + t770 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t747 + t763 * t612 + t768 * t611 - pkin(1) * (-m(3) * g(3) + t627) + pkin(7) * t780;
t1 = [-m(1) * g(1) + t781; -m(1) * g(2) + t786; (-m(1) - m(2) - m(3)) * g(3) + t627; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t786 - t764 * t609 + t769 * t610; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t781 + t769 * t609 + t764 * t610; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t773; t773; t775; t635; t771; t774;];
tauJB = t1;
