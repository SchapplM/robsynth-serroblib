% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRRPP4
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
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
% tauJB [(6+5)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRRPP4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP4_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP4_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP4_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP4_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP4_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP4_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP4_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:55:17
% EndTime: 2019-12-31 20:55:23
% DurationCPUTime: 6.15s
% Computational Cost: add. (63099->291), mult. (141397->358), div. (0->0), fcn. (95840->8), ass. (0->117)
t808 = Ifges(5,1) + Ifges(6,1);
t797 = Ifges(5,4) - Ifges(6,5);
t806 = Ifges(6,4) + Ifges(5,5);
t807 = Ifges(5,2) + Ifges(6,3);
t804 = Ifges(5,6) - Ifges(6,6);
t805 = -Ifges(6,2) - Ifges(5,3);
t768 = sin(qJ(3));
t769 = sin(qJ(2));
t771 = cos(qJ(3));
t772 = cos(qJ(2));
t740 = (-t768 * t769 + t771 * t772) * qJD(1);
t741 = (t768 * t772 + t769 * t771) * qJD(1);
t767 = sin(pkin(8));
t796 = cos(pkin(8));
t726 = -t796 * t740 + t767 * t741;
t727 = t767 * t740 + t796 * t741;
t764 = qJD(2) + qJD(3);
t803 = t807 * t726 - t797 * t727 - t804 * t764;
t802 = -t797 * t726 + t808 * t727 + t806 * t764;
t789 = qJD(1) * qJD(2);
t748 = t769 * qJDD(1) + t772 * t789;
t770 = sin(qJ(1));
t773 = cos(qJ(1));
t755 = -t773 * g(1) - t770 * g(2);
t774 = qJD(1) ^ 2;
t743 = -t774 * pkin(1) + qJDD(1) * pkin(6) + t755;
t795 = t769 * t743;
t799 = pkin(2) * t774;
t705 = qJDD(2) * pkin(2) - t748 * pkin(7) - t795 + (pkin(7) * t789 + t769 * t799 - g(3)) * t772;
t730 = -t769 * g(3) + t772 * t743;
t749 = t772 * qJDD(1) - t769 * t789;
t791 = qJD(1) * t769;
t753 = qJD(2) * pkin(2) - pkin(7) * t791;
t766 = t772 ^ 2;
t706 = t749 * pkin(7) - qJD(2) * t753 - t766 * t799 + t730;
t676 = t771 * t705 - t768 * t706;
t714 = t740 * qJD(3) + t771 * t748 + t768 * t749;
t763 = qJDD(2) + qJDD(3);
t668 = (t740 * t764 - t714) * qJ(4) + (t740 * t741 + t763) * pkin(3) + t676;
t677 = t768 * t705 + t771 * t706;
t713 = -t741 * qJD(3) - t768 * t748 + t771 * t749;
t732 = t764 * pkin(3) - t741 * qJ(4);
t736 = t740 ^ 2;
t670 = -t736 * pkin(3) + t713 * qJ(4) - t764 * t732 + t677;
t800 = -2 * qJD(4);
t666 = t767 * t668 + t796 * t670 + t726 * t800;
t686 = -t796 * t713 + t767 * t714;
t717 = t764 * mrSges(5,1) - t727 * mrSges(5,3);
t698 = t726 * pkin(4) - t727 * qJ(5);
t762 = t764 ^ 2;
t660 = -t762 * pkin(4) + t763 * qJ(5) + 0.2e1 * qJD(5) * t764 - t726 * t698 + t666;
t718 = -t764 * mrSges(6,1) + t727 * mrSges(6,2);
t788 = m(6) * t660 + t763 * mrSges(6,3) + t764 * t718;
t699 = t726 * mrSges(6,1) - t727 * mrSges(6,3);
t792 = -t726 * mrSges(5,1) - t727 * mrSges(5,2) - t699;
t798 = -mrSges(5,3) - mrSges(6,2);
t649 = m(5) * t666 - t763 * mrSges(5,2) + t798 * t686 - t764 * t717 + t792 * t726 + t788;
t780 = t796 * t668 - t767 * t670;
t665 = t727 * t800 + t780;
t687 = t767 * t713 + t796 * t714;
t716 = -t764 * mrSges(5,2) - t726 * mrSges(5,3);
t661 = -t763 * pkin(4) - t762 * qJ(5) + qJDD(5) + ((2 * qJD(4)) + t698) * t727 - t780;
t719 = -t726 * mrSges(6,2) + t764 * mrSges(6,3);
t783 = -m(6) * t661 + t763 * mrSges(6,1) + t764 * t719;
t651 = m(5) * t665 + t763 * mrSges(5,1) + t798 * t687 + t764 * t716 + t792 * t727 + t783;
t643 = t767 * t649 + t796 * t651;
t728 = -t740 * mrSges(4,1) + t741 * mrSges(4,2);
t731 = -t764 * mrSges(4,2) + t740 * mrSges(4,3);
t638 = m(4) * t676 + t763 * mrSges(4,1) - t714 * mrSges(4,3) - t741 * t728 + t764 * t731 + t643;
t733 = t764 * mrSges(4,1) - t741 * mrSges(4,3);
t784 = t796 * t649 - t767 * t651;
t639 = m(4) * t677 - t763 * mrSges(4,2) + t713 * mrSges(4,3) + t740 * t728 - t764 * t733 + t784;
t634 = t771 * t638 + t768 * t639;
t729 = -t772 * g(3) - t795;
t738 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t769 + Ifges(3,2) * t772) * qJD(1);
t739 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t769 + Ifges(3,4) * t772) * qJD(1);
t656 = t687 * mrSges(6,2) + t727 * t699 - t783;
t722 = Ifges(4,4) * t741 + Ifges(4,2) * t740 + Ifges(4,6) * t764;
t723 = Ifges(4,1) * t741 + Ifges(4,4) * t740 + Ifges(4,5) * t764;
t776 = -mrSges(4,1) * t676 - mrSges(5,1) * t665 + mrSges(6,1) * t661 + mrSges(4,2) * t677 + mrSges(5,2) * t666 - mrSges(6,3) * t660 - Ifges(4,5) * t714 - Ifges(4,6) * t713 - pkin(3) * t643 + pkin(4) * t656 - qJ(5) * t788 - t741 * t722 + t740 * t723 + t803 * t727 + (qJ(5) * t699 - t802) * t726 - t806 * t687 + (qJ(5) * mrSges(6,2) + t804) * t686 + (-Ifges(4,3) + t805) * t763;
t801 = mrSges(3,1) * t729 - mrSges(3,2) * t730 + Ifges(3,5) * t748 + Ifges(3,6) * t749 + Ifges(3,3) * qJDD(2) + pkin(2) * t634 + (t769 * t738 - t772 * t739) * qJD(1) - t776;
t747 = (-mrSges(3,1) * t772 + mrSges(3,2) * t769) * qJD(1);
t790 = qJD(1) * t772;
t752 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t790;
t632 = m(3) * t729 + qJDD(2) * mrSges(3,1) - t748 * mrSges(3,3) + qJD(2) * t752 - t747 * t791 + t634;
t751 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t791;
t785 = -t768 * t638 + t771 * t639;
t633 = m(3) * t730 - qJDD(2) * mrSges(3,2) + t749 * mrSges(3,3) - qJD(2) * t751 + t747 * t790 + t785;
t786 = -t769 * t632 + t772 * t633;
t624 = m(2) * t755 - t774 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t786;
t754 = t770 * g(1) - t773 * g(2);
t781 = -qJDD(1) * pkin(1) - t754;
t742 = -t774 * pkin(6) + t781;
t715 = -t749 * pkin(2) + t753 * t791 + (-pkin(7) * t766 - pkin(6)) * t774 + t781;
t672 = -t713 * pkin(3) - t736 * qJ(4) + t741 * t732 + qJDD(4) + t715;
t663 = -0.2e1 * qJD(5) * t727 + (t726 * t764 - t687) * qJ(5) + (t727 * t764 + t686) * pkin(4) + t672;
t657 = m(6) * t663 + t686 * mrSges(6,1) - t687 * mrSges(6,3) - t727 * t718 + t726 * t719;
t652 = m(5) * t672 + t686 * mrSges(5,1) + t687 * mrSges(5,2) + t726 * t716 + t727 * t717 + t657;
t778 = m(4) * t715 - t713 * mrSges(4,1) + t714 * mrSges(4,2) - t740 * t731 + t741 * t733 + t652;
t777 = -m(3) * t742 + t749 * mrSges(3,1) - t748 * mrSges(3,2) - t751 * t791 + t752 * t790 - t778;
t645 = m(2) * t754 + qJDD(1) * mrSges(2,1) - t774 * mrSges(2,2) + t777;
t794 = t770 * t624 + t773 * t645;
t626 = t772 * t632 + t769 * t633;
t793 = t804 * t726 - t806 * t727 + t805 * t764;
t787 = t773 * t624 - t770 * t645;
t640 = -mrSges(5,1) * t672 - mrSges(6,1) * t663 + mrSges(6,2) * t660 + mrSges(5,3) * t666 - pkin(4) * t657 - t807 * t686 + t797 * t687 + t793 * t727 + t804 * t763 + t802 * t764;
t641 = mrSges(5,2) * t672 + mrSges(6,2) * t661 - mrSges(5,3) * t665 - mrSges(6,3) * t663 - qJ(5) * t657 - t797 * t686 + t808 * t687 + t793 * t726 + t806 * t763 + t803 * t764;
t721 = Ifges(4,5) * t741 + Ifges(4,6) * t740 + Ifges(4,3) * t764;
t627 = -mrSges(4,1) * t715 + mrSges(4,3) * t677 + Ifges(4,4) * t714 + Ifges(4,2) * t713 + Ifges(4,6) * t763 - pkin(3) * t652 + qJ(4) * t784 + t796 * t640 + t767 * t641 - t741 * t721 + t764 * t723;
t628 = mrSges(4,2) * t715 - mrSges(4,3) * t676 + Ifges(4,1) * t714 + Ifges(4,4) * t713 + Ifges(4,5) * t763 - qJ(4) * t643 - t767 * t640 + t796 * t641 + t740 * t721 - t764 * t722;
t737 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t769 + Ifges(3,6) * t772) * qJD(1);
t618 = -mrSges(3,1) * t742 + mrSges(3,3) * t730 + Ifges(3,4) * t748 + Ifges(3,2) * t749 + Ifges(3,6) * qJDD(2) - pkin(2) * t778 + pkin(7) * t785 + qJD(2) * t739 + t771 * t627 + t768 * t628 - t737 * t791;
t621 = mrSges(3,2) * t742 - mrSges(3,3) * t729 + Ifges(3,1) * t748 + Ifges(3,4) * t749 + Ifges(3,5) * qJDD(2) - pkin(7) * t634 - qJD(2) * t738 - t768 * t627 + t771 * t628 + t737 * t790;
t779 = mrSges(2,1) * t754 - mrSges(2,2) * t755 + Ifges(2,3) * qJDD(1) + pkin(1) * t777 + pkin(6) * t786 + t772 * t618 + t769 * t621;
t619 = mrSges(2,1) * g(3) + mrSges(2,3) * t755 + t774 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t626 - t801;
t616 = -mrSges(2,2) * g(3) - mrSges(2,3) * t754 + Ifges(2,5) * qJDD(1) - t774 * Ifges(2,6) - pkin(6) * t626 - t769 * t618 + t772 * t621;
t1 = [-m(1) * g(1) + t787; -m(1) * g(2) + t794; (-m(1) - m(2)) * g(3) + t626; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t794 + t773 * t616 - t770 * t619; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t787 + t770 * t616 + t773 * t619; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t779; t779; t801; -t776; t652; t656;];
tauJB = t1;
