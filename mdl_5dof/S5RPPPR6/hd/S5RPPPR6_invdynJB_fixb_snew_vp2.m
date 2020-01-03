% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPPPR6
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
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
% Datum: 2019-12-31 17:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPPPR6_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR6_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR6_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR6_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR6_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR6_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR6_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR6_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR6_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:47:47
% EndTime: 2019-12-31 17:47:49
% DurationCPUTime: 2.32s
% Computational Cost: add. (17502->253), mult. (44883->327), div. (0->0), fcn. (25558->8), ass. (0->127)
t790 = -2 * qJD(4);
t730 = sin(qJ(1));
t732 = cos(qJ(1));
t708 = -t732 * g(1) - t730 * g(2);
t733 = qJD(1) ^ 2;
t789 = -t733 * pkin(1) + qJDD(1) * qJ(2) + (2 * qJD(1) * qJD(2)) + t708;
t788 = Ifges(4,5) - Ifges(3,6);
t726 = sin(pkin(7));
t766 = t726 * qJD(1);
t709 = -0.2e1 * qJD(3) * t766;
t722 = t726 ^ 2;
t728 = cos(pkin(7));
t723 = t728 ^ 2;
t707 = t730 * g(1) - t732 * g(2);
t752 = qJDD(2) - t707;
t778 = qJ(3) * t726;
t664 = t709 + (-qJ(2) + (-t722 - t723) * pkin(3)) * t733 + (-t778 - pkin(1) + (-pkin(2) - qJ(4)) * t728) * qJDD(1) + t752;
t765 = t728 * qJD(1);
t787 = t765 * t790 + t664;
t679 = -t726 * g(3) + t789 * t728;
t749 = -pkin(2) * t728 - t778;
t697 = t749 * qJD(1);
t670 = -t697 * t765 - t679;
t775 = t723 * t733;
t776 = t722 * t733;
t786 = -t775 - t776;
t785 = (Ifges(4,4) - Ifges(3,5)) * t726;
t678 = -t728 * g(3) - t789 * t726;
t784 = Ifges(3,1) + Ifges(4,2);
t783 = Ifges(3,4) + Ifges(4,6);
t782 = Ifges(3,2) + Ifges(4,3);
t781 = mrSges(3,2) * t726;
t780 = mrSges(4,2) * t728;
t725 = sin(pkin(8));
t779 = mrSges(5,2) * t725;
t774 = t725 * t728;
t747 = mrSges(5,1) * t726 + mrSges(5,3) * t774;
t694 = t747 * qJD(1);
t777 = t694 * t725;
t773 = t726 * t733;
t669 = t697 * t766 + qJDD(3) - t678;
t660 = (-qJ(4) * t728 * t733 + pkin(3) * qJDD(1)) * t726 + t669;
t727 = cos(pkin(8));
t772 = t727 * t660;
t771 = t727 * t728;
t698 = (-mrSges(4,3) * t726 + t780) * qJD(1);
t699 = (-mrSges(3,1) * t728 + t781) * qJD(1);
t657 = t725 * t660 + t727 * t787;
t688 = (mrSges(5,1) * t727 - t779) * t765;
t746 = -mrSges(5,2) * t726 - mrSges(5,3) * t771;
t689 = (pkin(4) * t727 + pkin(6) * t725) * t765;
t761 = t727 * t765;
t763 = qJDD(1) * t726;
t655 = -pkin(4) * t776 + pkin(6) * t763 - t689 * t761 + t657;
t762 = qJDD(1) * t728;
t662 = pkin(3) * t762 - qJ(4) * t775 + qJDD(4) - t670;
t658 = ((qJDD(1) * t725 + t727 * t773) * pkin(6) + (qJDD(1) * t727 - t725 * t773) * pkin(4)) * t728 + t662;
t729 = sin(qJ(5));
t731 = cos(qJ(5));
t652 = -t729 * t655 + t731 * t658;
t744 = t726 * t731 + t729 * t774;
t686 = t744 * qJD(1);
t745 = t726 * t729 - t731 * t774;
t687 = t745 * qJD(1);
t671 = -t686 * mrSges(6,1) + t687 * mrSges(6,2);
t674 = t686 * qJD(5) + t745 * qJDD(1);
t704 = qJD(5) + t761;
t676 = -t704 * mrSges(6,2) + t686 * mrSges(6,3);
t757 = t727 * t762;
t703 = qJDD(5) + t757;
t650 = m(6) * t652 + t703 * mrSges(6,1) - t674 * mrSges(6,3) - t687 * t671 + t704 * t676;
t653 = t731 * t655 + t729 * t658;
t673 = -t687 * qJD(5) + t744 * qJDD(1);
t677 = t704 * mrSges(6,1) - t687 * mrSges(6,3);
t651 = m(6) * t653 - t703 * mrSges(6,2) + t673 * mrSges(6,3) + t686 * t671 - t704 * t677;
t753 = -t729 * t650 + t731 * t651;
t640 = m(5) * t657 + t746 * qJDD(1) + (-t688 * t771 - t694 * t726) * qJD(1) + t753;
t656 = -t725 * t787 + t772;
t695 = t746 * qJD(1);
t654 = -pkin(4) * t763 - pkin(6) * t776 - t772 + (t664 + (t790 - t689) * t765) * t725;
t739 = -m(6) * t654 + t673 * mrSges(6,1) - t674 * mrSges(6,2) + t686 * t676 - t687 * t677;
t646 = m(5) * t656 + t747 * qJDD(1) + (t688 * t774 + t695 * t726) * qJD(1) + t739;
t635 = t725 * t640 + t727 * t646;
t741 = m(4) * t669 + t635;
t628 = m(3) * t678 + ((-mrSges(4,1) - mrSges(3,3)) * qJDD(1) + (-t698 - t699) * qJD(1)) * t726 - t741;
t642 = t731 * t650 + t729 * t651;
t748 = m(5) * t662 + mrSges(5,1) * t757 + t695 * t761 + t642;
t738 = -m(4) * t670 + mrSges(4,1) * t762 + t698 * t765 + t748;
t638 = m(3) * t679 + ((mrSges(3,3) - t779) * qJDD(1) + (t699 - t777) * qJD(1)) * t728 + t738;
t754 = -t726 * t628 + t728 * t638;
t622 = m(2) * t708 - t733 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t754;
t743 = -t733 * qJ(2) + t752;
t693 = -qJDD(1) * pkin(1) + t743;
t675 = t709 + (-pkin(1) + t749) * qJDD(1) + t743;
t769 = t727 * t640 - t725 * t646;
t742 = m(4) * t675 + mrSges(4,1) * t786 - mrSges(4,3) * t763 + t769;
t735 = m(3) * t693 - mrSges(3,1) * t762 + mrSges(3,3) * t786 + t742;
t751 = t780 + t781;
t630 = m(2) * t707 - t733 * mrSges(2,2) + (mrSges(2,1) - t751) * qJDD(1) - t735;
t770 = t730 * t622 + t732 * t630;
t624 = t728 * t628 + t726 * t638;
t767 = (t728 * t788 + t785) * qJD(1);
t755 = t732 * t622 - t730 * t630;
t750 = -Ifges(5,5) * t725 - Ifges(5,6) * t727;
t665 = Ifges(6,5) * t687 + Ifges(6,6) * t686 + Ifges(6,3) * t704;
t667 = Ifges(6,1) * t687 + Ifges(6,4) * t686 + Ifges(6,5) * t704;
t643 = -mrSges(6,1) * t654 + mrSges(6,3) * t653 + Ifges(6,4) * t674 + Ifges(6,2) * t673 + Ifges(6,6) * t703 - t687 * t665 + t704 * t667;
t666 = Ifges(6,4) * t687 + Ifges(6,2) * t686 + Ifges(6,6) * t704;
t644 = mrSges(6,2) * t654 - mrSges(6,3) * t652 + Ifges(6,1) * t674 + Ifges(6,4) * t673 + Ifges(6,5) * t703 + t686 * t665 - t704 * t666;
t680 = (Ifges(5,3) * t726 + t750 * t728) * qJD(1);
t736 = Ifges(5,6) * t726 + (-Ifges(5,4) * t725 - Ifges(5,2) * t727) * t728;
t681 = t736 * qJD(1);
t737 = Ifges(5,5) * t726 + (-Ifges(5,1) * t725 - Ifges(5,4) * t727) * t728;
t625 = mrSges(5,2) * t662 - mrSges(5,3) * t656 - pkin(6) * t642 - t729 * t643 + t731 * t644 + (-t680 * t771 - t681 * t726) * qJD(1) + t737 * qJDD(1);
t682 = t737 * qJD(1);
t734 = mrSges(6,1) * t652 - mrSges(6,2) * t653 + Ifges(6,5) * t674 + Ifges(6,6) * t673 + Ifges(6,3) * t703 + t687 * t666 - t686 * t667;
t626 = -mrSges(5,1) * t662 + mrSges(5,3) * t657 - pkin(4) * t642 + (t680 * t774 + t726 * t682) * qJD(1) + t736 * qJDD(1) - t734;
t634 = mrSges(4,2) * t762 + t742;
t641 = (-mrSges(5,2) * qJDD(1) - qJD(1) * t694) * t774 + t748;
t617 = -mrSges(3,1) * t693 + mrSges(3,3) * t679 - mrSges(4,1) * t670 + mrSges(4,2) * t675 - t725 * t625 - t727 * t626 + pkin(3) * t641 - qJ(4) * t769 - pkin(2) * t634 + t767 * t766 + (t783 * t726 + t782 * t728) * qJDD(1);
t619 = -qJ(3) * t634 - mrSges(3,3) * t678 + mrSges(3,2) * t693 + pkin(3) * t635 + mrSges(4,1) * t669 - mrSges(4,3) * t675 + pkin(4) * t739 + pkin(6) * t753 + mrSges(5,1) * t656 - mrSges(5,2) * t657 + t729 * t644 + t731 * t643 + (Ifges(5,3) + t784) * t763 + ((t750 + t783) * qJDD(1) + (-t681 * t725 + t682 * t727 - t767) * qJD(1)) * t728;
t632 = t751 * qJDD(1) + t735;
t740 = mrSges(2,1) * t707 - mrSges(2,2) * t708 + Ifges(2,3) * qJDD(1) - pkin(1) * t632 + qJ(2) * t754 + t728 * t617 + t726 * t619;
t633 = (qJDD(1) * mrSges(4,1) + qJD(1) * t698) * t726 + t741;
t615 = mrSges(2,1) * g(3) - pkin(1) * t624 + mrSges(2,3) * t708 + t733 * Ifges(2,5) + pkin(2) * t633 - qJ(3) * t738 + t725 * t626 + qJ(4) * t635 - mrSges(3,1) * t678 + mrSges(3,2) * t679 - t727 * t625 - mrSges(4,2) * t669 + mrSges(4,3) * t670 + (Ifges(2,6) + t785 + (qJ(3) * t779 + t788) * t728) * qJDD(1) + (-t783 * t722 * qJD(1) + (qJ(3) * t777 + t783 * t765 + (-t782 + t784) * t766) * t728) * qJD(1);
t614 = -mrSges(2,2) * g(3) - mrSges(2,3) * t707 + Ifges(2,5) * qJDD(1) - t733 * Ifges(2,6) - qJ(2) * t624 - t726 * t617 + t728 * t619;
t1 = [-m(1) * g(1) + t755; -m(1) * g(2) + t770; (-m(1) - m(2)) * g(3) + t624; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t770 + t732 * t614 - t730 * t615; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t755 + t730 * t614 + t732 * t615; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t740; t740; t632; t633; t641; t734;];
tauJB = t1;
