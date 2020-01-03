% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRRRP8
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 22:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRRRP8_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP8_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP8_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP8_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP8_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP8_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP8_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP8_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP8_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:00:14
% EndTime: 2019-12-31 22:00:23
% DurationCPUTime: 5.83s
% Computational Cost: add. (61774->290), mult. (123643->351), div. (0->0), fcn. (82853->8), ass. (0->116)
t775 = Ifges(5,4) + Ifges(6,4);
t784 = Ifges(5,2) + Ifges(6,2);
t780 = Ifges(5,6) + Ifges(6,6);
t746 = sin(qJ(3));
t750 = cos(qJ(3));
t747 = sin(qJ(2));
t771 = qJD(1) * t747;
t729 = t746 * qJD(2) + t750 * t771;
t751 = cos(qJ(2));
t769 = qJD(1) * qJD(2);
t766 = t751 * t769;
t732 = t747 * qJDD(1) + t766;
t701 = -t729 * qJD(3) + t750 * qJDD(2) - t746 * t732;
t728 = t750 * qJD(2) - t746 * t771;
t702 = t728 * qJD(3) + t746 * qJDD(2) + t750 * t732;
t745 = sin(qJ(4));
t749 = cos(qJ(4));
t704 = t749 * t728 - t745 * t729;
t668 = t704 * qJD(4) + t745 * t701 + t749 * t702;
t705 = t745 * t728 + t749 * t729;
t681 = -t704 * mrSges(6,1) + t705 * mrSges(6,2);
t748 = sin(qJ(1));
t752 = cos(qJ(1));
t737 = t748 * g(1) - t752 * g(2);
t754 = qJD(1) ^ 2;
t721 = -qJDD(1) * pkin(1) - t754 * pkin(6) - t737;
t741 = t747 * t769;
t733 = t751 * qJDD(1) - t741;
t685 = (-t732 - t766) * pkin(7) + (-t733 + t741) * pkin(2) + t721;
t738 = -t752 * g(1) - t748 * g(2);
t722 = -t754 * pkin(1) + qJDD(1) * pkin(6) + t738;
t710 = -t747 * g(3) + t751 * t722;
t731 = (-pkin(2) * t751 - pkin(7) * t747) * qJD(1);
t753 = qJD(2) ^ 2;
t770 = t751 * qJD(1);
t689 = -t753 * pkin(2) + qJDD(2) * pkin(7) + t731 * t770 + t710;
t669 = t750 * t685 - t746 * t689;
t727 = qJDD(3) - t733;
t740 = qJD(3) - t770;
t653 = (t728 * t740 - t702) * pkin(8) + (t728 * t729 + t727) * pkin(3) + t669;
t670 = t746 * t685 + t750 * t689;
t711 = t740 * pkin(3) - t729 * pkin(8);
t726 = t728 ^ 2;
t655 = -t726 * pkin(3) + t701 * pkin(8) - t740 * t711 + t670;
t647 = t749 * t653 - t745 * t655;
t723 = qJDD(4) + t727;
t739 = qJD(4) + t740;
t643 = -0.2e1 * qJD(5) * t705 + (t704 * t739 - t668) * qJ(5) + (t704 * t705 + t723) * pkin(4) + t647;
t690 = -t739 * mrSges(6,2) + t704 * mrSges(6,3);
t768 = m(6) * t643 + t723 * mrSges(6,1) + t739 * t690;
t639 = -t668 * mrSges(6,3) - t705 * t681 + t768;
t648 = t745 * t653 + t749 * t655;
t667 = -t705 * qJD(4) + t749 * t701 - t745 * t702;
t692 = t739 * pkin(4) - t705 * qJ(5);
t703 = t704 ^ 2;
t645 = -t703 * pkin(4) + t667 * qJ(5) + 0.2e1 * qJD(5) * t704 - t739 * t692 + t648;
t781 = Ifges(5,5) + Ifges(6,5);
t782 = Ifges(5,1) + Ifges(6,1);
t772 = -t775 * t704 - t782 * t705 - t781 * t739;
t778 = t784 * t704 + t775 * t705 + t780 * t739;
t779 = Ifges(5,3) + Ifges(6,3);
t783 = mrSges(5,1) * t647 + mrSges(6,1) * t643 - mrSges(5,2) * t648 - mrSges(6,2) * t645 + pkin(4) * t639 + t780 * t667 + t781 * t668 + t772 * t704 + t778 * t705 + t779 * t723;
t682 = -t704 * mrSges(5,1) + t705 * mrSges(5,2);
t691 = -t739 * mrSges(5,2) + t704 * mrSges(5,3);
t631 = m(5) * t647 + t723 * mrSges(5,1) + t739 * t691 + (-t681 - t682) * t705 + (-mrSges(5,3) - mrSges(6,3)) * t668 + t768;
t693 = t739 * mrSges(6,1) - t705 * mrSges(6,3);
t694 = t739 * mrSges(5,1) - t705 * mrSges(5,3);
t767 = m(6) * t645 + t667 * mrSges(6,3) + t704 * t681;
t634 = m(5) * t648 + t667 * mrSges(5,3) + t704 * t682 + (-t693 - t694) * t739 + (-mrSges(5,2) - mrSges(6,2)) * t723 + t767;
t629 = t749 * t631 + t745 * t634;
t696 = Ifges(4,4) * t729 + Ifges(4,2) * t728 + Ifges(4,6) * t740;
t697 = Ifges(4,1) * t729 + Ifges(4,4) * t728 + Ifges(4,5) * t740;
t777 = mrSges(4,1) * t669 - mrSges(4,2) * t670 + Ifges(4,5) * t702 + Ifges(4,6) * t701 + Ifges(4,3) * t727 + pkin(3) * t629 + t729 * t696 - t728 * t697 + t783;
t709 = -t751 * g(3) - t747 * t722;
t688 = -qJDD(2) * pkin(2) - t753 * pkin(7) + t731 * t771 - t709;
t656 = -t701 * pkin(3) - t726 * pkin(8) + t729 * t711 + t688;
t650 = -t667 * pkin(4) - t703 * qJ(5) + t705 * t692 + qJDD(5) + t656;
t640 = m(6) * t650 - t667 * mrSges(6,1) + t668 * mrSges(6,2) - t704 * t690 + t705 * t693;
t773 = -t780 * t704 - t781 * t705 - t779 * t739;
t624 = -mrSges(5,1) * t656 + mrSges(5,3) * t648 - mrSges(6,1) * t650 + mrSges(6,3) * t645 - pkin(4) * t640 + qJ(5) * t767 + (-qJ(5) * t693 - t772) * t739 + (-qJ(5) * mrSges(6,2) + t780) * t723 + t773 * t705 + t775 * t668 + t784 * t667;
t628 = mrSges(5,2) * t656 + mrSges(6,2) * t650 - mrSges(5,3) * t647 - mrSges(6,3) * t643 - qJ(5) * t639 + t775 * t667 + t782 * t668 - t773 * t704 + t781 * t723 - t778 * t739;
t695 = Ifges(4,5) * t729 + Ifges(4,6) * t728 + Ifges(4,3) * t740;
t759 = m(5) * t656 - t667 * mrSges(5,1) + t668 * mrSges(5,2) - t704 * t691 + t705 * t694 + t640;
t762 = -t745 * t631 + t749 * t634;
t609 = -mrSges(4,1) * t688 + mrSges(4,3) * t670 + Ifges(4,4) * t702 + Ifges(4,2) * t701 + Ifges(4,6) * t727 - pkin(3) * t759 + pkin(8) * t762 + t749 * t624 + t745 * t628 - t729 * t695 + t740 * t697;
t610 = mrSges(4,2) * t688 - mrSges(4,3) * t669 + Ifges(4,1) * t702 + Ifges(4,4) * t701 + Ifges(4,5) * t727 - pkin(8) * t629 - t745 * t624 + t749 * t628 + t728 * t695 - t740 * t696;
t706 = -t728 * mrSges(4,1) + t729 * mrSges(4,2);
t707 = -t740 * mrSges(4,2) + t728 * mrSges(4,3);
t626 = m(4) * t669 + t727 * mrSges(4,1) - t702 * mrSges(4,3) - t729 * t706 + t740 * t707 + t629;
t708 = t740 * mrSges(4,1) - t729 * mrSges(4,3);
t627 = m(4) * t670 - t727 * mrSges(4,2) + t701 * mrSges(4,3) + t728 * t706 - t740 * t708 + t762;
t623 = -t746 * t626 + t750 * t627;
t637 = -m(4) * t688 + t701 * mrSges(4,1) - t702 * mrSges(4,2) + t728 * t707 - t729 * t708 - t759;
t719 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t747 + Ifges(3,2) * t751) * qJD(1);
t720 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t747 + Ifges(3,4) * t751) * qJD(1);
t776 = mrSges(3,1) * t709 - mrSges(3,2) * t710 + Ifges(3,5) * t732 + Ifges(3,6) * t733 + Ifges(3,3) * qJDD(2) + pkin(2) * t637 + pkin(7) * t623 + t750 * t609 + t746 * t610 + (t747 * t719 - t751 * t720) * qJD(1);
t730 = (-mrSges(3,1) * t751 + mrSges(3,2) * t747) * qJD(1);
t735 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t771;
t621 = m(3) * t710 - qJDD(2) * mrSges(3,2) + t733 * mrSges(3,3) - qJD(2) * t735 + t730 * t770 + t623;
t736 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t770;
t636 = m(3) * t709 + qJDD(2) * mrSges(3,1) - t732 * mrSges(3,3) + qJD(2) * t736 - t730 * t771 + t637;
t763 = t751 * t621 - t747 * t636;
t613 = m(2) * t738 - t754 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t763;
t622 = t750 * t626 + t746 * t627;
t758 = -m(3) * t721 + t733 * mrSges(3,1) - t732 * mrSges(3,2) - t735 * t771 + t736 * t770 - t622;
t617 = m(2) * t737 + qJDD(1) * mrSges(2,1) - t754 * mrSges(2,2) + t758;
t774 = t748 * t613 + t752 * t617;
t615 = t747 * t621 + t751 * t636;
t764 = t752 * t613 - t748 * t617;
t718 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t747 + Ifges(3,6) * t751) * qJD(1);
t606 = mrSges(3,2) * t721 - mrSges(3,3) * t709 + Ifges(3,1) * t732 + Ifges(3,4) * t733 + Ifges(3,5) * qJDD(2) - pkin(7) * t622 - qJD(2) * t719 - t746 * t609 + t750 * t610 + t718 * t770;
t608 = -mrSges(3,1) * t721 + mrSges(3,3) * t710 + Ifges(3,4) * t732 + Ifges(3,2) * t733 + Ifges(3,6) * qJDD(2) - pkin(2) * t622 + qJD(2) * t720 - t718 * t771 - t777;
t760 = mrSges(2,1) * t737 - mrSges(2,2) * t738 + Ifges(2,3) * qJDD(1) + pkin(1) * t758 + pkin(6) * t763 + t747 * t606 + t751 * t608;
t604 = mrSges(2,1) * g(3) + mrSges(2,3) * t738 + t754 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t615 - t776;
t603 = -mrSges(2,2) * g(3) - mrSges(2,3) * t737 + Ifges(2,5) * qJDD(1) - t754 * Ifges(2,6) - pkin(6) * t615 + t751 * t606 - t747 * t608;
t1 = [-m(1) * g(1) + t764; -m(1) * g(2) + t774; (-m(1) - m(2)) * g(3) + t615; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t774 + t752 * t603 - t748 * t604; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t764 + t748 * t603 + t752 * t604; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t760; t760; t776; t777; t783; t640;];
tauJB = t1;
