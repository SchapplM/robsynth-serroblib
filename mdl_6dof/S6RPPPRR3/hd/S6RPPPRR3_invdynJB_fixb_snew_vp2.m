% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3,theta4]';
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
% tauJB [(6+6)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 13:42
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPPPRR3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR3_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR3_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPPRR3_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR3_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR3_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR3_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR3_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 13:41:03
% EndTime: 2019-05-05 13:41:08
% DurationCPUTime: 4.77s
% Computational Cost: add. (61809->275), mult. (124698->331), div. (0->0), fcn. (71359->10), ass. (0->124)
t722 = qJD(1) ^ 2;
t711 = sin(pkin(10));
t713 = cos(pkin(10));
t716 = sin(qJ(5));
t719 = cos(qJ(5));
t733 = t711 * t716 - t713 * t719;
t683 = t733 * qJD(1);
t734 = -t711 * t719 - t713 * t716;
t684 = t734 * qJD(1);
t743 = t684 * qJD(5);
t667 = t733 * qJDD(1) - t743;
t754 = -pkin(1) - pkin(2);
t753 = mrSges(2,1) + mrSges(3,1);
t752 = Ifges(3,4) + Ifges(2,5);
t751 = Ifges(2,6) - Ifges(3,6);
t750 = mrSges(5,1) * t713;
t749 = mrSges(5,2) * t711;
t701 = t713 ^ 2;
t748 = t701 * t722;
t717 = sin(qJ(1));
t720 = cos(qJ(1));
t689 = -t720 * g(1) - t717 * g(2);
t729 = qJDD(1) * qJ(2) + 0.2e1 * qJD(2) * qJD(1) + t689;
t680 = -pkin(1) * t722 + t729;
t674 = t754 * t722 + t729;
t688 = t717 * g(1) - t720 * g(2);
t728 = -t722 * qJ(2) + qJDD(2) - t688;
t679 = t754 * qJDD(1) + t728;
t712 = sin(pkin(9));
t714 = cos(pkin(9));
t656 = t714 * t674 + t712 * t679;
t652 = -pkin(3) * t722 - qJDD(1) * qJ(4) + t656;
t708 = g(3) + qJDD(3);
t742 = qJD(1) * qJD(4);
t746 = t713 * t708 + 0.2e1 * t711 * t742;
t638 = (pkin(4) * t713 * t722 + pkin(7) * qJDD(1) - t652) * t711 + t746;
t642 = t711 * t708 + (t652 - 0.2e1 * t742) * t713;
t741 = qJDD(1) * t713;
t639 = -pkin(4) * t748 - pkin(7) * t741 + t642;
t635 = t716 * t638 + t719 * t639;
t666 = -pkin(5) * t683 - pkin(8) * t684;
t721 = qJD(5) ^ 2;
t632 = -pkin(5) * t721 + qJDD(5) * pkin(8) + t666 * t683 + t635;
t700 = t711 ^ 2;
t655 = -t712 * t674 + t714 * t679;
t730 = qJDD(1) * pkin(3) + qJDD(4) - t655;
t640 = pkin(4) * t741 + (-qJ(4) + (-t700 - t701) * pkin(7)) * t722 + t730;
t744 = t683 * qJD(5);
t668 = t734 * qJDD(1) + t744;
t633 = (-t668 - t744) * pkin(8) + (-t667 + t743) * pkin(5) + t640;
t715 = sin(qJ(6));
t718 = cos(qJ(6));
t629 = -t632 * t715 + t633 * t718;
t671 = qJD(5) * t718 - t684 * t715;
t649 = qJD(6) * t671 + qJDD(5) * t715 + t668 * t718;
t672 = qJD(5) * t715 + t684 * t718;
t653 = -mrSges(7,1) * t671 + mrSges(7,2) * t672;
t681 = qJD(6) - t683;
t657 = -mrSges(7,2) * t681 + mrSges(7,3) * t671;
t665 = qJDD(6) - t667;
t626 = m(7) * t629 + mrSges(7,1) * t665 - mrSges(7,3) * t649 - t653 * t672 + t657 * t681;
t630 = t632 * t718 + t633 * t715;
t648 = -qJD(6) * t672 + qJDD(5) * t718 - t668 * t715;
t658 = mrSges(7,1) * t681 - mrSges(7,3) * t672;
t627 = m(7) * t630 - mrSges(7,2) * t665 + mrSges(7,3) * t648 + t653 * t671 - t658 * t681;
t620 = -t626 * t715 + t718 * t627;
t663 = -mrSges(6,1) * t683 + mrSges(6,2) * t684;
t676 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t684;
t618 = m(6) * t635 - qJDD(5) * mrSges(6,2) + mrSges(6,3) * t667 - qJD(5) * t676 + t663 * t683 + t620;
t634 = t638 * t719 - t639 * t716;
t631 = -qJDD(5) * pkin(5) - pkin(8) * t721 + t666 * t684 - t634;
t628 = -m(7) * t631 + t648 * mrSges(7,1) - mrSges(7,2) * t649 + t671 * t657 - t658 * t672;
t675 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t683;
t624 = m(6) * t634 + qJDD(5) * mrSges(6,1) - mrSges(6,3) * t668 + qJD(5) * t675 - t663 * t684 + t628;
t613 = t716 * t618 + t719 * t624;
t641 = -t652 * t711 + t746;
t732 = mrSges(5,3) * qJDD(1) + t722 * (-t749 + t750);
t611 = m(5) * t641 + t732 * t711 + t613;
t738 = t719 * t618 - t716 * t624;
t612 = m(5) * t642 - t732 * t713 + t738;
t607 = -t611 * t711 + t713 * t612;
t603 = m(4) * t656 - mrSges(4,1) * t722 + qJDD(1) * mrSges(4,2) + t607;
t651 = -t722 * qJ(4) + t730;
t619 = t718 * t626 + t715 * t627;
t727 = m(6) * t640 - t667 * mrSges(6,1) + t668 * mrSges(6,2) - t683 * t675 + t684 * t676 + t619;
t725 = -m(5) * t651 + qJDD(1) * t749 - t727 + (t700 * t722 + t748) * mrSges(5,3);
t614 = t725 - t722 * mrSges(4,2) + m(4) * t655 + (-mrSges(4,1) - t750) * qJDD(1);
t739 = t714 * t603 - t712 * t614;
t731 = m(3) * t680 + qJDD(1) * mrSges(3,3) + t739;
t596 = m(2) * t689 - qJDD(1) * mrSges(2,2) - t753 * t722 + t731;
t601 = t603 * t712 + t614 * t714;
t682 = -qJDD(1) * pkin(1) + t728;
t600 = m(3) * t682 - qJDD(1) * mrSges(3,1) - t722 * mrSges(3,3) + t601;
t597 = m(2) * t688 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t722 - t600;
t747 = t717 * t596 + t720 * t597;
t735 = -Ifges(5,5) * t711 - Ifges(5,6) * t713;
t745 = t722 * t735;
t740 = t720 * t596 - t597 * t717;
t737 = -Ifges(5,1) * t711 - Ifges(5,4) * t713;
t736 = -Ifges(5,4) * t711 - Ifges(5,2) * t713;
t606 = t713 * t611 + t711 * t612;
t605 = m(4) * t708 + t606;
t644 = Ifges(7,4) * t672 + Ifges(7,2) * t671 + Ifges(7,6) * t681;
t645 = Ifges(7,1) * t672 + Ifges(7,4) * t671 + Ifges(7,5) * t681;
t726 = mrSges(7,1) * t629 - mrSges(7,2) * t630 + Ifges(7,5) * t649 + Ifges(7,6) * t648 + Ifges(7,3) * t665 + t644 * t672 - t645 * t671;
t643 = Ifges(7,5) * t672 + Ifges(7,6) * t671 + Ifges(7,3) * t681;
t621 = -mrSges(7,1) * t631 + mrSges(7,3) * t630 + Ifges(7,4) * t649 + Ifges(7,2) * t648 + Ifges(7,6) * t665 - t643 * t672 + t645 * t681;
t622 = mrSges(7,2) * t631 - mrSges(7,3) * t629 + Ifges(7,1) * t649 + Ifges(7,4) * t648 + Ifges(7,5) * t665 + t643 * t671 - t644 * t681;
t660 = Ifges(6,4) * t684 + Ifges(6,2) * t683 + Ifges(6,6) * qJD(5);
t661 = Ifges(6,1) * t684 + Ifges(6,4) * t683 + Ifges(6,5) * qJD(5);
t724 = mrSges(6,1) * t634 - mrSges(6,2) * t635 + Ifges(6,5) * t668 + Ifges(6,6) * t667 + Ifges(6,3) * qJDD(5) + pkin(5) * t628 + pkin(8) * t620 + t718 * t621 + t715 * t622 + t684 * t660 - t683 * t661;
t659 = Ifges(6,5) * t684 + Ifges(6,6) * t683 + Ifges(6,3) * qJD(5);
t608 = mrSges(6,2) * t640 - mrSges(6,3) * t634 + Ifges(6,1) * t668 + Ifges(6,4) * t667 + Ifges(6,5) * qJDD(5) - pkin(8) * t619 - qJD(5) * t660 - t621 * t715 + t622 * t718 + t659 * t683;
t609 = -mrSges(6,1) * t640 + mrSges(6,3) * t635 + Ifges(6,4) * t668 + Ifges(6,2) * t667 + Ifges(6,6) * qJDD(5) - pkin(5) * t619 + qJD(5) * t661 - t659 * t684 - t726;
t591 = -mrSges(5,1) * t651 + mrSges(5,3) * t642 - pkin(4) * t727 + pkin(7) * t738 + t736 * qJDD(1) + t716 * t608 + t719 * t609 + t711 * t745;
t592 = mrSges(5,2) * t651 - mrSges(5,3) * t641 - pkin(7) * t613 + t737 * qJDD(1) + t719 * t608 - t716 * t609 - t713 * t745;
t615 = mrSges(5,1) * t741 - t725;
t723 = -mrSges(3,1) * t682 - mrSges(4,1) * t655 - mrSges(2,2) * t689 - pkin(2) * t601 + pkin(3) * t615 - qJ(4) * t607 - t591 * t713 - t592 * t711 + qJ(2) * (-mrSges(3,1) * t722 + t731) - pkin(1) * t600 + mrSges(4,2) * t656 + mrSges(3,3) * t680 + mrSges(2,1) * t688 + (Ifges(3,2) + Ifges(2,3) + Ifges(4,3)) * qJDD(1);
t604 = -m(3) * g(3) - t605;
t590 = -t724 - mrSges(4,1) * t708 + mrSges(4,3) * t656 - mrSges(5,1) * t641 + mrSges(5,2) * t642 - pkin(4) * t613 - pkin(3) * t606 + (-Ifges(4,6) - t735) * qJDD(1) + (t711 * t736 - t713 * t737 + Ifges(4,5)) * t722;
t589 = mrSges(4,2) * t708 - mrSges(4,3) * t655 - Ifges(4,5) * qJDD(1) - Ifges(4,6) * t722 - qJ(4) * t606 - t591 * t711 + t592 * t713;
t588 = mrSges(3,2) * t682 - mrSges(2,3) * t688 - qJ(2) * t604 - qJ(3) * t601 + t714 * t589 - t712 * t590 - t751 * t722 + t752 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t587 = mrSges(3,2) * t680 + mrSges(2,3) * t689 - pkin(1) * t604 + pkin(2) * t605 + t753 * g(3) - qJ(3) * t739 + t751 * qJDD(1) - t712 * t589 - t714 * t590 + t752 * t722;
t1 = [-m(1) * g(1) + t740; -m(1) * g(2) + t747; (-m(1) - m(2) - m(3)) * g(3) - t605; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t747 - t717 * t587 + t720 * t588; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t740 + t720 * t587 + t717 * t588; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t723; t723; t600; t605; t615; t724; t726;];
tauJB  = t1;
