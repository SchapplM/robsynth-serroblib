% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PRRRR11
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
%   pkin=[a2,a3,a4,a5,alpha3,d2,d3,d4,d5,theta1]';
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
% Datum: 2024-09-27 21:46
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PRRRR11_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR11_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR11_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR11_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR11_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR11_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR11_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR11_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR11_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 21:45:41
% EndTime: 2024-09-27 21:45:47
% DurationCPUTime: 3.11s
% Computational Cost: add. (145490->271), mult. (361971->361), div. (0->0), fcn. (274631->12), ass. (0->122)
t737 = sin(pkin(5));
t771 = pkin(7) * t737;
t748 = qJD(2) ^ 2;
t770 = t737 ^ 2 * t748;
t742 = sin(qJ(3));
t769 = t737 * t742;
t746 = cos(qJ(3));
t768 = t737 * t746;
t739 = cos(pkin(5));
t767 = t739 * t742;
t766 = t739 * t746;
t736 = sin(pkin(10));
t738 = cos(pkin(10));
t719 = t736 * g(1) - t738 * g(2);
t720 = -t738 * g(1) - t736 * g(2);
t743 = sin(qJ(2));
t747 = cos(qJ(2));
t700 = t747 * t719 - t743 * t720;
t691 = qJDD(2) * pkin(2) + t748 * t771 + t700;
t701 = t743 * t719 + t747 * t720;
t692 = -t748 * pkin(2) + qJDD(2) * t771 + t701;
t735 = -g(3) + qJDD(1);
t670 = t691 * t766 - t742 * t692 + t735 * t768;
t763 = qJD(2) * qJD(3);
t710 = (qJDD(2) * t742 + t746 * t763) * t737;
t725 = t739 * qJDD(2) + qJDD(3);
t726 = t739 * qJD(2) + qJD(3);
t764 = qJD(2) * t737;
t760 = t746 * t764;
t661 = (t726 * t760 - t710) * pkin(8) + (t742 * t746 * t770 + t725) * pkin(3) + t670;
t671 = t691 * t767 + t746 * t692 + t735 * t769;
t761 = t742 * t764;
t709 = t726 * pkin(3) - pkin(8) * t761;
t711 = (qJDD(2) * t746 - t742 * t763) * t737;
t762 = t746 ^ 2 * t770;
t662 = -pkin(3) * t762 + t711 * pkin(8) - t726 * t709 + t671;
t741 = sin(qJ(4));
t745 = cos(qJ(4));
t648 = t745 * t661 - t741 * t662;
t704 = (-t741 * t742 + t745 * t746) * t764;
t677 = t704 * qJD(4) + t745 * t710 + t741 * t711;
t705 = (t741 * t746 + t742 * t745) * t764;
t723 = qJDD(4) + t725;
t724 = qJD(4) + t726;
t645 = (t704 * t724 - t677) * pkin(9) + (t704 * t705 + t723) * pkin(4) + t648;
t649 = t741 * t661 + t745 * t662;
t676 = -t705 * qJD(4) - t741 * t710 + t745 * t711;
t695 = t724 * pkin(4) - t705 * pkin(9);
t703 = t704 ^ 2;
t646 = -t703 * pkin(4) + t676 * pkin(9) - t724 * t695 + t649;
t740 = sin(qJ(5));
t744 = cos(qJ(5));
t643 = t744 * t645 - t740 * t646;
t684 = t744 * t704 - t740 * t705;
t656 = t684 * qJD(5) + t740 * t676 + t744 * t677;
t685 = t740 * t704 + t744 * t705;
t667 = -t684 * mrSges(6,1) + t685 * mrSges(6,2);
t718 = qJD(5) + t724;
t678 = -t718 * mrSges(6,2) + t684 * mrSges(6,3);
t717 = qJDD(5) + t723;
t640 = m(6) * t643 + t717 * mrSges(6,1) - t656 * mrSges(6,3) - t685 * t667 + t718 * t678;
t644 = t740 * t645 + t744 * t646;
t655 = -t685 * qJD(5) + t744 * t676 - t740 * t677;
t679 = t718 * mrSges(6,1) - t685 * mrSges(6,3);
t641 = m(6) * t644 - t717 * mrSges(6,2) + t655 * mrSges(6,3) + t684 * t667 - t718 * t679;
t632 = t744 * t640 + t740 * t641;
t687 = -t704 * mrSges(5,1) + t705 * mrSges(5,2);
t693 = -t724 * mrSges(5,2) + t704 * mrSges(5,3);
t629 = m(5) * t648 + t723 * mrSges(5,1) - t677 * mrSges(5,3) - t705 * t687 + t724 * t693 + t632;
t694 = t724 * mrSges(5,1) - t705 * mrSges(5,3);
t756 = -t740 * t640 + t744 * t641;
t630 = m(5) * t649 - t723 * mrSges(5,2) + t676 * mrSges(5,3) + t704 * t687 - t724 * t694 + t756;
t625 = t745 * t629 + t741 * t630;
t707 = -t726 * mrSges(4,2) + mrSges(4,3) * t760;
t708 = (-mrSges(4,1) * t746 + mrSges(4,2) * t742) * t764;
t623 = m(4) * t670 + t725 * mrSges(4,1) - t710 * mrSges(4,3) + t726 * t707 - t708 * t761 + t625;
t706 = t726 * mrSges(4,1) - mrSges(4,3) * t761;
t757 = -t741 * t629 + t745 * t630;
t624 = m(4) * t671 - t725 * mrSges(4,2) + t711 * mrSges(4,3) - t726 * t706 + t708 * t760 + t757;
t686 = -t737 * t691 + t739 * t735;
t669 = -t711 * pkin(3) - pkin(8) * t762 + t709 * t761 + t686;
t651 = -t676 * pkin(4) - t703 * pkin(9) + t705 * t695 + t669;
t753 = m(6) * t651 - t655 * mrSges(6,1) + t656 * mrSges(6,2) - t684 * t678 + t685 * t679;
t750 = m(5) * t669 - t676 * mrSges(5,1) + t677 * mrSges(5,2) - t704 * t693 + t705 * t694 + t753;
t636 = t710 * mrSges(4,2) - t711 * mrSges(4,1) + m(4) * t686 + (t706 * t742 - t707 * t746) * t764 + t750;
t608 = t623 * t766 + t624 * t767 - t737 * t636;
t605 = m(3) * t700 + qJDD(2) * mrSges(3,1) - t748 * mrSges(3,2) + t608;
t613 = -t742 * t623 + t746 * t624;
t611 = m(3) * t701 - t748 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t613;
t599 = t747 * t605 + t743 * t611;
t597 = m(2) * t719 + t599;
t758 = -t743 * t605 + t747 * t611;
t598 = m(2) * t720 + t758;
t765 = t738 * t597 + t736 * t598;
t607 = t623 * t768 + t624 * t769 + t739 * t636;
t759 = -t736 * t597 + t738 * t598;
t755 = m(3) * t735 + t607;
t754 = m(2) * t735 + t755;
t663 = Ifges(6,5) * t685 + Ifges(6,6) * t684 + Ifges(6,3) * t718;
t665 = Ifges(6,1) * t685 + Ifges(6,4) * t684 + Ifges(6,5) * t718;
t633 = -mrSges(6,1) * t651 + mrSges(6,3) * t644 + Ifges(6,4) * t656 + Ifges(6,2) * t655 + Ifges(6,6) * t717 - t685 * t663 + t718 * t665;
t664 = Ifges(6,4) * t685 + Ifges(6,2) * t684 + Ifges(6,6) * t718;
t634 = mrSges(6,2) * t651 - mrSges(6,3) * t643 + Ifges(6,1) * t656 + Ifges(6,4) * t655 + Ifges(6,5) * t717 + t684 * t663 - t718 * t664;
t680 = Ifges(5,5) * t705 + Ifges(5,6) * t704 + Ifges(5,3) * t724;
t682 = Ifges(5,1) * t705 + Ifges(5,4) * t704 + Ifges(5,5) * t724;
t616 = -mrSges(5,1) * t669 + mrSges(5,3) * t649 + Ifges(5,4) * t677 + Ifges(5,2) * t676 + Ifges(5,6) * t723 - pkin(4) * t753 + pkin(9) * t756 + t744 * t633 + t740 * t634 - t705 * t680 + t724 * t682;
t681 = Ifges(5,4) * t705 + Ifges(5,2) * t704 + Ifges(5,6) * t724;
t617 = mrSges(5,2) * t669 - mrSges(5,3) * t648 + Ifges(5,1) * t677 + Ifges(5,4) * t676 + Ifges(5,5) * t723 - pkin(9) * t632 - t740 * t633 + t744 * t634 + t704 * t680 - t724 * t681;
t697 = Ifges(4,3) * t726 + (Ifges(4,5) * t742 + Ifges(4,6) * t746) * t764;
t699 = Ifges(4,5) * t726 + (Ifges(4,1) * t742 + Ifges(4,4) * t746) * t764;
t601 = -mrSges(4,1) * t686 + mrSges(4,3) * t671 + Ifges(4,4) * t710 + Ifges(4,2) * t711 + Ifges(4,6) * t725 - pkin(3) * t750 + pkin(8) * t757 + t745 * t616 + t741 * t617 - t697 * t761 + t726 * t699;
t698 = Ifges(4,6) * t726 + (Ifges(4,4) * t742 + Ifges(4,2) * t746) * t764;
t603 = mrSges(4,2) * t686 - mrSges(4,3) * t670 + Ifges(4,1) * t710 + Ifges(4,4) * t711 + Ifges(4,5) * t725 - pkin(8) * t625 - t741 * t616 + t745 * t617 + t697 * t760 - t726 * t698;
t751 = mrSges(6,1) * t643 - mrSges(6,2) * t644 + Ifges(6,5) * t656 + Ifges(6,6) * t655 + Ifges(6,3) * t717 + t685 * t664 - t684 * t665;
t749 = mrSges(5,1) * t648 - mrSges(5,2) * t649 + Ifges(5,5) * t677 + Ifges(5,6) * t676 + Ifges(5,3) * t723 + pkin(4) * t632 + t705 * t681 - t704 * t682 + t751;
t615 = (t698 * t742 - t699 * t746) * t764 + Ifges(4,3) * t725 + Ifges(4,5) * t710 + Ifges(4,6) * t711 + mrSges(4,1) * t670 - mrSges(4,2) * t671 + pkin(3) * t625 + t749;
t752 = mrSges(3,1) * t700 - mrSges(3,2) * t701 + Ifges(3,3) * qJDD(2) + pkin(2) * t608 + t601 * t768 + t603 * t769 + t613 * t771 + t739 * t615;
t593 = mrSges(3,2) * t735 - mrSges(3,3) * t700 + Ifges(3,5) * qJDD(2) - t748 * Ifges(3,6) - t742 * t601 + t746 * t603 + (-t607 * t737 - t608 * t739) * pkin(7);
t592 = -mrSges(3,1) * t735 + mrSges(3,3) * t701 + t748 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t607 - t737 * t615 + (pkin(7) * t613 + t601 * t746 + t603 * t742) * t739;
t591 = mrSges(2,2) * t735 - mrSges(2,3) * t719 - pkin(6) * t599 - t743 * t592 + t747 * t593;
t590 = -mrSges(2,1) * t735 + mrSges(2,3) * t720 - pkin(1) * t755 + pkin(6) * t758 + t747 * t592 + t743 * t593;
t1 = [-m(1) * g(1) + t759; -m(1) * g(2) + t765; -m(1) * g(3) + t754; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t765 - t736 * t590 + t738 * t591; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t759 + t738 * t590 + t736 * t591; -mrSges(1,1) * g(2) + mrSges(2,1) * t719 + mrSges(1,2) * g(1) - mrSges(2,2) * t720 + pkin(1) * t599 + t752; t754; t752; t615; t749; t751;];
tauJB = t1;
