% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PRRRR10
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PRRRR10_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR10_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR10_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR10_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR10_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR10_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR10_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR10_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR10_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:24:24
% EndTime: 2019-12-05 17:24:37
% DurationCPUTime: 11.51s
% Computational Cost: add. (193047->279), mult. (396496->377), div. (0->0), fcn. (310398->14), ass. (0->131)
t702 = sin(pkin(6));
t709 = sin(qJ(3));
t713 = cos(qJ(3));
t728 = qJD(2) * qJD(3);
t686 = (-qJDD(2) * t713 + t709 * t728) * t702;
t701 = sin(pkin(11));
t704 = cos(pkin(11));
t692 = t701 * g(1) - t704 * g(2);
t693 = -t704 * g(1) - t701 * g(2);
t700 = -g(3) + qJDD(1);
t710 = sin(qJ(2));
t706 = cos(pkin(5));
t714 = cos(qJ(2));
t731 = t706 * t714;
t703 = sin(pkin(5));
t734 = t703 * t714;
t664 = t692 * t731 - t710 * t693 + t700 * t734;
t715 = qJD(2) ^ 2;
t738 = pkin(8) * t702;
t660 = qJDD(2) * pkin(2) + t715 * t738 + t664;
t732 = t706 * t710;
t735 = t703 * t710;
t665 = t692 * t732 + t714 * t693 + t700 * t735;
t661 = -t715 * pkin(2) + qJDD(2) * t738 + t665;
t679 = -t703 * t692 + t706 * t700;
t705 = cos(pkin(6));
t636 = -t709 * t661 + (t660 * t705 + t679 * t702) * t713;
t698 = t705 * qJD(2) + qJD(3);
t729 = qJD(2) * t702;
t726 = t713 * t729;
t682 = -t698 * mrSges(4,2) + mrSges(4,3) * t726;
t683 = (-mrSges(4,1) * t713 + mrSges(4,2) * t709) * t729;
t685 = (qJDD(2) * t709 + t713 * t728) * t702;
t697 = t705 * qJDD(2) + qJDD(3);
t733 = t705 * t709;
t736 = t702 * t709;
t637 = t660 * t733 + t713 * t661 + t679 * t736;
t684 = (-pkin(3) * t713 - pkin(9) * t709) * t729;
t696 = t698 ^ 2;
t633 = -t696 * pkin(3) + t697 * pkin(9) + t684 * t726 + t637;
t675 = t705 * t679;
t635 = t686 * pkin(3) - t685 * pkin(9) + t675 + (-t660 + (pkin(3) * t709 - pkin(9) * t713) * t698 * qJD(2)) * t702;
t708 = sin(qJ(4));
t712 = cos(qJ(4));
t630 = t712 * t633 + t708 * t635;
t727 = t709 * t729;
t677 = t712 * t698 - t708 * t727;
t678 = t708 * t698 + t712 * t727;
t663 = -t677 * pkin(4) - t678 * pkin(10);
t680 = qJDD(4) + t686;
t691 = qJD(4) - t726;
t690 = t691 ^ 2;
t627 = -t690 * pkin(4) + t680 * pkin(10) + t677 * t663 + t630;
t632 = -t697 * pkin(3) - t696 * pkin(9) + t684 * t727 - t636;
t655 = -t678 * qJD(4) - t708 * t685 + t712 * t697;
t656 = t677 * qJD(4) + t712 * t685 + t708 * t697;
t628 = (-t677 * t691 - t656) * pkin(10) + (t678 * t691 - t655) * pkin(4) + t632;
t707 = sin(qJ(5));
t711 = cos(qJ(5));
t624 = -t707 * t627 + t711 * t628;
t666 = -t707 * t678 + t711 * t691;
t640 = t666 * qJD(5) + t711 * t656 + t707 * t680;
t667 = t711 * t678 + t707 * t691;
t645 = -t666 * mrSges(6,1) + t667 * mrSges(6,2);
t676 = qJD(5) - t677;
t647 = -t676 * mrSges(6,2) + t666 * mrSges(6,3);
t653 = qJDD(5) - t655;
t621 = m(6) * t624 + t653 * mrSges(6,1) - t640 * mrSges(6,3) - t667 * t645 + t676 * t647;
t625 = t711 * t627 + t707 * t628;
t639 = -t667 * qJD(5) - t707 * t656 + t711 * t680;
t648 = t676 * mrSges(6,1) - t667 * mrSges(6,3);
t622 = m(6) * t625 - t653 * mrSges(6,2) + t639 * mrSges(6,3) + t666 * t645 - t676 * t648;
t614 = t711 * t621 + t707 * t622;
t668 = -t691 * mrSges(5,2) + t677 * mrSges(5,3);
t669 = t691 * mrSges(5,1) - t678 * mrSges(5,3);
t718 = -m(5) * t632 + t655 * mrSges(5,1) - t656 * mrSges(5,2) + t677 * t668 - t678 * t669 - t614;
t610 = m(4) * t636 + t697 * mrSges(4,1) - t685 * mrSges(4,3) + t698 * t682 - t683 * t727 + t718;
t737 = t610 * t713;
t681 = t698 * mrSges(4,1) - mrSges(4,3) * t727;
t615 = -t707 * t621 + t711 * t622;
t662 = -t677 * mrSges(5,1) + t678 * mrSges(5,2);
t613 = m(5) * t630 - t680 * mrSges(5,2) + t655 * mrSges(5,3) + t677 * t662 - t691 * t669 + t615;
t629 = -t708 * t633 + t712 * t635;
t626 = -t680 * pkin(4) - t690 * pkin(10) + t678 * t663 - t629;
t623 = -m(6) * t626 + t639 * mrSges(6,1) - t640 * mrSges(6,2) + t666 * t647 - t667 * t648;
t619 = m(5) * t629 + t680 * mrSges(5,1) - t656 * mrSges(5,3) - t678 * t662 + t691 * t668 + t623;
t724 = t712 * t613 - t708 * t619;
t604 = m(4) * t637 - t697 * mrSges(4,2) - t686 * mrSges(4,3) - t698 * t681 + t683 * t726 + t724;
t607 = t708 * t613 + t712 * t619;
t646 = -t702 * t660 + t675;
t606 = m(4) * t646 + t686 * mrSges(4,1) + t685 * mrSges(4,2) + (t681 * t709 - t682 * t713) * t729 + t607;
t593 = t604 * t733 - t702 * t606 + t705 * t737;
t589 = m(3) * t664 + qJDD(2) * mrSges(3,1) - t715 * mrSges(3,2) + t593;
t592 = t604 * t736 + t705 * t606 + t702 * t737;
t591 = m(3) * t679 + t592;
t598 = t713 * t604 - t709 * t610;
t597 = m(3) * t665 - t715 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t598;
t579 = t589 * t731 - t703 * t591 + t597 * t732;
t577 = m(2) * t692 + t579;
t583 = -t710 * t589 + t714 * t597;
t582 = m(2) * t693 + t583;
t730 = t704 * t577 + t701 * t582;
t578 = t589 * t734 + t706 * t591 + t597 * t735;
t725 = -t701 * t577 + t704 * t582;
t723 = m(2) * t700 + t578;
t641 = Ifges(6,5) * t667 + Ifges(6,6) * t666 + Ifges(6,3) * t676;
t643 = Ifges(6,1) * t667 + Ifges(6,4) * t666 + Ifges(6,5) * t676;
t616 = -mrSges(6,1) * t626 + mrSges(6,3) * t625 + Ifges(6,4) * t640 + Ifges(6,2) * t639 + Ifges(6,6) * t653 - t667 * t641 + t676 * t643;
t642 = Ifges(6,4) * t667 + Ifges(6,2) * t666 + Ifges(6,6) * t676;
t617 = mrSges(6,2) * t626 - mrSges(6,3) * t624 + Ifges(6,1) * t640 + Ifges(6,4) * t639 + Ifges(6,5) * t653 + t666 * t641 - t676 * t642;
t649 = Ifges(5,5) * t678 + Ifges(5,6) * t677 + Ifges(5,3) * t691;
t650 = Ifges(5,4) * t678 + Ifges(5,2) * t677 + Ifges(5,6) * t691;
t599 = mrSges(5,2) * t632 - mrSges(5,3) * t629 + Ifges(5,1) * t656 + Ifges(5,4) * t655 + Ifges(5,5) * t680 - pkin(10) * t614 - t707 * t616 + t711 * t617 + t677 * t649 - t691 * t650;
t651 = Ifges(5,1) * t678 + Ifges(5,4) * t677 + Ifges(5,5) * t691;
t717 = mrSges(6,1) * t624 - mrSges(6,2) * t625 + Ifges(6,5) * t640 + Ifges(6,6) * t639 + Ifges(6,3) * t653 + t667 * t642 - t666 * t643;
t600 = -mrSges(5,1) * t632 + mrSges(5,3) * t630 + Ifges(5,4) * t656 + Ifges(5,2) * t655 + Ifges(5,6) * t680 - pkin(4) * t614 - t678 * t649 + t691 * t651 - t717;
t672 = Ifges(4,6) * t698 + (Ifges(4,4) * t709 + Ifges(4,2) * t713) * t729;
t673 = Ifges(4,5) * t698 + (Ifges(4,1) * t709 + Ifges(4,4) * t713) * t729;
t584 = Ifges(4,5) * t685 - Ifges(4,6) * t686 + Ifges(4,3) * t697 + mrSges(4,1) * t636 - mrSges(4,2) * t637 + t708 * t599 + t712 * t600 + pkin(3) * t718 + pkin(9) * t724 + (t672 * t709 - t673 * t713) * t729;
t671 = Ifges(4,3) * t698 + (Ifges(4,5) * t709 + Ifges(4,6) * t713) * t729;
t585 = mrSges(4,2) * t646 - mrSges(4,3) * t636 + Ifges(4,1) * t685 - Ifges(4,4) * t686 + Ifges(4,5) * t697 - pkin(9) * t607 + t712 * t599 - t708 * t600 + t671 * t726 - t698 * t672;
t716 = mrSges(5,1) * t629 - mrSges(5,2) * t630 + Ifges(5,5) * t656 + Ifges(5,6) * t655 + Ifges(5,3) * t680 + pkin(4) * t623 + pkin(10) * t615 + t711 * t616 + t707 * t617 + t678 * t650 - t677 * t651;
t586 = -mrSges(4,1) * t646 + mrSges(4,3) * t637 + Ifges(4,4) * t685 - Ifges(4,2) * t686 + Ifges(4,6) * t697 - pkin(3) * t607 - t671 * t727 + t698 * t673 - t716;
t719 = pkin(8) * t598 + t585 * t709 + t586 * t713;
t574 = -mrSges(3,1) * t679 + mrSges(3,3) * t665 + t715 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t592 - t702 * t584 + t719 * t705;
t575 = mrSges(3,2) * t679 - mrSges(3,3) * t664 + Ifges(3,5) * qJDD(2) - t715 * Ifges(3,6) + t713 * t585 - t709 * t586 + (-t592 * t702 - t593 * t705) * pkin(8);
t720 = pkin(7) * t583 + t574 * t714 + t575 * t710;
t573 = mrSges(3,1) * t664 - mrSges(3,2) * t665 + Ifges(3,3) * qJDD(2) + pkin(2) * t593 + t705 * t584 + t719 * t702;
t572 = mrSges(2,2) * t700 - mrSges(2,3) * t692 - t710 * t574 + t714 * t575 + (-t578 * t703 - t579 * t706) * pkin(7);
t571 = -mrSges(2,1) * t700 + mrSges(2,3) * t693 - pkin(1) * t578 - t703 * t573 + t720 * t706;
t1 = [-m(1) * g(1) + t725; -m(1) * g(2) + t730; -m(1) * g(3) + t723; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t730 - t701 * t571 + t704 * t572; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t725 + t704 * t571 + t701 * t572; -mrSges(1,1) * g(2) + mrSges(2,1) * t692 + mrSges(1,2) * g(1) - mrSges(2,2) * t693 + pkin(1) * t579 + t706 * t573 + t720 * t703; t723; t573; t584; t716; t717;];
tauJB = t1;
