% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRRPPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
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
% Datum: 2019-05-05 03:15
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRRPPR4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR4_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR4_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPPR4_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR4_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR4_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR4_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPPR4_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 03:08:54
% EndTime: 2019-05-05 03:09:05
% DurationCPUTime: 7.37s
% Computational Cost: add. (108639->321), mult. (226424->395), div. (0->0), fcn. (151427->12), ass. (0->135)
t729 = -2 * qJD(4);
t728 = Ifges(5,1) + Ifges(6,1);
t720 = Ifges(5,4) - Ifges(6,5);
t727 = Ifges(5,5) + Ifges(6,4);
t726 = -Ifges(5,2) - Ifges(6,3);
t725 = Ifges(6,2) + Ifges(5,3);
t724 = Ifges(5,6) - Ifges(6,6);
t676 = sin(pkin(10));
t678 = cos(pkin(10));
t663 = t676 * g(1) - t678 * g(2);
t664 = -t678 * g(1) - t676 * g(2);
t674 = -g(3) + qJDD(1);
t685 = cos(qJ(2));
t679 = cos(pkin(6));
t682 = sin(qJ(2));
t713 = t679 * t682;
t677 = sin(pkin(6));
t714 = t677 * t682;
t615 = t663 * t713 + t685 * t664 + t674 * t714;
t687 = qJD(2) ^ 2;
t609 = -t687 * pkin(2) + qJDD(2) * pkin(8) + t615;
t640 = -t677 * t663 + t679 * t674;
t681 = sin(qJ(3));
t684 = cos(qJ(3));
t604 = t684 * t609 + t681 * t640;
t659 = (-pkin(3) * t684 - qJ(4) * t681) * qJD(2);
t686 = qJD(3) ^ 2;
t705 = t684 * qJD(2);
t591 = -t686 * pkin(3) + qJDD(3) * qJ(4) + t659 * t705 + t604;
t614 = -t682 * t664 + (t663 * t679 + t674 * t677) * t685;
t608 = -qJDD(2) * pkin(2) - t687 * pkin(8) - t614;
t704 = qJD(2) * qJD(3);
t702 = t684 * t704;
t661 = t681 * qJDD(2) + t702;
t670 = t681 * t704;
t662 = t684 * qJDD(2) - t670;
t595 = (-t661 - t702) * qJ(4) + (-t662 + t670) * pkin(3) + t608;
t675 = sin(pkin(11));
t707 = qJD(2) * t681;
t717 = cos(pkin(11));
t652 = t675 * qJD(3) + t717 * t707;
t585 = -t675 * t591 + t717 * t595 + t652 * t729;
t638 = t675 * qJDD(3) + t717 * t661;
t603 = -t681 * t609 + t684 * t640;
t690 = qJDD(3) * pkin(3) + t686 * qJ(4) - t659 * t707 - qJDD(4) + t603;
t651 = -t717 * qJD(3) + t675 * t707;
t703 = t651 * t705;
t723 = -(t638 + t703) * qJ(5) - t690;
t722 = -2 * qJD(5);
t721 = -mrSges(5,3) - mrSges(6,2);
t586 = t717 * t591 + t675 * t595 + t651 * t729;
t635 = -mrSges(5,1) * t705 - t652 * mrSges(5,3);
t636 = mrSges(6,1) * t705 + t652 * mrSges(6,2);
t637 = -t717 * qJDD(3) + t675 * t661;
t626 = t651 * pkin(4) - t652 * qJ(5);
t715 = t684 ^ 2 * t687;
t581 = -pkin(4) * t715 - t662 * qJ(5) - t651 * t626 + t705 * t722 + t586;
t583 = t662 * pkin(4) - qJ(5) * t715 + t652 * t626 + qJDD(5) - t585;
t578 = (-t638 + t703) * pkin(9) + (t651 * t652 + t662) * pkin(5) + t583;
t639 = pkin(5) * t705 - t652 * pkin(9);
t649 = t651 ^ 2;
t579 = -t649 * pkin(5) + t637 * pkin(9) - t639 * t705 + t581;
t680 = sin(qJ(6));
t683 = cos(qJ(6));
t576 = t683 * t578 - t680 * t579;
t624 = t683 * t651 - t680 * t652;
t597 = t624 * qJD(6) + t680 * t637 + t683 * t638;
t625 = t680 * t651 + t683 * t652;
t605 = -t624 * mrSges(7,1) + t625 * mrSges(7,2);
t668 = qJD(6) + t705;
t610 = -t668 * mrSges(7,2) + t624 * mrSges(7,3);
t656 = qJDD(6) + t662;
t573 = m(7) * t576 + t656 * mrSges(7,1) - t597 * mrSges(7,3) - t625 * t605 + t668 * t610;
t577 = t680 * t578 + t683 * t579;
t596 = -t625 * qJD(6) + t683 * t637 - t680 * t638;
t611 = t668 * mrSges(7,1) - t625 * mrSges(7,3);
t574 = m(7) * t577 - t656 * mrSges(7,2) + t596 * mrSges(7,3) + t624 * t605 - t668 * t611;
t698 = -t680 * t573 + t683 * t574;
t693 = m(6) * t581 - t662 * mrSges(6,3) + t698;
t627 = t651 * mrSges(6,1) - t652 * mrSges(6,3);
t708 = -t651 * mrSges(5,1) - t652 * mrSges(5,2) - t627;
t565 = m(5) * t586 + t662 * mrSges(5,2) + t708 * t651 + t721 * t637 + (t635 - t636) * t705 + t693;
t633 = -t651 * mrSges(6,2) - mrSges(6,3) * t705;
t634 = mrSges(5,2) * t705 - t651 * mrSges(5,3);
t567 = t683 * t573 + t680 * t574;
t691 = -m(6) * t583 - t662 * mrSges(6,1) - t567;
t566 = m(5) * t585 - t662 * mrSges(5,1) + t708 * t652 + t721 * t638 + (-t633 - t634) * t705 + t691;
t563 = t675 * t565 + t717 * t566;
t665 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t707;
t666 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t705;
t689 = -m(4) * t608 + t662 * mrSges(4,1) - t661 * mrSges(4,2) - t665 * t707 + t666 * t705 - t563;
t559 = m(3) * t614 + qJDD(2) * mrSges(3,1) - t687 * mrSges(3,2) + t689;
t716 = t559 * t685;
t660 = (-mrSges(4,1) * t684 + mrSges(4,2) * t681) * qJD(2);
t699 = t717 * t565 - t675 * t566;
t562 = m(4) * t604 - qJDD(3) * mrSges(4,2) + t662 * mrSges(4,3) - qJD(3) * t665 + t660 * t705 + t699;
t587 = t652 * t722 + (-t652 * t705 + t637) * pkin(4) + t723;
t584 = -t649 * pkin(9) + (-pkin(4) - pkin(5)) * t637 + (pkin(4) * t705 + (2 * qJD(5)) + t639) * t652 - t723;
t695 = -m(7) * t584 + t596 * mrSges(7,1) - t597 * mrSges(7,2) + t624 * t610 - t625 * t611;
t575 = m(6) * t587 + t637 * mrSges(6,1) - t638 * mrSges(6,3) + t651 * t633 - t652 * t636 + t695;
t688 = m(5) * t690 - t637 * mrSges(5,1) - t638 * mrSges(5,2) - t651 * t634 - t652 * t635 - t575;
t571 = m(4) * t603 + qJDD(3) * mrSges(4,1) - t661 * mrSges(4,3) + qJD(3) * t666 - t660 * t707 + t688;
t700 = t684 * t562 - t681 * t571;
t552 = m(3) * t615 - t687 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t700;
t556 = t681 * t562 + t684 * t571;
t555 = m(3) * t640 + t556;
t542 = t552 * t713 - t677 * t555 + t679 * t716;
t540 = m(2) * t663 + t542;
t547 = t685 * t552 - t682 * t559;
t546 = m(2) * t664 + t547;
t712 = t678 * t540 + t676 * t546;
t711 = t726 * t651 + t720 * t652 - t724 * t705;
t710 = t724 * t651 - t727 * t652 + t725 * t705;
t709 = t720 * t651 - t728 * t652 + t727 * t705;
t541 = t552 * t714 + t679 * t555 + t677 * t716;
t701 = -t676 * t540 + t678 * t546;
t598 = Ifges(7,5) * t625 + Ifges(7,6) * t624 + Ifges(7,3) * t668;
t600 = Ifges(7,1) * t625 + Ifges(7,4) * t624 + Ifges(7,5) * t668;
t568 = -mrSges(7,1) * t584 + mrSges(7,3) * t577 + Ifges(7,4) * t597 + Ifges(7,2) * t596 + Ifges(7,6) * t656 - t625 * t598 + t668 * t600;
t599 = Ifges(7,4) * t625 + Ifges(7,2) * t624 + Ifges(7,6) * t668;
t569 = mrSges(7,2) * t584 - mrSges(7,3) * t576 + Ifges(7,1) * t597 + Ifges(7,4) * t596 + Ifges(7,5) * t656 + t624 * t598 - t668 * t599;
t548 = mrSges(5,1) * t690 - mrSges(6,1) * t587 + mrSges(6,2) * t581 + mrSges(5,3) * t586 - pkin(4) * t575 - pkin(5) * t695 - pkin(9) * t698 - t683 * t568 - t680 * t569 + t726 * t637 + t720 * t638 + t710 * t652 - t662 * t724 + t709 * t705;
t553 = -mrSges(5,2) * t690 + mrSges(6,2) * t583 - mrSges(5,3) * t585 - mrSges(6,3) * t587 - pkin(9) * t567 - qJ(5) * t575 - t680 * t568 + t683 * t569 - t720 * t637 + t728 * t638 + t710 * t651 - t662 * t727 + t711 * t705;
t645 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t681 + Ifges(4,6) * t684) * qJD(2);
t646 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t681 + Ifges(4,2) * t684) * qJD(2);
t538 = mrSges(4,2) * t608 - mrSges(4,3) * t603 + Ifges(4,1) * t661 + Ifges(4,4) * t662 + Ifges(4,5) * qJDD(3) - qJ(4) * t563 - qJD(3) * t646 - t675 * t548 + t717 * t553 + t645 * t705;
t647 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t681 + Ifges(4,4) * t684) * qJD(2);
t543 = Ifges(7,3) * t656 + Ifges(4,4) * t661 + Ifges(7,6) * t596 + Ifges(7,5) * t597 + mrSges(4,3) * t604 - mrSges(4,1) * t608 + qJD(3) * t647 - t624 * t600 + t625 * t599 + Ifges(4,6) * qJDD(3) - mrSges(7,2) * t577 - mrSges(6,3) * t581 + mrSges(6,1) * t583 - mrSges(5,1) * t585 + mrSges(5,2) * t586 + mrSges(7,1) * t576 - pkin(3) * t563 + pkin(5) * t567 - qJ(5) * (-t636 * t705 + t693) - pkin(4) * (-t633 * t705 + t691) - t645 * t707 + (qJ(5) * t627 + t709) * t651 + (pkin(4) * t627 - t711) * t652 + (qJ(5) * mrSges(6,2) + t724) * t637 + (pkin(4) * mrSges(6,2) - t727) * t638 + (Ifges(4,2) + t725) * t662;
t536 = mrSges(3,2) * t640 - mrSges(3,3) * t614 + Ifges(3,5) * qJDD(2) - t687 * Ifges(3,6) - pkin(8) * t556 + t684 * t538 - t681 * t543;
t537 = Ifges(3,6) * qJDD(2) + t687 * Ifges(3,5) - mrSges(3,1) * t640 + mrSges(3,3) * t615 - Ifges(4,5) * t661 - Ifges(4,6) * t662 - Ifges(4,3) * qJDD(3) - mrSges(4,1) * t603 + mrSges(4,2) * t604 - t675 * t553 - t717 * t548 - pkin(3) * t688 - qJ(4) * t699 - pkin(2) * t556 + (-t681 * t646 + t684 * t647) * qJD(2);
t692 = pkin(7) * t547 + t536 * t682 + t537 * t685;
t535 = mrSges(3,1) * t614 - mrSges(3,2) * t615 + Ifges(3,3) * qJDD(2) + pkin(2) * t689 + pkin(8) * t700 + t681 * t538 + t684 * t543;
t534 = mrSges(2,2) * t674 - mrSges(2,3) * t663 + t685 * t536 - t682 * t537 + (-t541 * t677 - t542 * t679) * pkin(7);
t533 = -mrSges(2,1) * t674 + mrSges(2,3) * t664 - pkin(1) * t541 - t677 * t535 + t692 * t679;
t1 = [-m(1) * g(1) + t701; -m(1) * g(2) + t712; -m(1) * g(3) + m(2) * t674 + t541; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t712 - t676 * t533 + t678 * t534; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t701 + t678 * t533 + t676 * t534; -mrSges(1,1) * g(2) + mrSges(2,1) * t663 + mrSges(1,2) * g(1) - mrSges(2,2) * t664 + pkin(1) * t542 + t679 * t535 + t692 * t677;];
tauB  = t1;
