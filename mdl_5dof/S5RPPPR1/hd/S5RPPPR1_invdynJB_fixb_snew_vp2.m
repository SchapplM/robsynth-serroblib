% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPPPR1
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% Datum: 2019-12-05 17:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPPPR1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR1_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR1_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR1_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR1_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR1_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:29:01
% EndTime: 2019-12-05 17:29:04
% DurationCPUTime: 3.61s
% Computational Cost: add. (35992->234), mult. (83592->324), div. (0->0), fcn. (50610->10), ass. (0->118)
t686 = sin(qJ(1));
t688 = cos(qJ(1));
t661 = t688 * g(2) + t686 * g(3);
t655 = qJDD(1) * pkin(1) + t661;
t660 = t686 * g(2) - t688 * g(3);
t689 = qJD(1) ^ 2;
t656 = -t689 * pkin(1) + t660;
t681 = sin(pkin(7));
t684 = cos(pkin(7));
t638 = t681 * t655 + t684 * t656;
t730 = -t689 * pkin(2) + qJDD(1) * qJ(3) + (2 * qJD(1) * qJD(3)) + t638;
t637 = t684 * t655 - t681 * t656;
t697 = -t689 * qJ(3) + qJDD(3) - t637;
t680 = sin(pkin(8));
t683 = cos(pkin(8));
t705 = -pkin(3) * t683 - qJ(4) * t680;
t723 = t680 * qJD(1);
t729 = (-pkin(2) + t705) * qJDD(1) + t697 - 0.2e1 * qJD(4) * t723;
t678 = -g(1) + qJDD(2);
t622 = t683 * t678 - t730 * t680;
t728 = mrSges(4,2) * t680;
t675 = t680 ^ 2;
t727 = t675 * t689;
t679 = sin(pkin(9));
t726 = t679 * t680;
t682 = cos(pkin(9));
t725 = t680 * t682;
t623 = t680 * t678 + t730 * t683;
t653 = (-mrSges(4,1) * t683 + t728) * qJD(1);
t652 = t705 * qJD(1);
t722 = t683 * qJD(1);
t616 = t652 * t722 + t623;
t702 = -pkin(4) * t683 - pkin(6) * t725;
t724 = t729 * t682;
t609 = t702 * qJDD(1) + (-t616 + (-pkin(4) * t675 * t682 + pkin(6) * t680 * t683) * t689) * t679 + t724;
t612 = t682 * t616 + t729 * t679;
t651 = t702 * qJD(1);
t718 = t679 ^ 2 * t727;
t720 = qJDD(1) * t680;
t610 = -t679 * pkin(6) * t720 - pkin(4) * t718 + t651 * t722 + t612;
t685 = sin(qJ(5));
t687 = cos(qJ(5));
t607 = t687 * t609 - t685 * t610;
t699 = (-t679 * t687 - t682 * t685) * t680;
t642 = qJD(1) * t699;
t698 = (-t679 * t685 + t682 * t687) * t680;
t643 = qJD(1) * t698;
t626 = -t642 * mrSges(6,1) + t643 * mrSges(6,2);
t629 = t642 * qJD(5) + qJDD(1) * t698;
t663 = qJD(5) - t722;
t635 = -t663 * mrSges(6,2) + t642 * mrSges(6,3);
t719 = t683 * qJDD(1);
t662 = qJDD(5) - t719;
t604 = m(6) * t607 + t662 * mrSges(6,1) - t629 * mrSges(6,3) - t643 * t626 + t663 * t635;
t608 = t685 * t609 + t687 * t610;
t628 = -t643 * qJD(5) + qJDD(1) * t699;
t636 = t663 * mrSges(6,1) - t643 * mrSges(6,3);
t605 = m(6) * t608 - t662 * mrSges(6,2) + t628 * mrSges(6,3) + t642 * t626 - t663 * t636;
t596 = t687 * t604 + t685 * t605;
t611 = -t679 * t616 + t724;
t709 = mrSges(5,1) * t679 + mrSges(5,2) * t682;
t644 = t709 * t723;
t700 = mrSges(5,2) * t683 - mrSges(5,3) * t726;
t646 = t700 * qJD(1);
t701 = -mrSges(5,1) * t683 - mrSges(5,3) * t725;
t594 = m(5) * t611 + t701 * qJDD(1) + (-t644 * t725 - t646 * t683) * qJD(1) + t596;
t647 = t701 * qJD(1);
t712 = -t685 * t604 + t687 * t605;
t595 = m(5) * t612 + t700 * qJDD(1) + (-t644 * t726 + t647 * t683) * qJD(1) + t712;
t713 = -t679 * t594 + t682 * t595;
t589 = m(4) * t623 + (qJDD(1) * mrSges(4,3) + qJD(1) * t653) * t683 + t713;
t615 = t652 * t723 + qJDD(4) - t622;
t613 = -pkin(6) * t718 + (pkin(4) * qJDD(1) * t679 + qJD(1) * t651 * t682) * t680 + t615;
t696 = m(6) * t613 - t628 * mrSges(6,1) + t629 * mrSges(6,2) - t642 * t635 + t643 * t636;
t692 = m(5) * t615 + t696;
t703 = t646 * t679 + t647 * t682;
t600 = m(4) * t622 + ((-mrSges(4,3) - t709) * qJDD(1) + (-t653 - t703) * qJD(1)) * t680 - t692;
t714 = t683 * t589 - t680 * t600;
t579 = m(3) * t638 - t689 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t714;
t592 = t682 * t594 + t679 * t595;
t632 = -qJDD(1) * pkin(2) + t697;
t693 = -m(4) * t632 + mrSges(4,1) * t719 - t592 + (t683 ^ 2 * t689 + t727) * mrSges(4,3);
t586 = m(3) * t637 - t689 * mrSges(3,2) + (mrSges(3,1) - t728) * qJDD(1) + t693;
t574 = t681 * t579 + t684 * t586;
t582 = t680 * t589 + t683 * t600;
t580 = m(3) * t678 + t582;
t715 = t684 * t579 - t681 * t586;
t571 = m(2) * t660 - t689 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t715;
t572 = m(2) * t661 + qJDD(1) * mrSges(2,1) - t689 * mrSges(2,2) + t574;
t716 = t688 * t571 - t686 * t572;
t708 = Ifges(4,1) * t680 + Ifges(4,4) * t683;
t707 = Ifges(4,5) * t680 + Ifges(4,6) * t683;
t706 = Ifges(5,5) * t682 - Ifges(5,6) * t679;
t704 = -t686 * t571 - t688 * t572;
t695 = -Ifges(5,5) * t683 + (Ifges(5,1) * t682 - Ifges(5,4) * t679) * t680;
t694 = -Ifges(5,6) * t683 + (Ifges(5,4) * t682 - Ifges(5,2) * t679) * t680;
t617 = Ifges(6,5) * t643 + Ifges(6,6) * t642 + Ifges(6,3) * t663;
t619 = Ifges(6,1) * t643 + Ifges(6,4) * t642 + Ifges(6,5) * t663;
t597 = -mrSges(6,1) * t613 + mrSges(6,3) * t608 + Ifges(6,4) * t629 + Ifges(6,2) * t628 + Ifges(6,6) * t662 - t643 * t617 + t663 * t619;
t618 = Ifges(6,4) * t643 + Ifges(6,2) * t642 + Ifges(6,6) * t663;
t598 = mrSges(6,2) * t613 - mrSges(6,3) * t607 + Ifges(6,1) * t629 + Ifges(6,4) * t628 + Ifges(6,5) * t662 + t642 * t617 - t663 * t618;
t639 = (-Ifges(5,3) * t683 + t706 * t680) * qJD(1);
t641 = t695 * qJD(1);
t583 = -mrSges(5,1) * t615 + mrSges(5,3) * t612 + t685 * t598 + t687 * t597 - pkin(4) * t696 + pkin(6) * t712 + (-t639 * t725 - t683 * t641) * qJD(1) + t694 * qJDD(1);
t640 = t694 * qJD(1);
t584 = mrSges(5,2) * t615 - mrSges(5,3) * t611 - pkin(6) * t596 - t685 * t597 + t687 * t598 + (-t639 * t726 + t640 * t683) * qJD(1) + t695 * qJDD(1);
t654 = t707 * qJD(1);
t569 = mrSges(4,2) * t632 - mrSges(4,3) * t622 - qJ(4) * t592 + t708 * qJDD(1) - t679 * t583 + t682 * t584 + t654 * t722;
t690 = mrSges(6,1) * t607 - mrSges(6,2) * t608 + Ifges(6,5) * t629 + Ifges(6,6) * t628 + Ifges(6,3) * t662 + t643 * t618 - t642 * t619;
t576 = -mrSges(4,1) * t632 - mrSges(5,1) * t611 + mrSges(5,2) * t612 + mrSges(4,3) * t623 - pkin(3) * t592 - pkin(4) * t596 + (Ifges(4,2) + Ifges(5,3)) * t719 + ((Ifges(4,4) - t706) * qJDD(1) + (-t640 * t682 - t641 * t679 - t654) * qJD(1)) * t680 - t690;
t591 = mrSges(4,2) * t720 - t693;
t691 = mrSges(2,1) * t661 + mrSges(3,1) * t637 - mrSges(2,2) * t660 - mrSges(3,2) * t638 + pkin(1) * t574 - pkin(2) * t591 + qJ(3) * t714 + t680 * t569 + t683 * t576 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t606 = (t703 * qJD(1) + t709 * qJDD(1)) * t680 + t692;
t567 = -mrSges(3,1) * t678 + mrSges(3,3) * t638 - mrSges(4,1) * t622 + mrSges(4,2) * t623 - t679 * t584 - t682 * t583 + pkin(3) * t606 - qJ(4) * t713 - pkin(2) * t582 + (Ifges(3,6) - t707) * qJDD(1) + (Ifges(3,5) - t680 * (Ifges(4,4) * t680 + Ifges(4,2) * t683) + t683 * t708) * t689;
t566 = mrSges(3,2) * t678 - mrSges(3,3) * t637 + Ifges(3,5) * qJDD(1) - t689 * Ifges(3,6) - qJ(3) * t582 + t683 * t569 - t680 * t576;
t565 = -mrSges(2,2) * g(1) - mrSges(2,3) * t661 + Ifges(2,5) * qJDD(1) - t689 * Ifges(2,6) - qJ(2) * t574 + t684 * t566 - t681 * t567;
t564 = mrSges(2,1) * g(1) + mrSges(2,3) * t660 + t689 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t580 + qJ(2) * t715 + t681 * t566 + t684 * t567;
t1 = [(-m(1) - m(2)) * g(1) + t580; -m(1) * g(2) + t704; -m(1) * g(3) + t716; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t691; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t716 - t688 * t564 - t686 * t565; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t704 - t686 * t564 + t688 * t565; t691; t580; t591; t606; t690;];
tauJB = t1;
