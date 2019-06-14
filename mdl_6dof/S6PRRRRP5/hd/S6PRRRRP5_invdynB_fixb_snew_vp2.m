% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRRRRP5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-05-05 10:09
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRRRRP5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP5_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP5_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRP5_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP5_invdynB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP5_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP5_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP5_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 10:01:56
% EndTime: 2019-05-05 10:02:18
% DurationCPUTime: 21.38s
% Computational Cost: add. (370754->331), mult. (761034->425), div. (0->0), fcn. (603963->14), ass. (0->147)
t730 = Ifges(6,1) + Ifges(7,1);
t724 = Ifges(6,4) + Ifges(7,4);
t723 = Ifges(6,5) + Ifges(7,5);
t729 = Ifges(6,2) + Ifges(7,2);
t728 = Ifges(6,6) + Ifges(7,6);
t727 = Ifges(6,3) + Ifges(7,3);
t680 = sin(pkin(7));
t687 = sin(qJ(3));
t691 = cos(qJ(3));
t708 = qJD(2) * qJD(3);
t665 = (-qJDD(2) * t691 + t687 * t708) * t680;
t679 = sin(pkin(12));
t682 = cos(pkin(12));
t671 = g(1) * t679 - g(2) * t682;
t672 = -g(1) * t682 - g(2) * t679;
t678 = -g(3) + qJDD(1);
t688 = sin(qJ(2));
t684 = cos(pkin(6));
t692 = cos(qJ(2));
t715 = t684 * t692;
t681 = sin(pkin(6));
t718 = t681 * t692;
t642 = t671 * t715 - t672 * t688 + t678 * t718;
t693 = qJD(2) ^ 2;
t726 = pkin(9) * t680;
t638 = qJDD(2) * pkin(2) + t693 * t726 + t642;
t716 = t684 * t688;
t719 = t681 * t688;
t643 = t671 * t716 + t692 * t672 + t678 * t719;
t639 = -pkin(2) * t693 + qJDD(2) * t726 + t643;
t658 = -t671 * t681 + t678 * t684;
t683 = cos(pkin(7));
t600 = -t687 * t639 + (t638 * t683 + t658 * t680) * t691;
t725 = -mrSges(6,2) - mrSges(7,2);
t677 = qJD(2) * t683 + qJD(3);
t709 = qJD(2) * t680;
t704 = t691 * t709;
t661 = -mrSges(4,2) * t677 + mrSges(4,3) * t704;
t662 = (-mrSges(4,1) * t691 + mrSges(4,2) * t687) * t709;
t664 = (qJDD(2) * t687 + t691 * t708) * t680;
t676 = qJDD(2) * t683 + qJDD(3);
t717 = t683 * t687;
t720 = t680 * t687;
t601 = t638 * t717 + t691 * t639 + t658 * t720;
t663 = (-pkin(3) * t691 - pkin(10) * t687) * t709;
t675 = t677 ^ 2;
t597 = -pkin(3) * t675 + pkin(10) * t676 + t663 * t704 + t601;
t654 = t683 * t658;
t599 = t665 * pkin(3) - t664 * pkin(10) + t654 + (-t638 + (pkin(3) * t687 - pkin(10) * t691) * t677 * qJD(2)) * t680;
t686 = sin(qJ(4));
t690 = cos(qJ(4));
t593 = t690 * t597 + t686 * t599;
t705 = t687 * t709;
t656 = t677 * t690 - t686 * t705;
t657 = t677 * t686 + t690 * t705;
t641 = -pkin(4) * t656 - pkin(11) * t657;
t659 = qJDD(4) + t665;
t670 = qJD(4) - t704;
t669 = t670 ^ 2;
t588 = -pkin(4) * t669 + pkin(11) * t659 + t641 * t656 + t593;
t596 = -t676 * pkin(3) - t675 * pkin(10) + t663 * t705 - t600;
t633 = -qJD(4) * t657 - t664 * t686 + t676 * t690;
t634 = qJD(4) * t656 + t664 * t690 + t676 * t686;
t591 = (-t656 * t670 - t634) * pkin(11) + (t657 * t670 - t633) * pkin(4) + t596;
t685 = sin(qJ(5));
t689 = cos(qJ(5));
t583 = -t685 * t588 + t689 * t591;
t645 = -t657 * t685 + t670 * t689;
t606 = qJD(5) * t645 + t634 * t689 + t659 * t685;
t646 = t657 * t689 + t670 * t685;
t617 = -mrSges(7,1) * t645 + mrSges(7,2) * t646;
t618 = -mrSges(6,1) * t645 + mrSges(6,2) * t646;
t655 = qJD(5) - t656;
t622 = -mrSges(6,2) * t655 + mrSges(6,3) * t645;
t631 = qJDD(5) - t633;
t580 = -0.2e1 * qJD(6) * t646 + (t645 * t655 - t606) * qJ(6) + (t645 * t646 + t631) * pkin(5) + t583;
t621 = -mrSges(7,2) * t655 + mrSges(7,3) * t645;
t707 = m(7) * t580 + t631 * mrSges(7,1) + t655 * t621;
t573 = m(6) * t583 + t631 * mrSges(6,1) + t655 * t622 + (-t617 - t618) * t646 + (-mrSges(6,3) - mrSges(7,3)) * t606 + t707;
t584 = t689 * t588 + t685 * t591;
t605 = -qJD(5) * t646 - t634 * t685 + t659 * t689;
t623 = pkin(5) * t655 - qJ(6) * t646;
t644 = t645 ^ 2;
t582 = -pkin(5) * t644 + qJ(6) * t605 + 0.2e1 * qJD(6) * t645 - t623 * t655 + t584;
t706 = m(7) * t582 + t605 * mrSges(7,3) + t645 * t617;
t624 = mrSges(7,1) * t655 - mrSges(7,3) * t646;
t710 = -mrSges(6,1) * t655 + mrSges(6,3) * t646 - t624;
t575 = m(6) * t584 + t605 * mrSges(6,3) + t645 * t618 + t725 * t631 + t710 * t655 + t706;
t572 = t573 * t689 + t575 * t685;
t647 = -mrSges(5,2) * t670 + mrSges(5,3) * t656;
t648 = mrSges(5,1) * t670 - mrSges(5,3) * t657;
t695 = -m(5) * t596 + t633 * mrSges(5,1) - mrSges(5,2) * t634 + t656 * t647 - t648 * t657 - t572;
t567 = m(4) * t600 + mrSges(4,1) * t676 - mrSges(4,3) * t664 + t661 * t677 - t662 * t705 + t695;
t721 = t567 * t691;
t660 = mrSges(4,1) * t677 - mrSges(4,3) * t705;
t640 = -mrSges(5,1) * t656 + mrSges(5,2) * t657;
t701 = -t573 * t685 + t689 * t575;
t570 = m(5) * t593 - mrSges(5,2) * t659 + mrSges(5,3) * t633 + t640 * t656 - t648 * t670 + t701;
t592 = -t686 * t597 + t599 * t690;
t587 = -pkin(4) * t659 - pkin(11) * t669 + t657 * t641 - t592;
t585 = -pkin(5) * t605 - qJ(6) * t644 + t623 * t646 + qJDD(6) + t587;
t700 = m(7) * t585 - t605 * mrSges(7,1) - t645 * t621;
t694 = -m(6) * t587 + t605 * mrSges(6,1) + t725 * t606 + t645 * t622 + t710 * t646 - t700;
t577 = m(5) * t592 + t659 * mrSges(5,1) - t634 * mrSges(5,3) - t657 * t640 + t670 * t647 + t694;
t702 = t690 * t570 - t577 * t686;
t560 = m(4) * t601 - mrSges(4,2) * t676 - mrSges(4,3) * t665 - t660 * t677 + t662 * t704 + t702;
t563 = t686 * t570 + t690 * t577;
t619 = -t680 * t638 + t654;
t562 = m(4) * t619 + t665 * mrSges(4,1) + t664 * mrSges(4,2) + (t660 * t687 - t661 * t691) * t709 + t563;
t549 = t560 * t717 - t562 * t680 + t683 * t721;
t545 = m(3) * t642 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t693 + t549;
t548 = t560 * t720 + t683 * t562 + t680 * t721;
t547 = m(3) * t658 + t548;
t555 = t691 * t560 - t567 * t687;
t554 = m(3) * t643 - mrSges(3,1) * t693 - qJDD(2) * mrSges(3,2) + t555;
t535 = t545 * t715 - t547 * t681 + t554 * t716;
t533 = m(2) * t671 + t535;
t541 = -t545 * t688 + t692 * t554;
t540 = m(2) * t672 + t541;
t714 = t682 * t533 + t679 * t540;
t713 = t645 * t728 + t646 * t723 + t655 * t727;
t712 = -t645 * t729 - t646 * t724 - t655 * t728;
t711 = t724 * t645 + t646 * t730 + t723 * t655;
t534 = t545 * t718 + t684 * t547 + t554 * t719;
t703 = -t533 * t679 + t682 * t540;
t564 = -mrSges(6,1) * t587 + mrSges(6,3) * t584 - mrSges(7,1) * t585 + mrSges(7,3) * t582 - pkin(5) * t700 + qJ(6) * t706 + (-qJ(6) * t624 + t711) * t655 + (-pkin(5) * t624 - t713) * t646 + (-mrSges(7,2) * qJ(6) + t728) * t631 + (-mrSges(7,2) * pkin(5) + t724) * t606 + t729 * t605;
t578 = -t606 * mrSges(7,3) - t646 * t617 + t707;
t571 = mrSges(6,2) * t587 + mrSges(7,2) * t585 - mrSges(6,3) * t583 - mrSges(7,3) * t580 - qJ(6) * t578 + t724 * t605 + t606 * t730 + t723 * t631 + t713 * t645 + t712 * t655;
t627 = Ifges(5,5) * t657 + Ifges(5,6) * t656 + Ifges(5,3) * t670;
t628 = Ifges(5,4) * t657 + Ifges(5,2) * t656 + Ifges(5,6) * t670;
t550 = mrSges(5,2) * t596 - mrSges(5,3) * t592 + Ifges(5,1) * t634 + Ifges(5,4) * t633 + Ifges(5,5) * t659 - pkin(11) * t572 - t564 * t685 + t571 * t689 + t627 * t656 - t628 * t670;
t629 = Ifges(5,1) * t657 + Ifges(5,4) * t656 + Ifges(5,5) * t670;
t556 = -mrSges(5,1) * t596 - mrSges(6,1) * t583 - mrSges(7,1) * t580 + mrSges(6,2) * t584 + mrSges(7,2) * t582 + mrSges(5,3) * t593 + Ifges(5,4) * t634 + Ifges(5,2) * t633 + Ifges(5,6) * t659 - pkin(4) * t572 - pkin(5) * t578 - t657 * t627 + t670 * t629 + t712 * t646 + t711 * t645 - t727 * t631 - t723 * t606 - t728 * t605;
t651 = Ifges(4,6) * t677 + (Ifges(4,4) * t687 + Ifges(4,2) * t691) * t709;
t652 = Ifges(4,5) * t677 + (Ifges(4,1) * t687 + Ifges(4,4) * t691) * t709;
t536 = Ifges(4,5) * t664 - Ifges(4,6) * t665 + Ifges(4,3) * t676 + mrSges(4,1) * t600 - mrSges(4,2) * t601 + t686 * t550 + t690 * t556 + pkin(3) * t695 + pkin(10) * t702 + (t651 * t687 - t652 * t691) * t709;
t650 = Ifges(4,3) * t677 + (Ifges(4,5) * t687 + Ifges(4,6) * t691) * t709;
t537 = mrSges(4,2) * t619 - mrSges(4,3) * t600 + Ifges(4,1) * t664 - Ifges(4,4) * t665 + Ifges(4,5) * t676 - pkin(10) * t563 + t550 * t690 - t556 * t686 + t650 * t704 - t651 * t677;
t542 = Ifges(4,4) * t664 - Ifges(4,2) * t665 + Ifges(4,6) * t676 - t650 * t705 + t677 * t652 - mrSges(4,1) * t619 + mrSges(4,3) * t601 - Ifges(5,5) * t634 - Ifges(5,6) * t633 - Ifges(5,3) * t659 - t657 * t628 + t656 * t629 - mrSges(5,1) * t592 + mrSges(5,2) * t593 - t685 * t571 - t689 * t564 - pkin(4) * t694 - pkin(11) * t701 - pkin(3) * t563;
t696 = pkin(9) * t555 + t537 * t687 + t542 * t691;
t530 = -mrSges(3,1) * t658 + mrSges(3,3) * t643 + t693 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t548 - t680 * t536 + t696 * t683;
t531 = mrSges(3,2) * t658 - mrSges(3,3) * t642 + Ifges(3,5) * qJDD(2) - t693 * Ifges(3,6) + t691 * t537 - t687 * t542 + (-t548 * t680 - t549 * t683) * pkin(9);
t697 = pkin(8) * t541 + t530 * t692 + t531 * t688;
t529 = mrSges(3,1) * t642 - mrSges(3,2) * t643 + Ifges(3,3) * qJDD(2) + pkin(2) * t549 + t683 * t536 + t696 * t680;
t528 = mrSges(2,2) * t678 - mrSges(2,3) * t671 - t688 * t530 + t692 * t531 + (-t534 * t681 - t535 * t684) * pkin(8);
t527 = -mrSges(2,1) * t678 + mrSges(2,3) * t672 - pkin(1) * t534 - t681 * t529 + t697 * t684;
t1 = [-m(1) * g(1) + t703; -m(1) * g(2) + t714; -m(1) * g(3) + m(2) * t678 + t534; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t714 - t679 * t527 + t682 * t528; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t703 + t682 * t527 + t679 * t528; -mrSges(1,1) * g(2) + mrSges(2,1) * t671 + mrSges(1,2) * g(1) - mrSges(2,2) * t672 + pkin(1) * t535 + t684 * t529 + t697 * t681;];
tauB  = t1;
