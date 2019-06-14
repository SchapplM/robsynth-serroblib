% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRRRPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
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
% Datum: 2019-05-05 07:33
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRRRPR3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR3_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR3_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR3_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR3_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR3_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR3_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR3_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 07:28:06
% EndTime: 2019-05-05 07:28:17
% DurationCPUTime: 8.14s
% Computational Cost: add. (127824->322), mult. (255244->398), div. (0->0), fcn. (176804->12), ass. (0->136)
t728 = Ifges(5,1) + Ifges(6,2);
t727 = -Ifges(6,1) - Ifges(5,3);
t723 = Ifges(5,4) + Ifges(6,6);
t722 = Ifges(5,5) - Ifges(6,4);
t726 = Ifges(5,2) + Ifges(6,3);
t721 = Ifges(5,6) - Ifges(6,5);
t683 = sin(pkin(11));
t685 = cos(pkin(11));
t670 = g(1) * t683 - g(2) * t685;
t671 = -g(1) * t685 - g(2) * t683;
t682 = -g(3) + qJDD(1);
t684 = sin(pkin(6));
t686 = cos(pkin(6));
t690 = sin(qJ(2));
t693 = cos(qJ(2));
t632 = -t690 * t671 + (t670 * t686 + t682 * t684) * t693;
t725 = -2 * qJD(5);
t724 = cos(qJ(4));
t694 = qJD(2) ^ 2;
t699 = -qJDD(2) * pkin(2) - t632;
t624 = -t694 * pkin(8) + t699;
t689 = sin(qJ(3));
t692 = cos(qJ(3));
t709 = qJD(2) * qJD(3);
t708 = t692 * t709;
t668 = qJDD(2) * t689 + t708;
t669 = qJDD(2) * t692 - t689 * t709;
t711 = qJD(2) * t689;
t672 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t711;
t710 = qJD(2) * t692;
t673 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t710;
t688 = sin(qJ(4));
t660 = (t688 * t692 + t724 * t689) * qJD(2);
t620 = t660 * qJD(4) + t688 * t668 - t724 * t669;
t675 = qJD(3) * pkin(3) - pkin(9) * t711;
t681 = t692 ^ 2;
t604 = -t669 * pkin(3) + t675 * t711 + (-pkin(9) * t681 - pkin(8)) * t694 + t699;
t659 = t688 * t711 - t724 * t710;
t621 = -t659 * qJD(4) + t724 * t668 + t688 * t669;
t680 = qJD(3) + qJD(4);
t719 = t659 * t680;
t695 = (-t621 + t719) * qJ(5) + t604 + (t680 * pkin(4) + t725) * t660;
t593 = t620 * pkin(4) + t695;
t647 = mrSges(6,1) * t659 - mrSges(6,3) * t680;
t648 = mrSges(6,1) * t660 + mrSges(6,2) * t680;
t717 = t686 * t690;
t718 = t684 * t690;
t633 = t670 * t717 + t693 * t671 + t682 * t718;
t625 = -pkin(2) * t694 + qJDD(2) * pkin(8) + t633;
t650 = -t670 * t684 + t682 * t686;
t609 = -t689 * t625 + t692 * t650;
t599 = (-t668 + t708) * pkin(9) + (t689 * t692 * t694 + qJDD(3)) * pkin(3) + t609;
t610 = t692 * t625 + t689 * t650;
t600 = -pkin(3) * t681 * t694 + pkin(9) * t669 - qJD(3) * t675 + t610;
t594 = t724 * t599 - t688 * t600;
t637 = pkin(4) * t659 - qJ(5) * t660;
t678 = t680 ^ 2;
t679 = qJDD(3) + qJDD(4);
t591 = -t679 * pkin(4) - t678 * qJ(5) + t660 * t637 + qJDD(5) - t594;
t586 = (t659 * t660 - t679) * pkin(10) + (t621 + t719) * pkin(5) + t591;
t649 = pkin(5) * t660 - pkin(10) * t680;
t655 = t659 ^ 2;
t589 = -t655 * pkin(5) - t660 * t649 + (pkin(4) + pkin(10)) * t620 + t695;
t687 = sin(qJ(6));
t691 = cos(qJ(6));
t584 = t586 * t691 - t589 * t687;
t641 = t659 * t691 - t680 * t687;
t603 = qJD(6) * t641 + t620 * t687 + t679 * t691;
t642 = t659 * t687 + t680 * t691;
t611 = -mrSges(7,1) * t641 + mrSges(7,2) * t642;
t618 = qJDD(6) + t621;
t653 = qJD(6) + t660;
t622 = -mrSges(7,2) * t653 + mrSges(7,3) * t641;
t582 = m(7) * t584 + mrSges(7,1) * t618 - mrSges(7,3) * t603 - t611 * t642 + t622 * t653;
t585 = t586 * t687 + t589 * t691;
t602 = -qJD(6) * t642 + t620 * t691 - t679 * t687;
t623 = mrSges(7,1) * t653 - mrSges(7,3) * t642;
t583 = m(7) * t585 - mrSges(7,2) * t618 + mrSges(7,3) * t602 + t611 * t641 - t623 * t653;
t715 = -t687 * t582 + t691 * t583;
t573 = m(6) * t593 - t620 * mrSges(6,2) - t621 * mrSges(6,3) - t659 * t647 - t660 * t648 + t715;
t645 = -mrSges(5,2) * t680 - mrSges(5,3) * t659;
t646 = mrSges(5,1) * t680 - mrSges(5,3) * t660;
t697 = m(5) * t604 + t620 * mrSges(5,1) + mrSges(5,2) * t621 + t659 * t645 + t646 * t660 + t573;
t696 = -m(4) * t624 + t669 * mrSges(4,1) - mrSges(4,2) * t668 - t672 * t711 + t673 * t710 - t697;
t570 = m(3) * t632 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t694 + t696;
t720 = t570 * t693;
t638 = mrSges(5,1) * t659 + mrSges(5,2) * t660;
t574 = t691 * t582 + t687 * t583;
t639 = -mrSges(6,2) * t659 - mrSges(6,3) * t660;
t701 = -m(6) * t591 - t621 * mrSges(6,1) - t660 * t639 - t574;
t572 = m(5) * t594 - t621 * mrSges(5,3) - t660 * t638 + (t645 - t647) * t680 + (mrSges(5,1) - mrSges(6,2)) * t679 + t701;
t595 = t688 * t599 + t724 * t600;
t700 = -t678 * pkin(4) + t679 * qJ(5) - t659 * t637 + t595;
t590 = t680 * t725 - t700;
t588 = -t620 * pkin(5) - t655 * pkin(10) + ((2 * qJD(5)) + t649) * t680 + t700;
t702 = -m(7) * t588 + t602 * mrSges(7,1) - t603 * mrSges(7,2) + t641 * t622 - t642 * t623;
t698 = -m(6) * t590 + t679 * mrSges(6,3) + t680 * t648 - t702;
t579 = m(5) * t595 - t679 * mrSges(5,2) - t680 * t646 + (-t638 - t639) * t659 + (-mrSges(5,3) - mrSges(6,1)) * t620 + t698;
t567 = t724 * t572 + t688 * t579;
t667 = (-mrSges(4,1) * t692 + mrSges(4,2) * t689) * qJD(2);
t565 = m(4) * t609 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t668 + qJD(3) * t673 - t667 * t711 + t567;
t705 = -t572 * t688 + t724 * t579;
t566 = m(4) * t610 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t669 - qJD(3) * t672 + t667 * t710 + t705;
t706 = -t565 * t689 + t692 * t566;
t557 = m(3) * t633 - mrSges(3,1) * t694 - qJDD(2) * mrSges(3,2) + t706;
t560 = t692 * t565 + t689 * t566;
t559 = m(3) * t650 + t560;
t548 = t557 * t717 - t559 * t684 + t686 * t720;
t546 = m(2) * t670 + t548;
t552 = t693 * t557 - t570 * t690;
t551 = m(2) * t671 + t552;
t716 = t685 * t546 + t683 * t551;
t714 = t659 * t721 - t660 * t722 + t680 * t727;
t713 = t659 * t726 - t660 * t723 - t680 * t721;
t712 = -t723 * t659 + t660 * t728 + t722 * t680;
t547 = t557 * t718 + t686 * t559 + t684 * t720;
t707 = -t546 * t683 + t685 * t551;
t605 = Ifges(7,5) * t642 + Ifges(7,6) * t641 + Ifges(7,3) * t653;
t607 = Ifges(7,1) * t642 + Ifges(7,4) * t641 + Ifges(7,5) * t653;
t575 = -mrSges(7,1) * t588 + mrSges(7,3) * t585 + Ifges(7,4) * t603 + Ifges(7,2) * t602 + Ifges(7,6) * t618 - t605 * t642 + t607 * t653;
t606 = Ifges(7,4) * t642 + Ifges(7,2) * t641 + Ifges(7,6) * t653;
t576 = mrSges(7,2) * t588 - mrSges(7,3) * t584 + Ifges(7,1) * t603 + Ifges(7,4) * t602 + Ifges(7,5) * t618 + t605 * t641 - t606 * t653;
t553 = -mrSges(5,1) * t604 - mrSges(6,1) * t590 + mrSges(6,2) * t593 + mrSges(5,3) * t595 - pkin(4) * t573 - pkin(5) * t702 - pkin(10) * t715 - t691 * t575 - t687 * t576 - t620 * t726 + t723 * t621 + t714 * t660 + t721 * t679 + t712 * t680;
t561 = mrSges(6,1) * t591 + mrSges(7,1) * t584 + mrSges(5,2) * t604 - mrSges(7,2) * t585 - mrSges(5,3) * t594 - mrSges(6,3) * t593 + Ifges(7,5) * t603 + Ifges(7,6) * t602 + Ifges(7,3) * t618 + pkin(5) * t574 - qJ(5) * t573 + t642 * t606 - t641 * t607 + t713 * t680 + t722 * t679 + t714 * t659 + t728 * t621 - t723 * t620;
t656 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t689 + Ifges(4,6) * t692) * qJD(2);
t658 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t689 + Ifges(4,4) * t692) * qJD(2);
t542 = -mrSges(4,1) * t624 + mrSges(4,3) * t610 + Ifges(4,4) * t668 + Ifges(4,2) * t669 + Ifges(4,6) * qJDD(3) - pkin(3) * t697 + pkin(9) * t705 + qJD(3) * t658 + t724 * t553 + t688 * t561 - t656 * t711;
t657 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t689 + Ifges(4,2) * t692) * qJD(2);
t544 = mrSges(4,2) * t624 - mrSges(4,3) * t609 + Ifges(4,1) * t668 + Ifges(4,4) * t669 + Ifges(4,5) * qJDD(3) - pkin(9) * t567 - qJD(3) * t657 - t688 * t553 + t724 * t561 + t656 * t710;
t541 = mrSges(3,2) * t650 - mrSges(3,3) * t632 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t694 - pkin(8) * t560 - t542 * t689 + t544 * t692;
t543 = (qJ(5) * t639 - t712) * t659 + t713 * t660 + (mrSges(6,1) * qJ(5) + t721) * t620 - t722 * t621 + Ifges(3,6) * qJDD(2) - Ifges(4,3) * qJDD(3) - t691 * t576 + t694 * Ifges(3,5) + t687 * t575 - Ifges(4,5) * t668 - Ifges(4,6) * t669 - mrSges(3,1) * t650 + mrSges(3,3) * t633 - mrSges(4,1) * t609 + mrSges(4,2) * t610 - mrSges(6,2) * t591 - mrSges(5,1) * t594 + mrSges(5,2) * t595 + mrSges(6,3) * t590 + pkin(10) * t574 - pkin(3) * t567 + (-t657 * t689 + t658 * t692) * qJD(2) - pkin(2) * t560 - qJ(5) * t698 - pkin(4) * (-t680 * t647 + t701) + (mrSges(6,2) * pkin(4) + t727) * t679;
t703 = pkin(7) * t552 + t541 * t690 + t543 * t693;
t540 = mrSges(3,1) * t632 - mrSges(3,2) * t633 + Ifges(3,3) * qJDD(2) + pkin(2) * t696 + pkin(8) * t706 + t692 * t542 + t689 * t544;
t539 = mrSges(2,2) * t682 - mrSges(2,3) * t670 + t693 * t541 - t690 * t543 + (-t547 * t684 - t548 * t686) * pkin(7);
t538 = -mrSges(2,1) * t682 + mrSges(2,3) * t671 - pkin(1) * t547 - t684 * t540 + t703 * t686;
t1 = [-m(1) * g(1) + t707; -m(1) * g(2) + t716; -m(1) * g(3) + m(2) * t682 + t547; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t716 - t683 * t538 + t685 * t539; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t707 + t685 * t538 + t683 * t539; -mrSges(1,1) * g(2) + mrSges(2,1) * t670 + mrSges(1,2) * g(1) - mrSges(2,2) * t671 + pkin(1) * t548 + t686 * t540 + t703 * t684;];
tauB  = t1;
