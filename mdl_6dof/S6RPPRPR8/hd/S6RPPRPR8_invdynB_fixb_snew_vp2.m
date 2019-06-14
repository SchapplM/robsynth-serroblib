% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPPRPR8
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3]';
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
% Datum: 2019-05-05 14:39
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPPRPR8_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR8_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR8_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR8_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR8_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR8_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR8_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR8_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR8_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:37:51
% EndTime: 2019-05-05 14:37:55
% DurationCPUTime: 3.42s
% Computational Cost: add. (31492->299), mult. (71769->346), div. (0->0), fcn. (46938->8), ass. (0->129)
t708 = Ifges(5,1) + Ifges(6,2);
t696 = Ifges(5,4) + Ifges(6,6);
t694 = Ifges(5,5) - Ifges(6,4);
t707 = Ifges(5,2) + Ifges(6,3);
t692 = Ifges(5,6) - Ifges(6,5);
t706 = -Ifges(5,3) - Ifges(6,1);
t649 = sin(qJ(1));
t652 = cos(qJ(1));
t627 = t649 * g(1) - t652 * g(2);
t654 = qJD(1) ^ 2;
t664 = -t654 * qJ(2) + qJDD(2) - t627;
t691 = -pkin(1) - qJ(3);
t705 = -(2 * qJD(1) * qJD(3)) + qJDD(1) * t691 + t664;
t646 = sin(pkin(9));
t639 = t646 ^ 2;
t647 = cos(pkin(9));
t685 = t647 ^ 2 + t639;
t676 = t685 * mrSges(4,3);
t628 = -t652 * g(1) - t649 * g(2);
t704 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t628;
t663 = qJDD(3) + t704;
t680 = qJDD(1) * t646;
t584 = pkin(3) * t680 + (-pkin(7) * t685 + t691) * t654 + t663;
t651 = cos(qJ(4));
t700 = sin(qJ(4));
t677 = t646 * t700;
t679 = qJDD(1) * t647;
t666 = t646 * t651 + t647 * t700;
t623 = t666 * qJD(1);
t683 = qJD(4) * t623;
t603 = -qJDD(1) * t677 + t651 * t679 - t683;
t684 = qJD(1) * t647;
t624 = -qJD(1) * t677 + t651 * t684;
t612 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t624;
t703 = m(5) * t584 + t603 * mrSges(5,2) + t624 * t612;
t609 = t654 * t691 + t663;
t702 = m(4) * t609 + mrSges(4,1) * t680 + mrSges(4,2) * t679;
t701 = -2 * qJD(5);
t699 = pkin(3) * t654;
t698 = mrSges(2,1) - mrSges(3,2);
t697 = mrSges(5,1) - mrSges(6,2);
t695 = Ifges(2,5) - Ifges(3,4);
t693 = -Ifges(2,6) + Ifges(3,5);
t600 = t646 * g(3) + t647 * t705;
t579 = (-pkin(7) * qJDD(1) - t646 * t699) * t647 + t600;
t601 = -g(3) * t647 + t646 * t705;
t580 = -pkin(7) * t680 - t639 * t699 + t601;
t564 = t651 * t579 - t700 * t580;
t594 = mrSges(5,1) * t623 + mrSges(5,2) * t624;
t682 = t624 * qJD(4);
t602 = qJDD(1) * t666 + t682;
t615 = pkin(5) * t624 - qJD(4) * pkin(8);
t622 = t623 ^ 2;
t656 = pkin(4) * t682 + t624 * t701 + (-t603 + t683) * qJ(5) + t584;
t556 = -t622 * pkin(5) - t624 * t615 + (pkin(4) + pkin(8)) * t602 + t656;
t593 = pkin(4) * t623 - qJ(5) * t624;
t653 = qJD(4) ^ 2;
t563 = -qJDD(4) * pkin(4) - t653 * qJ(5) + t624 * t593 + qJDD(5) - t564;
t557 = (t623 * t624 - qJDD(4)) * pkin(8) + (t603 + t683) * pkin(5) + t563;
t648 = sin(qJ(6));
t650 = cos(qJ(6));
t554 = -t556 * t648 + t557 * t650;
t605 = -qJD(4) * t648 + t623 * t650;
t572 = qJD(6) * t605 + qJDD(4) * t650 + t602 * t648;
t606 = qJD(4) * t650 + t623 * t648;
t574 = -mrSges(7,1) * t605 + mrSges(7,2) * t606;
t620 = qJD(6) + t624;
t581 = -mrSges(7,2) * t620 + mrSges(7,3) * t605;
t599 = qJDD(6) + t603;
t552 = m(7) * t554 + mrSges(7,1) * t599 - mrSges(7,3) * t572 - t574 * t606 + t581 * t620;
t555 = t556 * t650 + t557 * t648;
t571 = -qJD(6) * t606 - qJDD(4) * t648 + t602 * t650;
t582 = mrSges(7,1) * t620 - mrSges(7,3) * t606;
t553 = m(7) * t555 - mrSges(7,2) * t599 + mrSges(7,3) * t571 + t574 * t605 - t582 * t620;
t545 = t650 * t552 + t648 * t553;
t595 = -mrSges(6,2) * t623 - mrSges(6,3) * t624;
t660 = -m(6) * t563 - t603 * mrSges(6,1) - t624 * t595 - t545;
t611 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t623;
t613 = mrSges(6,1) * t623 - qJD(4) * mrSges(6,3);
t686 = t611 - t613;
t543 = m(5) * t564 - t603 * mrSges(5,3) + qJD(4) * t686 + qJDD(4) * t697 - t624 * t594 + t660;
t565 = t700 * t579 + t651 * t580;
t659 = -t653 * pkin(4) + qJDD(4) * qJ(5) - t623 * t593 + t565;
t562 = qJD(4) * t701 - t659;
t614 = mrSges(6,1) * t624 + qJD(4) * mrSges(6,2);
t559 = -t602 * pkin(5) - t622 * pkin(8) + ((2 * qJD(5)) + t615) * qJD(4) + t659;
t662 = -m(7) * t559 + t571 * mrSges(7,1) - t572 * mrSges(7,2) + t605 * t581 - t606 * t582;
t658 = -m(6) * t562 + qJDD(4) * mrSges(6,3) + qJD(4) * t614 - t662;
t550 = m(5) * t565 - qJDD(4) * mrSges(5,2) - qJD(4) * t612 + (-t594 - t595) * t623 + (-mrSges(5,3) - mrSges(6,1)) * t602 + t658;
t538 = t651 * t543 + t700 * t550;
t668 = -qJDD(1) * mrSges(4,3) - t654 * (mrSges(4,1) * t646 + mrSges(4,2) * t647);
t536 = m(4) * t600 + t647 * t668 + t538;
t672 = -t543 * t700 + t651 * t550;
t537 = m(4) * t601 + t646 * t668 + t672;
t533 = t647 * t536 + t646 * t537;
t621 = -qJDD(1) * pkin(1) + t664;
t661 = -m(3) * t621 + t654 * mrSges(3,3) - t533;
t531 = m(2) * t627 - t654 * mrSges(2,2) + qJDD(1) * t698 + t661;
t617 = t654 * pkin(1) - t704;
t561 = t602 * pkin(4) + t656;
t673 = -t648 * t552 + t650 * t553;
t665 = m(6) * t561 - t603 * mrSges(6,3) - t624 * t614 + t673;
t657 = t602 * t697 + t623 * t686 + t665 + t703;
t655 = -m(3) * t617 + t654 * mrSges(3,2) + qJDD(1) * mrSges(3,3) + t657 + t702;
t541 = (-mrSges(2,1) - t676) * t654 + m(2) * t628 - qJDD(1) * mrSges(2,2) + t655;
t690 = t652 * t531 + t649 * t541;
t689 = qJD(4) * t706 + t623 * t692 - t624 * t694;
t688 = -qJD(4) * t692 + t623 * t707 - t624 * t696;
t687 = t694 * qJD(4) - t696 * t623 + t624 * t708;
t675 = -t531 * t649 + t652 * t541;
t674 = -t646 * t536 + t647 * t537;
t671 = Ifges(4,1) * t647 - Ifges(4,4) * t646;
t670 = Ifges(4,4) * t647 - Ifges(4,2) * t646;
t669 = Ifges(4,5) * t647 - Ifges(4,6) * t646;
t544 = -t602 * mrSges(6,2) - t623 * t613 + t665;
t626 = t669 * qJD(1);
t568 = Ifges(7,1) * t606 + Ifges(7,4) * t605 + Ifges(7,5) * t620;
t567 = Ifges(7,4) * t606 + Ifges(7,2) * t605 + Ifges(7,6) * t620;
t566 = Ifges(7,5) * t606 + Ifges(7,6) * t605 + Ifges(7,3) * t620;
t547 = mrSges(7,2) * t559 - mrSges(7,3) * t554 + Ifges(7,1) * t572 + Ifges(7,4) * t571 + Ifges(7,5) * t599 + t566 * t605 - t567 * t620;
t546 = -mrSges(7,1) * t559 + mrSges(7,3) * t555 + Ifges(7,4) * t572 + Ifges(7,2) * t571 + Ifges(7,6) * t599 - t566 * t606 + t568 * t620;
t534 = mrSges(6,1) * t563 + mrSges(7,1) * t554 + mrSges(5,2) * t584 - mrSges(7,2) * t555 - mrSges(5,3) * t564 - mrSges(6,3) * t561 + Ifges(7,5) * t572 + Ifges(7,6) * t571 + Ifges(7,3) * t599 + pkin(5) * t545 - qJ(5) * t544 + t606 * t567 - t605 * t568 + t689 * t623 + t708 * t603 - t696 * t602 + t694 * qJDD(4) + t688 * qJD(4);
t532 = -m(3) * g(3) + t674;
t529 = -mrSges(5,1) * t584 - mrSges(6,1) * t562 + mrSges(6,2) * t561 + mrSges(5,3) * t565 - pkin(4) * t544 - pkin(5) * t662 - pkin(8) * t673 + t687 * qJD(4) + t692 * qJDD(4) - t650 * t546 - t648 * t547 - t602 * t707 + t696 * t603 + t689 * t624;
t528 = -t646 * qJD(1) * t626 + mrSges(4,2) * t609 - mrSges(4,3) * t600 - pkin(7) * t538 + qJDD(1) * t671 - t529 * t700 + t651 * t534;
t527 = -mrSges(4,1) * t609 + mrSges(4,3) * t601 - pkin(3) * t657 + pkin(7) * t672 + qJDD(1) * t670 + t651 * t529 + t534 * t700 - t626 * t684;
t526 = t650 * t547 - t648 * t546 - mrSges(2,3) * t627 + mrSges(3,1) * t621 + mrSges(4,1) * t600 - mrSges(4,2) * t601 - mrSges(6,3) * t562 + mrSges(6,2) * t563 + mrSges(5,1) * t564 - mrSges(5,2) * t565 - pkin(8) * t545 + pkin(3) * t538 + pkin(2) * t533 - qJ(2) * t532 + qJ(5) * t658 + pkin(4) * (-qJD(4) * t613 + t660) + (-qJ(5) * t595 + t687) * t623 - t688 * t624 + (-mrSges(6,2) * pkin(4) - t706) * qJDD(4) + (-mrSges(2,2) + mrSges(3,3)) * g(3) + (-mrSges(6,1) * qJ(5) - t692) * t602 + t694 * t603 + (t669 + t695) * qJDD(1) + (t646 * t671 + t647 * t670 + t693) * t654;
t525 = mrSges(2,3) * t628 - mrSges(3,1) * t617 - t646 * t528 - t647 * t527 - pkin(2) * (-t602 * mrSges(5,1) - t623 * t611 - t544 - t702 - t703) - qJ(3) * t674 - pkin(1) * t532 - t693 * qJDD(1) + t698 * g(3) + (-pkin(2) * t676 + t695) * t654;
t1 = [-m(1) * g(1) + t675; -m(1) * g(2) + t690; (-m(1) - m(2) - m(3)) * g(3) + t674; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t690 - t649 * t525 + t652 * t526; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t675 + t652 * t525 + t649 * t526; pkin(1) * t661 + qJ(2) * (-t654 * t676 + t655) - t646 * t527 - qJ(3) * t533 + mrSges(2,1) * t627 - mrSges(2,2) * t628 + t647 * t528 + mrSges(3,2) * t621 - mrSges(3,3) * t617 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-mrSges(3,2) * pkin(1) + Ifges(3,1) + Ifges(2,3)) * qJDD(1);];
tauB  = t1;
