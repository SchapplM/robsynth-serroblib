% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPPRRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
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
% Datum: 2019-05-05 15:42
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPPRRR4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR4_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR4_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR4_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR4_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR4_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR4_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR4_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:40:40
% EndTime: 2019-05-05 15:40:46
% DurationCPUTime: 5.62s
% Computational Cost: add. (83486->298), mult. (151823->359), div. (0->0), fcn. (83190->10), ass. (0->123)
t709 = qJD(1) ^ 2;
t703 = sin(qJ(1));
t707 = cos(qJ(1));
t681 = -t707 * g(1) - t703 * g(2);
t716 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t681;
t730 = -pkin(1) - pkin(2);
t656 = t730 * t709 + t716;
t680 = t703 * g(1) - t707 * g(2);
t715 = -t709 * qJ(2) + qJDD(2) - t680;
t659 = t730 * qJDD(1) + t715;
t698 = sin(pkin(10));
t699 = cos(pkin(10));
t636 = -t698 * t656 + t699 * t659;
t633 = qJDD(1) * pkin(3) - t709 * pkin(7) - t636;
t702 = sin(qJ(4));
t706 = cos(qJ(4));
t724 = qJD(1) * qJD(4);
t722 = t706 * t724;
t675 = -t702 * qJDD(1) - t722;
t723 = t702 * t724;
t676 = -t706 * qJDD(1) + t723;
t620 = (-t675 + t722) * pkin(8) + (-t676 - t723) * pkin(4) + t633;
t637 = t699 * t656 + t698 * t659;
t634 = -t709 * pkin(3) - qJDD(1) * pkin(7) + t637;
t695 = g(3) + qJDD(3);
t630 = t706 * t634 + t702 * t695;
t674 = (pkin(4) * t706 + pkin(8) * t702) * qJD(1);
t686 = t706 * qJD(1);
t708 = qJD(4) ^ 2;
t623 = -t708 * pkin(4) + qJDD(4) * pkin(8) - t674 * t686 + t630;
t701 = sin(qJ(5));
t705 = cos(qJ(5));
t609 = t705 * t620 - t701 * t623;
t725 = qJD(1) * t702;
t671 = t705 * qJD(4) + t701 * t725;
t646 = t671 * qJD(5) + t701 * qJDD(4) + t705 * t675;
t670 = qJDD(5) - t676;
t672 = t701 * qJD(4) - t705 * t725;
t683 = t686 + qJD(5);
t607 = (t671 * t683 - t646) * pkin(9) + (t671 * t672 + t670) * pkin(5) + t609;
t610 = t701 * t620 + t705 * t623;
t645 = -t672 * qJD(5) + t705 * qJDD(4) - t701 * t675;
t655 = t683 * pkin(5) - t672 * pkin(9);
t669 = t671 ^ 2;
t608 = -t669 * pkin(5) + t645 * pkin(9) - t683 * t655 + t610;
t700 = sin(qJ(6));
t704 = cos(qJ(6));
t606 = t700 * t607 + t704 * t608;
t629 = -t702 * t634 + t706 * t695;
t622 = -qJDD(4) * pkin(4) - t708 * pkin(8) - t674 * t725 - t629;
t611 = -t645 * pkin(5) - t669 * pkin(9) + t672 * t655 + t622;
t648 = t700 * t671 + t704 * t672;
t616 = -t648 * qJD(6) + t704 * t645 - t700 * t646;
t647 = t704 * t671 - t700 * t672;
t617 = t647 * qJD(6) + t700 * t645 + t704 * t646;
t682 = qJD(6) + t683;
t624 = Ifges(7,5) * t648 + Ifges(7,6) * t647 + Ifges(7,3) * t682;
t626 = Ifges(7,1) * t648 + Ifges(7,4) * t647 + Ifges(7,5) * t682;
t666 = qJDD(6) + t670;
t595 = -mrSges(7,1) * t611 + mrSges(7,3) * t606 + Ifges(7,4) * t617 + Ifges(7,2) * t616 + Ifges(7,6) * t666 - t648 * t624 + t682 * t626;
t605 = t704 * t607 - t700 * t608;
t625 = Ifges(7,4) * t648 + Ifges(7,2) * t647 + Ifges(7,6) * t682;
t596 = mrSges(7,2) * t611 - mrSges(7,3) * t605 + Ifges(7,1) * t617 + Ifges(7,4) * t616 + Ifges(7,5) * t666 + t647 * t624 - t682 * t625;
t640 = Ifges(6,5) * t672 + Ifges(6,6) * t671 + Ifges(6,3) * t683;
t642 = Ifges(6,1) * t672 + Ifges(6,4) * t671 + Ifges(6,5) * t683;
t638 = -t682 * mrSges(7,2) + t647 * mrSges(7,3);
t639 = t682 * mrSges(7,1) - t648 * mrSges(7,3);
t714 = m(7) * t611 - t616 * mrSges(7,1) + t617 * mrSges(7,2) - t647 * t638 + t648 * t639;
t628 = -t647 * mrSges(7,1) + t648 * mrSges(7,2);
t601 = m(7) * t605 + t666 * mrSges(7,1) - t617 * mrSges(7,3) - t648 * t628 + t682 * t638;
t602 = m(7) * t606 - t666 * mrSges(7,2) + t616 * mrSges(7,3) + t647 * t628 - t682 * t639;
t719 = -t700 * t601 + t704 * t602;
t579 = -mrSges(6,1) * t622 + mrSges(6,3) * t610 + Ifges(6,4) * t646 + Ifges(6,2) * t645 + Ifges(6,6) * t670 - pkin(5) * t714 + pkin(9) * t719 + t704 * t595 + t700 * t596 - t672 * t640 + t683 * t642;
t594 = t704 * t601 + t700 * t602;
t641 = Ifges(6,4) * t672 + Ifges(6,2) * t671 + Ifges(6,6) * t683;
t584 = mrSges(6,2) * t622 - mrSges(6,3) * t609 + Ifges(6,1) * t646 + Ifges(6,4) * t645 + Ifges(6,5) * t670 - pkin(9) * t594 - t700 * t595 + t704 * t596 + t671 * t640 - t683 * t641;
t649 = -t671 * mrSges(6,1) + t672 * mrSges(6,2);
t653 = -t683 * mrSges(6,2) + t671 * mrSges(6,3);
t592 = m(6) * t609 + t670 * mrSges(6,1) - t646 * mrSges(6,3) - t672 * t649 + t683 * t653 + t594;
t654 = t683 * mrSges(6,1) - t672 * mrSges(6,3);
t593 = m(6) * t610 - t670 * mrSges(6,2) + t645 * mrSges(6,3) + t671 * t649 - t683 * t654 + t719;
t590 = -t701 * t592 + t705 * t593;
t603 = -m(6) * t622 + t645 * mrSges(6,1) - t646 * mrSges(6,2) + t671 * t653 - t672 * t654 - t714;
t664 = (Ifges(5,6) * qJD(4)) + (-Ifges(5,4) * t702 - Ifges(5,2) * t706) * qJD(1);
t665 = (Ifges(5,5) * qJD(4)) + (-Ifges(5,1) * t702 - Ifges(5,4) * t706) * qJD(1);
t731 = mrSges(5,1) * t629 - mrSges(5,2) * t630 + Ifges(5,5) * t675 + Ifges(5,6) * t676 + Ifges(5,3) * qJDD(4) + pkin(4) * t603 + pkin(8) * t590 - (t702 * t664 - t706 * t665) * qJD(1) + t705 * t579 + t701 * t584;
t729 = mrSges(2,1) + mrSges(3,1);
t728 = Ifges(3,4) + Ifges(2,5);
t727 = Ifges(2,6) - Ifges(3,6);
t660 = -t709 * pkin(1) + t716;
t673 = (mrSges(5,1) * t706 - mrSges(5,2) * t702) * qJD(1);
t678 = qJD(4) * mrSges(5,1) + mrSges(5,3) * t725;
t587 = m(5) * t630 - qJDD(4) * mrSges(5,2) + t676 * mrSges(5,3) - qJD(4) * t678 - t673 * t686 + t590;
t679 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t686;
t597 = m(5) * t629 + qJDD(4) * mrSges(5,1) - t675 * mrSges(5,3) + qJD(4) * t679 + t673 * t725 + t603;
t583 = t706 * t587 - t702 * t597;
t578 = m(4) * t637 - t709 * mrSges(4,1) + qJDD(1) * mrSges(4,2) + t583;
t589 = t705 * t592 + t701 * t593;
t588 = -m(5) * t633 + t676 * mrSges(5,1) - t675 * mrSges(5,2) + t678 * t725 - t679 * t686 - t589;
t585 = m(4) * t636 - qJDD(1) * mrSges(4,1) - t709 * mrSges(4,2) + t588;
t720 = t699 * t578 - t698 * t585;
t717 = m(3) * t660 + qJDD(1) * mrSges(3,3) + t720;
t570 = m(2) * t681 - qJDD(1) * mrSges(2,2) - t729 * t709 + t717;
t575 = t698 * t578 + t699 * t585;
t661 = -qJDD(1) * pkin(1) + t715;
t574 = m(3) * t661 - qJDD(1) * mrSges(3,1) - t709 * mrSges(3,3) + t575;
t571 = m(2) * t680 + qJDD(1) * mrSges(2,1) - t709 * mrSges(2,2) - t574;
t726 = t703 * t570 + t707 * t571;
t721 = t707 * t570 - t703 * t571;
t582 = t702 * t587 + t706 * t597;
t581 = m(4) * t695 + t582;
t713 = -mrSges(7,1) * t605 + mrSges(7,2) * t606 - Ifges(7,5) * t617 - Ifges(7,6) * t616 - Ifges(7,3) * t666 - t648 * t625 + t647 * t626;
t663 = Ifges(5,3) * qJD(4) + (-Ifges(5,5) * t702 - Ifges(5,6) * t706) * qJD(1);
t566 = mrSges(5,2) * t633 - mrSges(5,3) * t629 + Ifges(5,1) * t675 + Ifges(5,4) * t676 + Ifges(5,5) * qJDD(4) - pkin(8) * t589 - qJD(4) * t664 - t701 * t579 + t705 * t584 - t663 * t686;
t710 = mrSges(6,1) * t609 - mrSges(6,2) * t610 + Ifges(6,5) * t646 + Ifges(6,6) * t645 + Ifges(6,3) * t670 + pkin(5) * t594 + t672 * t641 - t671 * t642 - t713;
t576 = -mrSges(5,1) * t633 + mrSges(5,3) * t630 + Ifges(5,4) * t675 + Ifges(5,2) * t676 + Ifges(5,6) * qJDD(4) - pkin(4) * t589 + qJD(4) * t665 + t663 * t725 - t710;
t711 = -mrSges(3,1) * t661 - mrSges(4,1) * t636 - mrSges(2,2) * t681 - pkin(2) * t575 - pkin(3) * t588 - pkin(7) * t583 - t702 * t566 - t706 * t576 + qJ(2) * (-t709 * mrSges(3,1) + t717) - pkin(1) * t574 + mrSges(4,2) * t637 + mrSges(3,3) * t660 + mrSges(2,1) * t680 + (Ifges(3,2) + Ifges(2,3) + Ifges(4,3)) * qJDD(1);
t580 = -m(3) * g(3) - t581;
t565 = -mrSges(4,1) * t695 + mrSges(4,3) * t637 + t709 * Ifges(4,5) - Ifges(4,6) * qJDD(1) - pkin(3) * t582 - t731;
t564 = mrSges(4,2) * t695 - mrSges(4,3) * t636 - Ifges(4,5) * qJDD(1) - t709 * Ifges(4,6) - pkin(7) * t582 + t706 * t566 - t702 * t576;
t563 = mrSges(3,2) * t661 - mrSges(2,3) * t680 - qJ(2) * t580 - qJ(3) * t575 + t699 * t564 - t698 * t565 - t727 * t709 + t728 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t562 = mrSges(3,2) * t660 + mrSges(2,3) * t681 - pkin(1) * t580 + pkin(2) * t581 + t729 * g(3) - qJ(3) * t720 + t727 * qJDD(1) - t698 * t564 - t699 * t565 + t728 * t709;
t1 = [-m(1) * g(1) + t721; -m(1) * g(2) + t726; (-m(1) - m(2) - m(3)) * g(3) - t581; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t726 - t703 * t562 + t707 * t563; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t721 + t707 * t562 + t703 * t563; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t711; t711; t574; t581; t731; t710; -t713;];
tauJB  = t1;
