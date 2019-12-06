% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRRPR3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-05 18:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRRPR3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR3_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR3_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR3_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR3_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR3_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:42:43
% EndTime: 2019-12-05 18:42:49
% DurationCPUTime: 6.23s
% Computational Cost: add. (107420->273), mult. (143754->347), div. (0->0), fcn. (91374->10), ass. (0->112)
t668 = qJDD(1) + qJDD(2);
t676 = sin(qJ(3));
t680 = cos(qJ(3));
t670 = qJD(1) + qJD(2);
t698 = qJD(3) * t670;
t648 = t676 * t668 + t680 * t698;
t678 = sin(qJ(1));
t682 = cos(qJ(1));
t660 = t682 * g(2) + t678 * g(3);
t653 = qJDD(1) * pkin(1) + t660;
t659 = t678 * g(2) - t682 * g(3);
t683 = qJD(1) ^ 2;
t654 = -t683 * pkin(1) + t659;
t677 = sin(qJ(2));
t681 = cos(qJ(2));
t632 = t677 * t653 + t681 * t654;
t666 = t670 ^ 2;
t626 = -t666 * pkin(2) + t668 * pkin(7) + t632;
t699 = t676 * t626;
t703 = pkin(3) * t666;
t610 = qJDD(3) * pkin(3) - t648 * qJ(4) - t699 + (qJ(4) * t698 + t676 * t703 - g(1)) * t680;
t616 = -t676 * g(1) + t680 * t626;
t649 = t680 * t668 - t676 * t698;
t701 = t670 * t676;
t655 = qJD(3) * pkin(3) - qJ(4) * t701;
t672 = t680 ^ 2;
t611 = t649 * qJ(4) - qJD(3) * t655 - t672 * t703 + t616;
t673 = sin(pkin(9));
t674 = cos(pkin(9));
t640 = (t673 * t680 + t674 * t676) * t670;
t590 = -0.2e1 * qJD(4) * t640 + t674 * t610 - t673 * t611;
t629 = t674 * t648 + t673 * t649;
t639 = (-t673 * t676 + t674 * t680) * t670;
t588 = (qJD(3) * t639 - t629) * pkin(8) + (t639 * t640 + qJDD(3)) * pkin(4) + t590;
t591 = 0.2e1 * qJD(4) * t639 + t673 * t610 + t674 * t611;
t628 = -t673 * t648 + t674 * t649;
t635 = qJD(3) * pkin(4) - t640 * pkin(8);
t638 = t639 ^ 2;
t589 = -t638 * pkin(4) + t628 * pkin(8) - qJD(3) * t635 + t591;
t675 = sin(qJ(5));
t679 = cos(qJ(5));
t586 = t679 * t588 - t675 * t589;
t620 = t679 * t639 - t675 * t640;
t600 = t620 * qJD(5) + t675 * t628 + t679 * t629;
t621 = t675 * t639 + t679 * t640;
t606 = -t620 * mrSges(6,1) + t621 * mrSges(6,2);
t669 = qJD(3) + qJD(5);
t613 = -t669 * mrSges(6,2) + t620 * mrSges(6,3);
t667 = qJDD(3) + qJDD(5);
t582 = m(6) * t586 + t667 * mrSges(6,1) - t600 * mrSges(6,3) - t621 * t606 + t669 * t613;
t587 = t675 * t588 + t679 * t589;
t599 = -t621 * qJD(5) + t679 * t628 - t675 * t629;
t614 = t669 * mrSges(6,1) - t621 * mrSges(6,3);
t583 = m(6) * t587 - t667 * mrSges(6,2) + t599 * mrSges(6,3) + t620 * t606 - t669 * t614;
t573 = t679 * t582 + t675 * t583;
t624 = -t639 * mrSges(5,1) + t640 * mrSges(5,2);
t633 = -qJD(3) * mrSges(5,2) + t639 * mrSges(5,3);
t571 = m(5) * t590 + qJDD(3) * mrSges(5,1) - t629 * mrSges(5,3) + qJD(3) * t633 - t640 * t624 + t573;
t634 = qJD(3) * mrSges(5,1) - t640 * mrSges(5,3);
t693 = -t675 * t582 + t679 * t583;
t572 = m(5) * t591 - qJDD(3) * mrSges(5,2) + t628 * mrSges(5,3) - qJD(3) * t634 + t639 * t624 + t693;
t567 = t674 * t571 + t673 * t572;
t615 = -t680 * g(1) - t699;
t618 = Ifges(5,4) * t640 + Ifges(5,2) * t639 + Ifges(5,6) * qJD(3);
t619 = Ifges(5,1) * t640 + Ifges(5,4) * t639 + Ifges(5,5) * qJD(3);
t642 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t676 + Ifges(4,2) * t680) * t670;
t643 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t676 + Ifges(4,4) * t680) * t670;
t602 = Ifges(6,4) * t621 + Ifges(6,2) * t620 + Ifges(6,6) * t669;
t603 = Ifges(6,1) * t621 + Ifges(6,4) * t620 + Ifges(6,5) * t669;
t687 = -mrSges(6,1) * t586 + mrSges(6,2) * t587 - Ifges(6,5) * t600 - Ifges(6,6) * t599 - Ifges(6,3) * t667 - t621 * t602 + t620 * t603;
t704 = mrSges(4,1) * t615 + mrSges(5,1) * t590 - mrSges(4,2) * t616 - mrSges(5,2) * t591 + Ifges(4,5) * t648 + Ifges(5,5) * t629 + Ifges(4,6) * t649 + Ifges(5,6) * t628 + pkin(3) * t567 + pkin(4) * t573 + (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + t640 * t618 - t639 * t619 + (t676 * t642 - t680 * t643) * t670 - t687;
t700 = t670 * t680;
t647 = (-mrSges(4,1) * t680 + mrSges(4,2) * t676) * t670;
t657 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t700;
t565 = m(4) * t615 + qJDD(3) * mrSges(4,1) - t648 * mrSges(4,3) + qJD(3) * t657 - t647 * t701 + t567;
t656 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t701;
t694 = -t673 * t571 + t674 * t572;
t566 = m(4) * t616 - qJDD(3) * mrSges(4,2) + t649 * mrSges(4,3) - qJD(3) * t656 + t647 * t700 + t694;
t695 = -t676 * t565 + t680 * t566;
t557 = m(3) * t632 - t666 * mrSges(3,1) - t668 * mrSges(3,2) + t695;
t631 = t681 * t653 - t677 * t654;
t689 = -t668 * pkin(2) - t631;
t612 = -t649 * pkin(3) + qJDD(4) + t655 * t701 + (-qJ(4) * t672 - pkin(7)) * t666 + t689;
t593 = -t628 * pkin(4) - t638 * pkin(8) + t640 * t635 + t612;
t692 = m(6) * t593 - t599 * mrSges(6,1) + t600 * mrSges(6,2) - t620 * t613 + t621 * t614;
t584 = m(5) * t612 - t628 * mrSges(5,1) + t629 * mrSges(5,2) - t639 * t633 + t640 * t634 + t692;
t625 = -t666 * pkin(7) + t689;
t685 = -m(4) * t625 + t649 * mrSges(4,1) - t648 * mrSges(4,2) - t656 * t701 + t657 * t700 - t584;
t577 = m(3) * t631 + t668 * mrSges(3,1) - t666 * mrSges(3,2) + t685;
t554 = t677 * t557 + t681 * t577;
t559 = t680 * t565 + t676 * t566;
t696 = t681 * t557 - t677 * t577;
t551 = m(2) * t659 - t683 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t696;
t552 = m(2) * t660 + qJDD(1) * mrSges(2,1) - t683 * mrSges(2,2) + t554;
t697 = t682 * t551 - t678 * t552;
t691 = -t678 * t551 - t682 * t552;
t601 = Ifges(6,5) * t621 + Ifges(6,6) * t620 + Ifges(6,3) * t669;
t574 = -mrSges(6,1) * t593 + mrSges(6,3) * t587 + Ifges(6,4) * t600 + Ifges(6,2) * t599 + Ifges(6,6) * t667 - t621 * t601 + t669 * t603;
t575 = mrSges(6,2) * t593 - mrSges(6,3) * t586 + Ifges(6,1) * t600 + Ifges(6,4) * t599 + Ifges(6,5) * t667 + t620 * t601 - t669 * t602;
t617 = Ifges(5,5) * t640 + Ifges(5,6) * t639 + Ifges(5,3) * qJD(3);
t560 = -mrSges(5,1) * t612 + mrSges(5,3) * t591 + Ifges(5,4) * t629 + Ifges(5,2) * t628 + Ifges(5,6) * qJDD(3) - pkin(4) * t692 + pkin(8) * t693 + qJD(3) * t619 + t679 * t574 + t675 * t575 - t640 * t617;
t561 = mrSges(5,2) * t612 - mrSges(5,3) * t590 + Ifges(5,1) * t629 + Ifges(5,4) * t628 + Ifges(5,5) * qJDD(3) - pkin(8) * t573 - qJD(3) * t618 - t675 * t574 + t679 * t575 + t639 * t617;
t641 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t676 + Ifges(4,6) * t680) * t670;
t547 = -mrSges(4,1) * t625 + mrSges(4,3) * t616 + Ifges(4,4) * t648 + Ifges(4,2) * t649 + Ifges(4,6) * qJDD(3) - pkin(3) * t584 + qJ(4) * t694 + qJD(3) * t643 + t674 * t560 + t673 * t561 - t641 * t701;
t549 = mrSges(4,2) * t625 - mrSges(4,3) * t615 + Ifges(4,1) * t648 + Ifges(4,4) * t649 + Ifges(4,5) * qJDD(3) - qJ(4) * t567 - qJD(3) * t642 - t673 * t560 + t674 * t561 + t641 * t700;
t688 = mrSges(3,1) * t631 - mrSges(3,2) * t632 + Ifges(3,3) * t668 + pkin(2) * t685 + pkin(7) * t695 + t680 * t547 + t676 * t549;
t686 = mrSges(2,1) * t660 - mrSges(2,2) * t659 + Ifges(2,3) * qJDD(1) + pkin(1) * t554 + t688;
t545 = mrSges(3,1) * g(1) + mrSges(3,3) * t632 + t666 * Ifges(3,5) + Ifges(3,6) * t668 - pkin(2) * t559 - t704;
t544 = -mrSges(3,2) * g(1) - mrSges(3,3) * t631 + Ifges(3,5) * t668 - t666 * Ifges(3,6) - pkin(7) * t559 - t676 * t547 + t680 * t549;
t543 = -mrSges(2,2) * g(1) - mrSges(2,3) * t660 + Ifges(2,5) * qJDD(1) - t683 * Ifges(2,6) - pkin(6) * t554 + t681 * t544 - t677 * t545;
t542 = Ifges(2,6) * qJDD(1) + t683 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t659 + t677 * t544 + t681 * t545 - pkin(1) * (-m(3) * g(1) + t559) + pkin(6) * t696;
t1 = [(-m(1) - m(2) - m(3)) * g(1) + t559; -m(1) * g(2) + t691; -m(1) * g(3) + t697; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t686; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t697 - t682 * t542 - t678 * t543; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t691 - t678 * t542 + t682 * t543; t686; t688; t704; t584; -t687;];
tauJB = t1;
