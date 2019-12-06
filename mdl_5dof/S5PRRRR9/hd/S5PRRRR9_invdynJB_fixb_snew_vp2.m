% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PRRRR9
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PRRRR9_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR9_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR9_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR9_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR9_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR9_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR9_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR9_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR9_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:19:51
% EndTime: 2019-12-05 17:20:01
% DurationCPUTime: 7.45s
% Computational Cost: add. (92748->268), mult. (178386->347), div. (0->0), fcn. (123743->12), ass. (0->119)
t660 = sin(pkin(10));
t662 = cos(pkin(10));
t650 = t660 * g(1) - t662 * g(2);
t651 = -t662 * g(1) - t660 * g(2);
t659 = -g(3) + qJDD(1);
t661 = sin(pkin(5));
t663 = cos(pkin(5));
t667 = sin(qJ(2));
t671 = cos(qJ(2));
t614 = -t667 * t651 + (t650 * t663 + t659 * t661) * t671;
t691 = t663 * t667;
t692 = t661 * t667;
t615 = t650 * t691 + t671 * t651 + t659 * t692;
t673 = qJD(2) ^ 2;
t611 = -t673 * pkin(2) + qJDD(2) * pkin(7) + t615;
t630 = -t661 * t650 + t663 * t659;
t666 = sin(qJ(3));
t670 = cos(qJ(3));
t606 = t670 * t611 + t666 * t630;
t647 = (-pkin(3) * t670 - pkin(8) * t666) * qJD(2);
t672 = qJD(3) ^ 2;
t688 = t670 * qJD(2);
t597 = -t672 * pkin(3) + qJDD(3) * pkin(8) + t647 * t688 + t606;
t610 = -qJDD(2) * pkin(2) - t673 * pkin(7) - t614;
t687 = qJD(2) * qJD(3);
t686 = t670 * t687;
t648 = t666 * qJDD(2) + t686;
t657 = t666 * t687;
t649 = t670 * qJDD(2) - t657;
t600 = (-t648 - t686) * pkin(8) + (-t649 + t657) * pkin(3) + t610;
t665 = sin(qJ(4));
t669 = cos(qJ(4));
t586 = -t665 * t597 + t669 * t600;
t689 = qJD(2) * t666;
t644 = t669 * qJD(3) - t665 * t689;
t622 = t644 * qJD(4) + t665 * qJDD(3) + t669 * t648;
t641 = qJDD(4) - t649;
t645 = t665 * qJD(3) + t669 * t689;
t656 = qJD(4) - t688;
t584 = (t644 * t656 - t622) * pkin(9) + (t644 * t645 + t641) * pkin(4) + t586;
t587 = t669 * t597 + t665 * t600;
t621 = -t645 * qJD(4) + t669 * qJDD(3) - t665 * t648;
t629 = t656 * pkin(4) - t645 * pkin(9);
t640 = t644 ^ 2;
t585 = -t640 * pkin(4) + t621 * pkin(9) - t656 * t629 + t587;
t664 = sin(qJ(5));
t668 = cos(qJ(5));
t583 = t664 * t584 + t668 * t585;
t605 = -t666 * t611 + t670 * t630;
t596 = -qJDD(3) * pkin(3) - t672 * pkin(8) + t647 * t689 - t605;
t588 = -t621 * pkin(4) - t640 * pkin(9) + t645 * t629 + t596;
t624 = t664 * t644 + t668 * t645;
t593 = -t624 * qJD(5) + t668 * t621 - t664 * t622;
t623 = t668 * t644 - t664 * t645;
t594 = t623 * qJD(5) + t664 * t621 + t668 * t622;
t655 = qJD(5) + t656;
t601 = Ifges(6,5) * t624 + Ifges(6,6) * t623 + Ifges(6,3) * t655;
t603 = Ifges(6,1) * t624 + Ifges(6,4) * t623 + Ifges(6,5) * t655;
t637 = qJDD(5) + t641;
t571 = -mrSges(6,1) * t588 + mrSges(6,3) * t583 + Ifges(6,4) * t594 + Ifges(6,2) * t593 + Ifges(6,6) * t637 - t624 * t601 + t655 * t603;
t582 = t668 * t584 - t664 * t585;
t602 = Ifges(6,4) * t624 + Ifges(6,2) * t623 + Ifges(6,6) * t655;
t572 = mrSges(6,2) * t588 - mrSges(6,3) * t582 + Ifges(6,1) * t594 + Ifges(6,4) * t593 + Ifges(6,5) * t637 + t623 * t601 - t655 * t602;
t616 = Ifges(5,5) * t645 + Ifges(5,6) * t644 + Ifges(5,3) * t656;
t618 = Ifges(5,1) * t645 + Ifges(5,4) * t644 + Ifges(5,5) * t656;
t612 = -t655 * mrSges(6,2) + t623 * mrSges(6,3);
t613 = t655 * mrSges(6,1) - t624 * mrSges(6,3);
t678 = m(6) * t588 - t593 * mrSges(6,1) + t594 * mrSges(6,2) - t623 * t612 + t624 * t613;
t607 = -t623 * mrSges(6,1) + t624 * mrSges(6,2);
t578 = m(6) * t582 + t637 * mrSges(6,1) - t594 * mrSges(6,3) - t624 * t607 + t655 * t612;
t579 = m(6) * t583 - t637 * mrSges(6,2) + t593 * mrSges(6,3) + t623 * t607 - t655 * t613;
t683 = -t664 * t578 + t668 * t579;
t554 = -mrSges(5,1) * t596 + mrSges(5,3) * t587 + Ifges(5,4) * t622 + Ifges(5,2) * t621 + Ifges(5,6) * t641 - pkin(4) * t678 + pkin(9) * t683 + t668 * t571 + t664 * t572 - t645 * t616 + t656 * t618;
t570 = t668 * t578 + t664 * t579;
t617 = Ifges(5,4) * t645 + Ifges(5,2) * t644 + Ifges(5,6) * t656;
t558 = mrSges(5,2) * t596 - mrSges(5,3) * t586 + Ifges(5,1) * t622 + Ifges(5,4) * t621 + Ifges(5,5) * t641 - pkin(9) * t570 - t664 * t571 + t668 * t572 + t644 * t616 - t656 * t617;
t625 = -t644 * mrSges(5,1) + t645 * mrSges(5,2);
t627 = -t656 * mrSges(5,2) + t644 * mrSges(5,3);
t568 = m(5) * t586 + t641 * mrSges(5,1) - t622 * mrSges(5,3) - t645 * t625 + t656 * t627 + t570;
t628 = t656 * mrSges(5,1) - t645 * mrSges(5,3);
t569 = m(5) * t587 - t641 * mrSges(5,2) + t621 * mrSges(5,3) + t644 * t625 - t656 * t628 + t683;
t566 = -t665 * t568 + t669 * t569;
t580 = -m(5) * t596 + t621 * mrSges(5,1) - t622 * mrSges(5,2) + t644 * t627 - t645 * t628 - t678;
t635 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t666 + Ifges(4,2) * t670) * qJD(2);
t636 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t666 + Ifges(4,4) * t670) * qJD(2);
t694 = mrSges(4,1) * t605 - mrSges(4,2) * t606 + Ifges(4,5) * t648 + Ifges(4,6) * t649 + Ifges(4,3) * qJDD(3) + pkin(3) * t580 + pkin(8) * t566 + t669 * t554 + t665 * t558 + (t666 * t635 - t670 * t636) * qJD(2);
t565 = t669 * t568 + t665 * t569;
t652 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t689;
t653 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t688;
t676 = -m(4) * t610 + t649 * mrSges(4,1) - t648 * mrSges(4,2) - t652 * t689 + t653 * t688 - t565;
t561 = m(3) * t614 + qJDD(2) * mrSges(3,1) - t673 * mrSges(3,2) + t676;
t693 = t561 * t671;
t646 = (-mrSges(4,1) * t670 + mrSges(4,2) * t666) * qJD(2);
t564 = m(4) * t606 - qJDD(3) * mrSges(4,2) + t649 * mrSges(4,3) - qJD(3) * t652 + t646 * t688 + t566;
t574 = m(4) * t605 + qJDD(3) * mrSges(4,1) - t648 * mrSges(4,3) + qJD(3) * t653 - t646 * t689 + t580;
t684 = t670 * t564 - t666 * t574;
t553 = m(3) * t615 - t673 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t684;
t557 = t666 * t564 + t670 * t574;
t556 = m(3) * t630 + t557;
t543 = t553 * t691 - t661 * t556 + t663 * t693;
t541 = m(2) * t650 + t543;
t548 = t671 * t553 - t667 * t561;
t547 = m(2) * t651 + t548;
t690 = t662 * t541 + t660 * t547;
t542 = t553 * t692 + t663 * t556 + t661 * t693;
t685 = -t660 * t541 + t662 * t547;
t682 = m(2) * t659 + t542;
t634 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t666 + Ifges(4,6) * t670) * qJD(2);
t544 = mrSges(4,2) * t610 - mrSges(4,3) * t605 + Ifges(4,1) * t648 + Ifges(4,4) * t649 + Ifges(4,5) * qJDD(3) - pkin(8) * t565 - qJD(3) * t635 - t665 * t554 + t669 * t558 + t634 * t688;
t677 = -mrSges(6,1) * t582 + mrSges(6,2) * t583 - Ifges(6,5) * t594 - Ifges(6,6) * t593 - Ifges(6,3) * t637 - t624 * t602 + t623 * t603;
t674 = mrSges(5,1) * t586 - mrSges(5,2) * t587 + Ifges(5,5) * t622 + Ifges(5,6) * t621 + Ifges(5,3) * t641 + pkin(4) * t570 + t645 * t617 - t644 * t618 - t677;
t549 = -mrSges(4,1) * t610 + mrSges(4,3) * t606 + Ifges(4,4) * t648 + Ifges(4,2) * t649 + Ifges(4,6) * qJDD(3) - pkin(3) * t565 + qJD(3) * t636 - t634 * t689 - t674;
t538 = mrSges(3,2) * t630 - mrSges(3,3) * t614 + Ifges(3,5) * qJDD(2) - t673 * Ifges(3,6) - pkin(7) * t557 + t670 * t544 - t666 * t549;
t539 = -mrSges(3,1) * t630 + mrSges(3,3) * t615 + t673 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t557 - t694;
t679 = pkin(6) * t548 + t538 * t667 + t539 * t671;
t537 = mrSges(3,1) * t614 - mrSges(3,2) * t615 + Ifges(3,3) * qJDD(2) + pkin(2) * t676 + pkin(7) * t684 + t666 * t544 + t670 * t549;
t536 = mrSges(2,2) * t659 - mrSges(2,3) * t650 + t671 * t538 - t667 * t539 + (-t542 * t661 - t543 * t663) * pkin(6);
t535 = -mrSges(2,1) * t659 + mrSges(2,3) * t651 - pkin(1) * t542 - t661 * t537 + t679 * t663;
t1 = [-m(1) * g(1) + t685; -m(1) * g(2) + t690; -m(1) * g(3) + t682; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t690 - t660 * t535 + t662 * t536; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t685 + t662 * t535 + t660 * t536; -mrSges(1,1) * g(2) + mrSges(2,1) * t650 + mrSges(1,2) * g(1) - mrSges(2,2) * t651 + pkin(1) * t543 + t663 * t537 + t679 * t661; t682; t537; t694; t674; -t677;];
tauJB = t1;
