% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPRRR6
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPRRR6_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR6_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR6_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR6_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR6_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR6_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR6_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR6_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR6_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:01:17
% EndTime: 2019-12-31 19:01:22
% DurationCPUTime: 4.50s
% Computational Cost: add. (54526->271), mult. (106635->343), div. (0->0), fcn. (67645->10), ass. (0->112)
t662 = sin(qJ(4));
t663 = sin(qJ(3));
t666 = cos(qJ(4));
t667 = cos(qJ(3));
t631 = (t662 * t663 - t666 * t667) * qJD(1);
t664 = sin(qJ(1));
t668 = cos(qJ(1));
t645 = t664 * g(1) - t668 * g(2);
t636 = qJDD(1) * pkin(1) + t645;
t646 = -t668 * g(1) - t664 * g(2);
t669 = qJD(1) ^ 2;
t638 = -t669 * pkin(1) + t646;
t659 = sin(pkin(9));
t660 = cos(pkin(9));
t619 = t659 * t636 + t660 * t638;
t614 = -t669 * pkin(2) + qJDD(1) * pkin(6) + t619;
t658 = -g(3) + qJDD(2);
t602 = -t663 * t614 + t667 * t658;
t686 = qJD(1) * qJD(3);
t685 = t667 * t686;
t639 = t663 * qJDD(1) + t685;
t589 = (-t639 + t685) * pkin(7) + (t663 * t667 * t669 + qJDD(3)) * pkin(3) + t602;
t603 = t667 * t614 + t663 * t658;
t640 = t667 * qJDD(1) - t663 * t686;
t688 = qJD(1) * t663;
t644 = qJD(3) * pkin(3) - pkin(7) * t688;
t657 = t667 ^ 2;
t590 = -t657 * t669 * pkin(3) + t640 * pkin(7) - qJD(3) * t644 + t603;
t583 = t662 * t589 + t666 * t590;
t632 = (t662 * t667 + t663 * t666) * qJD(1);
t604 = -t632 * qJD(4) - t662 * t639 + t666 * t640;
t615 = t631 * mrSges(5,1) + t632 * mrSges(5,2);
t654 = qJD(3) + qJD(4);
t623 = t654 * mrSges(5,1) - t632 * mrSges(5,3);
t653 = qJDD(3) + qJDD(4);
t616 = t631 * pkin(4) - t632 * pkin(8);
t652 = t654 ^ 2;
t579 = -t652 * pkin(4) + t653 * pkin(8) - t631 * t616 + t583;
t618 = t660 * t636 - t659 * t638;
t677 = -qJDD(1) * pkin(2) - t618;
t595 = -t640 * pkin(3) + t644 * t688 + (-pkin(7) * t657 - pkin(6)) * t669 + t677;
t605 = -t631 * qJD(4) + t666 * t639 + t662 * t640;
t580 = (t631 * t654 - t605) * pkin(8) + (t632 * t654 - t604) * pkin(4) + t595;
t661 = sin(qJ(5));
t665 = cos(qJ(5));
t576 = -t661 * t579 + t665 * t580;
t620 = -t661 * t632 + t665 * t654;
t586 = t620 * qJD(5) + t665 * t605 + t661 * t653;
t621 = t665 * t632 + t661 * t654;
t596 = -t620 * mrSges(6,1) + t621 * mrSges(6,2);
t601 = qJDD(5) - t604;
t624 = qJD(5) + t631;
t606 = -t624 * mrSges(6,2) + t620 * mrSges(6,3);
t572 = m(6) * t576 + t601 * mrSges(6,1) - t586 * mrSges(6,3) - t621 * t596 + t624 * t606;
t577 = t665 * t579 + t661 * t580;
t585 = -t621 * qJD(5) - t661 * t605 + t665 * t653;
t607 = t624 * mrSges(6,1) - t621 * mrSges(6,3);
t573 = m(6) * t577 - t601 * mrSges(6,2) + t585 * mrSges(6,3) + t620 * t596 - t624 * t607;
t680 = -t661 * t572 + t665 * t573;
t559 = m(5) * t583 - t653 * mrSges(5,2) + t604 * mrSges(5,3) - t631 * t615 - t654 * t623 + t680;
t582 = t666 * t589 - t662 * t590;
t622 = -t654 * mrSges(5,2) - t631 * mrSges(5,3);
t578 = -t653 * pkin(4) - t652 * pkin(8) + t632 * t616 - t582;
t676 = -m(6) * t578 + t585 * mrSges(6,1) - t586 * mrSges(6,2) + t620 * t606 - t621 * t607;
t568 = m(5) * t582 + t653 * mrSges(5,1) - t605 * mrSges(5,3) - t632 * t615 + t654 * t622 + t676;
t553 = t662 * t559 + t666 * t568;
t629 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t663 + Ifges(4,2) * t667) * qJD(1);
t630 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t663 + Ifges(4,4) * t667) * qJD(1);
t591 = Ifges(6,5) * t621 + Ifges(6,6) * t620 + Ifges(6,3) * t624;
t593 = Ifges(6,1) * t621 + Ifges(6,4) * t620 + Ifges(6,5) * t624;
t565 = -mrSges(6,1) * t578 + mrSges(6,3) * t577 + Ifges(6,4) * t586 + Ifges(6,2) * t585 + Ifges(6,6) * t601 - t621 * t591 + t624 * t593;
t592 = Ifges(6,4) * t621 + Ifges(6,2) * t620 + Ifges(6,6) * t624;
t566 = mrSges(6,2) * t578 - mrSges(6,3) * t576 + Ifges(6,1) * t586 + Ifges(6,4) * t585 + Ifges(6,5) * t601 + t620 * t591 - t624 * t592;
t610 = Ifges(5,4) * t632 - Ifges(5,2) * t631 + Ifges(5,6) * t654;
t611 = Ifges(5,1) * t632 - Ifges(5,4) * t631 + Ifges(5,5) * t654;
t673 = -mrSges(5,1) * t582 + mrSges(5,2) * t583 - Ifges(5,5) * t605 - Ifges(5,6) * t604 - Ifges(5,3) * t653 - pkin(4) * t676 - pkin(8) * t680 - t665 * t565 - t661 * t566 - t632 * t610 - t631 * t611;
t690 = mrSges(4,1) * t602 - mrSges(4,2) * t603 + Ifges(4,5) * t639 + Ifges(4,6) * t640 + Ifges(4,3) * qJDD(3) + pkin(3) * t553 + (t663 * t629 - t667 * t630) * qJD(1) - t673;
t637 = (-mrSges(4,1) * t667 + mrSges(4,2) * t663) * qJD(1);
t687 = qJD(1) * t667;
t643 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t687;
t551 = m(4) * t602 + qJDD(3) * mrSges(4,1) - t639 * mrSges(4,3) + qJD(3) * t643 - t637 * t688 + t553;
t642 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t688;
t681 = t666 * t559 - t662 * t568;
t552 = m(4) * t603 - qJDD(3) * mrSges(4,2) + t640 * mrSges(4,3) - qJD(3) * t642 + t637 * t687 + t681;
t682 = -t663 * t551 + t667 * t552;
t542 = m(3) * t619 - t669 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t682;
t613 = -t669 * pkin(6) + t677;
t561 = t665 * t572 + t661 * t573;
t675 = m(5) * t595 - t604 * mrSges(5,1) + t605 * mrSges(5,2) + t631 * t622 + t632 * t623 + t561;
t671 = -m(4) * t613 + t640 * mrSges(4,1) - t639 * mrSges(4,2) - t642 * t688 + t643 * t687 - t675;
t555 = m(3) * t618 + qJDD(1) * mrSges(3,1) - t669 * mrSges(3,2) + t671;
t539 = t659 * t542 + t660 * t555;
t536 = m(2) * t645 + qJDD(1) * mrSges(2,1) - t669 * mrSges(2,2) + t539;
t683 = t660 * t542 - t659 * t555;
t537 = m(2) * t646 - t669 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t683;
t689 = t668 * t536 + t664 * t537;
t545 = t667 * t551 + t663 * t552;
t543 = m(3) * t658 + t545;
t684 = -t664 * t536 + t668 * t537;
t609 = Ifges(5,5) * t632 - Ifges(5,6) * t631 + Ifges(5,3) * t654;
t546 = mrSges(5,2) * t595 - mrSges(5,3) * t582 + Ifges(5,1) * t605 + Ifges(5,4) * t604 + Ifges(5,5) * t653 - pkin(8) * t561 - t661 * t565 + t665 * t566 - t631 * t609 - t654 * t610;
t672 = mrSges(6,1) * t576 - mrSges(6,2) * t577 + Ifges(6,5) * t586 + Ifges(6,6) * t585 + Ifges(6,3) * t601 + t621 * t592 - t620 * t593;
t547 = -mrSges(5,1) * t595 + mrSges(5,3) * t583 + Ifges(5,4) * t605 + Ifges(5,2) * t604 + Ifges(5,6) * t653 - pkin(4) * t561 - t632 * t609 + t654 * t611 - t672;
t628 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t663 + Ifges(4,6) * t667) * qJD(1);
t530 = -mrSges(4,1) * t613 + mrSges(4,3) * t603 + Ifges(4,4) * t639 + Ifges(4,2) * t640 + Ifges(4,6) * qJDD(3) - pkin(3) * t675 + pkin(7) * t681 + qJD(3) * t630 + t662 * t546 + t666 * t547 - t628 * t688;
t532 = mrSges(4,2) * t613 - mrSges(4,3) * t602 + Ifges(4,1) * t639 + Ifges(4,4) * t640 + Ifges(4,5) * qJDD(3) - pkin(7) * t553 - qJD(3) * t629 + t666 * t546 - t662 * t547 + t628 * t687;
t674 = mrSges(2,1) * t645 + mrSges(3,1) * t618 - mrSges(2,2) * t646 - mrSges(3,2) * t619 + pkin(1) * t539 + pkin(2) * t671 + pkin(6) * t682 + t667 * t530 + t663 * t532 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t528 = -mrSges(3,1) * t658 + mrSges(3,3) * t619 + t669 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t545 - t690;
t527 = mrSges(3,2) * t658 - mrSges(3,3) * t618 + Ifges(3,5) * qJDD(1) - t669 * Ifges(3,6) - pkin(6) * t545 - t663 * t530 + t667 * t532;
t526 = -mrSges(2,2) * g(3) - mrSges(2,3) * t645 + Ifges(2,5) * qJDD(1) - t669 * Ifges(2,6) - qJ(2) * t539 + t660 * t527 - t659 * t528;
t525 = mrSges(2,1) * g(3) + mrSges(2,3) * t646 + t669 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t543 + qJ(2) * t683 + t659 * t527 + t660 * t528;
t1 = [-m(1) * g(1) + t684; -m(1) * g(2) + t689; (-m(1) - m(2)) * g(3) + t543; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t689 - t664 * t525 + t668 * t526; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t684 + t668 * t525 + t664 * t526; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t674; t674; t543; t690; -t673; t672;];
tauJB = t1;
