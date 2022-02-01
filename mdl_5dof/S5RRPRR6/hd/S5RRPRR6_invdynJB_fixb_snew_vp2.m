% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRPRR6
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% m [6x1]
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
% Datum: 2022-01-20 11:18
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRPRR6_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR6_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR6_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR6_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR6_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR6_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR6_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:17:12
% EndTime: 2022-01-20 11:17:17
% DurationCPUTime: 4.77s
% Computational Cost: add. (67358->251), mult. (92265->330), div. (0->0), fcn. (56620->10), ass. (0->117)
t668 = sin(qJ(1));
t672 = cos(qJ(1));
t649 = t668 * g(1) - g(2) * t672;
t640 = qJDD(1) * pkin(1) + t649;
t650 = -g(1) * t672 - g(2) * t668;
t673 = qJD(1) ^ 2;
t641 = -pkin(1) * t673 + t650;
t667 = sin(qJ(2));
t671 = cos(qJ(2));
t622 = t667 * t640 + t671 * t641;
t660 = (qJD(1) + qJD(2));
t657 = t660 ^ 2;
t658 = qJDD(1) + qJDD(2);
t707 = -pkin(2) * t657 + qJ(3) * t658 + (2 * qJD(3) * t660) + t622;
t663 = sin(pkin(9));
t664 = cos(pkin(9));
t608 = -t664 * g(3) - t663 * t707;
t706 = mrSges(4,2) * t663;
t705 = mrSges(4,3) * t658;
t704 = t657 * t663 ^ 2;
t703 = t658 * t664;
t702 = t660 * t663;
t701 = t660 * t664;
t666 = sin(qJ(4));
t700 = t663 * t666;
t670 = cos(qJ(4));
t699 = t663 * t670;
t609 = -g(3) * t663 + t664 * t707;
t633 = (-mrSges(4,1) * t664 + t706) * t660;
t685 = -pkin(3) * t664 - pkin(7) * t663;
t635 = t685 * t660;
t595 = t635 * t701 + t609;
t621 = t671 * t640 - t667 * t641;
t681 = -t657 * qJ(3) + qJDD(3) - t621;
t607 = (-pkin(2) + t685) * t658 + t681;
t606 = t670 * t607;
t695 = qJD(4) * t660;
t632 = (t658 * t670 - t666 * t695) * t663;
t645 = qJDD(4) - t703;
t646 = qJD(4) - t701;
t587 = t645 * pkin(4) - t632 * pkin(8) + t606 + (-pkin(4) * t670 * t704 - pkin(8) * t646 * t702 - t595) * t666;
t590 = t595 * t670 + t607 * t666;
t692 = t660 * t699;
t630 = pkin(4) * t646 - pkin(8) * t692;
t631 = (-t658 * t666 - t670 * t695) * t663;
t694 = t666 ^ 2 * t704;
t588 = -pkin(4) * t694 + pkin(8) * t631 - t630 * t646 + t590;
t665 = sin(qJ(5));
t669 = cos(qJ(5));
t585 = t587 * t669 - t588 * t665;
t623 = (-t665 * t670 - t666 * t669) * t702;
t600 = qJD(5) * t623 + t631 * t665 + t632 * t669;
t624 = (-t665 * t666 + t669 * t670) * t702;
t610 = -mrSges(6,1) * t623 + mrSges(6,2) * t624;
t644 = qJD(5) + t646;
t615 = -mrSges(6,2) * t644 + mrSges(6,3) * t623;
t642 = qJDD(5) + t645;
t582 = m(6) * t585 + mrSges(6,1) * t642 - t600 * mrSges(6,3) - t610 * t624 + t615 * t644;
t586 = t587 * t665 + t588 * t669;
t599 = -qJD(5) * t624 + t631 * t669 - t632 * t665;
t616 = mrSges(6,1) * t644 - mrSges(6,3) * t624;
t583 = m(6) * t586 - mrSges(6,2) * t642 + t599 * mrSges(6,3) + t610 * t623 - t616 * t644;
t574 = t582 * t669 + t583 * t665;
t589 = -t666 * t595 + t606;
t693 = t660 * t700;
t627 = -mrSges(5,2) * t646 - mrSges(5,3) * t693;
t629 = (mrSges(5,1) * t666 + mrSges(5,2) * t670) * t702;
t572 = m(5) * t589 + mrSges(5,1) * t645 - mrSges(5,3) * t632 + t627 * t646 - t629 * t692 + t574;
t628 = mrSges(5,1) * t646 - mrSges(5,3) * t692;
t686 = -t582 * t665 + t583 * t669;
t573 = m(5) * t590 - mrSges(5,2) * t645 + mrSges(5,3) * t631 - t628 * t646 - t629 * t693 + t686;
t687 = -t666 * t572 + t573 * t670;
t567 = m(4) * t609 + (t633 * t660 + t705) * t664 + t687;
t594 = t635 * t702 - t608;
t591 = -pkin(4) * t631 - pkin(8) * t694 + t630 * t692 + t594;
t679 = m(6) * t591 - t599 * mrSges(6,1) + t600 * mrSges(6,2) - t615 * t623 + t624 * t616;
t675 = -m(5) * t594 + t631 * mrSges(5,1) - t632 * mrSges(5,2) - t679;
t578 = m(4) * t608 + (-t705 + (-t627 * t666 - t628 * t670 - t633) * t660) * t663 + t675;
t688 = t567 * t664 - t578 * t663;
t558 = m(3) * t622 - mrSges(3,1) * t657 - mrSges(3,2) * t658 + t688;
t570 = t670 * t572 + t666 * t573;
t613 = -pkin(2) * t658 + t681;
t677 = -m(4) * t613 + mrSges(4,1) * t703 - t570 + (t657 * t664 ^ 2 + t704) * mrSges(4,3);
t564 = m(3) * t621 - t657 * mrSges(3,2) + (mrSges(3,1) - t706) * t658 + t677;
t553 = t558 * t667 + t564 * t671;
t550 = m(2) * t649 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t673 + t553;
t689 = t558 * t671 - t564 * t667;
t551 = m(2) * t650 - mrSges(2,1) * t673 - qJDD(1) * mrSges(2,2) + t689;
t698 = t550 * t672 + t551 * t668;
t561 = t567 * t663 + t578 * t664;
t690 = -t550 * t668 + t551 * t672;
t684 = Ifges(4,1) * t663 + Ifges(4,4) * t664;
t683 = Ifges(4,5) * t663 + Ifges(4,6) * t664;
t618 = Ifges(5,6) * t646 + (Ifges(5,4) * t670 - Ifges(5,2) * t666) * t702;
t619 = Ifges(5,5) * t646 + (Ifges(5,1) * t670 - Ifges(5,4) * t666) * t702;
t682 = t618 * t670 + t619 * t666;
t601 = Ifges(6,5) * t624 + Ifges(6,6) * t623 + Ifges(6,3) * t644;
t603 = Ifges(6,1) * t624 + Ifges(6,4) * t623 + Ifges(6,5) * t644;
t575 = -mrSges(6,1) * t591 + mrSges(6,3) * t586 + Ifges(6,4) * t600 + Ifges(6,2) * t599 + Ifges(6,6) * t642 - t601 * t624 + t603 * t644;
t602 = Ifges(6,4) * t624 + Ifges(6,2) * t623 + Ifges(6,6) * t644;
t576 = mrSges(6,2) * t591 - mrSges(6,3) * t585 + Ifges(6,1) * t600 + Ifges(6,4) * t599 + Ifges(6,5) * t642 + t601 * t623 - t602 * t644;
t617 = Ifges(5,3) * t646 + (Ifges(5,5) * t670 - Ifges(5,6) * t666) * t702;
t559 = -mrSges(5,1) * t594 + mrSges(5,3) * t590 + Ifges(5,4) * t632 + Ifges(5,2) * t631 + Ifges(5,6) * t645 - pkin(4) * t679 + pkin(8) * t686 + t669 * t575 + t665 * t576 - t617 * t692 + t646 * t619;
t562 = mrSges(5,2) * t594 - mrSges(5,3) * t589 + Ifges(5,1) * t632 + Ifges(5,4) * t631 + Ifges(5,5) * t645 - pkin(8) * t574 - t575 * t665 + t576 * t669 - t617 * t693 - t618 * t646;
t634 = t683 * t660;
t546 = mrSges(4,2) * t613 - mrSges(4,3) * t608 - pkin(7) * t570 - t666 * t559 + t670 * t562 + t634 * t701 + t658 * t684;
t678 = -mrSges(6,1) * t585 + mrSges(6,2) * t586 - Ifges(6,5) * t600 - Ifges(6,6) * t599 - Ifges(6,3) * t642 - t602 * t624 + t623 * t603;
t674 = mrSges(5,1) * t589 - mrSges(5,2) * t590 + Ifges(5,5) * t632 + Ifges(5,6) * t631 + Ifges(5,3) * t645 + pkin(4) * t574 - t678;
t555 = -pkin(3) * t570 + Ifges(4,2) * t703 - t674 + (Ifges(4,4) * t658 + (-t634 - t682) * t660) * t663 + mrSges(4,3) * t609 - mrSges(4,1) * t613;
t569 = t658 * t706 - t677;
t680 = mrSges(3,1) * t621 - mrSges(3,2) * t622 + Ifges(3,3) * t658 - pkin(2) * t569 + qJ(3) * t688 + t546 * t663 + t555 * t664;
t676 = mrSges(2,1) * t649 - mrSges(2,2) * t650 + Ifges(2,3) * qJDD(1) + pkin(1) * t553 + t680;
t544 = t657 * Ifges(3,5) + mrSges(3,1) * g(3) + mrSges(3,3) * t622 - mrSges(4,1) * t608 + mrSges(4,2) * t609 - t666 * t562 - t670 * t559 - pkin(3) * t675 - pkin(7) * t687 - pkin(2) * t561 + (Ifges(3,6) - t683) * t658 + (-pkin(3) * (-t627 * t700 - t628 * t699) + (-t663 * (Ifges(4,4) * t663 + Ifges(4,2) * t664) + t664 * t684) * t660) * t660;
t543 = -mrSges(3,2) * g(3) - mrSges(3,3) * t621 + Ifges(3,5) * t658 - Ifges(3,6) * t657 - qJ(3) * t561 + t546 * t664 - t555 * t663;
t542 = -mrSges(2,2) * g(3) - mrSges(2,3) * t649 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t673 - pkin(6) * t553 + t543 * t671 - t544 * t667;
t541 = Ifges(2,6) * qJDD(1) + t673 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t650 + t667 * t543 + t671 * t544 - pkin(1) * (-m(3) * g(3) + t561) + pkin(6) * t689;
t1 = [-m(1) * g(1) + t690; -m(1) * g(2) + t698; (-m(1) - m(2) - m(3)) * g(3) + t561; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t698 - t668 * t541 + t542 * t672; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t690 + t672 * t541 + t668 * t542; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t676; t676; t680; t569; t682 * t702 + t674; -t678;];
tauJB = t1;
