% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPPRP1
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
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
% Datum: 2022-01-23 09:13
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPPRP1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP1_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP1_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP1_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP1_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP1_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:12:24
% EndTime: 2022-01-23 09:12:27
% DurationCPUTime: 2.16s
% Computational Cost: add. (17596->229), mult. (36515->286), div. (0->0), fcn. (20211->8), ass. (0->107)
t653 = sin(qJ(1));
t655 = cos(qJ(1));
t631 = t653 * g(1) - t655 * g(2);
t627 = qJDD(1) * pkin(1) + t631;
t632 = -t655 * g(1) - t653 * g(2);
t656 = qJD(1) ^ 2;
t628 = -t656 * pkin(1) + t632;
t649 = sin(pkin(7));
t651 = cos(pkin(7));
t600 = t649 * t627 + t651 * t628;
t704 = -t656 * pkin(2) + qJDD(1) * qJ(3) + (2 * qJD(1) * qJD(3)) + t600;
t703 = Ifges(5,1) + Ifges(6,1);
t696 = Ifges(5,4) + Ifges(6,4);
t695 = Ifges(5,5) + Ifges(6,5);
t702 = Ifges(5,2) + Ifges(6,2);
t694 = Ifges(5,6) + Ifges(6,6);
t701 = Ifges(5,3) + Ifges(6,3);
t647 = -g(3) + qJDD(2);
t648 = sin(pkin(8));
t650 = cos(pkin(8));
t589 = t650 * t647 - t704 * t648;
t652 = sin(qJ(4));
t654 = cos(qJ(4));
t681 = t650 * qJD(1);
t634 = qJD(4) - t681;
t682 = t648 * qJD(1);
t685 = (-t696 * t652 + t703 * t654) * t682 + t695 * t634;
t686 = (-t702 * t652 + t696 * t654) * t682 + t694 * t634;
t700 = t685 * t652 + t686 * t654;
t616 = (mrSges(6,1) * t652 + mrSges(6,2) * t654) * t682;
t679 = qJD(1) * qJD(4);
t619 = (qJDD(1) * t654 - t652 * t679) * t648;
t674 = t654 * t682;
t590 = t648 * t647 + t704 * t650;
t665 = -pkin(3) * t650 - pkin(6) * t648;
t626 = t665 * qJD(1);
t588 = t626 * t681 + t590;
t599 = t651 * t627 - t649 * t628;
t661 = -t656 * qJ(3) + qJDD(3) - t599;
t593 = (-pkin(2) + t665) * qJDD(1) + t661;
t592 = t654 * t593;
t678 = t650 * qJDD(1);
t633 = qJDD(4) - t678;
t666 = -0.2e1 * qJD(5) * t682;
t690 = t648 ^ 2 * t656;
t580 = t654 * t666 + t633 * pkin(4) - t619 * qJ(5) + t592 + (-pkin(4) * t654 * t690 - qJ(5) * t634 * t682 - t588) * t652;
t675 = t652 * t682;
t611 = -t634 * mrSges(6,2) - mrSges(6,3) * t675;
t676 = m(6) * t580 + t633 * mrSges(6,1) + t634 * t611;
t577 = -t619 * mrSges(6,3) - t616 * t674 + t676;
t584 = t654 * t588 + t652 * t593;
t613 = t634 * pkin(4) - qJ(5) * t674;
t618 = (-qJDD(1) * t652 - t654 * t679) * t648;
t677 = t652 ^ 2 * t690;
t582 = -pkin(4) * t677 + t618 * qJ(5) - t634 * t613 + t652 * t666 + t584;
t583 = -t652 * t588 + t592;
t699 = mrSges(5,1) * t583 + mrSges(6,1) * t580 - mrSges(5,2) * t584 - mrSges(6,2) * t582 + pkin(4) * t577 + t694 * t618 + t695 * t619 + t701 * t633;
t698 = t650 ^ 2;
t697 = -mrSges(5,2) - mrSges(6,2);
t692 = mrSges(4,2) * t648;
t691 = t619 * mrSges(6,2);
t624 = (-mrSges(4,1) * t650 + t692) * qJD(1);
t612 = -t634 * mrSges(5,2) - mrSges(5,3) * t675;
t664 = (-t616 - (mrSges(5,1) * t652 + mrSges(5,2) * t654) * t682) * t682;
t573 = m(5) * t583 + t633 * mrSges(5,1) + t634 * t612 + (-mrSges(5,3) - mrSges(6,3)) * t619 + t654 * t664 + t676;
t614 = t634 * mrSges(6,1) - mrSges(6,3) * t674;
t684 = -t634 * mrSges(5,1) + mrSges(5,3) * t674 - t614;
t688 = m(6) * t582 + t618 * mrSges(6,3);
t574 = m(5) * t584 + t618 * mrSges(5,3) + t697 * t633 + t684 * t634 + t652 * t664 + t688;
t668 = -t652 * t573 + t654 * t574;
t683 = qJDD(1) * mrSges(4,3);
t567 = m(4) * t590 + (qJD(1) * t624 + t683) * t650 + t668;
t660 = t684 * t654 + (-t611 - t612) * t652;
t587 = t626 * t682 - t589;
t585 = -t618 * pkin(4) - qJ(5) * t677 + t613 * t674 + qJDD(5) + t587;
t672 = m(6) * t585 - t618 * mrSges(6,1);
t662 = -m(5) * t587 + t618 * mrSges(5,1) - t672;
t576 = m(4) * t589 + t697 * t619 + (-t683 + (-t624 + t660) * qJD(1)) * t648 + t662;
t669 = t650 * t567 - t648 * t576;
t558 = m(3) * t600 - t656 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t669;
t571 = t654 * t573 + t652 * t574;
t596 = -qJDD(1) * pkin(2) + t661;
t659 = -m(4) * t596 + mrSges(4,1) * t678 - t571 + (t656 * t698 + t690) * mrSges(4,3);
t564 = m(3) * t599 - t656 * mrSges(3,2) + (mrSges(3,1) - t692) * qJDD(1) + t659;
t553 = t649 * t558 + t651 * t564;
t550 = m(2) * t631 + qJDD(1) * mrSges(2,1) - t656 * mrSges(2,2) + t553;
t670 = t651 * t558 - t649 * t564;
t551 = m(2) * t632 - t656 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t670;
t689 = t655 * t550 + t653 * t551;
t561 = t648 * t567 + t650 * t576;
t687 = (t694 * t652 - t695 * t654) * t682 - t701 * t634;
t559 = m(3) * t647 + t561;
t671 = -t653 * t550 + t655 * t551;
t663 = Ifges(4,5) * t648 + Ifges(4,6) * t650;
t578 = t691 + (t611 * t652 + t614 * t654) * t682 + t672;
t562 = -mrSges(5,1) * t587 + mrSges(5,3) * t584 - mrSges(6,1) * t585 + mrSges(6,3) * t582 - pkin(4) * t578 + qJ(5) * t688 + (-qJ(5) * t614 + t685) * t634 + (-qJ(5) * mrSges(6,2) + t694) * t633 + t696 * t619 + t702 * t618 + (-qJ(5) * t616 * t652 + t687 * t654) * t682;
t570 = mrSges(5,2) * t587 + mrSges(6,2) * t585 - mrSges(5,3) * t583 - mrSges(6,3) * t580 - qJ(5) * t577 + t696 * t618 + t703 * t619 + t695 * t633 - t686 * t634 + t687 * t675;
t625 = t663 * qJD(1);
t546 = t625 * t681 + mrSges(4,2) * t596 - mrSges(4,3) * t589 - pkin(6) * t571 - t652 * t562 + t654 * t570 + (Ifges(4,1) * t648 + Ifges(4,4) * t650) * qJDD(1);
t555 = Ifges(4,2) * t678 - mrSges(4,1) * t596 + mrSges(4,3) * t590 - pkin(3) * t571 + (Ifges(4,4) * qJDD(1) + (-t625 - t700) * qJD(1)) * t648 - t699;
t569 = qJDD(1) * t692 - t659;
t657 = mrSges(2,1) * t631 + mrSges(3,1) * t599 - mrSges(2,2) * t632 - mrSges(3,2) * t600 + pkin(1) * t553 - pkin(2) * t569 + qJ(3) * t669 + t648 * t546 + t650 * t555 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t544 = t656 * Ifges(3,5) - mrSges(3,1) * t647 + mrSges(3,3) * t600 - mrSges(4,1) * t589 + mrSges(4,2) * t590 - t652 * t570 - t654 * t562 - pkin(3) * (-t619 * mrSges(5,2) + t662 - t691) - pkin(6) * t668 - pkin(2) * t561 + (Ifges(3,6) - t663) * qJDD(1) + (Ifges(4,4) * t698 * qJD(1) + (-pkin(3) * t660 + (-Ifges(4,4) * t648 + (Ifges(4,1) - Ifges(4,2)) * t650) * qJD(1)) * t648) * qJD(1);
t543 = mrSges(3,2) * t647 - mrSges(3,3) * t599 + Ifges(3,5) * qJDD(1) - t656 * Ifges(3,6) - qJ(3) * t561 + t650 * t546 - t648 * t555;
t542 = -mrSges(2,2) * g(3) - mrSges(2,3) * t631 + Ifges(2,5) * qJDD(1) - t656 * Ifges(2,6) - qJ(2) * t553 + t651 * t543 - t649 * t544;
t541 = mrSges(2,1) * g(3) + mrSges(2,3) * t632 + t656 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t559 + qJ(2) * t670 + t649 * t543 + t651 * t544;
t1 = [-m(1) * g(1) + t671; -m(1) * g(2) + t689; (-m(1) - m(2)) * g(3) + t559; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t689 - t653 * t541 + t655 * t542; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t671 + t655 * t541 + t653 * t542; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t657; t657; t559; t569; t700 * t682 + t699; t578;];
tauJB = t1;
