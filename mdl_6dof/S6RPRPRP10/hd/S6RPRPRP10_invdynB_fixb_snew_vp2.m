% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRPRP10
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5]';
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
% Datum: 2019-05-05 18:12
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRPRP10_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP10_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP10_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP10_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP10_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPRP10_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP10_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP10_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP10_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 18:09:49
% EndTime: 2019-05-05 18:09:53
% DurationCPUTime: 2.25s
% Computational Cost: add. (16329->297), mult. (31861->335), div. (0->0), fcn. (15796->6), ass. (0->122)
t718 = -2 * qJD(4);
t717 = Ifges(6,1) + Ifges(7,1);
t701 = Ifges(6,4) - Ifges(7,5);
t700 = Ifges(7,4) + Ifges(6,5);
t698 = (Ifges(4,5) - Ifges(5,4));
t716 = Ifges(4,2) + Ifges(5,3);
t715 = Ifges(5,2) + Ifges(4,1);
t714 = Ifges(6,2) + Ifges(7,3);
t696 = (Ifges(4,6) - Ifges(5,5));
t695 = -Ifges(5,6) - Ifges(4,4);
t694 = Ifges(6,6) - Ifges(7,6);
t713 = (-Ifges(4,3) - Ifges(5,1));
t712 = Ifges(6,3) + Ifges(7,2);
t653 = sin(qJ(1));
t655 = cos(qJ(1));
t633 = -t655 * g(1) - t653 * g(2);
t668 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t633;
t657 = qJD(1) ^ 2;
t709 = (-pkin(1) - pkin(7));
t680 = t709 * t657;
t598 = t680 + t668;
t652 = sin(qJ(3));
t654 = cos(qJ(3));
t683 = qJD(1) * qJD(3);
t678 = t654 * t683;
t623 = qJDD(1) * t652 + t678;
t640 = t652 * t683;
t624 = qJDD(1) * t654 - t640;
t685 = qJD(1) * t652;
t627 = -(qJD(3) * mrSges(4,2)) - mrSges(4,3) * t685;
t684 = qJD(1) * t654;
t628 = (qJD(3) * mrSges(4,1)) - mrSges(4,3) * t684;
t662 = pkin(3) * t678 + t684 * t718 + t668 + (-t624 + t640) * qJ(4);
t570 = t623 * pkin(3) + t662 + t680;
t629 = mrSges(5,1) * t685 - (qJD(3) * mrSges(5,3));
t630 = mrSges(5,1) * t684 + (qJD(3) * mrSges(5,2));
t631 = pkin(4) * t684 - (qJD(3) * pkin(8));
t648 = t652 ^ 2;
t708 = pkin(3) + pkin(8);
t564 = -t631 * t684 + t708 * t623 + (-pkin(4) * t648 + t709) * t657 + t662;
t620 = (pkin(3) * t652 - qJ(4) * t654) * qJD(1);
t656 = qJD(3) ^ 2;
t632 = t653 * g(1) - t655 * g(2);
t667 = -t657 * qJ(2) + qJDD(2) - t632;
t599 = t709 * qJDD(1) + t667;
t693 = t654 * t599;
t666 = -t656 * qJ(4) + t620 * t684 + qJDD(4) - t693;
t706 = pkin(8) * t657;
t568 = t624 * pkin(4) - t708 * qJDD(3) + (pkin(4) * t683 + t654 * t706 - g(3)) * t652 + t666;
t651 = sin(qJ(5));
t707 = cos(qJ(5));
t562 = t707 * t564 + t651 * t568;
t619 = t707 * qJD(3) + t651 * t685;
t583 = qJD(5) * t619 + qJDD(3) * t651 - t707 * t623;
t638 = qJD(5) + t684;
t593 = mrSges(6,1) * t638 - mrSges(6,3) * t619;
t617 = qJDD(5) + t624;
t618 = qJD(3) * t651 - t707 * t685;
t587 = pkin(5) * t618 - qJ(6) * t619;
t634 = t638 ^ 2;
t557 = -pkin(5) * t634 + qJ(6) * t617 + 0.2e1 * qJD(6) * t638 - t587 * t618 + t562;
t594 = -mrSges(7,1) * t638 + mrSges(7,2) * t619;
t679 = m(7) * t557 + t617 * mrSges(7,3) + t638 * t594;
t588 = mrSges(7,1) * t618 - mrSges(7,3) * t619;
t688 = -mrSges(6,1) * t618 - mrSges(6,2) * t619 - t588;
t702 = -mrSges(6,3) - mrSges(7,2);
t553 = m(6) * t562 - t617 * mrSges(6,2) + t702 * t583 - t638 * t593 + t688 * t618 + t679;
t561 = -t651 * t564 + t707 * t568;
t584 = -t618 * qJD(5) + t707 * qJDD(3) + t651 * t623;
t592 = -mrSges(6,2) * t638 - mrSges(6,3) * t618;
t558 = -t617 * pkin(5) - t634 * qJ(6) + t619 * t587 + qJDD(6) - t561;
t595 = -mrSges(7,2) * t618 + mrSges(7,3) * t638;
t672 = -m(7) * t558 + t617 * mrSges(7,1) + t638 * t595;
t554 = m(6) * t561 + t617 * mrSges(6,1) + t702 * t584 + t638 * t592 + t688 * t619 + t672;
t675 = t707 * t553 - t651 * t554;
t661 = m(5) * t570 - t624 * mrSges(5,3) - (t629 * t652 + t630 * t654) * qJD(1) + t675;
t703 = mrSges(4,1) - mrSges(5,2);
t711 = -m(4) * t598 - t624 * mrSges(4,2) - t703 * t623 - t627 * t685 - t628 * t684 - t661;
t591 = -t654 * g(3) + t652 * t599;
t571 = t656 * pkin(3) - qJDD(3) * qJ(4) + (qJD(3) * t718) + t620 * t685 - t591;
t705 = t652 * g(3);
t704 = mrSges(2,1) - mrSges(3,2);
t699 = Ifges(2,5) - Ifges(3,4);
t697 = (-Ifges(2,6) + Ifges(3,5));
t590 = t693 + t705;
t549 = t651 * t553 + t707 * t554;
t572 = -qJDD(3) * pkin(3) + t666 - t705;
t664 = -m(5) * t572 - t624 * mrSges(5,1) - t549;
t621 = (-mrSges(5,2) * t652 - mrSges(5,3) * t654) * qJD(1);
t673 = qJD(1) * (-t621 - (mrSges(4,1) * t652 + mrSges(4,2) * t654) * qJD(1));
t545 = m(4) * t590 - t624 * mrSges(4,3) + t703 * qJDD(3) + (t627 - t629) * qJD(3) + t654 * t673 + t664;
t567 = -t623 * pkin(4) + qJD(3) * t631 - t648 * t706 - t571;
t560 = -0.2e1 * qJD(6) * t619 + (t618 * t638 - t584) * qJ(6) + (t619 * t638 + t583) * pkin(5) + t567;
t555 = m(7) * t560 + t583 * mrSges(7,1) - t584 * mrSges(7,3) - t619 * t594 + t618 * t595;
t660 = m(6) * t567 + t583 * mrSges(6,1) + t584 * mrSges(6,2) + t618 * t592 + t619 * t593 + t555;
t659 = -m(5) * t571 + qJDD(3) * mrSges(5,3) + qJD(3) * t630 + t660;
t551 = t659 + (-mrSges(4,3) - mrSges(5,1)) * t623 + m(4) * t591 - qJDD(3) * mrSges(4,2) - qJD(3) * t628 + t652 * t673;
t541 = t654 * t545 + t652 * t551;
t601 = -qJDD(1) * pkin(1) + t667;
t665 = -m(3) * t601 + (t657 * mrSges(3,3)) - t541;
t539 = m(2) * t632 - (t657 * mrSges(2,2)) + t704 * qJDD(1) + t665;
t600 = t657 * pkin(1) - t668;
t658 = -m(3) * t600 + (t657 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) - t711;
t544 = m(2) * t633 - (t657 * mrSges(2,1)) - qJDD(1) * mrSges(2,2) + t658;
t692 = t655 * t539 + t653 * t544;
t691 = t714 * t618 - t701 * t619 - t694 * t638;
t690 = t694 * t618 - t700 * t619 - t712 * t638;
t689 = -t701 * t618 + t717 * t619 + t700 * t638;
t687 = -(t696 * qJD(3)) + (t716 * t652 + t695 * t654) * qJD(1);
t686 = (t698 * qJD(3)) + (t695 * t652 + t715 * t654) * qJD(1);
t677 = -t539 * t653 + t655 * t544;
t676 = -t652 * t545 + t654 * t551;
t674 = qJD(1) * ((t713 * qJD(3)) + (t696 * t652 - t698 * t654) * qJD(1));
t548 = mrSges(6,2) * t567 + mrSges(7,2) * t558 - mrSges(6,3) * t561 - mrSges(7,3) * t560 - qJ(6) * t555 - t701 * t583 + t717 * t584 + t700 * t617 + t690 * t618 + t691 * t638;
t547 = -mrSges(6,1) * t567 - mrSges(7,1) * t560 + mrSges(7,2) * t557 + mrSges(6,3) * t562 - pkin(5) * t555 - t714 * t583 + t701 * t584 + t694 * t617 + t690 * t619 + t689 * t638;
t546 = -t623 * mrSges(5,2) + t661;
t540 = -m(3) * g(3) + t676;
t537 = -mrSges(4,3) * t590 - qJ(4) * t546 + mrSges(7,3) * t557 - mrSges(7,1) * t558 - mrSges(5,3) * t570 + mrSges(5,1) * t572 + qJ(6) * t679 + pkin(5) * t672 + mrSges(6,1) * t561 - mrSges(6,2) * t562 + pkin(4) * t549 + mrSges(4,2) * t598 + t715 * t624 + t695 * t623 + (-pkin(5) * t588 - t691) * t619 + (-qJ(6) * t588 + t689) * t618 + t712 * t617 + (-mrSges(7,2) * pkin(5) + t700) * t584 + (-mrSges(7,2) * qJ(6) - t694) * t583 + t698 * qJDD(3) + t687 * qJD(3) + t652 * t674;
t536 = -mrSges(4,1) * t598 - mrSges(5,1) * t571 + mrSges(5,2) * t570 + mrSges(4,3) * t591 - pkin(3) * t546 + pkin(4) * t660 - pkin(8) * t675 + t686 * qJD(3) + t696 * qJDD(3) - t707 * t547 - t651 * t548 - t716 * t623 - t695 * t624 + t654 * t674;
t535 = mrSges(4,1) * t590 - mrSges(4,2) * t591 - t651 * t547 - mrSges(5,3) * t571 + mrSges(5,2) * t572 - qJ(2) * t540 - pkin(8) * t549 - mrSges(2,3) * t632 + mrSges(3,1) * t601 + pkin(2) * t541 + qJ(4) * t659 + pkin(3) * (-qJD(3) * t629 + t664) + t707 * t548 + (t697 * t657) + t698 * t624 + (-mrSges(5,1) * qJ(4) - t696) * t623 + (-mrSges(5,2) * pkin(3) - t713) * qJDD(3) + t699 * qJDD(1) + (mrSges(3,3) - mrSges(2,2)) * g(3) + ((-pkin(3) * t621 - t687) * t654 + (-qJ(4) * t621 + t686) * t652) * qJD(1);
t534 = -mrSges(3,1) * t600 + mrSges(2,3) * t633 - pkin(1) * t540 - pkin(2) * t711 - pkin(7) * t676 + t704 * g(3) - t697 * qJDD(1) - t654 * t536 - t652 * t537 + t699 * t657;
t1 = [-m(1) * g(1) + t677; -m(1) * g(2) + t692; (-m(1) - m(2) - m(3)) * g(3) + t676; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t692 - t653 * t534 + t655 * t535; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t677 + t655 * t534 + t653 * t535; pkin(1) * t665 + qJ(2) * t658 + t654 * t537 - t652 * t536 - pkin(7) * t541 + mrSges(2,1) * t632 - mrSges(2,2) * t633 + mrSges(3,2) * t601 - mrSges(3,3) * t600 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-mrSges(3,2) * pkin(1) + Ifges(3,1) + Ifges(2,3)) * qJDD(1);];
tauB  = t1;
