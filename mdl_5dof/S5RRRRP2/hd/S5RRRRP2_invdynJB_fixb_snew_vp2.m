% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRRRP2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-05 18:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRRRP2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP2_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP2_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP2_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP2_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP2_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:47:50
% EndTime: 2019-12-05 18:47:54
% DurationCPUTime: 3.29s
% Computational Cost: add. (46882->251), mult. (59433->304), div. (0->0), fcn. (35123->8), ass. (0->104)
t713 = Ifges(5,4) + Ifges(6,4);
t723 = Ifges(5,2) + Ifges(6,2);
t719 = Ifges(5,6) + Ifges(6,6);
t679 = qJD(1) + qJD(2);
t682 = sin(qJ(4));
t683 = sin(qJ(3));
t686 = cos(qJ(4));
t687 = cos(qJ(3));
t649 = (-t682 * t683 + t686 * t687) * t679;
t677 = qJDD(1) + qJDD(2);
t707 = qJD(3) * t679;
t655 = t683 * t677 + t687 * t707;
t656 = t687 * t677 - t683 * t707;
t619 = t649 * qJD(4) + t686 * t655 + t682 * t656;
t650 = (t682 * t687 + t683 * t686) * t679;
t632 = -t649 * mrSges(6,1) + t650 * mrSges(6,2);
t685 = sin(qJ(1));
t689 = cos(qJ(1));
t667 = t689 * g(2) + t685 * g(3);
t660 = qJDD(1) * pkin(1) + t667;
t666 = t685 * g(2) - t689 * g(3);
t690 = qJD(1) ^ 2;
t661 = -t690 * pkin(1) + t666;
t684 = sin(qJ(2));
t688 = cos(qJ(2));
t638 = t684 * t660 + t688 * t661;
t675 = t679 ^ 2;
t635 = -t675 * pkin(2) + t677 * pkin(7) + t638;
t710 = t683 * t635;
t714 = pkin(3) * t675;
t604 = qJDD(3) * pkin(3) - t655 * pkin(8) - t710 + (pkin(8) * t707 + t683 * t714 - g(1)) * t687;
t621 = -t683 * g(1) + t687 * t635;
t712 = t679 * t683;
t664 = qJD(3) * pkin(3) - pkin(8) * t712;
t681 = t687 ^ 2;
t605 = t656 * pkin(8) - qJD(3) * t664 - t681 * t714 + t621;
t599 = t686 * t604 - t682 * t605;
t676 = qJDD(3) + qJDD(4);
t678 = qJD(3) + qJD(4);
t593 = -0.2e1 * qJD(5) * t650 + (t649 * t678 - t619) * qJ(5) + (t649 * t650 + t676) * pkin(4) + t599;
t640 = -t678 * mrSges(6,2) + t649 * mrSges(6,3);
t706 = m(6) * t593 + t676 * mrSges(6,1) + t678 * t640;
t589 = -t619 * mrSges(6,3) - t650 * t632 + t706;
t600 = t682 * t604 + t686 * t605;
t618 = -t650 * qJD(4) - t682 * t655 + t686 * t656;
t642 = t678 * pkin(4) - t650 * qJ(5);
t645 = t649 ^ 2;
t595 = -t645 * pkin(4) + t618 * qJ(5) + 0.2e1 * qJD(5) * t649 - t678 * t642 + t600;
t720 = Ifges(5,5) + Ifges(6,5);
t721 = Ifges(5,1) + Ifges(6,1);
t708 = -t713 * t649 - t721 * t650 - t720 * t678;
t717 = t723 * t649 + t713 * t650 + t719 * t678;
t718 = Ifges(5,3) + Ifges(6,3);
t722 = mrSges(5,1) * t599 + mrSges(6,1) * t593 - mrSges(5,2) * t600 - mrSges(6,2) * t595 + pkin(4) * t589 + t719 * t618 + t720 * t619 + t708 * t649 + t717 * t650 + t718 * t676;
t633 = -t649 * mrSges(5,1) + t650 * mrSges(5,2);
t641 = -t678 * mrSges(5,2) + t649 * mrSges(5,3);
t584 = m(5) * t599 + t676 * mrSges(5,1) + t678 * t641 + (-t632 - t633) * t650 + (-mrSges(5,3) - mrSges(6,3)) * t619 + t706;
t643 = t678 * mrSges(6,1) - t650 * mrSges(6,3);
t644 = t678 * mrSges(5,1) - t650 * mrSges(5,3);
t705 = m(6) * t595 + t618 * mrSges(6,3) + t649 * t632;
t587 = m(5) * t600 + t618 * mrSges(5,3) + t649 * t633 + (-t643 - t644) * t678 + (-mrSges(5,2) - mrSges(6,2)) * t676 + t705;
t579 = t686 * t584 + t682 * t587;
t620 = -t687 * g(1) - t710;
t647 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t683 + Ifges(4,2) * t687) * t679;
t648 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t683 + Ifges(4,4) * t687) * t679;
t716 = mrSges(4,1) * t620 - mrSges(4,2) * t621 + Ifges(4,5) * t655 + Ifges(4,6) * t656 + Ifges(4,3) * qJDD(3) + pkin(3) * t579 + (t683 * t647 - t687 * t648) * t679 + t722;
t711 = t679 * t687;
t654 = (-mrSges(4,1) * t687 + mrSges(4,2) * t683) * t679;
t663 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t711;
t577 = m(4) * t620 + qJDD(3) * mrSges(4,1) - t655 * mrSges(4,3) + qJD(3) * t663 - t654 * t712 + t579;
t662 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t712;
t700 = -t682 * t584 + t686 * t587;
t578 = m(4) * t621 - qJDD(3) * mrSges(4,2) + t656 * mrSges(4,3) - qJD(3) * t662 + t654 * t711 + t700;
t701 = -t683 * t577 + t687 * t578;
t569 = m(3) * t638 - t675 * mrSges(3,1) - t677 * mrSges(3,2) + t701;
t637 = t688 * t660 - t684 * t661;
t697 = -t677 * pkin(2) - t637;
t634 = -t675 * pkin(7) + t697;
t606 = -t656 * pkin(3) + t664 * t712 + (-pkin(8) * t681 - pkin(7)) * t675 + t697;
t597 = -t618 * pkin(4) - t645 * qJ(5) + t650 * t642 + qJDD(5) + t606;
t590 = m(6) * t597 - t618 * mrSges(6,1) + t619 * mrSges(6,2) - t649 * t640 + t650 * t643;
t694 = m(5) * t606 - t618 * mrSges(5,1) + t619 * mrSges(5,2) - t649 * t641 + t650 * t644 + t590;
t692 = -m(4) * t634 + t656 * mrSges(4,1) - t655 * mrSges(4,2) - t662 * t712 + t663 * t711 - t694;
t581 = m(3) * t637 + t677 * mrSges(3,1) - t675 * mrSges(3,2) + t692;
t566 = t684 * t569 + t688 * t581;
t571 = t687 * t577 + t683 * t578;
t709 = -t719 * t649 - t720 * t650 - t718 * t678;
t702 = t688 * t569 - t684 * t581;
t563 = m(2) * t666 - t690 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t702;
t564 = m(2) * t667 + qJDD(1) * mrSges(2,1) - t690 * mrSges(2,2) + t566;
t703 = t689 * t563 - t685 * t564;
t699 = -t685 * t563 - t689 * t564;
t572 = -mrSges(5,1) * t606 + mrSges(5,3) * t600 - mrSges(6,1) * t597 + mrSges(6,3) * t595 - pkin(4) * t590 + qJ(5) * t705 + (-qJ(5) * t643 - t708) * t678 + (-qJ(5) * mrSges(6,2) + t719) * t676 + t709 * t650 + t713 * t619 + t723 * t618;
t573 = mrSges(5,2) * t606 + mrSges(6,2) * t597 - mrSges(5,3) * t599 - mrSges(6,3) * t593 - qJ(5) * t589 + t713 * t618 + t721 * t619 - t709 * t649 + t720 * t676 - t717 * t678;
t646 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t683 + Ifges(4,6) * t687) * t679;
t559 = -mrSges(4,1) * t634 + mrSges(4,3) * t621 + Ifges(4,4) * t655 + Ifges(4,2) * t656 + Ifges(4,6) * qJDD(3) - pkin(3) * t694 + pkin(8) * t700 + qJD(3) * t648 + t686 * t572 + t682 * t573 - t646 * t712;
t561 = mrSges(4,2) * t634 - mrSges(4,3) * t620 + Ifges(4,1) * t655 + Ifges(4,4) * t656 + Ifges(4,5) * qJDD(3) - pkin(8) * t579 - qJD(3) * t647 - t682 * t572 + t686 * t573 + t646 * t711;
t696 = mrSges(3,1) * t637 - mrSges(3,2) * t638 + Ifges(3,3) * t677 + pkin(2) * t692 + pkin(7) * t701 + t687 * t559 + t683 * t561;
t695 = mrSges(2,1) * t667 - mrSges(2,2) * t666 + Ifges(2,3) * qJDD(1) + pkin(1) * t566 + t696;
t557 = mrSges(3,1) * g(1) + mrSges(3,3) * t638 + t675 * Ifges(3,5) + Ifges(3,6) * t677 - pkin(2) * t571 - t716;
t556 = -mrSges(3,2) * g(1) - mrSges(3,3) * t637 + Ifges(3,5) * t677 - t675 * Ifges(3,6) - pkin(7) * t571 - t683 * t559 + t687 * t561;
t555 = -mrSges(2,2) * g(1) - mrSges(2,3) * t667 + Ifges(2,5) * qJDD(1) - t690 * Ifges(2,6) - pkin(6) * t566 + t688 * t556 - t684 * t557;
t554 = Ifges(2,6) * qJDD(1) + t690 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t666 + t684 * t556 + t688 * t557 - pkin(1) * (-m(3) * g(1) + t571) + pkin(6) * t702;
t1 = [(-m(1) - m(2) - m(3)) * g(1) + t571; -m(1) * g(2) + t699; -m(1) * g(3) + t703; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t695; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t703 - t689 * t554 - t685 * t555; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t699 - t685 * t554 + t689 * t555; t695; t696; t716; t722; t590;];
tauJB = t1;
