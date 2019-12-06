% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPPRR2
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-05 17:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPPRR2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR2_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR2_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR2_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR2_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR2_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR2_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR2_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:39:51
% EndTime: 2019-12-05 17:39:55
% DurationCPUTime: 3.28s
% Computational Cost: add. (31571->240), mult. (70919->292), div. (0->0), fcn. (47562->8), ass. (0->108)
t669 = sin(qJ(1));
t672 = cos(qJ(1));
t643 = t669 * g(1) - t672 * g(2);
t673 = qJD(1) ^ 2;
t682 = -t673 * qJ(2) + qJDD(2) - t643;
t703 = -pkin(1) - qJ(3);
t711 = -(2 * qJD(1) * qJD(3)) + t703 * qJDD(1) + t682;
t644 = -t672 * g(1) - t669 * g(2);
t710 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t644;
t634 = t673 * pkin(1) - t710;
t709 = -m(3) * t634 + t673 * mrSges(3,2) + qJDD(1) * mrSges(3,3);
t707 = pkin(3) * t673;
t706 = mrSges(2,1) - mrSges(3,2);
t705 = -Ifges(3,4) + Ifges(2,5);
t704 = -Ifges(2,6) + Ifges(3,5);
t666 = cos(pkin(8));
t702 = mrSges(4,2) * t666;
t665 = sin(pkin(8));
t622 = t665 * g(3) + t711 * t666;
t606 = (-pkin(6) * qJDD(1) - t665 * t707) * t666 + t622;
t623 = -t666 * g(3) + t711 * t665;
t655 = t665 ^ 2;
t696 = qJDD(1) * t665;
t607 = -pkin(6) * t696 - t655 * t707 + t623;
t668 = sin(qJ(4));
t671 = cos(qJ(4));
t595 = t671 * t606 - t668 * t607;
t685 = -t665 * t668 + t666 * t671;
t686 = -t665 * t671 - t666 * t668;
t638 = t686 * qJD(1);
t698 = t638 * qJD(4);
t625 = t685 * qJDD(1) + t698;
t639 = t685 * qJD(1);
t584 = (-t625 + t698) * pkin(7) + (t638 * t639 + qJDD(4)) * pkin(4) + t595;
t596 = t668 * t606 + t671 * t607;
t624 = -t639 * qJD(4) + t686 * qJDD(1);
t632 = qJD(4) * pkin(4) - t639 * pkin(7);
t637 = t638 ^ 2;
t585 = -t637 * pkin(4) + t624 * pkin(7) - qJD(4) * t632 + t596;
t667 = sin(qJ(5));
t670 = cos(qJ(5));
t582 = t670 * t584 - t667 * t585;
t615 = t670 * t638 - t667 * t639;
t593 = t615 * qJD(5) + t667 * t624 + t670 * t625;
t616 = t667 * t638 + t670 * t639;
t601 = -t615 * mrSges(6,1) + t616 * mrSges(6,2);
t657 = qJD(4) + qJD(5);
t608 = -t657 * mrSges(6,2) + t615 * mrSges(6,3);
t654 = qJDD(4) + qJDD(5);
t579 = m(6) * t582 + t654 * mrSges(6,1) - t593 * mrSges(6,3) - t616 * t601 + t657 * t608;
t583 = t667 * t584 + t670 * t585;
t592 = -t616 * qJD(5) + t670 * t624 - t667 * t625;
t609 = t657 * mrSges(6,1) - t616 * mrSges(6,3);
t580 = m(6) * t583 - t654 * mrSges(6,2) + t592 * mrSges(6,3) + t615 * t601 - t657 * t609;
t568 = t670 * t579 + t667 * t580;
t618 = -t638 * mrSges(5,1) + t639 * mrSges(5,2);
t630 = -qJD(4) * mrSges(5,2) + t638 * mrSges(5,3);
t565 = m(5) * t595 + qJDD(4) * mrSges(5,1) - t625 * mrSges(5,3) + qJD(4) * t630 - t639 * t618 + t568;
t631 = qJD(4) * mrSges(5,1) - t639 * mrSges(5,3);
t690 = -t667 * t579 + t670 * t580;
t566 = m(5) * t596 - qJDD(4) * mrSges(5,2) + t624 * mrSges(5,3) - qJD(4) * t631 + t638 * t618 + t690;
t561 = t671 * t565 + t668 * t566;
t684 = -qJDD(1) * mrSges(4,3) - t673 * (mrSges(4,1) * t665 + t702);
t559 = m(4) * t622 + t684 * t666 + t561;
t691 = -t668 * t565 + t671 * t566;
t560 = m(4) * t623 + t684 * t665 + t691;
t555 = t666 * t559 + t665 * t560;
t636 = -qJDD(1) * pkin(1) + t682;
t679 = -m(3) * t636 + t673 * mrSges(3,3) - t555;
t551 = m(2) * t643 - t673 * mrSges(2,2) + t706 * qJDD(1) + t679;
t681 = qJDD(3) + t710;
t629 = t703 * t673 + t681;
t700 = -t666 ^ 2 - t655;
t611 = pkin(3) * t696 + (t700 * pkin(6) + t703) * t673 + t681;
t587 = -t624 * pkin(4) - t637 * pkin(7) + t639 * t632 + t611;
t680 = m(6) * t587 - t592 * mrSges(6,1) + t593 * mrSges(6,2) - t615 * t608 + t616 * t609;
t677 = m(5) * t611 - t624 * mrSges(5,1) + t625 * mrSges(5,2) - t638 * t630 + t639 * t631 + t680;
t675 = m(4) * t629 + mrSges(4,1) * t696 + qJDD(1) * t702 + t677;
t694 = t700 * mrSges(4,3);
t573 = (-mrSges(2,1) + t694) * t673 + m(2) * t644 + t675 - qJDD(1) * mrSges(2,2) + t709;
t701 = t672 * t551 + t669 * t573;
t687 = Ifges(4,5) * t666 - Ifges(4,6) * t665;
t699 = t673 * t687;
t693 = -t669 * t551 + t672 * t573;
t692 = -t665 * t559 + t666 * t560;
t689 = Ifges(4,1) * t666 - Ifges(4,4) * t665;
t688 = Ifges(4,4) * t666 - Ifges(4,2) * t665;
t598 = Ifges(6,4) * t616 + Ifges(6,2) * t615 + Ifges(6,6) * t657;
t599 = Ifges(6,1) * t616 + Ifges(6,4) * t615 + Ifges(6,5) * t657;
t678 = mrSges(6,1) * t582 - mrSges(6,2) * t583 + Ifges(6,5) * t593 + Ifges(6,6) * t592 + Ifges(6,3) * t654 + t616 * t598 - t615 * t599;
t597 = Ifges(6,5) * t616 + Ifges(6,6) * t615 + Ifges(6,3) * t657;
t569 = -mrSges(6,1) * t587 + mrSges(6,3) * t583 + Ifges(6,4) * t593 + Ifges(6,2) * t592 + Ifges(6,6) * t654 - t616 * t597 + t657 * t599;
t570 = mrSges(6,2) * t587 - mrSges(6,3) * t582 + Ifges(6,1) * t593 + Ifges(6,4) * t592 + Ifges(6,5) * t654 + t615 * t597 - t657 * t598;
t612 = Ifges(5,5) * t639 + Ifges(5,6) * t638 + Ifges(5,3) * qJD(4);
t614 = Ifges(5,1) * t639 + Ifges(5,4) * t638 + Ifges(5,5) * qJD(4);
t556 = -mrSges(5,1) * t611 + mrSges(5,3) * t596 + Ifges(5,4) * t625 + Ifges(5,2) * t624 + Ifges(5,6) * qJDD(4) - pkin(4) * t680 + pkin(7) * t690 + qJD(4) * t614 + t670 * t569 + t667 * t570 - t639 * t612;
t613 = Ifges(5,4) * t639 + Ifges(5,2) * t638 + Ifges(5,6) * qJD(4);
t557 = mrSges(5,2) * t611 - mrSges(5,3) * t595 + Ifges(5,1) * t625 + Ifges(5,4) * t624 + Ifges(5,5) * qJDD(4) - pkin(7) * t568 - qJD(4) * t613 - t667 * t569 + t670 * t570 + t638 * t612;
t547 = -mrSges(4,1) * t629 + mrSges(4,3) * t623 - pkin(3) * t677 + pkin(6) * t691 + t688 * qJDD(1) + t671 * t556 + t668 * t557 - t666 * t699;
t549 = mrSges(4,2) * t629 - mrSges(4,3) * t622 - pkin(6) * t561 + t689 * qJDD(1) - t668 * t556 + t671 * t557 - t665 * t699;
t553 = qJDD(1) * mrSges(3,2) - t679;
t575 = t673 * t694 + t675;
t676 = -mrSges(2,2) * t644 - mrSges(3,3) * t634 - pkin(1) * t553 - qJ(3) * t555 - t665 * t547 + t666 * t549 + qJ(2) * (t575 + t709) + mrSges(3,2) * t636 + mrSges(2,1) * t643 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1);
t674 = mrSges(5,1) * t595 - mrSges(5,2) * t596 + Ifges(5,5) * t625 + Ifges(5,6) * t624 + Ifges(5,3) * qJDD(4) + pkin(4) * t568 + t639 * t613 - t638 * t614 + t678;
t554 = -m(3) * g(3) + t692;
t546 = (t687 + t705) * qJDD(1) - mrSges(2,3) * t643 + mrSges(3,1) * t636 + mrSges(4,1) * t622 - mrSges(4,2) * t623 + pkin(3) * t561 + pkin(2) * t555 - qJ(2) * t554 + t674 + (-mrSges(2,2) + mrSges(3,3)) * g(3) + (t665 * t689 + t666 * t688 + t704) * t673;
t545 = -mrSges(3,1) * t634 + mrSges(2,3) * t644 - pkin(1) * t554 + pkin(2) * t575 + t706 * g(3) - qJ(3) * t692 - t704 * qJDD(1) - t666 * t547 - t665 * t549 + t705 * t673;
t1 = [-m(1) * g(1) + t693; -m(1) * g(2) + t701; (-m(1) - m(2) - m(3)) * g(3) + t692; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t701 - t669 * t545 + t672 * t546; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t693 + t672 * t545 + t669 * t546; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t676; t676; t553; t575; t674; t678;];
tauJB = t1;
