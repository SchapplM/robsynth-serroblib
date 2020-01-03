% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRRRP4
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
% Datum: 2019-12-31 21:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRRRP4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP4_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP4_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP4_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP4_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP4_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP4_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP4_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:50:55
% EndTime: 2019-12-31 21:50:58
% DurationCPUTime: 3.19s
% Computational Cost: add. (44430->250), mult. (56088->303), div. (0->0), fcn. (32755->8), ass. (0->105)
t713 = Ifges(5,1) + Ifges(6,1);
t702 = Ifges(5,4) - Ifges(6,5);
t711 = Ifges(6,4) + Ifges(5,5);
t712 = Ifges(5,2) + Ifges(6,3);
t709 = Ifges(5,6) - Ifges(6,6);
t710 = -Ifges(6,2) - Ifges(5,3);
t673 = sin(qJ(4));
t670 = qJD(1) + qJD(2);
t677 = cos(qJ(3));
t700 = t670 * t677;
t674 = sin(qJ(3));
t701 = t670 * t674;
t705 = cos(qJ(4));
t638 = t673 * t701 - t705 * t700;
t639 = (t673 * t677 + t705 * t674) * t670;
t669 = qJD(3) + qJD(4);
t708 = t712 * t638 - t702 * t639 - t709 * t669;
t707 = -t702 * t638 + t713 * t639 + t711 * t669;
t668 = qJDD(1) + qJDD(2);
t695 = qJD(3) * t670;
t645 = t674 * t668 + t677 * t695;
t676 = sin(qJ(1));
t679 = cos(qJ(1));
t657 = t676 * g(1) - t679 * g(2);
t650 = qJDD(1) * pkin(1) + t657;
t658 = -t679 * g(1) - t676 * g(2);
t680 = qJD(1) ^ 2;
t651 = -t680 * pkin(1) + t658;
t675 = sin(qJ(2));
t678 = cos(qJ(2));
t628 = t675 * t650 + t678 * t651;
t666 = t670 ^ 2;
t625 = -pkin(2) * t666 + pkin(7) * t668 + t628;
t699 = t674 * t625;
t704 = pkin(3) * t666;
t595 = qJDD(3) * pkin(3) - t645 * pkin(8) - t699 + (pkin(8) * t695 + t674 * t704 - g(3)) * t677;
t610 = -g(3) * t674 + t677 * t625;
t646 = t677 * t668 - t674 * t695;
t654 = qJD(3) * pkin(3) - pkin(8) * t701;
t672 = t677 ^ 2;
t596 = pkin(8) * t646 - qJD(3) * t654 - t672 * t704 + t610;
t592 = t673 * t595 + t705 * t596;
t607 = t639 * qJD(4) + t673 * t645 - t705 * t646;
t632 = t669 * mrSges(5,1) - t639 * mrSges(5,3);
t667 = qJDD(3) + qJDD(4);
t621 = pkin(4) * t638 - qJ(5) * t639;
t665 = t669 ^ 2;
t588 = -pkin(4) * t665 + qJ(5) * t667 + 0.2e1 * qJD(5) * t669 - t621 * t638 + t592;
t633 = -t669 * mrSges(6,1) + t639 * mrSges(6,2);
t694 = m(6) * t588 + t667 * mrSges(6,3) + t669 * t633;
t622 = mrSges(6,1) * t638 - mrSges(6,3) * t639;
t696 = -mrSges(5,1) * t638 - mrSges(5,2) * t639 - t622;
t703 = -mrSges(5,3) - mrSges(6,2);
t576 = m(5) * t592 - t667 * mrSges(5,2) + t703 * t607 - t669 * t632 + t696 * t638 + t694;
t591 = t705 * t595 - t673 * t596;
t608 = -t638 * qJD(4) + t705 * t645 + t673 * t646;
t631 = -mrSges(5,2) * t669 - mrSges(5,3) * t638;
t589 = -t667 * pkin(4) - t665 * qJ(5) + t639 * t621 + qJDD(5) - t591;
t634 = -t638 * mrSges(6,2) + t669 * mrSges(6,3);
t689 = -m(6) * t589 + t667 * mrSges(6,1) + t669 * t634;
t578 = m(5) * t591 + t667 * mrSges(5,1) + t703 * t608 + t669 * t631 + t696 * t639 + t689;
t570 = t673 * t576 + t705 * t578;
t609 = -t677 * g(3) - t699;
t636 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t674 + Ifges(4,2) * t677) * t670;
t637 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t674 + Ifges(4,4) * t677) * t670;
t583 = t608 * mrSges(6,2) + t639 * t622 - t689;
t683 = -mrSges(5,1) * t591 + mrSges(6,1) * t589 + mrSges(5,2) * t592 - mrSges(6,3) * t588 + pkin(4) * t583 - qJ(5) * t694 + t710 * t667 + t708 * t639 + (qJ(5) * t622 - t707) * t638 - t711 * t608 + (qJ(5) * mrSges(6,2) + t709) * t607;
t706 = mrSges(4,1) * t609 - mrSges(4,2) * t610 + Ifges(4,5) * t645 + Ifges(4,6) * t646 + Ifges(4,3) * qJDD(3) + pkin(3) * t570 + (t636 * t674 - t637 * t677) * t670 - t683;
t644 = (-mrSges(4,1) * t677 + mrSges(4,2) * t674) * t670;
t653 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t700;
t568 = m(4) * t609 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t645 + qJD(3) * t653 - t644 * t701 + t570;
t652 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t701;
t690 = t705 * t576 - t578 * t673;
t569 = m(4) * t610 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t646 - qJD(3) * t652 + t644 * t700 + t690;
t691 = -t568 * t674 + t677 * t569;
t560 = m(3) * t628 - mrSges(3,1) * t666 - mrSges(3,2) * t668 + t691;
t627 = t678 * t650 - t675 * t651;
t687 = -t668 * pkin(2) - t627;
t624 = -t666 * pkin(7) + t687;
t597 = -t646 * pkin(3) + t654 * t701 + (-pkin(8) * t672 - pkin(7)) * t666 + t687;
t585 = -0.2e1 * qJD(5) * t639 + (t638 * t669 - t608) * qJ(5) + (t639 * t669 + t607) * pkin(4) + t597;
t579 = m(6) * t585 + t607 * mrSges(6,1) - t608 * mrSges(6,3) - t639 * t633 + t638 * t634;
t684 = m(5) * t597 + t607 * mrSges(5,1) + mrSges(5,2) * t608 + t638 * t631 + t632 * t639 + t579;
t682 = -m(4) * t624 + t646 * mrSges(4,1) - mrSges(4,2) * t645 - t652 * t701 + t653 * t700 - t684;
t572 = m(3) * t627 + mrSges(3,1) * t668 - mrSges(3,2) * t666 + t682;
t557 = t675 * t560 + t678 * t572;
t554 = m(2) * t657 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t680 + t557;
t692 = t678 * t560 - t572 * t675;
t555 = m(2) * t658 - mrSges(2,1) * t680 - qJDD(1) * mrSges(2,2) + t692;
t698 = t679 * t554 + t676 * t555;
t562 = t677 * t568 + t674 * t569;
t697 = t709 * t638 - t711 * t639 + t710 * t669;
t693 = -t554 * t676 + t679 * t555;
t563 = -mrSges(5,1) * t597 - mrSges(6,1) * t585 + mrSges(6,2) * t588 + mrSges(5,3) * t592 - pkin(4) * t579 - t712 * t607 + t702 * t608 + t697 * t639 + t709 * t667 + t707 * t669;
t564 = mrSges(5,2) * t597 + mrSges(6,2) * t589 - mrSges(5,3) * t591 - mrSges(6,3) * t585 - qJ(5) * t579 - t702 * t607 + t713 * t608 + t697 * t638 + t711 * t667 + t708 * t669;
t635 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t674 + Ifges(4,6) * t677) * t670;
t548 = -mrSges(4,1) * t624 + mrSges(4,3) * t610 + Ifges(4,4) * t645 + Ifges(4,2) * t646 + Ifges(4,6) * qJDD(3) - pkin(3) * t684 + pkin(8) * t690 + qJD(3) * t637 + t705 * t563 + t673 * t564 - t635 * t701;
t550 = mrSges(4,2) * t624 - mrSges(4,3) * t609 + Ifges(4,1) * t645 + Ifges(4,4) * t646 + Ifges(4,5) * qJDD(3) - pkin(8) * t570 - qJD(3) * t636 - t673 * t563 + t705 * t564 + t635 * t700;
t686 = mrSges(3,1) * t627 - mrSges(3,2) * t628 + Ifges(3,3) * t668 + pkin(2) * t682 + pkin(7) * t691 + t677 * t548 + t674 * t550;
t685 = mrSges(2,1) * t657 - mrSges(2,2) * t658 + Ifges(2,3) * qJDD(1) + pkin(1) * t557 + t686;
t546 = mrSges(3,1) * g(3) + mrSges(3,3) * t628 + t666 * Ifges(3,5) + Ifges(3,6) * t668 - pkin(2) * t562 - t706;
t545 = -mrSges(3,2) * g(3) - mrSges(3,3) * t627 + Ifges(3,5) * t668 - Ifges(3,6) * t666 - pkin(7) * t562 - t548 * t674 + t550 * t677;
t544 = -mrSges(2,2) * g(3) - mrSges(2,3) * t657 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t680 - pkin(6) * t557 + t545 * t678 - t546 * t675;
t543 = Ifges(2,6) * qJDD(1) + t680 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t658 + t675 * t545 + t678 * t546 - pkin(1) * (-m(3) * g(3) + t562) + pkin(6) * t692;
t1 = [-m(1) * g(1) + t693; -m(1) * g(2) + t698; (-m(1) - m(2) - m(3)) * g(3) + t562; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t698 - t676 * t543 + t679 * t544; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t693 + t679 * t543 + t676 * t544; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t685; t685; t686; t706; -t683; t583;];
tauJB = t1;
