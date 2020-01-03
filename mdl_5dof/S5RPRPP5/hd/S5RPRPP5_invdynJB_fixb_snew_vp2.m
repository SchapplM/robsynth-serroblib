% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPRPP5
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3]';
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
% Datum: 2019-12-31 18:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPRPP5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP5_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP5_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP5_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RPRPP5_invdynJB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP5_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP5_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPP5_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:16:15
% EndTime: 2019-12-31 18:16:17
% DurationCPUTime: 1.44s
% Computational Cost: add. (5082->223), mult. (9759->254), div. (0->0), fcn. (3699->4), ass. (0->90)
t684 = Ifges(4,4) - Ifges(5,5) - Ifges(6,4);
t714 = -Ifges(4,2) - Ifges(5,3) - Ifges(6,2);
t706 = Ifges(4,6) - Ifges(5,6) + Ifges(6,6);
t713 = Ifges(4,1) + Ifges(5,1) + Ifges(6,1);
t707 = Ifges(4,5) + Ifges(5,4) - Ifges(6,5);
t712 = -2 * qJD(1);
t661 = sin(qJ(3));
t663 = cos(qJ(3));
t708 = t706 * qJD(3) + (t714 * t661 + t684 * t663) * qJD(1);
t705 = Ifges(4,3) + Ifges(5,2) + Ifges(6,3);
t704 = t707 * qJD(3) + (-t684 * t661 + t713 * t663) * qJD(1);
t662 = sin(qJ(1));
t664 = cos(qJ(1));
t638 = -t664 * g(1) - t662 * g(2);
t666 = qJD(1) ^ 2;
t591 = (t666 * pkin(1)) - qJDD(1) * qJ(2) + (qJD(2) * t712) - t638;
t585 = -(t666 * pkin(6)) - t591;
t686 = qJD(1) * qJD(3);
t623 = t661 * qJDD(1) + t663 * t686;
t624 = t663 * qJDD(1) - t661 * t686;
t688 = qJD(1) * t661;
t631 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t688;
t687 = qJD(1) * t663;
t634 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t687;
t703 = -m(4) * t585 - t623 * mrSges(4,1) - t624 * mrSges(4,2) - t631 * t688 - t634 * t687;
t702 = 2 * qJD(4);
t701 = -pkin(3) - pkin(4);
t700 = pkin(4) * t666;
t699 = mrSges(2,1) - mrSges(3,2);
t698 = -mrSges(4,3) - mrSges(5,2);
t697 = Ifges(2,5) - Ifges(3,4);
t696 = (-Ifges(2,6) + Ifges(3,5));
t694 = t624 * mrSges(6,2);
t693 = qJ(4) * t661;
t692 = t624 * qJ(4);
t637 = t662 * g(1) - t664 * g(2);
t675 = -t666 * qJ(2) + qJDD(2) - t637;
t586 = (-pkin(1) - pkin(6)) * qJDD(1) + t675;
t582 = -t663 * g(3) + t661 * t586;
t619 = (pkin(3) * t661 - qJ(4) * t663) * qJD(1);
t665 = qJD(3) ^ 2;
t673 = -t665 * pkin(3) + qJDD(3) * qJ(4) + qJD(3) * t702 + t582;
t578 = -t619 * t688 + t673;
t635 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t687;
t632 = -qJD(3) * pkin(4) - qJ(5) * t687;
t657 = t661 ^ 2;
t572 = -t657 * t700 + t623 * qJ(5) + qJD(3) * t632 + ((2 * qJD(5)) - t619) * t688 + t673;
t621 = (-mrSges(6,1) * t661 + mrSges(6,2) * t663) * qJD(1);
t633 = -qJD(3) * mrSges(6,1) - mrSges(6,3) * t687;
t676 = m(6) * t572 + qJDD(3) * mrSges(6,2) + t623 * mrSges(6,3) + qJD(3) * t633 + t621 * t688;
t672 = m(5) * t578 + qJDD(3) * mrSges(5,3) + qJD(3) * t635 + t676;
t620 = (mrSges(5,1) * t661 - mrSges(5,3) * t663) * qJD(1);
t681 = qJD(1) * (-t620 - (mrSges(4,1) * t661 + mrSges(4,2) * t663) * qJD(1));
t561 = m(4) * t582 - qJDD(3) * mrSges(4,2) - qJD(3) * t634 + t698 * t623 + t661 * t681 + t672;
t581 = t661 * g(3) + t663 * t586;
t677 = -t665 * qJ(4) + t619 * t687 + qJDD(4);
t573 = -t624 * qJ(5) + ((qJD(5) * t712) - t586) * t663 + t701 * qJDD(3) + (-qJ(5) * t686 + t663 * t700 - g(3)) * t661 + t677;
t630 = qJD(3) * mrSges(6,2) + mrSges(6,3) * t688;
t568 = m(6) * t573 - qJDD(3) * mrSges(6,1) - t624 * mrSges(6,3) - qJD(3) * t630 - t621 * t687;
t579 = -qJDD(3) * pkin(3) - t581 + t677;
t636 = -mrSges(5,2) * t688 + qJD(3) * mrSges(5,3);
t670 = -m(5) * t579 + qJDD(3) * mrSges(5,1) + qJD(3) * t636 - t568;
t562 = m(4) * t581 + qJDD(3) * mrSges(4,1) + qJD(3) * t631 + t698 * t624 + t663 * t681 + t670;
t555 = t661 * t561 + t663 * t562;
t593 = -qJDD(1) * pkin(1) + t675;
t671 = -m(3) * t593 + (t666 * mrSges(3,3)) - t555;
t551 = m(2) * t637 - (t666 * mrSges(2,2)) + t699 * qJDD(1) + t671;
t575 = t623 * pkin(3) - t692 + (-0.2e1 * qJD(4) * t663 + (pkin(3) * t663 + t693) * qJD(3)) * qJD(1) + t585;
t569 = t692 + qJDD(5) + (-qJ(5) * t657 + pkin(6)) * t666 + t701 * t623 + (-qJD(3) * t693 + (-pkin(3) * qJD(3) + t632 + t702) * t663) * qJD(1) + t591;
t679 = -m(6) * t569 + t623 * mrSges(6,1) + t630 * t688;
t674 = m(5) * t575 + t623 * mrSges(5,1) + t636 * t688 + t679;
t689 = t633 + t635;
t563 = (-mrSges(6,2) - mrSges(5,3)) * t624 - t689 * t687 + t674;
t668 = -m(3) * t591 + (t666 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) + t563 - t703;
t558 = m(2) * t638 - (t666 * mrSges(2,1)) - qJDD(1) * mrSges(2,2) + t668;
t691 = t664 * t551 + t662 * t558;
t683 = -t662 * t551 + t664 * t558;
t682 = t663 * t561 - t661 * t562;
t678 = qJD(1) * (-t705 * qJD(3) + (t706 * t661 - t707 * t663) * qJD(1));
t567 = t633 * t687 - t679 + t694;
t547 = -mrSges(4,1) * t585 - mrSges(5,1) * t575 + mrSges(6,1) * t569 + mrSges(5,2) * t578 + mrSges(4,3) * t582 - mrSges(6,3) * t572 - pkin(3) * t563 + pkin(4) * t567 - qJ(5) * t676 + t704 * qJD(3) + t706 * qJDD(3) + t714 * t623 + t684 * t624 + t663 * t678;
t549 = mrSges(4,2) * t585 + mrSges(5,2) * t579 + mrSges(6,2) * t569 - mrSges(4,3) * t581 - mrSges(5,3) * t575 - mrSges(6,3) * t573 - qJ(4) * t563 - qJ(5) * t568 - t708 * qJD(3) + t707 * qJDD(3) - t684 * t623 + t713 * t624 + t661 * t678;
t553 = qJDD(1) * mrSges(3,2) - t671;
t669 = mrSges(2,1) * t637 - mrSges(2,2) * t638 + mrSges(3,2) * t593 - mrSges(3,3) * t591 - pkin(1) * t553 - pkin(6) * t555 + qJ(2) * t668 - t661 * t547 + t663 * t549 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1);
t566 = t624 * mrSges(5,2) + t620 * t687 - t670;
t667 = mrSges(4,1) * t581 - mrSges(5,1) * t579 - mrSges(6,1) * t573 - mrSges(4,2) * t582 + mrSges(6,2) * t572 + mrSges(5,3) * t578 - pkin(3) * t566 - pkin(4) * t568 + qJ(4) * t672 + t708 * t687 + (-qJ(4) * t620 + t704) * t688 + t707 * t624 + (-qJ(4) * mrSges(5,2) - t706) * t623 + t705 * qJDD(3);
t554 = -m(3) * g(3) + t682;
t546 = -qJ(2) * t554 + pkin(2) * t555 + t697 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3) + t667 + (t696 * t666) + mrSges(3,1) * t593 - mrSges(2,3) * t637;
t545 = mrSges(2,3) * t638 - mrSges(3,1) * t591 - t661 * t549 - pkin(2) * (t624 * mrSges(5,3) - t674 + t694 + t703) - pkin(6) * t682 - pkin(1) * t554 + t697 * t666 + (-pkin(2) * t689 * qJD(1) - t547) * t663 - t696 * qJDD(1) + t699 * g(3);
t1 = [-m(1) * g(1) + t683; -m(1) * g(2) + t691; (-m(1) - m(2) - m(3)) * g(3) + t682; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t691 - t662 * t545 + t664 * t546; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t683 + t664 * t545 + t662 * t546; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t669; t669; t553; t667; t566; t567;];
tauJB = t1;
