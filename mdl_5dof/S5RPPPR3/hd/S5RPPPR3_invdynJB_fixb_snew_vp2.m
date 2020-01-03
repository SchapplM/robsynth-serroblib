% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPPPR3
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3]';
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
% Datum: 2019-12-31 17:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPPPR3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR3_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR3_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR3_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR3_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR3_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR3_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR3_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:43:58
% EndTime: 2019-12-31 17:44:00
% DurationCPUTime: 1.67s
% Computational Cost: add. (14039->213), mult. (29464->260), div. (0->0), fcn. (16529->8), ass. (0->98)
t629 = sin(qJ(1));
t631 = cos(qJ(1));
t603 = t629 * g(1) - t631 * g(2);
t600 = qJDD(1) * pkin(1) + t603;
t604 = -t631 * g(1) - t629 * g(2);
t632 = qJD(1) ^ 2;
t601 = -t632 * pkin(1) + t604;
t625 = sin(pkin(7));
t627 = cos(pkin(7));
t583 = t625 * t600 + t627 * t601;
t670 = -t632 * pkin(2) + qJDD(1) * qJ(3) + (2 * qJD(1) * qJD(3)) + t583;
t624 = sin(pkin(8));
t626 = cos(pkin(8));
t669 = (Ifges(4,6) - Ifges(5,6)) * t626 + (Ifges(5,4) + Ifges(4,5)) * t624;
t618 = t624 ^ 2;
t619 = t626 ^ 2;
t658 = t619 * t632;
t666 = t618 * t632 + t658;
t663 = Ifges(4,4) - Ifges(5,5);
t665 = t663 * t624;
t623 = -g(3) + qJDD(2);
t562 = t626 * t623 - t670 * t624;
t664 = Ifges(4,1) + Ifges(5,1);
t662 = Ifges(4,2) + Ifges(5,3);
t661 = mrSges(4,2) * t624;
t660 = qJ(4) * t624;
t657 = t632 * qJ(3);
t596 = (-mrSges(5,1) * t626 - mrSges(5,3) * t624) * qJD(1);
t597 = (-mrSges(4,1) * t626 + t661) * qJD(1);
t642 = -pkin(3) * t626 - t660;
t595 = t642 * qJD(1);
t654 = t624 * qJD(1);
t558 = t595 * t654 + qJDD(4) - t562;
t554 = (-pkin(4) * t626 * t632 - pkin(6) * qJDD(1)) * t624 + t558;
t563 = t624 * t623 + t670 * t626;
t653 = t626 * qJD(1);
t560 = t595 * t653 + t563;
t650 = qJDD(1) * t626;
t555 = -pkin(4) * t658 - pkin(6) * t650 + t560;
t628 = sin(qJ(5));
t630 = cos(qJ(5));
t552 = t630 * t554 - t628 * t555;
t639 = -t624 * t628 - t626 * t630;
t589 = t639 * qJD(1);
t640 = t624 * t630 - t626 * t628;
t590 = t640 * qJD(1);
t574 = -t589 * mrSges(6,1) + t590 * mrSges(6,2);
t580 = t589 * qJD(5) + t640 * qJDD(1);
t584 = -qJD(5) * mrSges(6,2) + t589 * mrSges(6,3);
t549 = m(6) * t552 + qJDD(5) * mrSges(6,1) - t580 * mrSges(6,3) + qJD(5) * t584 - t590 * t574;
t553 = t628 * t554 + t630 * t555;
t579 = -t590 * qJD(5) + t639 * qJDD(1);
t585 = qJD(5) * mrSges(6,1) - t590 * mrSges(6,3);
t550 = m(6) * t553 - qJDD(5) * mrSges(6,2) + t579 * mrSges(6,3) - qJD(5) * t585 + t589 * t574;
t540 = t630 * t549 + t628 * t550;
t636 = m(5) * t558 + t540;
t537 = m(4) * t562 + ((-mrSges(5,2) - mrSges(4,3)) * qJDD(1) + (-t596 - t597) * qJD(1)) * t624 - t636;
t644 = -t628 * t549 + t630 * t550;
t638 = m(5) * t560 + mrSges(5,2) * t650 + t596 * t653 + t644;
t538 = m(4) * t563 + (qJDD(1) * mrSges(4,3) + qJD(1) * t597) * t626 + t638;
t645 = -t624 * t537 + t626 * t538;
t530 = m(3) * t583 - t632 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t645;
t582 = t627 * t600 - t625 * t601;
t648 = -qJDD(3) + t582;
t637 = -0.2e1 * qJD(4) * t654 - t648;
t561 = -t657 + (-pkin(2) + t642) * qJDD(1) + t637;
t557 = (qJ(3) + (-t618 - t619) * pkin(6)) * t632 + (t660 + pkin(2) + (pkin(3) + pkin(4)) * t626) * qJDD(1) - t637;
t641 = -m(6) * t557 + t579 * mrSges(6,1) - t580 * mrSges(6,2) + t589 * t584 - t590 * t585;
t651 = qJDD(1) * t624;
t547 = m(5) * t561 - mrSges(5,1) * t650 - t666 * mrSges(5,2) - mrSges(5,3) * t651 + t641;
t567 = -qJDD(1) * pkin(2) - t648 - t657;
t633 = -m(4) * t567 + mrSges(4,1) * t650 + t666 * mrSges(4,3) - t547;
t544 = (mrSges(3,1) - t661) * qJDD(1) - t632 * mrSges(3,2) + m(3) * t582 + t633;
t527 = t625 * t530 + t627 * t544;
t524 = m(2) * t603 + qJDD(1) * mrSges(2,1) - t632 * mrSges(2,2) + t527;
t646 = t627 * t530 - t625 * t544;
t525 = m(2) * t604 - t632 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t646;
t656 = t631 * t524 + t629 * t525;
t533 = t626 * t537 + t624 * t538;
t655 = t669 * qJD(1);
t531 = m(3) * t623 + t533;
t647 = -t629 * t524 + t631 * t525;
t569 = Ifges(6,4) * t590 + Ifges(6,2) * t589 + Ifges(6,6) * qJD(5);
t570 = Ifges(6,1) * t590 + Ifges(6,4) * t589 + Ifges(6,5) * qJD(5);
t635 = mrSges(6,1) * t552 - mrSges(6,2) * t553 + Ifges(6,5) * t580 + Ifges(6,6) * t579 + Ifges(6,3) * qJDD(5) + t590 * t569 - t589 * t570;
t568 = Ifges(6,5) * t590 + Ifges(6,6) * t589 + Ifges(6,3) * qJD(5);
t541 = -mrSges(6,1) * t557 + mrSges(6,3) * t553 + Ifges(6,4) * t580 + Ifges(6,2) * t579 + Ifges(6,6) * qJDD(5) + qJD(5) * t570 - t590 * t568;
t542 = mrSges(6,2) * t557 - mrSges(6,3) * t552 + Ifges(6,1) * t580 + Ifges(6,4) * t579 + Ifges(6,5) * qJDD(5) - qJD(5) * t569 + t589 * t568;
t518 = -mrSges(4,1) * t567 + mrSges(4,3) * t563 - mrSges(5,1) * t561 + mrSges(5,2) * t560 - t628 * t542 - t630 * t541 - pkin(4) * t641 - pkin(6) * t644 - pkin(3) * t547 - t655 * t654 + (t662 * t626 + t665) * qJDD(1);
t520 = mrSges(4,2) * t567 + mrSges(5,2) * t558 - mrSges(4,3) * t562 - mrSges(5,3) * t561 - pkin(6) * t540 - qJ(4) * t547 - t628 * t541 + t630 * t542 + t655 * t653 + (t664 * t624 + t663 * t626) * qJDD(1);
t546 = mrSges(4,2) * t651 - t633;
t634 = mrSges(2,1) * t603 + mrSges(3,1) * t582 - mrSges(2,2) * t604 - mrSges(3,2) * t583 + pkin(1) * t527 - pkin(2) * t546 + qJ(3) * t645 + t626 * t518 + t624 * t520 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t539 = (qJDD(1) * mrSges(5,2) + qJD(1) * t596) * t624 + t636;
t516 = t635 + (Ifges(3,6) - t669) * qJDD(1) - qJ(4) * t638 - mrSges(3,1) * t623 + mrSges(3,3) * t583 - mrSges(5,3) * t560 - mrSges(4,1) * t562 + mrSges(4,2) * t563 + mrSges(5,1) * t558 + pkin(4) * t540 + pkin(3) * t539 - pkin(2) * t533 + (t663 * t619 + (-t665 + (-t662 + t664) * t626) * t624 + Ifges(3,5)) * t632;
t515 = mrSges(3,2) * t623 - mrSges(3,3) * t582 + Ifges(3,5) * qJDD(1) - t632 * Ifges(3,6) - qJ(3) * t533 - t624 * t518 + t626 * t520;
t514 = -mrSges(2,2) * g(3) - mrSges(2,3) * t603 + Ifges(2,5) * qJDD(1) - t632 * Ifges(2,6) - qJ(2) * t527 + t627 * t515 - t625 * t516;
t513 = mrSges(2,1) * g(3) + mrSges(2,3) * t604 + t632 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t531 + qJ(2) * t646 + t625 * t515 + t627 * t516;
t1 = [-m(1) * g(1) + t647; -m(1) * g(2) + t656; (-m(1) - m(2)) * g(3) + t531; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t656 - t629 * t513 + t631 * t514; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t647 + t631 * t513 + t629 * t514; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t634; t634; t531; t546; t539; t635;];
tauJB = t1;
