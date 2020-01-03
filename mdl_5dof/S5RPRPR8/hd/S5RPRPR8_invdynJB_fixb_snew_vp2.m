% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPRPR8
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-31 18:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPRPR8_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR8_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR8_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR8_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR8_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR8_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR8_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR8_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR8_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:21:39
% EndTime: 2019-12-31 18:21:44
% DurationCPUTime: 4.58s
% Computational Cost: add. (53711->268), mult. (109769->339), div. (0->0), fcn. (67801->10), ass. (0->109)
t629 = sin(qJ(1));
t632 = cos(qJ(1));
t613 = t629 * g(1) - t632 * g(2);
t604 = qJDD(1) * pkin(1) + t613;
t614 = -t632 * g(1) - t629 * g(2);
t634 = qJD(1) ^ 2;
t607 = -t634 * pkin(1) + t614;
t624 = sin(pkin(8));
t626 = cos(pkin(8));
t581 = t626 * t604 - t624 * t607;
t572 = -qJDD(1) * pkin(2) - t634 * pkin(6) - t581;
t628 = sin(qJ(3));
t631 = cos(qJ(3));
t647 = qJD(1) * qJD(3);
t646 = t631 * t647;
t608 = t628 * qJDD(1) + t646;
t617 = t628 * t647;
t609 = t631 * qJDD(1) - t617;
t561 = (-t608 - t646) * qJ(4) + (-t609 + t617) * pkin(3) + t572;
t582 = t624 * t604 + t626 * t607;
t573 = -t634 * pkin(2) + qJDD(1) * pkin(6) + t582;
t622 = -g(3) + qJDD(2);
t567 = t631 * t573 + t628 * t622;
t605 = (-pkin(3) * t631 - qJ(4) * t628) * qJD(1);
t633 = qJD(3) ^ 2;
t648 = t631 * qJD(1);
t565 = -t633 * pkin(3) + qJDD(3) * qJ(4) + t605 * t648 + t567;
t623 = sin(pkin(9));
t625 = cos(pkin(9));
t649 = qJD(1) * t628;
t601 = t623 * qJD(3) + t625 * t649;
t549 = -0.2e1 * qJD(4) * t601 + t625 * t561 - t623 * t565;
t587 = t623 * qJDD(3) + t625 * t608;
t600 = t625 * qJD(3) - t623 * t649;
t547 = (-t600 * t648 - t587) * pkin(7) + (t600 * t601 - t609) * pkin(4) + t549;
t550 = 0.2e1 * qJD(4) * t600 + t623 * t561 + t625 * t565;
t586 = t625 * qJDD(3) - t623 * t608;
t588 = -pkin(4) * t648 - t601 * pkin(7);
t599 = t600 ^ 2;
t548 = -t599 * pkin(4) + t586 * pkin(7) + t588 * t648 + t550;
t627 = sin(qJ(5));
t630 = cos(qJ(5));
t546 = t627 * t547 + t630 * t548;
t566 = -t628 * t573 + t631 * t622;
t563 = -qJDD(3) * pkin(3) - t633 * qJ(4) + t605 * t649 + qJDD(4) - t566;
t551 = -t586 * pkin(4) - t599 * pkin(7) + t601 * t588 + t563;
t579 = t627 * t600 + t630 * t601;
t553 = -t579 * qJD(5) + t630 * t586 - t627 * t587;
t578 = t630 * t600 - t627 * t601;
t554 = t578 * qJD(5) + t627 * t586 + t630 * t587;
t615 = qJD(5) - t648;
t557 = Ifges(6,5) * t579 + Ifges(6,6) * t578 + Ifges(6,3) * t615;
t559 = Ifges(6,1) * t579 + Ifges(6,4) * t578 + Ifges(6,5) * t615;
t603 = qJDD(5) - t609;
t535 = -mrSges(6,1) * t551 + mrSges(6,3) * t546 + Ifges(6,4) * t554 + Ifges(6,2) * t553 + Ifges(6,6) * t603 - t579 * t557 + t615 * t559;
t545 = t630 * t547 - t627 * t548;
t558 = Ifges(6,4) * t579 + Ifges(6,2) * t578 + Ifges(6,6) * t615;
t536 = mrSges(6,2) * t551 - mrSges(6,3) * t545 + Ifges(6,1) * t554 + Ifges(6,4) * t553 + Ifges(6,5) * t603 + t578 * t557 - t615 * t558;
t574 = Ifges(5,5) * t601 + Ifges(5,6) * t600 - Ifges(5,3) * t648;
t576 = Ifges(5,1) * t601 + Ifges(5,4) * t600 - Ifges(5,5) * t648;
t568 = -t615 * mrSges(6,2) + t578 * mrSges(6,3);
t569 = t615 * mrSges(6,1) - t579 * mrSges(6,3);
t639 = m(6) * t551 - t553 * mrSges(6,1) + t554 * mrSges(6,2) - t578 * t568 + t579 * t569;
t564 = -t578 * mrSges(6,1) + t579 * mrSges(6,2);
t542 = m(6) * t545 + t603 * mrSges(6,1) - t554 * mrSges(6,3) - t579 * t564 + t615 * t568;
t543 = m(6) * t546 - t603 * mrSges(6,2) + t553 * mrSges(6,3) + t578 * t564 - t615 * t569;
t642 = -t627 * t542 + t630 * t543;
t520 = -mrSges(5,1) * t563 + mrSges(5,3) * t550 + Ifges(5,4) * t587 + Ifges(5,2) * t586 - Ifges(5,6) * t609 - pkin(4) * t639 + pkin(7) * t642 + t630 * t535 + t627 * t536 - t601 * t574 - t576 * t648;
t534 = t630 * t542 + t627 * t543;
t575 = Ifges(5,4) * t601 + Ifges(5,2) * t600 - Ifges(5,6) * t648;
t522 = mrSges(5,2) * t563 - mrSges(5,3) * t549 + Ifges(5,1) * t587 + Ifges(5,4) * t586 - Ifges(5,5) * t609 - pkin(7) * t534 - t627 * t535 + t630 * t536 + t600 * t574 + t575 * t648;
t583 = -t600 * mrSges(5,1) + t601 * mrSges(5,2);
t640 = mrSges(5,2) * t648 + t600 * mrSges(5,3);
t532 = m(5) * t549 - t609 * mrSges(5,1) - t587 * mrSges(5,3) - t601 * t583 - t640 * t648 + t534;
t585 = -mrSges(5,1) * t648 - t601 * mrSges(5,3);
t533 = m(5) * t550 + t609 * mrSges(5,2) + t586 * mrSges(5,3) + t600 * t583 + t585 * t648 + t642;
t530 = -t623 * t532 + t625 * t533;
t544 = m(5) * t563 - t586 * mrSges(5,1) + t587 * mrSges(5,2) + t601 * t585 - t600 * t640 + t639;
t596 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t628 + Ifges(4,2) * t631) * qJD(1);
t597 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t628 + Ifges(4,4) * t631) * qJD(1);
t651 = mrSges(4,1) * t566 - mrSges(4,2) * t567 + Ifges(4,5) * t608 + Ifges(4,6) * t609 + Ifges(4,3) * qJDD(3) - pkin(3) * t544 + qJ(4) * t530 + t625 * t520 + t623 * t522 + (t628 * t596 - t631 * t597) * qJD(1);
t606 = (-mrSges(4,1) * t631 + mrSges(4,2) * t628) * qJD(1);
t611 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t649;
t528 = m(4) * t567 - qJDD(3) * mrSges(4,2) + t609 * mrSges(4,3) - qJD(3) * t611 + t606 * t648 + t530;
t612 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t648;
t538 = m(4) * t566 + qJDD(3) * mrSges(4,1) - t608 * mrSges(4,3) + qJD(3) * t612 - t606 * t649 - t544;
t643 = t631 * t528 - t628 * t538;
t517 = m(3) * t582 - t634 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t643;
t529 = t625 * t532 + t623 * t533;
t637 = -m(4) * t572 + t609 * mrSges(4,1) - t608 * mrSges(4,2) - t611 * t649 + t612 * t648 - t529;
t524 = m(3) * t581 + qJDD(1) * mrSges(3,1) - t634 * mrSges(3,2) + t637;
t512 = t624 * t517 + t626 * t524;
t509 = m(2) * t613 + qJDD(1) * mrSges(2,1) - t634 * mrSges(2,2) + t512;
t644 = t626 * t517 - t624 * t524;
t510 = m(2) * t614 - t634 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t644;
t650 = t632 * t509 + t629 * t510;
t521 = t628 * t528 + t631 * t538;
t518 = m(3) * t622 + t521;
t645 = -t629 * t509 + t632 * t510;
t595 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t628 + Ifges(4,6) * t631) * qJD(1);
t505 = mrSges(4,2) * t572 - mrSges(4,3) * t566 + Ifges(4,1) * t608 + Ifges(4,4) * t609 + Ifges(4,5) * qJDD(3) - qJ(4) * t529 - qJD(3) * t596 - t623 * t520 + t625 * t522 + t595 * t648;
t636 = mrSges(6,1) * t545 - mrSges(6,2) * t546 + Ifges(6,5) * t554 + Ifges(6,6) * t553 + Ifges(6,3) * t603 + t579 * t558 - t578 * t559;
t514 = Ifges(4,6) * qJDD(3) - t636 - t595 * t649 + (Ifges(4,2) + Ifges(5,3)) * t609 - pkin(3) * t529 - pkin(4) * t534 - mrSges(5,1) * t549 + mrSges(5,2) * t550 + mrSges(4,3) * t567 - mrSges(4,1) * t572 - Ifges(5,6) * t586 - Ifges(5,5) * t587 + qJD(3) * t597 + t600 * t576 - t601 * t575 + Ifges(4,4) * t608;
t638 = mrSges(2,1) * t613 + mrSges(3,1) * t581 - mrSges(2,2) * t614 - mrSges(3,2) * t582 + pkin(1) * t512 + pkin(2) * t637 + pkin(6) * t643 + t628 * t505 + t631 * t514 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t503 = -mrSges(3,1) * t622 + mrSges(3,3) * t582 + t634 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t521 - t651;
t502 = mrSges(3,2) * t622 - mrSges(3,3) * t581 + Ifges(3,5) * qJDD(1) - t634 * Ifges(3,6) - pkin(6) * t521 + t631 * t505 - t628 * t514;
t501 = -mrSges(2,2) * g(3) - mrSges(2,3) * t613 + Ifges(2,5) * qJDD(1) - t634 * Ifges(2,6) - qJ(2) * t512 + t626 * t502 - t624 * t503;
t500 = mrSges(2,1) * g(3) + mrSges(2,3) * t614 + t634 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t518 + qJ(2) * t644 + t624 * t502 + t626 * t503;
t1 = [-m(1) * g(1) + t645; -m(1) * g(2) + t650; (-m(1) - m(2)) * g(3) + t518; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t650 - t629 * t500 + t632 * t501; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t645 + t632 * t500 + t629 * t501; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t638; t638; t518; t651; t544; t636;];
tauJB = t1;
