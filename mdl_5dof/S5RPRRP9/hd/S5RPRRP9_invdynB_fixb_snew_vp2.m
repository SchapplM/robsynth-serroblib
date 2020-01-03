% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPRRP9
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPRRP9_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP9_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP9_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP9_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP9_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP9_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP9_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP9_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP9_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:48:54
% EndTime: 2019-12-31 18:48:58
% DurationCPUTime: 4.16s
% Computational Cost: add. (41076->267), mult. (99452->324), div. (0->0), fcn. (71380->8), ass. (0->112)
t612 = Ifges(5,1) + Ifges(6,1);
t606 = Ifges(5,4) - Ifges(6,5);
t605 = Ifges(6,4) + Ifges(5,5);
t611 = Ifges(5,2) + Ifges(6,3);
t610 = -Ifges(6,2) - Ifges(5,3);
t604 = Ifges(5,6) - Ifges(6,6);
t575 = qJD(1) ^ 2;
t609 = cos(qJ(4));
t569 = cos(pkin(8));
t608 = pkin(2) * t569;
t607 = -mrSges(5,3) - mrSges(6,2);
t568 = sin(pkin(8));
t603 = mrSges(3,2) * t568;
t566 = t569 ^ 2;
t602 = t566 * t575;
t572 = sin(qJ(1));
t574 = cos(qJ(1));
t554 = -t574 * g(1) - t572 * g(2);
t550 = -t575 * pkin(1) + qJDD(1) * qJ(2) + t554;
t594 = qJD(1) * qJD(2);
t591 = -t569 * g(3) - 0.2e1 * t568 * t594;
t524 = (-pkin(6) * qJDD(1) + t575 * t608 - t550) * t568 + t591;
t540 = -t568 * g(3) + (t550 + 0.2e1 * t594) * t569;
t593 = qJDD(1) * t569;
t525 = -pkin(2) * t602 + pkin(6) * t593 + t540;
t571 = sin(qJ(3));
t573 = cos(qJ(3));
t500 = t573 * t524 - t571 * t525;
t581 = t568 * t573 + t569 * t571;
t580 = -t568 * t571 + t569 * t573;
t548 = t580 * qJD(1);
t595 = t548 * qJD(3);
t538 = t581 * qJDD(1) + t595;
t549 = t581 * qJD(1);
t489 = (-t538 + t595) * pkin(7) + (t548 * t549 + qJDD(3)) * pkin(3) + t500;
t501 = t571 * t524 + t573 * t525;
t537 = -t549 * qJD(3) + t580 * qJDD(1);
t543 = qJD(3) * pkin(3) - t549 * pkin(7);
t547 = t548 ^ 2;
t491 = -t547 * pkin(3) + t537 * pkin(7) - qJD(3) * t543 + t501;
t570 = sin(qJ(4));
t487 = t570 * t489 + t609 * t491;
t530 = t570 * t548 + t609 * t549;
t498 = t530 * qJD(4) - t609 * t537 + t570 * t538;
t567 = qJD(3) + qJD(4);
t520 = t567 * mrSges(5,1) - t530 * mrSges(5,3);
t529 = -t609 * t548 + t570 * t549;
t564 = qJDD(3) + qJDD(4);
t512 = t529 * pkin(4) - t530 * qJ(5);
t563 = t567 ^ 2;
t482 = -t563 * pkin(4) + t564 * qJ(5) + 0.2e1 * qJD(5) * t567 - t529 * t512 + t487;
t521 = -t567 * mrSges(6,1) + t530 * mrSges(6,2);
t592 = m(6) * t482 + t564 * mrSges(6,3) + t567 * t521;
t513 = t529 * mrSges(6,1) - t530 * mrSges(6,3);
t597 = -t529 * mrSges(5,1) - t530 * mrSges(5,2) - t513;
t477 = m(5) * t487 - t564 * mrSges(5,2) + t607 * t498 - t567 * t520 + t597 * t529 + t592;
t486 = t609 * t489 - t570 * t491;
t499 = -t529 * qJD(4) + t570 * t537 + t609 * t538;
t519 = -t567 * mrSges(5,2) - t529 * mrSges(5,3);
t483 = -t564 * pkin(4) - t563 * qJ(5) + t530 * t512 + qJDD(5) - t486;
t522 = -t529 * mrSges(6,2) + t567 * mrSges(6,3);
t586 = -m(6) * t483 + t564 * mrSges(6,1) + t567 * t522;
t479 = m(5) * t486 + t564 * mrSges(5,1) + t607 * t499 + t567 * t519 + t597 * t530 + t586;
t472 = t570 * t477 + t609 * t479;
t534 = -t548 * mrSges(4,1) + t549 * mrSges(4,2);
t541 = -qJD(3) * mrSges(4,2) + t548 * mrSges(4,3);
t470 = m(4) * t500 + qJDD(3) * mrSges(4,1) - t538 * mrSges(4,3) + qJD(3) * t541 - t549 * t534 + t472;
t542 = qJD(3) * mrSges(4,1) - t549 * mrSges(4,3);
t587 = t609 * t477 - t570 * t479;
t471 = m(4) * t501 - qJDD(3) * mrSges(4,2) + t537 * mrSges(4,3) - qJD(3) * t542 + t548 * t534 + t587;
t464 = t573 * t470 + t571 * t471;
t539 = -t568 * t550 + t591;
t579 = mrSges(3,3) * qJDD(1) + t575 * (-mrSges(3,1) * t569 + t603);
t462 = m(3) * t539 - t579 * t568 + t464;
t588 = -t571 * t470 + t573 * t471;
t463 = m(3) * t540 + t579 * t569 + t588;
t589 = -t568 * t462 + t569 * t463;
t455 = m(2) * t554 - t575 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t589;
t553 = t572 * g(1) - t574 * g(2);
t585 = qJDD(2) - t553;
t546 = -qJDD(1) * pkin(1) - t575 * qJ(2) + t585;
t565 = t568 ^ 2;
t536 = (-pkin(1) - t608) * qJDD(1) + (-qJ(2) + (-t565 - t566) * pkin(6)) * t575 + t585;
t493 = -t537 * pkin(3) - t547 * pkin(7) + t549 * t543 + t536;
t485 = -0.2e1 * qJD(5) * t530 + (t529 * t567 - t499) * qJ(5) + (t530 * t567 + t498) * pkin(4) + t493;
t480 = m(6) * t485 + t498 * mrSges(6,1) - t499 * mrSges(6,3) - t530 * t521 + t529 * t522;
t578 = m(5) * t493 + t498 * mrSges(5,1) + t499 * mrSges(5,2) + t529 * t519 + t530 * t520 + t480;
t577 = m(4) * t536 - t537 * mrSges(4,1) + t538 * mrSges(4,2) - t548 * t541 + t549 * t542 + t578;
t576 = -m(3) * t546 + mrSges(3,1) * t593 - t577 + (t565 * t575 + t602) * mrSges(3,3);
t474 = t576 + (mrSges(2,1) - t603) * qJDD(1) - t575 * mrSges(2,2) + m(2) * t553;
t601 = t572 * t455 + t574 * t474;
t456 = t569 * t462 + t568 * t463;
t600 = t611 * t529 - t606 * t530 - t604 * t567;
t599 = t604 * t529 - t605 * t530 + t610 * t567;
t598 = -t606 * t529 + t612 * t530 + t605 * t567;
t582 = Ifges(3,5) * t568 + Ifges(3,6) * t569;
t596 = t575 * t582;
t590 = t574 * t455 - t572 * t474;
t584 = Ifges(3,1) * t568 + Ifges(3,4) * t569;
t583 = Ifges(3,4) * t568 + Ifges(3,2) * t569;
t528 = Ifges(4,1) * t549 + Ifges(4,4) * t548 + Ifges(4,5) * qJD(3);
t527 = Ifges(4,4) * t549 + Ifges(4,2) * t548 + Ifges(4,6) * qJD(3);
t526 = Ifges(4,5) * t549 + Ifges(4,6) * t548 + Ifges(4,3) * qJD(3);
t466 = mrSges(5,2) * t493 + mrSges(6,2) * t483 - mrSges(5,3) * t486 - mrSges(6,3) * t485 - qJ(5) * t480 - t606 * t498 + t612 * t499 + t599 * t529 + t605 * t564 + t600 * t567;
t465 = -mrSges(5,1) * t493 - mrSges(6,1) * t485 + mrSges(6,2) * t482 + mrSges(5,3) * t487 - pkin(4) * t480 - t611 * t498 + t606 * t499 + t599 * t530 + t604 * t564 + t598 * t567;
t458 = mrSges(4,2) * t536 - mrSges(4,3) * t500 + Ifges(4,1) * t538 + Ifges(4,4) * t537 + Ifges(4,5) * qJDD(3) - pkin(7) * t472 - qJD(3) * t527 - t570 * t465 + t609 * t466 + t548 * t526;
t457 = -mrSges(4,1) * t536 + mrSges(4,3) * t501 + Ifges(4,4) * t538 + Ifges(4,2) * t537 + Ifges(4,6) * qJDD(3) - pkin(3) * t578 + pkin(7) * t587 + qJD(3) * t528 + t609 * t465 + t570 * t466 - t549 * t526;
t452 = (qJ(5) * t513 - t598) * t529 + (pkin(4) * t513 + t600) * t530 - qJ(5) * t592 - pkin(4) * t586 + (Ifges(2,6) - t582) * qJDD(1) - pkin(1) * t456 + t610 * t564 - Ifges(4,3) * qJDD(3) + mrSges(2,1) * g(3) + (qJ(5) * mrSges(6,2) + t604) * t498 + (pkin(4) * mrSges(6,2) - t605) * t499 + t548 * t528 - t549 * t527 + mrSges(2,3) * t554 - Ifges(4,6) * t537 - Ifges(4,5) * t538 - mrSges(3,1) * t539 + mrSges(3,2) * t540 - pkin(2) * t464 + (-t568 * t583 + t569 * t584 + Ifges(2,5)) * t575 - pkin(3) * t472 - mrSges(6,3) * t482 + mrSges(6,1) * t483 - mrSges(5,1) * t486 + mrSges(5,2) * t487 - mrSges(4,1) * t500 + mrSges(4,2) * t501;
t451 = mrSges(3,2) * t546 - mrSges(3,3) * t539 - pkin(6) * t464 + t584 * qJDD(1) - t571 * t457 + t573 * t458 + t569 * t596;
t450 = -mrSges(3,1) * t546 + mrSges(3,3) * t540 - pkin(2) * t577 + pkin(6) * t588 + t583 * qJDD(1) + t573 * t457 + t571 * t458 - t568 * t596;
t449 = -mrSges(2,2) * g(3) - mrSges(2,3) * t553 + Ifges(2,5) * qJDD(1) - t575 * Ifges(2,6) - qJ(2) * t456 - t568 * t450 + t569 * t451;
t1 = [-m(1) * g(1) + t590; -m(1) * g(2) + t601; (-m(1) - m(2)) * g(3) + t456; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t601 + t574 * t449 - t572 * t452; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t590 + t572 * t449 + t574 * t452; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(2,3) * qJDD(1) + mrSges(2,1) * t553 - mrSges(2,2) * t554 + t568 * t451 + t569 * t450 + pkin(1) * (-qJDD(1) * t603 + t576) + qJ(2) * t589;];
tauB = t1;
