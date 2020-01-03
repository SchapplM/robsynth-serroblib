% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPRRP4
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
% Datum: 2020-01-03 11:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPRRP4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP4_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP4_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP4_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP4_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP4_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:49:37
% EndTime: 2020-01-03 11:49:43
% DurationCPUTime: 3.53s
% Computational Cost: add. (28740->267), mult. (69192->337), div. (0->0), fcn. (45196->8), ass. (0->113)
t588 = sin(qJ(1));
t591 = cos(qJ(1));
t568 = -t588 * g(2) + t591 * g(3);
t592 = qJD(1) ^ 2;
t633 = -t592 * pkin(1) + qJDD(1) * qJ(2) + (2 * qJD(1) * qJD(2)) + t568;
t632 = Ifges(5,1) + Ifges(6,1);
t628 = Ifges(5,4) + Ifges(6,4);
t627 = Ifges(5,5) + Ifges(6,5);
t631 = Ifges(5,2) + Ifges(6,2);
t626 = -Ifges(5,6) - Ifges(6,6);
t630 = -Ifges(5,3) - Ifges(6,3);
t584 = sin(pkin(8));
t585 = cos(pkin(8));
t541 = -t585 * g(1) - t633 * t584;
t614 = t585 * qJD(1);
t572 = qJD(3) - t614;
t587 = sin(qJ(3));
t615 = t584 * qJD(1);
t607 = t587 * t615;
t554 = -t572 * mrSges(4,2) - mrSges(4,3) * t607;
t590 = cos(qJ(3));
t606 = t590 * t615;
t555 = t572 * mrSges(4,1) - mrSges(4,3) * t606;
t629 = -t554 * t587 - t555 * t590;
t625 = mrSges(3,2) * t584;
t622 = t584 ^ 2 * t592;
t542 = -t584 * g(1) + t633 * t585;
t563 = (-mrSges(3,1) * t585 + t625) * qJD(1);
t600 = -pkin(2) * t585 - pkin(6) * t584;
t565 = t600 * qJD(1);
t530 = t565 * t614 + t542;
t569 = -t591 * g(2) - t588 * g(3);
t596 = -t592 * qJ(2) + qJDD(2) - t569;
t543 = (-pkin(1) + t600) * qJDD(1) + t596;
t540 = t590 * t543;
t612 = qJD(1) * qJD(3);
t559 = (qJDD(1) * t590 - t587 * t612) * t584;
t611 = t585 * qJDD(1);
t571 = qJDD(3) - t611;
t507 = t571 * pkin(3) - t559 * pkin(7) + t540 + (-pkin(3) * t590 * t622 - pkin(7) * t572 * t615 - t530) * t587;
t511 = t590 * t530 + t587 * t543;
t557 = t572 * pkin(3) - pkin(7) * t606;
t558 = (-qJDD(1) * t587 - t590 * t612) * t584;
t610 = t587 ^ 2 * t622;
t508 = -pkin(3) * t610 + t558 * pkin(7) - t572 * t557 + t511;
t586 = sin(qJ(4));
t589 = cos(qJ(4));
t500 = t589 * t507 - t586 * t508;
t549 = (-t586 * t590 - t587 * t589) * t615;
t518 = t549 * qJD(4) + t586 * t558 + t589 * t559;
t550 = (-t586 * t587 + t589 * t590) * t615;
t531 = -t549 * mrSges(6,1) + t550 * mrSges(6,2);
t532 = -t549 * mrSges(5,1) + t550 * mrSges(5,2);
t570 = qJD(4) + t572;
t535 = -t570 * mrSges(5,2) + t549 * mrSges(5,3);
t567 = qJDD(4) + t571;
t497 = -0.2e1 * qJD(5) * t550 + (t549 * t570 - t518) * qJ(5) + (t549 * t550 + t567) * pkin(4) + t500;
t534 = -t570 * mrSges(6,2) + t549 * mrSges(6,3);
t609 = m(6) * t497 + t567 * mrSges(6,1) + t570 * t534;
t489 = m(5) * t500 + t567 * mrSges(5,1) + t570 * t535 + (-t531 - t532) * t550 + (-mrSges(5,3) - mrSges(6,3)) * t518 + t609;
t501 = t586 * t507 + t589 * t508;
t517 = -t550 * qJD(4) + t589 * t558 - t586 * t559;
t537 = t570 * mrSges(6,1) - t550 * mrSges(6,3);
t538 = t570 * mrSges(5,1) - t550 * mrSges(5,3);
t536 = t570 * pkin(4) - t550 * qJ(5);
t548 = t549 ^ 2;
t499 = -t548 * pkin(4) + t517 * qJ(5) + 0.2e1 * qJD(5) * t549 - t570 * t536 + t501;
t608 = m(6) * t499 + t517 * mrSges(6,3) + t549 * t531;
t492 = m(5) * t501 + t517 * mrSges(5,3) + t549 * t532 + (-t537 - t538) * t570 + (-mrSges(5,2) - mrSges(6,2)) * t567 + t608;
t487 = t589 * t489 + t586 * t492;
t510 = -t587 * t530 + t540;
t556 = (mrSges(4,1) * t587 + mrSges(4,2) * t590) * t615;
t484 = m(4) * t510 + t571 * mrSges(4,1) - t559 * mrSges(4,3) + t572 * t554 - t556 * t606 + t487;
t601 = -t586 * t489 + t589 * t492;
t485 = m(4) * t511 - t571 * mrSges(4,2) + t558 * mrSges(4,3) - t572 * t555 - t556 * t607 + t601;
t602 = -t587 * t484 + t590 * t485;
t616 = qJDD(1) * mrSges(3,3);
t480 = m(3) * t542 + (qJD(1) * t563 + t616) * t585 + t602;
t529 = t565 * t615 - t541;
t509 = -t558 * pkin(3) - pkin(7) * t610 + t557 * t606 + t529;
t503 = -t517 * pkin(4) - t548 * qJ(5) + t550 * t536 + qJDD(5) + t509;
t597 = m(6) * t503 - t517 * mrSges(6,1) + t518 * mrSges(6,2) - t549 * t534 + t550 * t537;
t594 = m(5) * t509 - t517 * mrSges(5,1) + t518 * mrSges(5,2) - t549 * t535 + t550 * t538 + t597;
t593 = -m(4) * t529 + t558 * mrSges(4,1) - t559 * mrSges(4,2) - t594;
t494 = t593 + m(3) * t541 + (-t616 + (-t563 + t629) * qJD(1)) * t584;
t603 = t585 * t480 - t584 * t494;
t473 = m(2) * t568 - t592 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t603;
t481 = t590 * t484 + t587 * t485;
t561 = -qJDD(1) * pkin(1) + t596;
t595 = -m(3) * t561 + mrSges(3,1) * t611 - t481 + (t585 ^ 2 * t592 + t622) * mrSges(3,3);
t477 = m(2) * t569 - t592 * mrSges(2,2) + (mrSges(2,1) - t625) * qJDD(1) + t595;
t621 = t588 * t473 + t591 * t477;
t474 = t584 * t480 + t585 * t494;
t620 = t626 * t549 - t627 * t550 + t630 * t570;
t619 = -t631 * t549 - t628 * t550 + t626 * t570;
t618 = t628 * t549 + t632 * t550 + t627 * t570;
t604 = -t591 * t473 + t588 * t477;
t599 = Ifges(3,1) * t584 + Ifges(3,4) * t585;
t598 = Ifges(3,5) * t584 + Ifges(3,6) * t585;
t564 = t598 * qJD(1);
t546 = Ifges(4,5) * t572 + (Ifges(4,1) * t590 - Ifges(4,4) * t587) * t615;
t545 = Ifges(4,6) * t572 + (Ifges(4,4) * t590 - Ifges(4,2) * t587) * t615;
t544 = Ifges(4,3) * t572 + (Ifges(4,5) * t590 - Ifges(4,6) * t587) * t615;
t495 = -t518 * mrSges(6,3) - t550 * t531 + t609;
t486 = mrSges(5,2) * t509 + mrSges(6,2) * t503 - mrSges(5,3) * t500 - mrSges(6,3) * t497 - qJ(5) * t495 + t628 * t517 + t632 * t518 - t620 * t549 + t627 * t567 + t619 * t570;
t482 = -mrSges(5,1) * t509 + mrSges(5,3) * t501 - mrSges(6,1) * t503 + mrSges(6,3) * t499 - pkin(4) * t597 + qJ(5) * t608 + (-qJ(5) * t537 + t618) * t570 + (-qJ(5) * mrSges(6,2) - t626) * t567 + t620 * t550 + t628 * t518 + t631 * t517;
t471 = mrSges(4,2) * t529 - mrSges(4,3) * t510 + Ifges(4,1) * t559 + Ifges(4,4) * t558 + Ifges(4,5) * t571 - pkin(7) * t487 - t586 * t482 + t589 * t486 - t544 * t607 - t572 * t545;
t470 = -mrSges(4,1) * t529 + mrSges(4,3) * t511 + Ifges(4,4) * t559 + Ifges(4,2) * t558 + Ifges(4,6) * t571 - pkin(3) * t594 + pkin(7) * t601 + t589 * t482 + t586 * t486 - t544 * t606 + t572 * t546;
t469 = Ifges(3,2) * t611 - mrSges(3,1) * t561 - mrSges(4,1) * t510 - mrSges(5,1) * t500 - mrSges(6,1) * t497 + mrSges(4,2) * t511 + mrSges(5,2) * t501 + mrSges(6,2) * t499 + mrSges(3,3) * t542 - Ifges(4,5) * t559 - Ifges(4,6) * t558 - Ifges(4,3) * t571 - pkin(2) * t481 - pkin(3) * t487 - pkin(4) * t495 + t630 * t567 + t619 * t550 + t618 * t549 - t627 * t518 + t626 * t517 + (Ifges(3,4) * qJDD(1) + (-t545 * t590 - t546 * t587 - t564) * qJD(1)) * t584;
t468 = mrSges(3,2) * t561 - mrSges(3,3) * t541 - pkin(6) * t481 + t599 * qJDD(1) - t587 * t470 + t590 * t471 + t564 * t614;
t467 = t592 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t568 - mrSges(3,1) * t541 + mrSges(3,2) * t542 - t587 * t471 - t590 * t470 - pkin(2) * t593 - pkin(6) * t602 - pkin(1) * t474 + (Ifges(2,6) - t598) * qJDD(1) + (-pkin(2) * t629 * t584 + (-t584 * (Ifges(3,4) * t584 + Ifges(3,2) * t585) + t585 * t599) * qJD(1)) * qJD(1);
t466 = -mrSges(2,2) * g(1) - mrSges(2,3) * t569 + Ifges(2,5) * qJDD(1) - t592 * Ifges(2,6) - qJ(2) * t474 + t585 * t468 - t584 * t469;
t1 = [(-m(1) - m(2)) * g(1) + t474; -m(1) * g(2) + t621; -m(1) * g(3) + t604; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + Ifges(2,3) * qJDD(1) + mrSges(2,1) * t569 - mrSges(2,2) * t568 + t584 * t468 + t585 * t469 + pkin(1) * (-qJDD(1) * t625 + t595) + qJ(2) * t603; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t604 + t588 * t466 + t591 * t467; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t621 - t591 * t466 + t588 * t467;];
tauB = t1;
