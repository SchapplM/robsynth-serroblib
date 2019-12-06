% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPRPR2
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
% Datum: 2019-12-05 17:50
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPRPR2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR2_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR2_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR2_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR2_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR2_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:49:33
% EndTime: 2019-12-05 17:49:36
% DurationCPUTime: 2.93s
% Computational Cost: add. (42687->206), mult. (61438->259), div. (0->0), fcn. (35545->10), ass. (0->99)
t580 = qJD(1) + qJD(3);
t576 = t580 ^ 2;
t586 = cos(pkin(9));
t620 = pkin(4) * t586;
t584 = sin(pkin(9));
t619 = mrSges(5,2) * t584;
t579 = t586 ^ 2;
t618 = t576 * t579;
t577 = qJDD(1) + qJDD(3);
t617 = t577 * t586;
t605 = Ifges(5,5) * t584 + Ifges(5,6) * t586;
t616 = t576 * t605;
t590 = sin(qJ(1));
t593 = cos(qJ(1));
t564 = t593 * g(2) + t590 * g(3);
t558 = qJDD(1) * pkin(1) + t564;
t563 = t590 * g(2) - t593 * g(3);
t594 = qJD(1) ^ 2;
t559 = -t594 * pkin(1) + t563;
t585 = sin(pkin(8));
t587 = cos(pkin(8));
t547 = t587 * t558 - t585 * t559;
t544 = qJDD(1) * pkin(2) + t547;
t548 = t585 * t558 + t587 * t559;
t545 = -t594 * pkin(2) + t548;
t589 = sin(qJ(3));
t592 = cos(qJ(3));
t532 = t589 * t544 + t592 * t545;
t529 = -t576 * pkin(3) + t577 * qJ(4) + t532;
t583 = -g(1) + qJDD(2);
t614 = qJD(4) * t580;
t615 = t586 * t583 - 0.2e1 * t584 * t614;
t522 = (-pkin(7) * t577 + t576 * t620 - t529) * t584 + t615;
t526 = t584 * t583 + (t529 + 0.2e1 * t614) * t586;
t523 = -pkin(4) * t618 + pkin(7) * t617 + t526;
t588 = sin(qJ(5));
t591 = cos(qJ(5));
t520 = t591 * t522 - t588 * t523;
t600 = -t584 * t588 + t586 * t591;
t551 = t600 * t580;
t601 = t584 * t591 + t586 * t588;
t552 = t601 * t580;
t538 = -t551 * mrSges(6,1) + t552 * mrSges(6,2);
t540 = t551 * qJD(5) + t601 * t577;
t549 = -qJD(5) * mrSges(6,2) + t551 * mrSges(6,3);
t518 = m(6) * t520 + qJDD(5) * mrSges(6,1) - t540 * mrSges(6,3) + qJD(5) * t549 - t552 * t538;
t521 = t588 * t522 + t591 * t523;
t539 = -t552 * qJD(5) + t600 * t577;
t550 = qJD(5) * mrSges(6,1) - t552 * mrSges(6,3);
t519 = m(6) * t521 - qJDD(5) * mrSges(6,2) + t539 * mrSges(6,3) - qJD(5) * t550 + t551 * t538;
t508 = t591 * t518 + t588 * t519;
t525 = -t584 * t529 + t615;
t603 = mrSges(5,3) * t577 + (-mrSges(5,1) * t586 + t619) * t576;
t506 = m(5) * t525 - t603 * t584 + t508;
t608 = -t588 * t518 + t591 * t519;
t507 = m(5) * t526 + t603 * t586 + t608;
t609 = -t584 * t506 + t586 * t507;
t499 = m(4) * t532 - t576 * mrSges(4,1) - t577 * mrSges(4,2) + t609;
t531 = t592 * t544 - t589 * t545;
t604 = qJDD(4) - t531;
t528 = -t577 * pkin(3) - t576 * qJ(4) + t604;
t578 = t584 ^ 2;
t524 = (-pkin(3) - t620) * t577 + (-qJ(4) + (-t578 - t579) * pkin(7)) * t576 + t604;
t598 = m(6) * t524 - t539 * mrSges(6,1) + t540 * mrSges(6,2) - t551 * t549 + t552 * t550;
t597 = -m(5) * t528 + mrSges(5,1) * t617 - t598 + (t576 * t578 + t618) * mrSges(5,3);
t512 = m(4) * t531 - t576 * mrSges(4,2) + (mrSges(4,1) - t619) * t577 + t597;
t494 = t589 * t499 + t592 * t512;
t489 = m(3) * t547 + qJDD(1) * mrSges(3,1) - t594 * mrSges(3,2) + t494;
t610 = t592 * t499 - t589 * t512;
t490 = m(3) * t548 - t594 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t610;
t484 = t587 * t489 + t585 * t490;
t502 = t586 * t506 + t584 * t507;
t613 = m(4) * t583 + t502;
t611 = -t585 * t489 + t587 * t490;
t481 = m(2) * t563 - t594 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t611;
t482 = m(2) * t564 + qJDD(1) * mrSges(2,1) - t594 * mrSges(2,2) + t484;
t612 = t593 * t481 - t590 * t482;
t500 = m(3) * t583 + t613;
t607 = Ifges(5,1) * t584 + Ifges(5,4) * t586;
t606 = Ifges(5,4) * t584 + Ifges(5,2) * t586;
t602 = -t590 * t481 - t593 * t482;
t533 = Ifges(6,5) * t552 + Ifges(6,6) * t551 + Ifges(6,3) * qJD(5);
t535 = Ifges(6,1) * t552 + Ifges(6,4) * t551 + Ifges(6,5) * qJD(5);
t509 = -mrSges(6,1) * t524 + mrSges(6,3) * t521 + Ifges(6,4) * t540 + Ifges(6,2) * t539 + Ifges(6,6) * qJDD(5) + qJD(5) * t535 - t552 * t533;
t534 = Ifges(6,4) * t552 + Ifges(6,2) * t551 + Ifges(6,6) * qJD(5);
t510 = mrSges(6,2) * t524 - mrSges(6,3) * t520 + Ifges(6,1) * t540 + Ifges(6,4) * t539 + Ifges(6,5) * qJDD(5) - qJD(5) * t534 + t551 * t533;
t492 = -mrSges(5,1) * t528 + mrSges(5,3) * t526 - pkin(4) * t598 + pkin(7) * t608 + t591 * t509 + t588 * t510 + t606 * t577 - t584 * t616;
t496 = mrSges(5,2) * t528 - mrSges(5,3) * t525 - pkin(7) * t508 - t588 * t509 + t591 * t510 + t607 * t577 + t586 * t616;
t514 = t577 * t619 - t597;
t599 = mrSges(4,1) * t531 - mrSges(4,2) * t532 + Ifges(4,3) * t577 - pkin(3) * t514 + qJ(4) * t609 + t586 * t492 + t584 * t496;
t596 = mrSges(6,1) * t520 - mrSges(6,2) * t521 + Ifges(6,5) * t540 + Ifges(6,6) * t539 + Ifges(6,3) * qJDD(5) + t552 * t534 - t551 * t535;
t595 = mrSges(2,1) * t564 + mrSges(3,1) * t547 - mrSges(2,2) * t563 - mrSges(3,2) * t548 + pkin(1) * t484 + pkin(2) * t494 + t599 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t485 = -mrSges(4,1) * t583 - mrSges(5,1) * t525 + mrSges(5,2) * t526 + mrSges(4,3) * t532 - pkin(3) * t502 - pkin(4) * t508 + (Ifges(4,6) - t605) * t577 - t596 + (-t584 * t606 + t586 * t607 + Ifges(4,5)) * t576;
t479 = mrSges(4,2) * t583 - mrSges(4,3) * t531 + Ifges(4,5) * t577 - t576 * Ifges(4,6) - qJ(4) * t502 - t584 * t492 + t586 * t496;
t478 = mrSges(3,2) * t583 - mrSges(3,3) * t547 + Ifges(3,5) * qJDD(1) - t594 * Ifges(3,6) - pkin(6) * t494 + t592 * t479 - t589 * t485;
t477 = -mrSges(3,1) * t583 + mrSges(3,3) * t548 + t594 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t613 + pkin(6) * t610 + t589 * t479 + t592 * t485;
t476 = -mrSges(2,2) * g(1) - mrSges(2,3) * t564 + Ifges(2,5) * qJDD(1) - t594 * Ifges(2,6) - qJ(2) * t484 - t585 * t477 + t587 * t478;
t475 = mrSges(2,1) * g(1) + mrSges(2,3) * t563 + t594 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t500 + qJ(2) * t611 + t587 * t477 + t585 * t478;
t1 = [(-m(1) - m(2)) * g(1) + t500; -m(1) * g(2) + t602; -m(1) * g(3) + t612; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t595; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t612 - t593 * t475 - t590 * t476; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t602 - t590 * t475 + t593 * t476; t595; t500; t599; t514; t596;];
tauJB = t1;
