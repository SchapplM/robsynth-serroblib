% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PRPRR4
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PRPRR4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR4_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR4_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR4_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR4_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR4_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR4_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR4_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:50:14
% EndTime: 2019-12-05 15:50:19
% DurationCPUTime: 3.98s
% Computational Cost: add. (48537->225), mult. (87704->293), div. (0->0), fcn. (59923->12), ass. (0->106)
t577 = sin(pkin(9));
t580 = cos(pkin(9));
t566 = t577 * g(1) - t580 * g(2);
t567 = -t580 * g(1) - t577 * g(2);
t575 = -g(3) + qJDD(1);
t584 = sin(qJ(2));
t581 = cos(pkin(5));
t587 = cos(qJ(2));
t606 = t581 * t587;
t578 = sin(pkin(5));
t608 = t578 * t587;
t535 = t566 * t606 - t584 * t567 + t575 * t608;
t533 = qJDD(2) * pkin(2) + t535;
t607 = t581 * t584;
t609 = t578 * t584;
t536 = t566 * t607 + t587 * t567 + t575 * t609;
t589 = qJD(2) ^ 2;
t534 = -t589 * pkin(2) + t536;
t576 = sin(pkin(10));
t579 = cos(pkin(10));
t529 = t576 * t533 + t579 * t534;
t527 = -t589 * pkin(3) + qJDD(2) * pkin(7) + t529;
t550 = -t578 * t566 + t581 * t575;
t549 = qJDD(3) + t550;
t583 = sin(qJ(4));
t586 = cos(qJ(4));
t524 = t586 * t527 + t583 * t549;
t563 = (-pkin(4) * t586 - pkin(8) * t583) * qJD(2);
t588 = qJD(4) ^ 2;
t602 = t586 * qJD(2);
t521 = -t588 * pkin(4) + qJDD(4) * pkin(8) + t563 * t602 + t524;
t528 = t579 * t533 - t576 * t534;
t526 = -qJDD(2) * pkin(3) - t589 * pkin(7) - t528;
t601 = qJD(2) * qJD(4);
t599 = t586 * t601;
t564 = t583 * qJDD(2) + t599;
t600 = t583 * t601;
t565 = t586 * qJDD(2) - t600;
t522 = (-t564 - t599) * pkin(8) + (-t565 + t600) * pkin(4) + t526;
t582 = sin(qJ(5));
t585 = cos(qJ(5));
t518 = -t582 * t521 + t585 * t522;
t603 = qJD(2) * t583;
t560 = t585 * qJD(4) - t582 * t603;
t543 = t560 * qJD(5) + t582 * qJDD(4) + t585 * t564;
t561 = t582 * qJD(4) + t585 * t603;
t544 = -t560 * mrSges(6,1) + t561 * mrSges(6,2);
t572 = qJD(5) - t602;
t547 = -t572 * mrSges(6,2) + t560 * mrSges(6,3);
t558 = qJDD(5) - t565;
t515 = m(6) * t518 + t558 * mrSges(6,1) - t543 * mrSges(6,3) - t561 * t544 + t572 * t547;
t519 = t585 * t521 + t582 * t522;
t542 = -t561 * qJD(5) + t585 * qJDD(4) - t582 * t564;
t548 = t572 * mrSges(6,1) - t561 * mrSges(6,3);
t516 = m(6) * t519 - t558 * mrSges(6,2) + t542 * mrSges(6,3) + t560 * t544 - t572 * t548;
t509 = -t582 * t515 + t585 * t516;
t605 = t586 * t549;
t520 = -qJDD(4) * pkin(4) - t588 * pkin(8) - t605 + (qJD(2) * t563 + t527) * t583;
t537 = Ifges(6,5) * t561 + Ifges(6,6) * t560 + Ifges(6,3) * t572;
t539 = Ifges(6,1) * t561 + Ifges(6,4) * t560 + Ifges(6,5) * t572;
t510 = -mrSges(6,1) * t520 + mrSges(6,3) * t519 + Ifges(6,4) * t543 + Ifges(6,2) * t542 + Ifges(6,6) * t558 - t561 * t537 + t572 * t539;
t538 = Ifges(6,4) * t561 + Ifges(6,2) * t560 + Ifges(6,6) * t572;
t511 = mrSges(6,2) * t520 - mrSges(6,3) * t518 + Ifges(6,1) * t543 + Ifges(6,4) * t542 + Ifges(6,5) * t558 + t560 * t537 - t572 * t538;
t517 = -m(6) * t520 + t542 * mrSges(6,1) - t543 * mrSges(6,2) + t560 * t547 - t561 * t548;
t523 = -t583 * t527 + t605;
t554 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t583 + Ifges(5,2) * t586) * qJD(2);
t555 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t583 + Ifges(5,4) * t586) * qJD(2);
t610 = mrSges(5,1) * t523 - mrSges(5,2) * t524 + Ifges(5,5) * t564 + Ifges(5,6) * t565 + Ifges(5,3) * qJDD(4) + pkin(4) * t517 + pkin(8) * t509 + t585 * t510 + t582 * t511 + (t583 * t554 - t586 * t555) * qJD(2);
t562 = (-mrSges(5,1) * t586 + mrSges(5,2) * t583) * qJD(2);
t568 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t603;
t507 = m(5) * t524 - qJDD(4) * mrSges(5,2) + t565 * mrSges(5,3) - qJD(4) * t568 + t562 * t602 + t509;
t569 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t602;
t513 = m(5) * t523 + qJDD(4) * mrSges(5,1) - t564 * mrSges(5,3) + qJD(4) * t569 - t562 * t603 + t517;
t596 = t586 * t507 - t583 * t513;
t498 = m(4) * t529 - t589 * mrSges(4,1) - qJDD(2) * mrSges(4,2) + t596;
t508 = t585 * t515 + t582 * t516;
t592 = -m(5) * t526 + t565 * mrSges(5,1) - t564 * mrSges(5,2) - t568 * t603 + t569 * t602 - t508;
t504 = m(4) * t528 + qJDD(2) * mrSges(4,1) - t589 * mrSges(4,2) + t592;
t493 = t576 * t498 + t579 * t504;
t491 = m(3) * t535 + qJDD(2) * mrSges(3,1) - t589 * mrSges(3,2) + t493;
t597 = t579 * t498 - t576 * t504;
t492 = m(3) * t536 - t589 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t597;
t502 = t583 * t507 + t586 * t513;
t501 = m(4) * t549 + t502;
t500 = m(3) * t550 + t501;
t479 = t491 * t606 + t492 * t607 - t578 * t500;
t477 = m(2) * t566 + t479;
t484 = -t584 * t491 + t587 * t492;
t483 = m(2) * t567 + t484;
t604 = t580 * t477 + t577 * t483;
t478 = t491 * t608 + t492 * t609 + t581 * t500;
t598 = -t577 * t477 + t580 * t483;
t595 = m(2) * t575 + t478;
t553 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t583 + Ifges(5,6) * t586) * qJD(2);
t494 = mrSges(5,2) * t526 - mrSges(5,3) * t523 + Ifges(5,1) * t564 + Ifges(5,4) * t565 + Ifges(5,5) * qJDD(4) - pkin(8) * t508 - qJD(4) * t554 - t582 * t510 + t585 * t511 + t553 * t602;
t591 = mrSges(6,1) * t518 - mrSges(6,2) * t519 + Ifges(6,5) * t543 + Ifges(6,6) * t542 + Ifges(6,3) * t558 + t561 * t538 - t560 * t539;
t495 = -mrSges(5,1) * t526 + mrSges(5,3) * t524 + Ifges(5,4) * t564 + Ifges(5,2) * t565 + Ifges(5,6) * qJDD(4) - pkin(4) * t508 + qJD(4) * t555 - t553 * t603 - t591;
t480 = mrSges(4,2) * t549 - mrSges(4,3) * t528 + Ifges(4,5) * qJDD(2) - t589 * Ifges(4,6) - pkin(7) * t502 + t586 * t494 - t583 * t495;
t485 = -mrSges(4,1) * t549 + mrSges(4,3) * t529 + t589 * Ifges(4,5) + Ifges(4,6) * qJDD(2) - pkin(3) * t502 - t610;
t473 = -mrSges(3,1) * t550 + mrSges(3,3) * t536 + t589 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t501 + qJ(3) * t597 + t576 * t480 + t579 * t485;
t474 = mrSges(3,2) * t550 - mrSges(3,3) * t535 + Ifges(3,5) * qJDD(2) - t589 * Ifges(3,6) - qJ(3) * t493 + t579 * t480 - t576 * t485;
t593 = pkin(6) * t484 + t473 * t587 + t474 * t584;
t475 = mrSges(3,1) * t535 - mrSges(3,2) * t536 + mrSges(4,1) * t528 - mrSges(4,2) * t529 + t583 * t494 + t586 * t495 + pkin(3) * t592 + pkin(7) * t596 + pkin(2) * t493 + (Ifges(3,3) + Ifges(4,3)) * qJDD(2);
t472 = mrSges(2,2) * t575 - mrSges(2,3) * t566 - t584 * t473 + t587 * t474 + (-t478 * t578 - t479 * t581) * pkin(6);
t471 = -mrSges(2,1) * t575 + mrSges(2,3) * t567 - pkin(1) * t478 - t578 * t475 + t593 * t581;
t1 = [-m(1) * g(1) + t598; -m(1) * g(2) + t604; -m(1) * g(3) + t595; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t604 - t577 * t471 + t580 * t472; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t598 + t580 * t471 + t577 * t472; -mrSges(1,1) * g(2) + mrSges(2,1) * t566 + mrSges(1,2) * g(1) - mrSges(2,2) * t567 + pkin(1) * t479 + t581 * t475 + t593 * t578; t595; t475; t501; t610; t591;];
tauJB = t1;
