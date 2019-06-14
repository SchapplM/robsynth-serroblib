% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJ [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 03:01
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRRPPR3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR3_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPPR3_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR3_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR3_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPPR3_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 02:57:49
% EndTime: 2019-05-05 02:57:51
% DurationCPUTime: 1.58s
% Computational Cost: add. (4700->257), mult. (9354->298), div. (0->0), fcn. (5104->10), ass. (0->112)
t538 = sin(qJ(3));
t541 = cos(qJ(3));
t564 = Ifges(4,4) - Ifges(5,5) + Ifges(6,4);
t593 = t538 * (Ifges(4,1) + Ifges(5,1) + Ifges(6,2)) + t541 * t564;
t592 = -t541 * (Ifges(4,2) + Ifges(5,3) + Ifges(6,1)) - t538 * t564;
t532 = sin(pkin(10));
t534 = cos(pkin(10));
t507 = g(1) * t532 - g(2) * t534;
t530 = -g(3) + qJDD(1);
t533 = sin(pkin(6));
t535 = cos(pkin(6));
t589 = t507 * t535 + t530 * t533;
t563 = Ifges(4,5) + Ifges(5,4) + Ifges(6,6);
t562 = Ifges(4,6) - Ifges(5,6) + Ifges(6,5);
t586 = t538 * (t592 * qJD(2) - t562 * qJD(3));
t508 = -g(1) * t534 - g(2) * t532;
t539 = sin(qJ(2));
t542 = cos(qJ(2));
t448 = t542 * t508 + t539 * t589;
t544 = qJD(2) ^ 2;
t446 = -pkin(2) * t544 + qJDD(2) * pkin(8) + t448;
t462 = -t507 * t533 + t530 * t535;
t441 = t541 * t446 + t538 * t462;
t497 = (-pkin(3) * t541 - qJ(4) * t538) * qJD(2);
t567 = qJD(2) * t541;
t585 = qJDD(3) * qJ(4) + t497 * t567 + t441;
t566 = qJD(2) * qJD(3);
t559 = t541 * t566;
t502 = qJDD(2) * t538 + t559;
t565 = qJD(2) * qJD(5);
t584 = -0.2e1 * t538 * t565 + (-t502 + t559) * qJ(5);
t543 = qJD(3) ^ 2;
t440 = -t538 * t446 + t541 * t462;
t568 = qJD(2) * t538;
t555 = t497 * t568 + qJDD(4) - t440;
t438 = -qJDD(3) * pkin(3) - t543 * qJ(4) + t555;
t515 = mrSges(5,2) * t567 + qJD(3) * mrSges(5,3);
t583 = m(5) * t438 - qJDD(3) * mrSges(5,1) - qJD(3) * t515;
t582 = -2 * qJD(4);
t581 = 2 * qJD(4);
t580 = -pkin(3) - pkin(9);
t579 = pkin(4) + pkin(9);
t577 = t543 * pkin(3);
t576 = mrSges(4,3) + mrSges(5,2);
t558 = t538 * t566;
t503 = qJDD(2) * t541 - t558;
t575 = mrSges(5,1) * t503;
t572 = t541 ^ 2 * t544;
t570 = t541 * t544;
t509 = -qJD(3) * pkin(4) - qJ(5) * t568;
t447 = -t539 * t508 + t542 * t589;
t445 = -qJDD(2) * pkin(2) - t544 * pkin(8) - t447;
t554 = -t503 * pkin(3) + t445 + (-t502 - t559) * qJ(4);
t547 = -qJ(5) * t572 + qJDD(5) - t554 + (t509 + t581) * t568;
t428 = t579 * t503 + t547 + (pkin(5) * t541 + t538 * t580) * t566 + t502 * pkin(5);
t501 = (pkin(5) * t538 + pkin(9) * t541) * qJD(2);
t431 = (-pkin(5) - qJ(4)) * t543 + (-pkin(4) * t570 - qJD(2) * t501) * t538 + (-pkin(3) - t579) * qJDD(3) + t555 + t584;
t537 = sin(qJ(6));
t540 = cos(qJ(6));
t426 = t428 * t540 - t431 * t537;
t495 = -qJD(3) * t540 + t537 * t567;
t457 = qJD(6) * t495 - qJDD(3) * t537 - t503 * t540;
t496 = -qJD(3) * t537 - t540 * t567;
t458 = -mrSges(7,1) * t495 + mrSges(7,2) * t496;
t520 = qJD(6) + t568;
t460 = -mrSges(7,2) * t520 + mrSges(7,3) * t495;
t491 = qJDD(6) + t502;
t423 = m(7) * t426 + mrSges(7,1) * t491 - mrSges(7,3) * t457 - t458 * t496 + t460 * t520;
t427 = t428 * t537 + t431 * t540;
t456 = -qJD(6) * t496 - qJDD(3) * t540 + t503 * t537;
t461 = mrSges(7,1) * t520 - mrSges(7,3) * t496;
t424 = m(7) * t427 - mrSges(7,2) * t491 + mrSges(7,3) * t456 + t458 * t495 - t461 * t520;
t414 = t540 * t423 + t537 * t424;
t569 = -t537 * t423 + t540 * t424;
t560 = t593 * qJD(2) + t563 * qJD(3);
t498 = (-mrSges(5,1) * t541 - mrSges(5,3) * t538) * qJD(2);
t499 = (-mrSges(4,1) * t541 + mrSges(4,2) * t538) * qJD(2);
t513 = -qJD(3) * mrSges(6,1) + mrSges(6,3) * t567;
t514 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t567;
t433 = (-t538 * t570 - qJDD(3)) * pkin(4) + t438 + t584;
t500 = (mrSges(6,1) * t538 - mrSges(6,2) * t541) * qJD(2);
t556 = -m(6) * t433 + t500 * t568 - t569;
t409 = m(4) * t440 + (mrSges(4,1) - mrSges(6,2)) * qJDD(3) + (-t513 + t514) * qJD(3) + (-t498 - t499) * t568 + (mrSges(6,3) - t576) * t502 + t556 - t583;
t511 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t568;
t524 = qJD(3) * t581;
t437 = t524 - t577 + t585;
t512 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t568;
t552 = pkin(4) * t572 + t503 * qJ(5) - t585;
t432 = 0.2e1 * t541 * t565 + t577 + (t582 - t509) * qJD(3) + t552;
t510 = qJD(3) * mrSges(6,2) - mrSges(6,3) * t568;
t430 = qJDD(3) * pkin(5) + qJD(3) * t509 + t524 + t580 * t543 + (-0.2e1 * qJD(5) - t501) * t567 - t552;
t553 = -m(7) * t430 + t456 * mrSges(7,1) - t457 * mrSges(7,2) + t495 * t460 - t496 * t461;
t548 = -m(6) * t432 + qJDD(3) * mrSges(6,1) - t503 * mrSges(6,3) + qJD(3) * t510 - t553;
t546 = m(5) * t437 + qJDD(3) * mrSges(5,3) + qJD(3) * t512 + t498 * t567 + t548;
t416 = t546 + (t499 - t500) * t567 - qJDD(3) * mrSges(4,2) + t576 * t503 + m(4) * t441 - qJD(3) * t511;
t557 = -t409 * t538 + t541 * t416;
t435 = -pkin(3) * t558 + t503 * pkin(4) + t547;
t412 = m(6) * t435 + t502 * mrSges(6,1) - t503 * mrSges(6,2) + t510 * t568 - t513 * t567 + t414;
t450 = Ifges(7,4) * t496 + Ifges(7,2) * t495 + Ifges(7,6) * t520;
t451 = Ifges(7,1) * t496 + Ifges(7,4) * t495 + Ifges(7,5) * t520;
t551 = mrSges(7,1) * t426 - mrSges(7,2) * t427 + Ifges(7,5) * t457 + Ifges(7,6) * t456 + Ifges(7,3) * t491 + t496 * t450 - t495 * t451;
t550 = qJDD(3) * mrSges(6,2) + qJD(3) * t513 - t556;
t439 = (pkin(3) * qJD(3) + t582) * t568 + t554;
t549 = m(5) * t439 - t502 * mrSges(5,3) - t512 * t568 - t515 * t567 - t412;
t545 = -m(4) * t445 + t503 * mrSges(4,1) - t511 * t568 + t514 * t567 - t549;
t449 = Ifges(7,5) * t496 + Ifges(7,6) * t495 + Ifges(7,3) * t520;
t418 = mrSges(7,2) * t430 - mrSges(7,3) * t426 + Ifges(7,1) * t457 + Ifges(7,4) * t456 + Ifges(7,5) * t491 + t449 * t495 - t450 * t520;
t417 = -mrSges(7,1) * t430 + mrSges(7,3) * t427 + Ifges(7,4) * t457 + Ifges(7,2) * t456 + Ifges(7,6) * t491 - t449 * t496 + t451 * t520;
t413 = -t502 * mrSges(6,3) + t550;
t411 = t498 * t568 + (mrSges(5,2) - mrSges(6,3)) * t502 + t550 + t583;
t410 = t549 - t575;
t1 = [m(2) * t530 + t535 * (m(3) * t462 + t409 * t541 + t416 * t538) + (t539 * (m(3) * t448 - mrSges(3,1) * t544 - qJDD(2) * mrSges(3,2) + t557) + t542 * (m(3) * t447 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t544 - mrSges(4,2) * t502 + t545 + t575)) * t533; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t447 - mrSges(3,2) * t448 + t538 * (mrSges(6,1) * t435 + mrSges(4,2) * t445 + mrSges(5,2) * t438 - mrSges(4,3) * t440 - mrSges(5,3) * t439 - mrSges(6,3) * t433 + pkin(5) * t414 - qJ(4) * t410 - qJ(5) * t413 + t551) + t541 * (-pkin(3) * t410 + pkin(4) * t412 + pkin(9) * t414 + mrSges(6,3) * t432 - mrSges(6,2) * t435 + mrSges(5,2) * t437 - mrSges(5,1) * t439 + mrSges(4,3) * t441 - mrSges(4,1) * t445 + t537 * t417 - t540 * t418 - qJ(5) * (-t500 * t567 + t548)) + pkin(2) * t545 + pkin(8) * t557 + (pkin(2) * mrSges(5,1) - t592) * t503 + (-pkin(2) * mrSges(4,2) + t593) * t502 + (t538 * t563 + t541 * t562) * qJDD(3) + (t541 * t560 + t586) * qJD(3); -pkin(3) * t411 - pkin(4) * t413 - pkin(9) * t569 - mrSges(6,1) * t432 + mrSges(6,2) * t433 + mrSges(5,3) * t437 - mrSges(5,1) * t438 + mrSges(4,1) * t440 - mrSges(4,2) * t441 - pkin(5) * t553 - t537 * t418 - t540 * t417 + qJ(4) * t546 + (mrSges(5,2) * qJ(4) + t562) * t503 + t563 * t502 + (Ifges(4,3) + Ifges(5,2) + Ifges(6,3)) * qJDD(3) + (-t586 + (-qJ(4) * t500 - t560) * t541) * qJD(2); t411; t412; t551;];
tauJ  = t1;
