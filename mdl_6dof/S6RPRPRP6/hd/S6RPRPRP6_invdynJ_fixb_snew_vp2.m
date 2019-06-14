% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRPRP6
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
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
% Datum: 2019-05-05 17:55
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRPRP6_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP6_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP6_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP6_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP6_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP6_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP6_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP6_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP6_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:52:46
% EndTime: 2019-05-05 17:52:50
% DurationCPUTime: 2.56s
% Computational Cost: add. (14937->286), mult. (35878->329), div. (0->0), fcn. (25175->8), ass. (0->124)
t580 = Ifges(6,4) + Ifges(7,4);
t596 = Ifges(6,2) + Ifges(7,2);
t590 = Ifges(6,6) + Ifges(7,6);
t595 = -2 * qJD(4);
t594 = Ifges(4,1) + Ifges(5,2);
t593 = Ifges(6,1) + Ifges(7,1);
t579 = Ifges(4,5) - Ifges(5,4);
t592 = Ifges(6,5) + Ifges(7,5);
t591 = -Ifges(4,2) - Ifges(5,3);
t578 = Ifges(4,6) - Ifges(5,5);
t577 = -Ifges(5,6) - Ifges(4,4);
t589 = Ifges(4,3) + Ifges(5,1);
t588 = Ifges(6,3) + Ifges(7,3);
t533 = sin(qJ(3));
t531 = cos(pkin(9));
t584 = cos(qJ(3));
t556 = t531 * t584;
t530 = sin(pkin(9));
t564 = qJD(1) * t530;
t515 = -qJD(1) * t556 + t533 * t564;
t532 = sin(qJ(5));
t535 = cos(qJ(5));
t500 = -qJD(3) * t532 + t515 * t535;
t501 = qJD(3) * t535 + t515 * t532;
t546 = t584 * t530 + t531 * t533;
t516 = t546 * qJD(1);
t512 = qJD(5) + t516;
t587 = t596 * t500 + t580 * t501 + t590 * t512;
t538 = qJD(1) ^ 2;
t534 = sin(qJ(1));
t536 = cos(qJ(1));
t551 = -g(1) * t536 - g(2) * t534;
t517 = -pkin(1) * t538 + qJDD(1) * qJ(2) + t551;
t561 = qJD(1) * qJD(2);
t554 = -t531 * g(3) - 0.2e1 * t530 * t561;
t575 = pkin(7) * qJDD(1);
t583 = pkin(2) * t538;
t476 = (t531 * t583 - t517 - t575) * t530 + t554;
t498 = -g(3) * t530 + (t517 + 0.2e1 * t561) * t531;
t527 = t531 ^ 2;
t477 = -t527 * t583 + t531 * t575 + t498;
t443 = t533 * t476 + t584 * t477;
t488 = pkin(3) * t515 - qJ(4) * t516;
t537 = qJD(3) ^ 2;
t438 = pkin(3) * t537 - qJDD(3) * qJ(4) + qJD(3) * t595 + t515 * t488 - t443;
t563 = qJD(3) * t516;
t495 = t563 + (t530 * t533 - t556) * qJDD(1);
t460 = -qJD(5) * t501 - qJDD(3) * t532 + t495 * t535;
t469 = -mrSges(7,2) * t512 + mrSges(7,3) * t500;
t586 = -t460 * mrSges(7,1) - t500 * t469;
t470 = -mrSges(6,2) * t512 + mrSges(6,3) * t500;
t585 = -(mrSges(6,1) + mrSges(7,1)) * t460 - (t469 + t470) * t500;
t582 = mrSges(4,1) - mrSges(5,2);
t442 = t584 * t476 - t533 * t477;
t489 = mrSges(4,1) * t515 + mrSges(4,2) * t516;
t562 = t515 * qJD(3);
t496 = t546 * qJDD(1) - t562;
t508 = pkin(4) * t516 - qJD(3) * pkin(8);
t514 = t515 ^ 2;
t555 = t534 * g(1) - t536 * g(2);
t550 = qJDD(2) - t555;
t565 = -t530 ^ 2 - t527;
t494 = (-pkin(2) * t531 - pkin(1)) * qJDD(1) + (t565 * pkin(7) - qJ(2)) * t538 + t550;
t539 = pkin(3) * t563 + t516 * t595 + (-t496 + t562) * qJ(4) + t494;
t430 = -pkin(4) * t514 - t508 * t516 + (pkin(3) + pkin(8)) * t495 + t539;
t439 = -qJDD(3) * pkin(3) - t537 * qJ(4) + t516 * t488 + qJDD(4) - t442;
t433 = (t515 * t516 - qJDD(3)) * pkin(8) + (t496 + t562) * pkin(4) + t439;
t425 = -t532 * t430 + t535 * t433;
t461 = qJD(5) * t500 + qJDD(3) * t535 + t495 * t532;
t463 = -mrSges(7,1) * t500 + mrSges(7,2) * t501;
t464 = -mrSges(6,1) * t500 + mrSges(6,2) * t501;
t493 = qJDD(5) + t496;
t420 = -0.2e1 * qJD(6) * t501 + (t500 * t512 - t461) * qJ(6) + (t500 * t501 + t493) * pkin(5) + t425;
t559 = m(7) * t420 + t493 * mrSges(7,1) + t512 * t469;
t410 = m(6) * t425 + mrSges(6,1) * t493 + t470 * t512 + (-t463 - t464) * t501 + (-mrSges(6,3) - mrSges(7,3)) * t461 + t559;
t426 = t535 * t430 + t532 * t433;
t472 = mrSges(7,1) * t512 - mrSges(7,3) * t501;
t473 = mrSges(6,1) * t512 - mrSges(6,3) * t501;
t471 = pkin(5) * t512 - qJ(6) * t501;
t499 = t500 ^ 2;
t422 = -pkin(5) * t499 + qJ(6) * t460 + 0.2e1 * qJD(6) * t500 - t471 * t512 + t426;
t558 = m(7) * t422 + t460 * mrSges(7,3) + t500 * t463;
t412 = m(6) * t426 + mrSges(6,3) * t460 + t464 * t500 + (-t472 - t473) * t512 + (-mrSges(6,2) - mrSges(7,2)) * t493 + t558;
t409 = t410 * t535 + t412 * t532;
t490 = -mrSges(5,2) * t515 - mrSges(5,3) * t516;
t543 = -m(5) * t439 - t496 * mrSges(5,1) - t516 * t490 - t409;
t506 = mrSges(5,1) * t515 - qJD(3) * mrSges(5,3);
t566 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t515 - t506;
t405 = m(4) * t442 - mrSges(4,3) * t496 + t566 * qJD(3) + t582 * qJDD(3) - t489 * t516 + t543;
t505 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t516;
t507 = mrSges(5,1) * t516 + qJD(3) * mrSges(5,2);
t435 = -pkin(4) * t495 - pkin(8) * t514 + qJD(3) * t508 - t438;
t428 = -pkin(5) * t460 - qJ(6) * t499 + t471 * t501 + qJDD(6) + t435;
t557 = m(7) * t428 + t461 * mrSges(7,2) + t501 * t472;
t547 = -m(6) * t435 - t461 * mrSges(6,2) - t501 * t473 - t557;
t542 = -m(5) * t438 + qJDD(3) * mrSges(5,3) + qJD(3) * t507 - t547;
t415 = t542 - qJDD(3) * mrSges(4,2) + (-t489 - t490) * t515 + (-mrSges(4,3) - mrSges(5,1)) * t495 + m(4) * t443 - qJD(3) * t505 + t585;
t573 = t584 * t405 + t533 * t415;
t572 = -t590 * t500 - t592 * t501 - t588 * t512;
t571 = -t580 * t500 - t593 * t501 - t592 * t512;
t569 = -t589 * qJD(3) + t578 * t515 - t579 * t516;
t568 = t578 * qJD(3) + t591 * t515 - t577 * t516;
t567 = t579 * qJD(3) + t577 * t515 + t594 * t516;
t553 = -t533 * t405 + t584 * t415;
t552 = -t410 * t532 + t535 * t412;
t549 = -mrSges(3,1) * t531 + mrSges(3,2) * t530;
t548 = mrSges(3,3) * qJDD(1) + t538 * t549;
t437 = pkin(3) * t495 + t539;
t544 = m(5) * t437 - t496 * mrSges(5,3) - t516 * t507 + t552;
t417 = -mrSges(7,3) * t461 - t463 * t501 + t559;
t541 = mrSges(6,1) * t425 + mrSges(7,1) * t420 - mrSges(6,2) * t426 - mrSges(7,2) * t422 + pkin(5) * t417 + t590 * t460 + t592 * t461 + t588 * t493 + t571 * t500 + t587 * t501;
t540 = m(4) * t494 + mrSges(4,2) * t496 + t582 * t495 + t505 * t516 + t566 * t515 + t544;
t519 = (Ifges(3,5) * t530 + Ifges(3,6) * t531) * qJD(1);
t513 = -qJDD(1) * pkin(1) - t538 * qJ(2) + t550;
t497 = -t530 * t517 + t554;
t423 = t557 + t586;
t408 = mrSges(6,2) * t435 + mrSges(7,2) * t428 - mrSges(6,3) * t425 - mrSges(7,3) * t420 - qJ(6) * t417 + t580 * t460 + t593 * t461 + t592 * t493 - t572 * t500 - t587 * t512;
t407 = qJDD(3) * mrSges(5,2) + qJD(3) * t506 - t543;
t406 = -mrSges(5,2) * t495 - t506 * t515 + t544;
t403 = t565 * t538 * mrSges(3,3) + m(3) * t513 + t549 * qJDD(1) + t540;
t402 = -mrSges(6,1) * t435 + mrSges(6,3) * t426 - mrSges(7,1) * t428 + mrSges(7,3) * t422 - pkin(5) * t423 + qJ(6) * t558 + (-qJ(6) * t472 - t571) * t512 + t572 * t501 + (-mrSges(7,2) * qJ(6) + t590) * t493 + t580 * t461 + t596 * t460;
t401 = mrSges(5,1) * t439 + mrSges(4,2) * t494 - mrSges(4,3) * t442 - mrSges(5,3) * t437 + pkin(4) * t409 - qJ(4) * t406 - t568 * qJD(3) + t579 * qJDD(3) + t577 * t495 + t594 * t496 + t569 * t515 + t541;
t400 = -mrSges(4,1) * t494 + mrSges(4,3) * t443 - mrSges(5,1) * t438 + mrSges(5,2) * t437 - t532 * t408 - t535 * t402 - pkin(4) * (t547 - t585) - pkin(8) * t552 - pkin(3) * t406 + t569 * t516 - t577 * t496 + t591 * t495 + t578 * qJDD(3) + t567 * qJD(3);
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t555 - mrSges(2,2) * t551 + t530 * (t531 * qJD(1) * t519 + mrSges(3,2) * t513 - mrSges(3,3) * t497 + t584 * t401 - t533 * t400 - pkin(7) * t573 + (Ifges(3,1) * t530 + Ifges(3,4) * t531) * qJDD(1)) + t531 * (-t519 * t564 - mrSges(3,1) * t513 + mrSges(3,3) * t498 + t533 * t401 + t584 * t400 - pkin(2) * t540 + pkin(7) * t553 + (Ifges(3,4) * t530 + Ifges(3,2) * t531) * qJDD(1)) - pkin(1) * t403 + qJ(2) * ((m(3) * t498 + t548 * t531 + t553) * t531 + (-m(3) * t497 + t548 * t530 - t573) * t530); t403; mrSges(4,1) * t442 - mrSges(4,2) * t443 + mrSges(5,2) * t439 - mrSges(5,3) * t438 + t535 * t408 - t532 * t402 - pkin(8) * t409 - pkin(3) * t407 + qJ(4) * (-t460 * mrSges(6,1) - t500 * t470 + t542 + t586) + t568 * t516 + (-qJ(4) * t490 + t567) * t515 + t579 * t496 + (-mrSges(5,1) * qJ(4) - t578) * t495 + t589 * qJDD(3); t407; t541; t423;];
tauJ  = t1;
