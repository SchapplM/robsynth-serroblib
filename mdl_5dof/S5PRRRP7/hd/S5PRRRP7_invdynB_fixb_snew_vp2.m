% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PRRRP7
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PRRRP7_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP7_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP7_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP7_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP7_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP7_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP7_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP7_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP7_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:54:27
% EndTime: 2019-12-05 16:54:33
% DurationCPUTime: 3.47s
% Computational Cost: add. (35669->245), mult. (67578->305), div. (0->0), fcn. (43962->10), ass. (0->106)
t541 = Ifges(5,1) + Ifges(6,1);
t536 = Ifges(5,4) + Ifges(6,4);
t535 = Ifges(5,5) + Ifges(6,5);
t540 = Ifges(5,2) + Ifges(6,2);
t539 = Ifges(5,6) + Ifges(6,6);
t538 = Ifges(5,3) + Ifges(6,3);
t499 = sin(pkin(9));
t501 = cos(pkin(9));
t492 = t499 * g(1) - t501 * g(2);
t493 = -t501 * g(1) - t499 * g(2);
t498 = -g(3) + qJDD(1);
t500 = sin(pkin(5));
t502 = cos(pkin(5));
t505 = sin(qJ(2));
t508 = cos(qJ(2));
t449 = -t505 * t493 + (t492 * t502 + t498 * t500) * t508;
t537 = -mrSges(5,2) - mrSges(6,2);
t510 = qJD(2) ^ 2;
t531 = t502 * t505;
t532 = t500 * t505;
t450 = t492 * t531 + t508 * t493 + t498 * t532;
t447 = -t510 * pkin(2) + qJDD(2) * pkin(7) + t450;
t473 = -t500 * t492 + t502 * t498;
t504 = sin(qJ(3));
t507 = cos(qJ(3));
t443 = t507 * t447 + t504 * t473;
t489 = (-pkin(3) * t507 - pkin(8) * t504) * qJD(2);
t509 = qJD(3) ^ 2;
t524 = t507 * qJD(2);
t438 = -t509 * pkin(3) + qJDD(3) * pkin(8) + t489 * t524 + t443;
t446 = -qJDD(2) * pkin(2) - t510 * pkin(7) - t449;
t523 = qJD(2) * qJD(3);
t519 = t507 * t523;
t490 = t504 * qJDD(2) + t519;
t520 = t504 * t523;
t491 = t507 * qJDD(2) - t520;
t441 = (-t490 - t519) * pkin(8) + (-t491 + t520) * pkin(3) + t446;
t503 = sin(qJ(4));
t506 = cos(qJ(4));
t433 = -t503 * t438 + t506 * t441;
t525 = qJD(2) * t504;
t486 = t506 * qJD(3) - t503 * t525;
t463 = t486 * qJD(4) + t503 * qJDD(3) + t506 * t490;
t487 = t503 * qJD(3) + t506 * t525;
t465 = -t486 * mrSges(6,1) + t487 * mrSges(6,2);
t466 = -t486 * mrSges(5,1) + t487 * mrSges(5,2);
t497 = qJD(4) - t524;
t469 = -t497 * mrSges(5,2) + t486 * mrSges(5,3);
t483 = qJDD(4) - t491;
t430 = -0.2e1 * qJD(5) * t487 + (t486 * t497 - t463) * qJ(5) + (t486 * t487 + t483) * pkin(4) + t433;
t468 = -t497 * mrSges(6,2) + t486 * mrSges(6,3);
t522 = m(6) * t430 + t483 * mrSges(6,1) + t497 * t468;
t423 = m(5) * t433 + t483 * mrSges(5,1) + t497 * t469 + (-t465 - t466) * t487 + (-mrSges(5,3) - mrSges(6,3)) * t463 + t522;
t434 = t506 * t438 + t503 * t441;
t462 = -t487 * qJD(4) + t506 * qJDD(3) - t503 * t490;
t470 = t497 * pkin(4) - t487 * qJ(5);
t482 = t486 ^ 2;
t432 = -t482 * pkin(4) + t462 * qJ(5) + 0.2e1 * qJD(5) * t486 - t497 * t470 + t434;
t521 = m(6) * t432 + t462 * mrSges(6,3) + t486 * t465;
t471 = t497 * mrSges(6,1) - t487 * mrSges(6,3);
t526 = -t497 * mrSges(5,1) + t487 * mrSges(5,3) - t471;
t425 = m(5) * t434 + t462 * mrSges(5,3) + t486 * t466 + t537 * t483 + t526 * t497 + t521;
t422 = t506 * t423 + t503 * t425;
t494 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t525;
t495 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t524;
t512 = -m(4) * t446 + t491 * mrSges(4,1) - t490 * mrSges(4,2) - t494 * t525 + t495 * t524 - t422;
t417 = m(3) * t449 + qJDD(2) * mrSges(3,1) - t510 * mrSges(3,2) + t512;
t533 = t417 * t508;
t488 = (-mrSges(4,1) * t507 + mrSges(4,2) * t504) * qJD(2);
t516 = -t503 * t423 + t506 * t425;
t420 = m(4) * t443 - qJDD(3) * mrSges(4,2) + t491 * mrSges(4,3) - qJD(3) * t494 + t488 * t524 + t516;
t442 = -t504 * t447 + t507 * t473;
t437 = -qJDD(3) * pkin(3) - t509 * pkin(8) + t489 * t525 - t442;
t435 = -t462 * pkin(4) - t482 * qJ(5) + t487 * t470 + qJDD(5) + t437;
t515 = m(6) * t435 - t462 * mrSges(6,1) - t486 * t468;
t511 = -m(5) * t437 + t462 * mrSges(5,1) + t537 * t463 + t486 * t469 + t526 * t487 - t515;
t427 = m(4) * t442 + qJDD(3) * mrSges(4,1) - t490 * mrSges(4,3) + qJD(3) * t495 - t488 * t525 + t511;
t517 = t507 * t420 - t504 * t427;
t410 = m(3) * t450 - t510 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t517;
t413 = t504 * t420 + t507 * t427;
t412 = m(3) * t473 + t413;
t400 = t410 * t531 - t500 * t412 + t502 * t533;
t398 = m(2) * t492 + t400;
t405 = t508 * t410 - t505 * t417;
t404 = m(2) * t493 + t405;
t530 = t501 * t398 + t499 * t404;
t529 = t539 * t486 + t535 * t487 + t538 * t497;
t528 = -t540 * t486 - t536 * t487 - t539 * t497;
t527 = t536 * t486 + t541 * t487 + t535 * t497;
t399 = t410 * t532 + t502 * t412 + t500 * t533;
t518 = -t499 * t398 + t501 * t404;
t414 = -mrSges(5,1) * t437 + mrSges(5,3) * t434 - mrSges(6,1) * t435 + mrSges(6,3) * t432 - pkin(4) * t515 + qJ(5) * t521 + (-qJ(5) * t471 + t527) * t497 + (-pkin(4) * t471 - t529) * t487 + (-qJ(5) * mrSges(6,2) + t539) * t483 + (-pkin(4) * mrSges(6,2) + t536) * t463 + t540 * t462;
t428 = -t463 * mrSges(6,3) - t487 * t465 + t522;
t421 = mrSges(5,2) * t437 + mrSges(6,2) * t435 - mrSges(5,3) * t433 - mrSges(6,3) * t430 - qJ(5) * t428 + t536 * t462 + t541 * t463 + t535 * t483 + t529 * t486 + t528 * t497;
t476 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t504 + Ifges(4,6) * t507) * qJD(2);
t477 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t504 + Ifges(4,2) * t507) * qJD(2);
t401 = mrSges(4,2) * t446 - mrSges(4,3) * t442 + Ifges(4,1) * t490 + Ifges(4,4) * t491 + Ifges(4,5) * qJDD(3) - pkin(8) * t422 - qJD(3) * t477 - t503 * t414 + t506 * t421 + t476 * t524;
t478 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t504 + Ifges(4,4) * t507) * qJD(2);
t406 = -t476 * t525 - mrSges(4,1) * t446 - mrSges(5,1) * t433 - mrSges(6,1) * t430 + mrSges(5,2) * t434 + mrSges(6,2) * t432 + mrSges(4,3) * t443 + Ifges(4,4) * t490 + Ifges(4,2) * t491 + Ifges(4,6) * qJDD(3) - pkin(3) * t422 - pkin(4) * t428 + qJD(3) * t478 + t528 * t487 + t527 * t486 - t538 * t483 - t535 * t463 - t539 * t462;
t395 = mrSges(3,2) * t473 - mrSges(3,3) * t449 + Ifges(3,5) * qJDD(2) - t510 * Ifges(3,6) - pkin(7) * t413 + t507 * t401 - t504 * t406;
t396 = Ifges(3,6) * qJDD(2) + t510 * Ifges(3,5) - mrSges(3,1) * t473 + mrSges(3,3) * t450 - Ifges(4,5) * t490 - Ifges(4,6) * t491 - Ifges(4,3) * qJDD(3) - mrSges(4,1) * t442 + mrSges(4,2) * t443 - t503 * t421 - t506 * t414 - pkin(3) * t511 - pkin(8) * t516 - pkin(2) * t413 + (-t504 * t477 + t507 * t478) * qJD(2);
t513 = pkin(6) * t405 + t395 * t505 + t396 * t508;
t394 = mrSges(3,1) * t449 - mrSges(3,2) * t450 + Ifges(3,3) * qJDD(2) + pkin(2) * t512 + pkin(7) * t517 + t504 * t401 + t507 * t406;
t393 = mrSges(2,2) * t498 - mrSges(2,3) * t492 + t508 * t395 - t505 * t396 + (-t399 * t500 - t400 * t502) * pkin(6);
t392 = -mrSges(2,1) * t498 + mrSges(2,3) * t493 - pkin(1) * t399 - t500 * t394 + t513 * t502;
t1 = [-m(1) * g(1) + t518; -m(1) * g(2) + t530; -m(1) * g(3) + m(2) * t498 + t399; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t530 - t499 * t392 + t501 * t393; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t518 + t501 * t392 + t499 * t393; -mrSges(1,1) * g(2) + mrSges(2,1) * t492 + mrSges(1,2) * g(1) - mrSges(2,2) * t493 + pkin(1) * t400 + t502 * t394 + t513 * t500;];
tauB = t1;
