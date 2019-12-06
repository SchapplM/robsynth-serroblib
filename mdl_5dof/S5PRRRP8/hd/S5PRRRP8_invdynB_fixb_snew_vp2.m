% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PRRRP8
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
% Datum: 2019-12-05 17:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PRRRP8_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP8_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP8_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP8_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP8_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP8_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP8_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP8_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP8_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:58:36
% EndTime: 2019-12-05 16:58:42
% DurationCPUTime: 3.40s
% Computational Cost: add. (34857->243), mult. (65485->305), div. (0->0), fcn. (42429->10), ass. (0->105)
t539 = Ifges(5,1) + Ifges(6,1);
t533 = Ifges(5,4) - Ifges(6,5);
t538 = -Ifges(5,5) - Ifges(6,4);
t537 = Ifges(5,2) + Ifges(6,3);
t531 = Ifges(5,6) - Ifges(6,6);
t536 = -Ifges(5,3) - Ifges(6,2);
t498 = sin(pkin(9));
t500 = cos(pkin(9));
t488 = t498 * g(1) - t500 * g(2);
t489 = -t500 * g(1) - t498 * g(2);
t497 = -g(3) + qJDD(1);
t499 = sin(pkin(5));
t501 = cos(pkin(5));
t504 = sin(qJ(2));
t506 = cos(qJ(2));
t445 = -t504 * t489 + (t488 * t501 + t497 * t499) * t506;
t535 = cos(qJ(4));
t534 = -mrSges(5,3) - mrSges(6,2);
t508 = qJD(2) ^ 2;
t528 = t501 * t504;
t529 = t499 * t504;
t446 = t488 * t528 + t506 * t489 + t497 * t529;
t444 = -t508 * pkin(2) + qJDD(2) * pkin(7) + t446;
t469 = -t499 * t488 + t501 * t497;
t503 = sin(qJ(3));
t505 = cos(qJ(3));
t440 = t505 * t444 + t503 * t469;
t485 = (-pkin(3) * t505 - pkin(8) * t503) * qJD(2);
t507 = qJD(3) ^ 2;
t521 = t505 * qJD(2);
t436 = -t507 * pkin(3) + qJDD(3) * pkin(8) + t485 * t521 + t440;
t443 = -qJDD(2) * pkin(2) - t508 * pkin(7) - t445;
t520 = qJD(2) * qJD(3);
t517 = t505 * t520;
t486 = t503 * qJDD(2) + t517;
t518 = t503 * t520;
t487 = t505 * qJDD(2) - t518;
t438 = (-t486 - t517) * pkin(8) + (-t487 + t518) * pkin(3) + t443;
t502 = sin(qJ(4));
t433 = t535 * t436 + t502 * t438;
t522 = qJD(2) * t503;
t483 = t502 * qJD(3) + t535 * t522;
t457 = t483 * qJD(4) - t535 * qJDD(3) + t502 * t486;
t494 = qJD(4) - t521;
t466 = t494 * mrSges(5,1) - t483 * mrSges(5,3);
t479 = qJDD(4) - t487;
t482 = -t535 * qJD(3) + t502 * t522;
t461 = t482 * pkin(4) - t483 * qJ(5);
t493 = t494 ^ 2;
t429 = -t493 * pkin(4) + t479 * qJ(5) + 0.2e1 * qJD(5) * t494 - t482 * t461 + t433;
t467 = -t494 * mrSges(6,1) + t483 * mrSges(6,2);
t519 = m(6) * t429 + t479 * mrSges(6,3) + t494 * t467;
t462 = t482 * mrSges(6,1) - t483 * mrSges(6,3);
t523 = -t482 * mrSges(5,1) - t483 * mrSges(5,2) - t462;
t425 = m(5) * t433 - t479 * mrSges(5,2) + t534 * t457 - t494 * t466 + t523 * t482 + t519;
t432 = -t502 * t436 + t535 * t438;
t458 = -t482 * qJD(4) + t502 * qJDD(3) + t535 * t486;
t465 = -t494 * mrSges(5,2) - t482 * mrSges(5,3);
t430 = -t479 * pkin(4) - t493 * qJ(5) + t483 * t461 + qJDD(5) - t432;
t468 = -t482 * mrSges(6,2) + t494 * mrSges(6,3);
t513 = -m(6) * t430 + t479 * mrSges(6,1) + t494 * t468;
t426 = m(5) * t432 + t479 * mrSges(5,1) + t534 * t458 + t494 * t465 + t523 * t483 + t513;
t421 = t502 * t425 + t535 * t426;
t490 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t522;
t491 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t521;
t510 = -m(4) * t443 + t487 * mrSges(4,1) - t486 * mrSges(4,2) - t490 * t522 + t491 * t521 - t421;
t415 = m(3) * t445 + qJDD(2) * mrSges(3,1) - t508 * mrSges(3,2) + t510;
t530 = t415 * t506;
t484 = (-mrSges(4,1) * t505 + mrSges(4,2) * t503) * qJD(2);
t514 = t535 * t425 - t502 * t426;
t418 = m(4) * t440 - qJDD(3) * mrSges(4,2) + t487 * mrSges(4,3) - qJD(3) * t490 + t484 * t521 + t514;
t439 = -t503 * t444 + t505 * t469;
t435 = -qJDD(3) * pkin(3) - t507 * pkin(8) + t485 * t522 - t439;
t431 = -0.2e1 * qJD(5) * t483 + (t482 * t494 - t458) * qJ(5) + (t483 * t494 + t457) * pkin(4) + t435;
t427 = m(6) * t431 + t457 * mrSges(6,1) - t458 * mrSges(6,3) - t483 * t467 + t482 * t468;
t509 = -m(5) * t435 - t457 * mrSges(5,1) - t458 * mrSges(5,2) - t482 * t465 - t483 * t466 - t427;
t423 = m(4) * t439 + qJDD(3) * mrSges(4,1) - t486 * mrSges(4,3) + qJD(3) * t491 - t484 * t522 + t509;
t515 = t505 * t418 - t503 * t423;
t409 = m(3) * t446 - t508 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t515;
t412 = t503 * t418 + t505 * t423;
t411 = m(3) * t469 + t412;
t399 = t409 * t528 - t499 * t411 + t501 * t530;
t397 = m(2) * t488 + t399;
t404 = t506 * t409 - t504 * t415;
t403 = m(2) * t489 + t404;
t527 = t500 * t397 + t498 * t403;
t526 = t537 * t482 - t533 * t483 - t531 * t494;
t525 = t531 * t482 + t538 * t483 + t536 * t494;
t524 = -t533 * t482 + t539 * t483 - t538 * t494;
t398 = t409 * t529 + t501 * t411 + t499 * t530;
t516 = -t498 * t397 + t500 * t403;
t419 = -mrSges(5,1) * t435 - mrSges(6,1) * t431 + mrSges(6,2) * t429 + mrSges(5,3) * t433 - pkin(4) * t427 - t537 * t457 + t533 * t458 + t531 * t479 + t525 * t483 + t524 * t494;
t420 = mrSges(5,2) * t435 + mrSges(6,2) * t430 - mrSges(5,3) * t432 - mrSges(6,3) * t431 - qJ(5) * t427 - t533 * t457 + t539 * t458 - t479 * t538 + t525 * t482 + t526 * t494;
t472 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t503 + Ifges(4,6) * t505) * qJD(2);
t473 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t503 + Ifges(4,2) * t505) * qJD(2);
t400 = mrSges(4,2) * t443 - mrSges(4,3) * t439 + Ifges(4,1) * t486 + Ifges(4,4) * t487 + Ifges(4,5) * qJDD(3) - pkin(8) * t421 - qJD(3) * t473 - t502 * t419 + t535 * t420 + t472 * t521;
t474 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t503 + Ifges(4,4) * t505) * qJD(2);
t405 = Ifges(4,4) * t486 + Ifges(4,2) * t487 + Ifges(4,6) * qJDD(3) - t472 * t522 + qJD(3) * t474 - mrSges(4,1) * t443 + mrSges(4,3) * t440 - mrSges(5,1) * t432 + mrSges(5,2) * t433 + mrSges(6,1) * t430 - mrSges(6,3) * t429 - pkin(4) * t513 - qJ(5) * t519 - pkin(3) * t421 + (pkin(4) * t462 + t526) * t483 + (qJ(5) * t462 - t524) * t482 + t536 * t479 + (pkin(4) * mrSges(6,2) + t538) * t458 + (qJ(5) * mrSges(6,2) + t531) * t457;
t394 = mrSges(3,2) * t469 - mrSges(3,3) * t445 + Ifges(3,5) * qJDD(2) - t508 * Ifges(3,6) - pkin(7) * t412 + t505 * t400 - t503 * t405;
t395 = Ifges(3,6) * qJDD(2) + t508 * Ifges(3,5) - mrSges(3,1) * t469 + mrSges(3,3) * t446 - Ifges(4,5) * t486 - Ifges(4,6) * t487 - Ifges(4,3) * qJDD(3) - mrSges(4,1) * t439 + mrSges(4,2) * t440 - t502 * t420 - t535 * t419 - pkin(3) * t509 - pkin(8) * t514 - pkin(2) * t412 + (-t503 * t473 + t505 * t474) * qJD(2);
t511 = pkin(6) * t404 + t394 * t504 + t395 * t506;
t393 = mrSges(3,1) * t445 - mrSges(3,2) * t446 + Ifges(3,3) * qJDD(2) + pkin(2) * t510 + pkin(7) * t515 + t503 * t400 + t505 * t405;
t392 = mrSges(2,2) * t497 - mrSges(2,3) * t488 + t506 * t394 - t504 * t395 + (-t398 * t499 - t399 * t501) * pkin(6);
t391 = -mrSges(2,1) * t497 + mrSges(2,3) * t489 - pkin(1) * t398 - t499 * t393 + t511 * t501;
t1 = [-m(1) * g(1) + t516; -m(1) * g(2) + t527; -m(1) * g(3) + m(2) * t497 + t398; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t527 - t498 * t391 + t500 * t392; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t516 + t500 * t391 + t498 * t392; -mrSges(1,1) * g(2) + mrSges(2,1) * t488 + mrSges(1,2) * g(1) - mrSges(2,2) * t489 + pkin(1) * t399 + t501 * t393 + t511 * t499;];
tauB = t1;
