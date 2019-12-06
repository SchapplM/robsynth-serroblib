% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRRPR3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-05 18:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRRPR3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR3_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR3_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR3_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR3_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR3_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:42:38
% EndTime: 2019-12-05 18:42:43
% DurationCPUTime: 5.06s
% Computational Cost: add. (86594->272), mult. (115951->347), div. (0->0), fcn. (73622->10), ass. (0->108)
t514 = qJD(1) + qJD(2);
t510 = t514 ^ 2;
t541 = pkin(3) * t510;
t519 = sin(qJ(3));
t540 = t514 * t519;
t523 = cos(qJ(3));
t539 = t514 * t523;
t521 = sin(qJ(1));
t525 = cos(qJ(1));
t506 = t525 * g(2) + t521 * g(3);
t500 = qJDD(1) * pkin(1) + t506;
t505 = t521 * g(2) - t525 * g(3);
t526 = qJD(1) ^ 2;
t501 = -t526 * pkin(1) + t505;
t520 = sin(qJ(2));
t524 = cos(qJ(2));
t479 = t520 * t500 + t524 * t501;
t512 = qJDD(1) + qJDD(2);
t474 = -t510 * pkin(2) + t512 * pkin(7) + t479;
t538 = t519 * t474;
t537 = qJD(3) * t514;
t495 = t519 * t512 + t523 * t537;
t458 = qJDD(3) * pkin(3) - t495 * qJ(4) - t538 + (qJ(4) * t537 + t519 * t541 - g(1)) * t523;
t464 = -t519 * g(1) + t523 * t474;
t496 = t523 * t512 - t519 * t537;
t502 = qJD(3) * pkin(3) - qJ(4) * t540;
t515 = t523 ^ 2;
t459 = t496 * qJ(4) - qJD(3) * t502 - t515 * t541 + t464;
t516 = sin(pkin(9));
t517 = cos(pkin(9));
t487 = (t516 * t523 + t517 * t519) * t514;
t441 = -0.2e1 * qJD(4) * t487 + t517 * t458 - t516 * t459;
t477 = t517 * t495 + t516 * t496;
t486 = (-t516 * t519 + t517 * t523) * t514;
t439 = (qJD(3) * t486 - t477) * pkin(8) + (t486 * t487 + qJDD(3)) * pkin(4) + t441;
t442 = 0.2e1 * qJD(4) * t486 + t516 * t458 + t517 * t459;
t476 = -t516 * t495 + t517 * t496;
t482 = qJD(3) * pkin(4) - t487 * pkin(8);
t485 = t486 ^ 2;
t440 = -t485 * pkin(4) + t476 * pkin(8) - qJD(3) * t482 + t442;
t518 = sin(qJ(5));
t522 = cos(qJ(5));
t437 = t522 * t439 - t518 * t440;
t468 = t522 * t486 - t518 * t487;
t448 = t468 * qJD(5) + t518 * t476 + t522 * t477;
t469 = t518 * t486 + t522 * t487;
t454 = -t468 * mrSges(6,1) + t469 * mrSges(6,2);
t513 = qJD(3) + qJD(5);
t461 = -t513 * mrSges(6,2) + t468 * mrSges(6,3);
t511 = qJDD(3) + qJDD(5);
t435 = m(6) * t437 + t511 * mrSges(6,1) - t448 * mrSges(6,3) - t469 * t454 + t513 * t461;
t438 = t518 * t439 + t522 * t440;
t447 = -t469 * qJD(5) + t522 * t476 - t518 * t477;
t462 = t513 * mrSges(6,1) - t469 * mrSges(6,3);
t436 = m(6) * t438 - t511 * mrSges(6,2) + t447 * mrSges(6,3) + t468 * t454 - t513 * t462;
t427 = t522 * t435 + t518 * t436;
t472 = -t486 * mrSges(5,1) + t487 * mrSges(5,2);
t480 = -qJD(3) * mrSges(5,2) + t486 * mrSges(5,3);
t425 = m(5) * t441 + qJDD(3) * mrSges(5,1) - t477 * mrSges(5,3) + qJD(3) * t480 - t487 * t472 + t427;
t481 = qJD(3) * mrSges(5,1) - t487 * mrSges(5,3);
t532 = -t518 * t435 + t522 * t436;
t426 = m(5) * t442 - qJDD(3) * mrSges(5,2) + t476 * mrSges(5,3) - qJD(3) * t481 + t486 * t472 + t532;
t421 = t517 * t425 + t516 * t426;
t463 = -t523 * g(1) - t538;
t494 = (-mrSges(4,1) * t523 + mrSges(4,2) * t519) * t514;
t504 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t539;
t419 = m(4) * t463 + qJDD(3) * mrSges(4,1) - t495 * mrSges(4,3) + qJD(3) * t504 - t494 * t540 + t421;
t503 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t540;
t533 = -t516 * t425 + t517 * t426;
t420 = m(4) * t464 - qJDD(3) * mrSges(4,2) + t496 * mrSges(4,3) - qJD(3) * t503 + t494 * t539 + t533;
t534 = -t519 * t419 + t523 * t420;
t412 = m(3) * t479 - t510 * mrSges(3,1) - t512 * mrSges(3,2) + t534;
t478 = t524 * t500 - t520 * t501;
t529 = -t512 * pkin(2) - t478;
t473 = -t510 * pkin(7) + t529;
t460 = -t496 * pkin(3) + qJDD(4) + t502 * t540 + (-qJ(4) * t515 - pkin(7)) * t510 + t529;
t444 = -t476 * pkin(4) - t485 * pkin(8) + t487 * t482 + t460;
t531 = m(6) * t444 - t447 * mrSges(6,1) + t448 * mrSges(6,2) - t468 * t461 + t469 * t462;
t528 = m(5) * t460 - t476 * mrSges(5,1) + t477 * mrSges(5,2) - t486 * t480 + t487 * t481 + t531;
t527 = -m(4) * t473 + t496 * mrSges(4,1) - t495 * mrSges(4,2) - t503 * t540 + t504 * t539 - t528;
t431 = m(3) * t478 + t512 * mrSges(3,1) - t510 * mrSges(3,2) + t527;
t409 = t520 * t412 + t524 * t431;
t413 = t523 * t419 + t519 * t420;
t535 = t524 * t412 - t520 * t431;
t407 = m(2) * t505 - t526 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t535;
t408 = m(2) * t506 + qJDD(1) * mrSges(2,1) - t526 * mrSges(2,2) + t409;
t536 = t525 * t407 - t521 * t408;
t530 = -t521 * t407 - t525 * t408;
t490 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t519 + Ifges(4,4) * t523) * t514;
t489 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t519 + Ifges(4,2) * t523) * t514;
t488 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t519 + Ifges(4,6) * t523) * t514;
t467 = Ifges(5,1) * t487 + Ifges(5,4) * t486 + Ifges(5,5) * qJD(3);
t466 = Ifges(5,4) * t487 + Ifges(5,2) * t486 + Ifges(5,6) * qJD(3);
t465 = Ifges(5,5) * t487 + Ifges(5,6) * t486 + Ifges(5,3) * qJD(3);
t451 = Ifges(6,1) * t469 + Ifges(6,4) * t468 + Ifges(6,5) * t513;
t450 = Ifges(6,4) * t469 + Ifges(6,2) * t468 + Ifges(6,6) * t513;
t449 = Ifges(6,5) * t469 + Ifges(6,6) * t468 + Ifges(6,3) * t513;
t429 = mrSges(6,2) * t444 - mrSges(6,3) * t437 + Ifges(6,1) * t448 + Ifges(6,4) * t447 + Ifges(6,5) * t511 + t468 * t449 - t513 * t450;
t428 = -mrSges(6,1) * t444 + mrSges(6,3) * t438 + Ifges(6,4) * t448 + Ifges(6,2) * t447 + Ifges(6,6) * t511 - t469 * t449 + t513 * t451;
t415 = mrSges(5,2) * t460 - mrSges(5,3) * t441 + Ifges(5,1) * t477 + Ifges(5,4) * t476 + Ifges(5,5) * qJDD(3) - pkin(8) * t427 - qJD(3) * t466 - t518 * t428 + t522 * t429 + t486 * t465;
t414 = -mrSges(5,1) * t460 + mrSges(5,3) * t442 + Ifges(5,4) * t477 + Ifges(5,2) * t476 + Ifges(5,6) * qJDD(3) - pkin(4) * t531 + pkin(8) * t532 + qJD(3) * t467 + t522 * t428 + t518 * t429 - t487 * t465;
t405 = mrSges(4,2) * t473 - mrSges(4,3) * t463 + Ifges(4,1) * t495 + Ifges(4,4) * t496 + Ifges(4,5) * qJDD(3) - qJ(4) * t421 - qJD(3) * t489 - t516 * t414 + t517 * t415 + t488 * t539;
t404 = -mrSges(4,1) * t473 + mrSges(4,3) * t464 + Ifges(4,4) * t495 + Ifges(4,2) * t496 + Ifges(4,6) * qJDD(3) - pkin(3) * t528 + qJ(4) * t533 + qJD(3) * t490 + t517 * t414 + t516 * t415 - t488 * t540;
t403 = mrSges(3,1) * g(1) + (-t519 * t489 + t523 * t490) * t514 + (-Ifges(4,3) - Ifges(5,3)) * qJDD(3) + t510 * Ifges(3,5) - Ifges(6,3) * t511 + Ifges(3,6) * t512 - Ifges(4,6) * t496 + t486 * t467 - t487 * t466 - Ifges(4,5) * t495 - Ifges(5,5) * t477 + mrSges(3,3) * t479 - pkin(2) * t413 - pkin(3) * t421 - pkin(4) * t427 - mrSges(6,1) * t437 + mrSges(6,2) * t438 - mrSges(5,1) * t441 + mrSges(5,2) * t442 - Ifges(6,6) * t447 - Ifges(6,5) * t448 - mrSges(4,1) * t463 + mrSges(4,2) * t464 + t468 * t451 - t469 * t450 - Ifges(5,6) * t476;
t402 = -mrSges(3,2) * g(1) - mrSges(3,3) * t478 + Ifges(3,5) * t512 - t510 * Ifges(3,6) - pkin(7) * t413 - t519 * t404 + t523 * t405;
t401 = -mrSges(2,2) * g(1) - mrSges(2,3) * t506 + Ifges(2,5) * qJDD(1) - t526 * Ifges(2,6) - pkin(6) * t409 + t524 * t402 - t520 * t403;
t400 = Ifges(2,6) * qJDD(1) + t526 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t505 + t520 * t402 + t524 * t403 - pkin(1) * (-m(3) * g(1) + t413) + pkin(6) * t535;
t1 = [(-m(1) - m(2) - m(3)) * g(1) + t413; -m(1) * g(2) + t530; -m(1) * g(3) + t536; mrSges(2,1) * t506 + mrSges(3,1) * t478 - mrSges(1,2) * g(3) - mrSges(2,2) * t505 - mrSges(3,2) * t479 + mrSges(1,3) * g(2) + Ifges(2,3) * qJDD(1) + Ifges(3,3) * t512 + pkin(1) * t409 + pkin(2) * t527 + pkin(7) * t534 + t523 * t404 + t519 * t405; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t536 - t525 * t400 - t521 * t401; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t530 - t521 * t400 + t525 * t401;];
tauB = t1;
