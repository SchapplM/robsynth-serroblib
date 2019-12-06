% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PRRPR6
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:33
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PRRPR6_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR6_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR6_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR6_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR6_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR6_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR6_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR6_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR6_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:30:44
% EndTime: 2019-12-05 16:30:55
% DurationCPUTime: 6.67s
% Computational Cost: add. (75776->266), mult. (154779->348), div. (0->0), fcn. (106192->12), ass. (0->113)
t499 = sin(pkin(9));
t502 = cos(pkin(9));
t490 = t499 * g(1) - t502 * g(2);
t491 = -t502 * g(1) - t499 * g(2);
t497 = -g(3) + qJDD(1);
t500 = sin(pkin(5));
t503 = cos(pkin(5));
t506 = sin(qJ(2));
t509 = cos(qJ(2));
t454 = -t506 * t491 + (t490 * t503 + t497 * t500) * t509;
t511 = qJD(2) ^ 2;
t526 = t503 * t506;
t527 = t500 * t506;
t455 = t490 * t526 + t509 * t491 + t497 * t527;
t451 = -t511 * pkin(2) + qJDD(2) * pkin(7) + t455;
t470 = -t500 * t490 + t503 * t497;
t505 = sin(qJ(3));
t508 = cos(qJ(3));
t446 = t508 * t451 + t505 * t470;
t486 = (-pkin(3) * t508 - qJ(4) * t505) * qJD(2);
t510 = qJD(3) ^ 2;
t523 = t508 * qJD(2);
t434 = -t510 * pkin(3) + qJDD(3) * qJ(4) + t486 * t523 + t446;
t450 = -qJDD(2) * pkin(2) - t511 * pkin(7) - t454;
t522 = qJD(2) * qJD(3);
t521 = t508 * t522;
t488 = t505 * qJDD(2) + t521;
t496 = t505 * t522;
t489 = t508 * qJDD(2) - t496;
t438 = (-t488 - t521) * qJ(4) + (-t489 + t496) * pkin(3) + t450;
t498 = sin(pkin(10));
t501 = cos(pkin(10));
t524 = qJD(2) * t505;
t481 = t498 * qJD(3) + t501 * t524;
t429 = -0.2e1 * qJD(4) * t481 - t498 * t434 + t501 * t438;
t468 = t498 * qJDD(3) + t501 * t488;
t480 = t501 * qJD(3) - t498 * t524;
t427 = (-t480 * t523 - t468) * pkin(8) + (t480 * t481 - t489) * pkin(4) + t429;
t430 = 0.2e1 * qJD(4) * t480 + t501 * t434 + t498 * t438;
t467 = t501 * qJDD(3) - t498 * t488;
t469 = -pkin(4) * t523 - t481 * pkin(8);
t479 = t480 ^ 2;
t428 = -t479 * pkin(4) + t467 * pkin(8) + t469 * t523 + t430;
t504 = sin(qJ(5));
t507 = cos(qJ(5));
t425 = t507 * t427 - t504 * t428;
t460 = t507 * t480 - t504 * t481;
t440 = t460 * qJD(5) + t504 * t467 + t507 * t468;
t461 = t504 * t480 + t507 * t481;
t447 = -t460 * mrSges(6,1) + t461 * mrSges(6,2);
t495 = qJD(5) - t523;
t452 = -t495 * mrSges(6,2) + t460 * mrSges(6,3);
t483 = qJDD(5) - t489;
t423 = m(6) * t425 + t483 * mrSges(6,1) - t440 * mrSges(6,3) - t461 * t447 + t495 * t452;
t426 = t504 * t427 + t507 * t428;
t439 = -t461 * qJD(5) + t507 * t467 - t504 * t468;
t453 = t495 * mrSges(6,1) - t461 * mrSges(6,3);
t424 = m(6) * t426 - t483 * mrSges(6,2) + t439 * mrSges(6,3) + t460 * t447 - t495 * t453;
t415 = t507 * t423 + t504 * t424;
t462 = -t480 * mrSges(5,1) + t481 * mrSges(5,2);
t465 = mrSges(5,2) * t523 + t480 * mrSges(5,3);
t413 = m(5) * t429 - t489 * mrSges(5,1) - t468 * mrSges(5,3) - t481 * t462 - t465 * t523 + t415;
t466 = -mrSges(5,1) * t523 - t481 * mrSges(5,3);
t517 = -t504 * t423 + t507 * t424;
t414 = m(5) * t430 + t489 * mrSges(5,2) + t467 * mrSges(5,3) + t480 * t462 + t466 * t523 + t517;
t411 = t501 * t413 + t498 * t414;
t492 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t524;
t493 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t523;
t513 = -m(4) * t450 + t489 * mrSges(4,1) - t488 * mrSges(4,2) - t492 * t524 + t493 * t523 - t411;
t407 = m(3) * t454 + qJDD(2) * mrSges(3,1) - t511 * mrSges(3,2) + t513;
t528 = t407 * t509;
t487 = (-mrSges(4,1) * t508 + mrSges(4,2) * t505) * qJD(2);
t518 = -t498 * t413 + t501 * t414;
t410 = m(4) * t446 - qJDD(3) * mrSges(4,2) + t489 * mrSges(4,3) - qJD(3) * t492 + t487 * t523 + t518;
t445 = -t505 * t451 + t508 * t470;
t433 = -qJDD(3) * pkin(3) - t510 * qJ(4) + t486 * t524 + qJDD(4) - t445;
t431 = -t467 * pkin(4) - t479 * pkin(8) + t481 * t469 + t433;
t514 = m(6) * t431 - t439 * mrSges(6,1) + t440 * mrSges(6,2) - t460 * t452 + t461 * t453;
t512 = -m(5) * t433 + t467 * mrSges(5,1) - t468 * mrSges(5,2) + t480 * t465 - t481 * t466 - t514;
t419 = m(4) * t445 + qJDD(3) * mrSges(4,1) - t488 * mrSges(4,3) + qJD(3) * t493 - t487 * t524 + t512;
t519 = t508 * t410 - t505 * t419;
t399 = m(3) * t455 - t511 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t519;
t402 = t505 * t410 + t508 * t419;
t401 = m(3) * t470 + t402;
t389 = t399 * t526 - t500 * t401 + t503 * t528;
t387 = m(2) * t490 + t389;
t394 = t509 * t399 - t506 * t407;
t393 = m(2) * t491 + t394;
t525 = t502 * t387 + t499 * t393;
t388 = t399 * t527 + t503 * t401 + t500 * t528;
t520 = -t499 * t387 + t502 * t393;
t441 = Ifges(6,5) * t461 + Ifges(6,6) * t460 + Ifges(6,3) * t495;
t443 = Ifges(6,1) * t461 + Ifges(6,4) * t460 + Ifges(6,5) * t495;
t416 = -mrSges(6,1) * t431 + mrSges(6,3) * t426 + Ifges(6,4) * t440 + Ifges(6,2) * t439 + Ifges(6,6) * t483 - t461 * t441 + t495 * t443;
t442 = Ifges(6,4) * t461 + Ifges(6,2) * t460 + Ifges(6,6) * t495;
t417 = mrSges(6,2) * t431 - mrSges(6,3) * t425 + Ifges(6,1) * t440 + Ifges(6,4) * t439 + Ifges(6,5) * t483 + t460 * t441 - t495 * t442;
t456 = Ifges(5,5) * t481 + Ifges(5,6) * t480 - Ifges(5,3) * t523;
t458 = Ifges(5,1) * t481 + Ifges(5,4) * t480 - Ifges(5,5) * t523;
t403 = -mrSges(5,1) * t433 + mrSges(5,3) * t430 + Ifges(5,4) * t468 + Ifges(5,2) * t467 - Ifges(5,6) * t489 - pkin(4) * t514 + pkin(8) * t517 + t507 * t416 + t504 * t417 - t481 * t456 - t458 * t523;
t457 = Ifges(5,4) * t481 + Ifges(5,2) * t480 - Ifges(5,6) * t523;
t404 = mrSges(5,2) * t433 - mrSges(5,3) * t429 + Ifges(5,1) * t468 + Ifges(5,4) * t467 - Ifges(5,5) * t489 - pkin(8) * t415 - t504 * t416 + t507 * t417 + t480 * t456 + t457 * t523;
t475 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t505 + Ifges(4,6) * t508) * qJD(2);
t476 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t505 + Ifges(4,2) * t508) * qJD(2);
t390 = mrSges(4,2) * t450 - mrSges(4,3) * t445 + Ifges(4,1) * t488 + Ifges(4,4) * t489 + Ifges(4,5) * qJDD(3) - qJ(4) * t411 - qJD(3) * t476 - t498 * t403 + t501 * t404 + t475 * t523;
t477 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t505 + Ifges(4,4) * t508) * qJD(2);
t395 = Ifges(4,4) * t488 + Ifges(4,6) * qJDD(3) - t475 * t524 + qJD(3) * t477 - mrSges(4,1) * t450 + mrSges(4,3) * t446 - Ifges(5,5) * t468 - Ifges(5,6) * t467 - t481 * t457 + t480 * t458 - mrSges(5,1) * t429 + mrSges(5,2) * t430 - Ifges(6,5) * t440 - Ifges(6,6) * t439 - Ifges(6,3) * t483 - t461 * t442 + t460 * t443 - mrSges(6,1) * t425 + mrSges(6,2) * t426 - pkin(4) * t415 - pkin(3) * t411 + (Ifges(4,2) + Ifges(5,3)) * t489;
t384 = mrSges(3,2) * t470 - mrSges(3,3) * t454 + Ifges(3,5) * qJDD(2) - t511 * Ifges(3,6) - pkin(7) * t402 + t508 * t390 - t505 * t395;
t385 = Ifges(3,6) * qJDD(2) + t511 * Ifges(3,5) - mrSges(3,1) * t470 + mrSges(3,3) * t455 - Ifges(4,5) * t488 - Ifges(4,6) * t489 - Ifges(4,3) * qJDD(3) - mrSges(4,1) * t445 + mrSges(4,2) * t446 - t498 * t404 - t501 * t403 - pkin(3) * t512 - qJ(4) * t518 - pkin(2) * t402 + (-t505 * t476 + t508 * t477) * qJD(2);
t515 = pkin(6) * t394 + t384 * t506 + t385 * t509;
t383 = mrSges(3,1) * t454 - mrSges(3,2) * t455 + Ifges(3,3) * qJDD(2) + pkin(2) * t513 + pkin(7) * t519 + t505 * t390 + t508 * t395;
t382 = mrSges(2,2) * t497 - mrSges(2,3) * t490 + t509 * t384 - t506 * t385 + (-t388 * t500 - t389 * t503) * pkin(6);
t381 = -mrSges(2,1) * t497 + mrSges(2,3) * t491 - pkin(1) * t388 - t500 * t383 + t503 * t515;
t1 = [-m(1) * g(1) + t520; -m(1) * g(2) + t525; -m(1) * g(3) + m(2) * t497 + t388; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t525 - t499 * t381 + t502 * t382; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t520 + t502 * t381 + t499 * t382; -mrSges(1,1) * g(2) + mrSges(2,1) * t490 + mrSges(1,2) * g(1) - mrSges(2,2) * t491 + pkin(1) * t389 + t503 * t383 + t515 * t500;];
tauB = t1;
