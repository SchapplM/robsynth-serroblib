% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PRRRR7
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PRRRR7_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR7_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR7_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR7_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR7_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR7_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR7_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR7_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR7_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:12:10
% EndTime: 2019-12-05 17:12:16
% DurationCPUTime: 4.40s
% Computational Cost: add. (51310->260), mult. (106352->334), div. (0->0), fcn. (72242->10), ass. (0->105)
t515 = cos(pkin(9));
t491 = sin(pkin(9));
t479 = -t515 * g(1) - t491 * g(2);
t490 = -g(3) + qJDD(1);
t495 = sin(qJ(2));
t499 = cos(qJ(2));
t461 = t499 * t479 + t495 * t490;
t500 = qJD(2) ^ 2;
t456 = -t500 * pkin(2) + qJDD(2) * pkin(6) + t461;
t478 = t491 * g(1) - t515 * g(2);
t494 = sin(qJ(3));
t498 = cos(qJ(3));
t445 = -t494 * t456 - t498 * t478;
t511 = qJD(2) * qJD(3);
t510 = t498 * t511;
t476 = t494 * qJDD(2) + t510;
t437 = (-t476 + t510) * pkin(7) + (t494 * t498 * t500 + qJDD(3)) * pkin(3) + t445;
t446 = t498 * t456 - t494 * t478;
t477 = t498 * qJDD(2) - t494 * t511;
t513 = qJD(2) * t494;
t482 = qJD(3) * pkin(3) - pkin(7) * t513;
t489 = t498 ^ 2;
t438 = -t489 * t500 * pkin(3) + t477 * pkin(7) - qJD(3) * t482 + t446;
t493 = sin(qJ(4));
t497 = cos(qJ(4));
t422 = t497 * t437 - t493 * t438;
t466 = (-t493 * t494 + t497 * t498) * qJD(2);
t442 = t466 * qJD(4) + t497 * t476 + t493 * t477;
t467 = (t493 * t498 + t494 * t497) * qJD(2);
t487 = qJDD(3) + qJDD(4);
t488 = qJD(3) + qJD(4);
t418 = (t466 * t488 - t442) * pkin(8) + (t466 * t467 + t487) * pkin(4) + t422;
t423 = t493 * t437 + t497 * t438;
t441 = -t467 * qJD(4) - t493 * t476 + t497 * t477;
t459 = t488 * pkin(4) - t467 * pkin(8);
t462 = t466 ^ 2;
t419 = -t462 * pkin(4) + t441 * pkin(8) - t488 * t459 + t423;
t492 = sin(qJ(5));
t496 = cos(qJ(5));
t416 = t496 * t418 - t492 * t419;
t451 = t496 * t466 - t492 * t467;
t427 = t451 * qJD(5) + t492 * t441 + t496 * t442;
t452 = t492 * t466 + t496 * t467;
t433 = -t451 * mrSges(6,1) + t452 * mrSges(6,2);
t485 = qJD(5) + t488;
t443 = -t485 * mrSges(6,2) + t451 * mrSges(6,3);
t484 = qJDD(5) + t487;
t414 = m(6) * t416 + t484 * mrSges(6,1) - t427 * mrSges(6,3) - t452 * t433 + t485 * t443;
t417 = t492 * t418 + t496 * t419;
t426 = -t452 * qJD(5) + t496 * t441 - t492 * t442;
t444 = t485 * mrSges(6,1) - t452 * mrSges(6,3);
t415 = m(6) * t417 - t484 * mrSges(6,2) + t426 * mrSges(6,3) + t451 * t433 - t485 * t444;
t406 = t496 * t414 + t492 * t415;
t453 = -t466 * mrSges(5,1) + t467 * mrSges(5,2);
t457 = -t488 * mrSges(5,2) + t466 * mrSges(5,3);
t404 = m(5) * t422 + t487 * mrSges(5,1) - t442 * mrSges(5,3) - t467 * t453 + t488 * t457 + t406;
t458 = t488 * mrSges(5,1) - t467 * mrSges(5,3);
t505 = -t492 * t414 + t496 * t415;
t405 = m(5) * t423 - t487 * mrSges(5,2) + t441 * mrSges(5,3) + t466 * t453 - t488 * t458 + t505;
t400 = t497 * t404 + t493 * t405;
t475 = (-mrSges(4,1) * t498 + mrSges(4,2) * t494) * qJD(2);
t512 = qJD(2) * t498;
t481 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t512;
t398 = m(4) * t445 + qJDD(3) * mrSges(4,1) - t476 * mrSges(4,3) + qJD(3) * t481 - t475 * t513 + t400;
t480 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t513;
t506 = -t493 * t404 + t497 * t405;
t399 = m(4) * t446 - qJDD(3) * mrSges(4,2) + t477 * mrSges(4,3) - qJD(3) * t480 + t475 * t512 + t506;
t507 = -t494 * t398 + t498 * t399;
t391 = m(3) * t461 - t500 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t507;
t460 = -t495 * t479 + t499 * t490;
t503 = -qJDD(2) * pkin(2) - t460;
t455 = -t500 * pkin(6) + t503;
t439 = -t477 * pkin(3) + t482 * t513 + (-pkin(7) * t489 - pkin(6)) * t500 + t503;
t421 = -t441 * pkin(4) - t462 * pkin(8) + t467 * t459 + t439;
t504 = m(6) * t421 - t426 * mrSges(6,1) + t427 * mrSges(6,2) - t451 * t443 + t452 * t444;
t502 = m(5) * t439 - t441 * mrSges(5,1) + t442 * mrSges(5,2) - t466 * t457 + t467 * t458 + t504;
t501 = -m(4) * t455 + t477 * mrSges(4,1) - t476 * mrSges(4,2) - t480 * t513 + t481 * t512 - t502;
t410 = m(3) * t460 + qJDD(2) * mrSges(3,1) - t500 * mrSges(3,2) + t501;
t508 = t499 * t391 - t495 * t410;
t387 = m(2) * t479 + t508;
t394 = t498 * t398 + t494 * t399;
t393 = (m(2) + m(3)) * t478 - t394;
t514 = t491 * t387 + t515 * t393;
t388 = t495 * t391 + t499 * t410;
t509 = t515 * t387 - t491 * t393;
t465 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t494 + Ifges(4,4) * t498) * qJD(2);
t464 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t494 + Ifges(4,2) * t498) * qJD(2);
t463 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t494 + Ifges(4,6) * t498) * qJD(2);
t449 = Ifges(5,1) * t467 + Ifges(5,4) * t466 + Ifges(5,5) * t488;
t448 = Ifges(5,4) * t467 + Ifges(5,2) * t466 + Ifges(5,6) * t488;
t447 = Ifges(5,5) * t467 + Ifges(5,6) * t466 + Ifges(5,3) * t488;
t430 = Ifges(6,1) * t452 + Ifges(6,4) * t451 + Ifges(6,5) * t485;
t429 = Ifges(6,4) * t452 + Ifges(6,2) * t451 + Ifges(6,6) * t485;
t428 = Ifges(6,5) * t452 + Ifges(6,6) * t451 + Ifges(6,3) * t485;
t408 = mrSges(6,2) * t421 - mrSges(6,3) * t416 + Ifges(6,1) * t427 + Ifges(6,4) * t426 + Ifges(6,5) * t484 + t451 * t428 - t485 * t429;
t407 = -mrSges(6,1) * t421 + mrSges(6,3) * t417 + Ifges(6,4) * t427 + Ifges(6,2) * t426 + Ifges(6,6) * t484 - t452 * t428 + t485 * t430;
t396 = mrSges(5,2) * t439 - mrSges(5,3) * t422 + Ifges(5,1) * t442 + Ifges(5,4) * t441 + Ifges(5,5) * t487 - pkin(8) * t406 - t492 * t407 + t496 * t408 + t466 * t447 - t488 * t448;
t395 = -mrSges(5,1) * t439 + mrSges(5,3) * t423 + Ifges(5,4) * t442 + Ifges(5,2) * t441 + Ifges(5,6) * t487 - pkin(4) * t504 + pkin(8) * t505 + t496 * t407 + t492 * t408 - t467 * t447 + t488 * t449;
t384 = mrSges(4,2) * t455 - mrSges(4,3) * t445 + Ifges(4,1) * t476 + Ifges(4,4) * t477 + Ifges(4,5) * qJDD(3) - pkin(7) * t400 - qJD(3) * t464 - t493 * t395 + t497 * t396 + t463 * t512;
t383 = -mrSges(4,1) * t455 + mrSges(4,3) * t446 + Ifges(4,4) * t476 + Ifges(4,2) * t477 + Ifges(4,6) * qJDD(3) - pkin(3) * t502 + pkin(7) * t506 + qJD(3) * t465 + t497 * t395 + t493 * t396 - t463 * t513;
t382 = -Ifges(4,3) * qJDD(3) + Ifges(3,6) * qJDD(2) + (-t494 * t464 + t498 * t465) * qJD(2) + t500 * Ifges(3,5) - Ifges(6,3) * t484 - Ifges(5,3) * t487 - Ifges(4,5) * t476 - Ifges(4,6) * t477 + mrSges(3,1) * t478 - pkin(2) * t394 - pkin(3) * t400 - pkin(4) * t406 - mrSges(6,1) * t416 + mrSges(6,2) * t417 - mrSges(5,1) * t422 + mrSges(5,2) * t423 - Ifges(6,6) * t426 - Ifges(6,5) * t427 - Ifges(5,6) * t441 - Ifges(5,5) * t442 - mrSges(4,1) * t445 + mrSges(4,2) * t446 + t451 * t430 - t452 * t429 + mrSges(3,3) * t461 + t466 * t449 - t467 * t448;
t381 = -mrSges(3,2) * t478 - mrSges(3,3) * t460 + Ifges(3,5) * qJDD(2) - t500 * Ifges(3,6) - pkin(6) * t394 - t494 * t383 + t498 * t384;
t380 = -mrSges(2,1) * t490 - mrSges(3,1) * t460 + mrSges(3,2) * t461 + mrSges(2,3) * t479 - Ifges(3,3) * qJDD(2) - pkin(1) * t388 - pkin(2) * t501 - pkin(6) * t507 - t498 * t383 - t494 * t384;
t379 = mrSges(2,2) * t490 - mrSges(2,3) * t478 - pkin(5) * t388 + t499 * t381 - t495 * t382;
t1 = [-m(1) * g(1) + t509; -m(1) * g(2) + t514; -m(1) * g(3) + m(2) * t490 + t388; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t514 + t515 * t379 - t491 * t380; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t509 + t491 * t379 + t515 * t380; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + mrSges(2,1) * t478 - mrSges(2,2) * t479 + t495 * t381 + t499 * t382 + pkin(1) * (m(3) * t478 - t394) + pkin(5) * t508;];
tauB = t1;
