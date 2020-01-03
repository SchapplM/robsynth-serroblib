% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPRPR7
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPRPR7_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR7_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR7_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR7_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR7_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR7_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR7_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR7_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR7_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:19:10
% EndTime: 2019-12-31 18:19:14
% DurationCPUTime: 3.68s
% Computational Cost: add. (40636->270), mult. (86093->344), div. (0->0), fcn. (53523->10), ass. (0->108)
t488 = sin(pkin(9));
t490 = cos(pkin(9));
t493 = sin(qJ(3));
t496 = cos(qJ(3));
t461 = (t488 * t493 - t490 * t496) * qJD(1);
t517 = 2 * qJD(4);
t494 = sin(qJ(1));
t497 = cos(qJ(1));
t480 = t494 * g(1) - t497 * g(2);
t472 = qJDD(1) * pkin(1) + t480;
t481 = -t497 * g(1) - t494 * g(2);
t499 = qJD(1) ^ 2;
t474 = -t499 * pkin(1) + t481;
t489 = sin(pkin(8));
t491 = cos(pkin(8));
t452 = t489 * t472 + t491 * t474;
t445 = -t499 * pkin(2) + qJDD(1) * pkin(6) + t452;
t487 = -g(3) + qJDD(2);
t436 = -t493 * t445 + t496 * t487;
t513 = qJD(1) * qJD(3);
t511 = t496 * t513;
t475 = t493 * qJDD(1) + t511;
t425 = (-t475 + t511) * qJ(4) + (t493 * t496 * t499 + qJDD(3)) * pkin(3) + t436;
t437 = t496 * t445 + t493 * t487;
t476 = t496 * qJDD(1) - t493 * t513;
t515 = qJD(1) * t493;
t477 = qJD(3) * pkin(3) - qJ(4) * t515;
t486 = t496 ^ 2;
t426 = -t486 * t499 * pkin(3) + t476 * qJ(4) - qJD(3) * t477 + t437;
t421 = t488 * t425 + t490 * t426 - t461 * t517;
t462 = (t488 * t496 + t490 * t493) * qJD(1);
t447 = t461 * mrSges(5,1) + t462 * mrSges(5,2);
t453 = -t488 * t475 + t490 * t476;
t458 = qJD(3) * mrSges(5,1) - t462 * mrSges(5,3);
t448 = t461 * pkin(4) - t462 * pkin(7);
t498 = qJD(3) ^ 2;
t419 = -t498 * pkin(4) + qJDD(3) * pkin(7) - t461 * t448 + t421;
t451 = t491 * t472 - t489 * t474;
t503 = -qJDD(1) * pkin(2) - t451;
t427 = -t476 * pkin(3) + qJDD(4) + t477 * t515 + (-qJ(4) * t486 - pkin(6)) * t499 + t503;
t454 = t490 * t475 + t488 * t476;
t422 = (qJD(3) * t461 - t454) * pkin(7) + (qJD(3) * t462 - t453) * pkin(4) + t427;
t492 = sin(qJ(5));
t495 = cos(qJ(5));
t416 = -t492 * t419 + t495 * t422;
t455 = t495 * qJD(3) - t492 * t462;
t434 = t455 * qJD(5) + t492 * qJDD(3) + t495 * t454;
t456 = t492 * qJD(3) + t495 * t462;
t435 = -t455 * mrSges(6,1) + t456 * mrSges(6,2);
t460 = qJD(5) + t461;
t438 = -t460 * mrSges(6,2) + t455 * mrSges(6,3);
t450 = qJDD(5) - t453;
t414 = m(6) * t416 + t450 * mrSges(6,1) - t434 * mrSges(6,3) - t456 * t435 + t460 * t438;
t417 = t495 * t419 + t492 * t422;
t433 = -t456 * qJD(5) + t495 * qJDD(3) - t492 * t454;
t439 = t460 * mrSges(6,1) - t456 * mrSges(6,3);
t415 = m(6) * t417 - t450 * mrSges(6,2) + t433 * mrSges(6,3) + t455 * t435 - t460 * t439;
t506 = -t492 * t414 + t495 * t415;
t405 = m(5) * t421 - qJDD(3) * mrSges(5,2) + t453 * mrSges(5,3) - qJD(3) * t458 - t461 * t447 + t506;
t505 = -t490 * t425 + t488 * t426;
t420 = -0.2e1 * qJD(4) * t462 - t505;
t457 = -qJD(3) * mrSges(5,2) - t461 * mrSges(5,3);
t418 = -qJDD(3) * pkin(4) - t498 * pkin(7) + (t517 + t448) * t462 + t505;
t502 = -m(6) * t418 + t433 * mrSges(6,1) - t434 * mrSges(6,2) + t455 * t438 - t456 * t439;
t410 = m(5) * t420 + qJDD(3) * mrSges(5,1) - t454 * mrSges(5,3) + qJD(3) * t457 - t462 * t447 + t502;
t400 = t488 * t405 + t490 * t410;
t473 = (-mrSges(4,1) * t496 + mrSges(4,2) * t493) * qJD(1);
t514 = qJD(1) * t496;
t479 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t514;
t398 = m(4) * t436 + qJDD(3) * mrSges(4,1) - t475 * mrSges(4,3) + qJD(3) * t479 - t473 * t515 + t400;
t478 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t515;
t507 = t490 * t405 - t488 * t410;
t399 = m(4) * t437 - qJDD(3) * mrSges(4,2) + t476 * mrSges(4,3) - qJD(3) * t478 + t473 * t514 + t507;
t508 = -t493 * t398 + t496 * t399;
t391 = m(3) * t452 - t499 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t508;
t444 = -t499 * pkin(6) + t503;
t406 = t495 * t414 + t492 * t415;
t501 = m(5) * t427 - t453 * mrSges(5,1) + t454 * mrSges(5,2) + t461 * t457 + t462 * t458 + t406;
t500 = -m(4) * t444 + t476 * mrSges(4,1) - t475 * mrSges(4,2) - t478 * t515 + t479 * t514 - t501;
t402 = m(3) * t451 + qJDD(1) * mrSges(3,1) - t499 * mrSges(3,2) + t500;
t388 = t489 * t391 + t491 * t402;
t386 = m(2) * t480 + qJDD(1) * mrSges(2,1) - t499 * mrSges(2,2) + t388;
t509 = t491 * t391 - t489 * t402;
t387 = m(2) * t481 - t499 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t509;
t516 = t497 * t386 + t494 * t387;
t392 = t496 * t398 + t493 * t399;
t512 = m(3) * t487 + t392;
t510 = -t494 * t386 + t497 * t387;
t468 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t493 + Ifges(4,4) * t496) * qJD(1);
t467 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t493 + Ifges(4,2) * t496) * qJD(1);
t466 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t493 + Ifges(4,6) * t496) * qJD(1);
t443 = Ifges(5,1) * t462 - Ifges(5,4) * t461 + Ifges(5,5) * qJD(3);
t442 = Ifges(5,4) * t462 - Ifges(5,2) * t461 + Ifges(5,6) * qJD(3);
t441 = Ifges(5,5) * t462 - Ifges(5,6) * t461 + Ifges(5,3) * qJD(3);
t430 = Ifges(6,1) * t456 + Ifges(6,4) * t455 + Ifges(6,5) * t460;
t429 = Ifges(6,4) * t456 + Ifges(6,2) * t455 + Ifges(6,6) * t460;
t428 = Ifges(6,5) * t456 + Ifges(6,6) * t455 + Ifges(6,3) * t460;
t408 = mrSges(6,2) * t418 - mrSges(6,3) * t416 + Ifges(6,1) * t434 + Ifges(6,4) * t433 + Ifges(6,5) * t450 + t455 * t428 - t460 * t429;
t407 = -mrSges(6,1) * t418 + mrSges(6,3) * t417 + Ifges(6,4) * t434 + Ifges(6,2) * t433 + Ifges(6,6) * t450 - t456 * t428 + t460 * t430;
t394 = -mrSges(5,1) * t427 - mrSges(6,1) * t416 + mrSges(6,2) * t417 + mrSges(5,3) * t421 + Ifges(5,4) * t454 - Ifges(6,5) * t434 + Ifges(5,2) * t453 + Ifges(5,6) * qJDD(3) - Ifges(6,6) * t433 - Ifges(6,3) * t450 - pkin(4) * t406 + qJD(3) * t443 - t456 * t429 + t455 * t430 - t462 * t441;
t393 = mrSges(5,2) * t427 - mrSges(5,3) * t420 + Ifges(5,1) * t454 + Ifges(5,4) * t453 + Ifges(5,5) * qJDD(3) - pkin(7) * t406 - qJD(3) * t442 - t492 * t407 + t495 * t408 - t461 * t441;
t382 = mrSges(4,2) * t444 - mrSges(4,3) * t436 + Ifges(4,1) * t475 + Ifges(4,4) * t476 + Ifges(4,5) * qJDD(3) - qJ(4) * t400 - qJD(3) * t467 + t490 * t393 - t488 * t394 + t466 * t514;
t381 = -mrSges(4,1) * t444 + mrSges(4,3) * t437 + Ifges(4,4) * t475 + Ifges(4,2) * t476 + Ifges(4,6) * qJDD(3) - pkin(3) * t501 + qJ(4) * t507 + qJD(3) * t468 + t488 * t393 + t490 * t394 - t466 * t515;
t380 = -pkin(2) * t392 - mrSges(3,1) * t487 + mrSges(3,3) * t452 - pkin(3) * t400 - Ifges(4,5) * t475 - Ifges(4,6) * t476 - mrSges(4,1) * t436 + mrSges(4,2) * t437 - pkin(7) * t506 - Ifges(5,5) * t454 - Ifges(5,6) * t453 - mrSges(5,1) * t420 + mrSges(5,2) * t421 - t492 * t408 - t495 * t407 - pkin(4) * t502 + t499 * Ifges(3,5) - t462 * t442 - t461 * t443 + Ifges(3,6) * qJDD(1) + (-Ifges(4,3) - Ifges(5,3)) * qJDD(3) + (-t493 * t467 + t496 * t468) * qJD(1);
t379 = mrSges(3,2) * t487 - mrSges(3,3) * t451 + Ifges(3,5) * qJDD(1) - t499 * Ifges(3,6) - pkin(6) * t392 - t493 * t381 + t496 * t382;
t378 = -mrSges(2,2) * g(3) - mrSges(2,3) * t480 + Ifges(2,5) * qJDD(1) - t499 * Ifges(2,6) - qJ(2) * t388 + t491 * t379 - t489 * t380;
t377 = mrSges(2,1) * g(3) + mrSges(2,3) * t481 + t499 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t512 + qJ(2) * t509 + t489 * t379 + t491 * t380;
t1 = [-m(1) * g(1) + t510; -m(1) * g(2) + t516; (-m(1) - m(2)) * g(3) + t512; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t516 - t494 * t377 + t497 * t378; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t510 + t497 * t377 + t494 * t378; pkin(1) * t388 + mrSges(2,1) * t480 - mrSges(2,2) * t481 + t493 * t382 + t496 * t381 + pkin(2) * t500 + pkin(6) * t508 + mrSges(3,1) * t451 - mrSges(3,2) * t452 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);];
tauB = t1;
