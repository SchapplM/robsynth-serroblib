% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PPRRR3
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
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PPRRR3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR3_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR3_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRR3_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR3_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR3_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR3_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRR3_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:16:36
% EndTime: 2019-12-05 15:16:39
% DurationCPUTime: 2.14s
% Computational Cost: add. (22983->207), mult. (41706->269), div. (0->0), fcn. (27328->10), ass. (0->93)
t490 = sin(pkin(8));
t492 = cos(pkin(8));
t478 = -t492 * g(1) - t490 * g(2);
t488 = -g(3) + qJDD(1);
t489 = sin(pkin(9));
t491 = cos(pkin(9));
t464 = t491 * t478 + t489 * t488;
t477 = t490 * g(1) - t492 * g(2);
t476 = qJDD(2) - t477;
t495 = sin(qJ(3));
t498 = cos(qJ(3));
t454 = t498 * t464 + t495 * t476;
t499 = qJD(3) ^ 2;
t449 = -t499 * pkin(3) + qJDD(3) * pkin(6) + t454;
t463 = t489 * t478 - t491 * t488;
t494 = sin(qJ(4));
t497 = cos(qJ(4));
t439 = -t494 * t449 + t497 * t463;
t512 = qJD(3) * qJD(4);
t510 = t497 * t512;
t474 = t494 * qJDD(3) + t510;
t436 = (-t474 + t510) * pkin(7) + (t494 * t497 * t499 + qJDD(4)) * pkin(4) + t439;
t440 = t497 * t449 + t494 * t463;
t475 = t497 * qJDD(3) - t494 * t512;
t514 = qJD(3) * t494;
t481 = qJD(4) * pkin(4) - pkin(7) * t514;
t487 = t497 ^ 2;
t437 = -t487 * t499 * pkin(4) + t475 * pkin(7) - qJD(4) * t481 + t440;
t493 = sin(qJ(5));
t496 = cos(qJ(5));
t434 = t496 * t436 - t493 * t437;
t468 = (-t493 * t494 + t496 * t497) * qJD(3);
t447 = t468 * qJD(5) + t496 * t474 + t493 * t475;
t469 = (t493 * t497 + t494 * t496) * qJD(3);
t456 = -t468 * mrSges(6,1) + t469 * mrSges(6,2);
t486 = qJD(4) + qJD(5);
t461 = -t486 * mrSges(6,2) + t468 * mrSges(6,3);
t485 = qJDD(4) + qJDD(5);
t431 = m(6) * t434 + t485 * mrSges(6,1) - t447 * mrSges(6,3) - t469 * t456 + t486 * t461;
t435 = t493 * t436 + t496 * t437;
t446 = -t469 * qJD(5) - t493 * t474 + t496 * t475;
t462 = t486 * mrSges(6,1) - t469 * mrSges(6,3);
t432 = m(6) * t435 - t485 * mrSges(6,2) + t446 * mrSges(6,3) + t468 * t456 - t486 * t462;
t423 = t496 * t431 + t493 * t432;
t466 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t494 + Ifges(5,2) * t497) * qJD(3);
t467 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t494 + Ifges(5,4) * t497) * qJD(3);
t451 = Ifges(6,4) * t469 + Ifges(6,2) * t468 + Ifges(6,6) * t486;
t452 = Ifges(6,1) * t469 + Ifges(6,4) * t468 + Ifges(6,5) * t486;
t502 = -mrSges(6,1) * t434 + mrSges(6,2) * t435 - Ifges(6,5) * t447 - Ifges(6,6) * t446 - Ifges(6,3) * t485 - t469 * t451 + t468 * t452;
t516 = mrSges(5,1) * t439 - mrSges(5,2) * t440 + Ifges(5,5) * t474 + Ifges(5,6) * t475 + Ifges(5,3) * qJDD(4) + pkin(4) * t423 + (t494 * t466 - t497 * t467) * qJD(3) - t502;
t473 = (-mrSges(5,1) * t497 + mrSges(5,2) * t494) * qJD(3);
t513 = qJD(3) * t497;
t480 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t513;
t421 = m(5) * t439 + qJDD(4) * mrSges(5,1) - t474 * mrSges(5,3) + qJD(4) * t480 - t473 * t514 + t423;
t479 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t514;
t506 = -t493 * t431 + t496 * t432;
t422 = m(5) * t440 - qJDD(4) * mrSges(5,2) + t475 * mrSges(5,3) - qJD(4) * t479 + t473 * t513 + t506;
t419 = -t494 * t421 + t497 * t422;
t416 = m(4) * t454 - t499 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t419;
t453 = -t495 * t464 + t498 * t476;
t504 = -qJDD(3) * pkin(3) - t453;
t448 = -t499 * pkin(6) + t504;
t438 = t481 * t514 - t475 * pkin(4) + (-pkin(7) * t487 - pkin(6)) * t499 + t504;
t503 = m(6) * t438 - t446 * mrSges(6,1) + t447 * mrSges(6,2) - t468 * t461 + t469 * t462;
t427 = -m(5) * t448 + t475 * mrSges(5,1) - t474 * mrSges(5,2) - t479 * t514 + t480 * t513 - t503;
t426 = m(4) * t453 + qJDD(3) * mrSges(4,1) - t499 * mrSges(4,2) + t427;
t507 = t498 * t416 - t495 * t426;
t410 = m(3) * t464 + t507;
t418 = t497 * t421 + t494 * t422;
t417 = (-m(3) - m(4)) * t463 - t418;
t508 = t491 * t410 - t489 * t417;
t402 = m(2) * t478 + t508;
t412 = t495 * t416 + t498 * t426;
t411 = m(3) * t476 + t412;
t409 = m(2) * t477 - t411;
t515 = t490 * t402 + t492 * t409;
t403 = t489 * t410 + t491 * t417;
t511 = m(2) * t488 + t403;
t509 = t492 * t402 - t490 * t409;
t450 = Ifges(6,5) * t469 + Ifges(6,6) * t468 + Ifges(6,3) * t486;
t424 = -mrSges(6,1) * t438 + mrSges(6,3) * t435 + Ifges(6,4) * t447 + Ifges(6,2) * t446 + Ifges(6,6) * t485 - t469 * t450 + t486 * t452;
t425 = mrSges(6,2) * t438 - mrSges(6,3) * t434 + Ifges(6,1) * t447 + Ifges(6,4) * t446 + Ifges(6,5) * t485 + t468 * t450 - t486 * t451;
t465 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t494 + Ifges(5,6) * t497) * qJD(3);
t405 = -mrSges(5,1) * t448 + mrSges(5,3) * t440 + Ifges(5,4) * t474 + Ifges(5,2) * t475 + Ifges(5,6) * qJDD(4) - pkin(4) * t503 + pkin(7) * t506 + qJD(4) * t467 + t496 * t424 + t493 * t425 - t465 * t514;
t413 = mrSges(5,2) * t448 - mrSges(5,3) * t439 + Ifges(5,1) * t474 + Ifges(5,4) * t475 + Ifges(5,5) * qJDD(4) - pkin(7) * t423 - qJD(4) * t466 - t493 * t424 + t496 * t425 + t465 * t513;
t501 = mrSges(4,1) * t453 - mrSges(4,2) * t454 + Ifges(4,3) * qJDD(3) + pkin(3) * t427 + pkin(6) * t419 + t497 * t405 + t494 * t413;
t404 = -mrSges(4,1) * t463 + mrSges(4,3) * t454 + t499 * Ifges(4,5) + Ifges(4,6) * qJDD(3) - pkin(3) * t418 - t516;
t399 = mrSges(4,2) * t463 - mrSges(4,3) * t453 + Ifges(4,5) * qJDD(3) - t499 * Ifges(4,6) - pkin(6) * t418 - t494 * t405 + t497 * t413;
t398 = -mrSges(3,1) * t476 + mrSges(3,3) * t464 - pkin(2) * t412 - t501;
t397 = mrSges(3,2) * t476 + mrSges(3,3) * t463 - pkin(5) * t412 + t498 * t399 - t495 * t404;
t396 = -mrSges(2,1) * t488 + mrSges(2,3) * t478 + mrSges(3,1) * t463 + mrSges(3,2) * t464 - t495 * t399 - t498 * t404 - pkin(2) * (-m(4) * t463 - t418) - pkin(5) * t507 - pkin(1) * t403;
t395 = mrSges(2,2) * t488 - mrSges(2,3) * t477 - qJ(2) * t403 + t491 * t397 - t489 * t398;
t1 = [-m(1) * g(1) + t509; -m(1) * g(2) + t515; -m(1) * g(3) + t511; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t515 + t492 * t395 - t490 * t396; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t509 + t490 * t395 + t492 * t396; -mrSges(1,1) * g(2) + mrSges(2,1) * t477 + mrSges(1,2) * g(1) - mrSges(2,2) * t478 - pkin(1) * t411 + qJ(2) * t508 + t489 * t397 + t491 * t398; t511; t411; t501; t516; -t502;];
tauJB = t1;
