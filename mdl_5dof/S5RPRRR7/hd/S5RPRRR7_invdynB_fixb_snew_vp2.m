% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPRRR7
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPRRR7_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR7_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR7_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR7_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR7_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR7_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR7_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR7_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR7_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:03:32
% EndTime: 2019-12-31 19:03:37
% DurationCPUTime: 4.06s
% Computational Cost: add. (50279->269), mult. (96518->338), div. (0->0), fcn. (60567->10), ass. (0->108)
t485 = sin(qJ(1));
t489 = cos(qJ(1));
t471 = t485 * g(1) - t489 * g(2);
t463 = qJDD(1) * pkin(1) + t471;
t472 = -t489 * g(1) - t485 * g(2);
t491 = qJD(1) ^ 2;
t465 = -t491 * pkin(1) + t472;
t480 = sin(pkin(9));
t481 = cos(pkin(9));
t444 = t480 * t463 + t481 * t465;
t433 = -t491 * pkin(2) + qJDD(1) * pkin(6) + t444;
t479 = -g(3) + qJDD(2);
t484 = sin(qJ(3));
t488 = cos(qJ(3));
t427 = t488 * t433 + t484 * t479;
t464 = (-mrSges(4,1) * t488 + mrSges(4,2) * t484) * qJD(1);
t502 = qJD(1) * qJD(3);
t476 = t484 * t502;
t468 = t488 * qJDD(1) - t476;
t504 = qJD(1) * t484;
t469 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t504;
t443 = t481 * t463 - t480 * t465;
t432 = -qJDD(1) * pkin(2) - t491 * pkin(6) - t443;
t500 = t488 * t502;
t467 = t484 * qJDD(1) + t500;
t420 = (-t467 - t500) * pkin(7) + (-t468 + t476) * pkin(3) + t432;
t466 = (-pkin(3) * t488 - pkin(7) * t484) * qJD(1);
t490 = qJD(3) ^ 2;
t503 = t488 * qJD(1);
t424 = -t490 * pkin(3) + qJDD(3) * pkin(7) + t466 * t503 + t427;
t483 = sin(qJ(4));
t487 = cos(qJ(4));
t410 = t487 * t420 - t483 * t424;
t461 = t487 * qJD(3) - t483 * t504;
t440 = t461 * qJD(4) + t483 * qJDD(3) + t487 * t467;
t460 = qJDD(4) - t468;
t462 = t483 * qJD(3) + t487 * t504;
t474 = qJD(4) - t503;
t407 = (t461 * t474 - t440) * pkin(8) + (t461 * t462 + t460) * pkin(4) + t410;
t411 = t483 * t420 + t487 * t424;
t439 = -t462 * qJD(4) + t487 * qJDD(3) - t483 * t467;
t448 = t474 * pkin(4) - t462 * pkin(8);
t459 = t461 ^ 2;
t408 = -t459 * pkin(4) + t439 * pkin(8) - t474 * t448 + t411;
t482 = sin(qJ(5));
t486 = cos(qJ(5));
t405 = t486 * t407 - t482 * t408;
t441 = t486 * t461 - t482 * t462;
t414 = t441 * qJD(5) + t482 * t439 + t486 * t440;
t442 = t482 * t461 + t486 * t462;
t425 = -t441 * mrSges(6,1) + t442 * mrSges(6,2);
t473 = qJD(5) + t474;
t428 = -t473 * mrSges(6,2) + t441 * mrSges(6,3);
t456 = qJDD(5) + t460;
t403 = m(6) * t405 + t456 * mrSges(6,1) - t414 * mrSges(6,3) - t442 * t425 + t473 * t428;
t406 = t482 * t407 + t486 * t408;
t413 = -t442 * qJD(5) + t486 * t439 - t482 * t440;
t429 = t473 * mrSges(6,1) - t442 * mrSges(6,3);
t404 = m(6) * t406 - t456 * mrSges(6,2) + t413 * mrSges(6,3) + t441 * t425 - t473 * t429;
t395 = t486 * t403 + t482 * t404;
t445 = -t461 * mrSges(5,1) + t462 * mrSges(5,2);
t446 = -t474 * mrSges(5,2) + t461 * mrSges(5,3);
t393 = m(5) * t410 + t460 * mrSges(5,1) - t440 * mrSges(5,3) - t462 * t445 + t474 * t446 + t395;
t447 = t474 * mrSges(5,1) - t462 * mrSges(5,3);
t495 = -t482 * t403 + t486 * t404;
t394 = m(5) * t411 - t460 * mrSges(5,2) + t439 * mrSges(5,3) + t461 * t445 - t474 * t447 + t495;
t496 = -t483 * t393 + t487 * t394;
t390 = m(4) * t427 - qJDD(3) * mrSges(4,2) + t468 * mrSges(4,3) - qJD(3) * t469 + t464 * t503 + t496;
t426 = -t484 * t433 + t488 * t479;
t470 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t503;
t423 = -qJDD(3) * pkin(3) - t490 * pkin(7) + t466 * t504 - t426;
t409 = -t439 * pkin(4) - t459 * pkin(8) + t462 * t448 + t423;
t494 = m(6) * t409 - t413 * mrSges(6,1) + t414 * mrSges(6,2) - t441 * t428 + t442 * t429;
t492 = -m(5) * t423 + t439 * mrSges(5,1) - t440 * mrSges(5,2) + t461 * t446 - t462 * t447 - t494;
t399 = m(4) * t426 + qJDD(3) * mrSges(4,1) - t467 * mrSges(4,3) + qJD(3) * t470 - t464 * t504 + t492;
t497 = t488 * t390 - t484 * t399;
t383 = m(3) * t444 - t491 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t497;
t391 = t487 * t393 + t483 * t394;
t493 = -m(4) * t432 + t468 * mrSges(4,1) - t467 * mrSges(4,2) - t469 * t504 + t470 * t503 - t391;
t387 = m(3) * t443 + qJDD(1) * mrSges(3,1) - t491 * mrSges(3,2) + t493;
t378 = t480 * t383 + t481 * t387;
t376 = m(2) * t471 + qJDD(1) * mrSges(2,1) - t491 * mrSges(2,2) + t378;
t498 = t481 * t383 - t480 * t387;
t377 = m(2) * t472 - t491 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t498;
t505 = t489 * t376 + t485 * t377;
t384 = t484 * t390 + t488 * t399;
t501 = m(3) * t479 + t384;
t499 = -t485 * t376 + t489 * t377;
t455 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t484 + Ifges(4,4) * t488) * qJD(1);
t454 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t484 + Ifges(4,2) * t488) * qJD(1);
t453 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t484 + Ifges(4,6) * t488) * qJD(1);
t436 = Ifges(5,1) * t462 + Ifges(5,4) * t461 + Ifges(5,5) * t474;
t435 = Ifges(5,4) * t462 + Ifges(5,2) * t461 + Ifges(5,6) * t474;
t434 = Ifges(5,5) * t462 + Ifges(5,6) * t461 + Ifges(5,3) * t474;
t419 = Ifges(6,1) * t442 + Ifges(6,4) * t441 + Ifges(6,5) * t473;
t418 = Ifges(6,4) * t442 + Ifges(6,2) * t441 + Ifges(6,6) * t473;
t417 = Ifges(6,5) * t442 + Ifges(6,6) * t441 + Ifges(6,3) * t473;
t397 = mrSges(6,2) * t409 - mrSges(6,3) * t405 + Ifges(6,1) * t414 + Ifges(6,4) * t413 + Ifges(6,5) * t456 + t441 * t417 - t473 * t418;
t396 = -mrSges(6,1) * t409 + mrSges(6,3) * t406 + Ifges(6,4) * t414 + Ifges(6,2) * t413 + Ifges(6,6) * t456 - t442 * t417 + t473 * t419;
t385 = mrSges(5,2) * t423 - mrSges(5,3) * t410 + Ifges(5,1) * t440 + Ifges(5,4) * t439 + Ifges(5,5) * t460 - pkin(8) * t395 - t482 * t396 + t486 * t397 + t461 * t434 - t474 * t435;
t380 = -mrSges(5,1) * t423 + mrSges(5,3) * t411 + Ifges(5,4) * t440 + Ifges(5,2) * t439 + Ifges(5,6) * t460 - pkin(4) * t494 + pkin(8) * t495 + t486 * t396 + t482 * t397 - t462 * t434 + t474 * t436;
t379 = Ifges(4,4) * t467 + Ifges(4,2) * t468 + Ifges(4,6) * qJDD(3) - t453 * t504 + qJD(3) * t455 - mrSges(4,1) * t432 + mrSges(4,3) * t427 - Ifges(5,5) * t440 - Ifges(5,6) * t439 - Ifges(5,3) * t460 - t462 * t435 + t461 * t436 - mrSges(5,1) * t410 + mrSges(5,2) * t411 - Ifges(6,5) * t414 - Ifges(6,6) * t413 - Ifges(6,3) * t456 - t442 * t418 + t441 * t419 - mrSges(6,1) * t405 + mrSges(6,2) * t406 - pkin(4) * t395 - pkin(3) * t391;
t372 = mrSges(4,2) * t432 - mrSges(4,3) * t426 + Ifges(4,1) * t467 + Ifges(4,4) * t468 + Ifges(4,5) * qJDD(3) - pkin(7) * t391 - qJD(3) * t454 - t483 * t380 + t487 * t385 + t453 * t503;
t371 = Ifges(3,6) * qJDD(1) + t491 * Ifges(3,5) - mrSges(3,1) * t479 + mrSges(3,3) * t444 - Ifges(4,5) * t467 - Ifges(4,6) * t468 - Ifges(4,3) * qJDD(3) - mrSges(4,1) * t426 + mrSges(4,2) * t427 - t483 * t385 - t487 * t380 - pkin(3) * t492 - pkin(7) * t496 - pkin(2) * t384 + (-t484 * t454 + t488 * t455) * qJD(1);
t370 = mrSges(3,2) * t479 - mrSges(3,3) * t443 + Ifges(3,5) * qJDD(1) - t491 * Ifges(3,6) - pkin(6) * t384 + t488 * t372 - t484 * t379;
t369 = -mrSges(2,2) * g(3) - mrSges(2,3) * t471 + Ifges(2,5) * qJDD(1) - t491 * Ifges(2,6) - qJ(2) * t378 + t481 * t370 - t480 * t371;
t368 = mrSges(2,1) * g(3) + mrSges(2,3) * t472 + t491 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t501 + qJ(2) * t498 + t480 * t370 + t481 * t371;
t1 = [-m(1) * g(1) + t499; -m(1) * g(2) + t505; (-m(1) - m(2)) * g(3) + t501; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t505 - t485 * t368 + t489 * t369; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t499 + t489 * t368 + t485 * t369; pkin(1) * t378 + mrSges(2,1) * t471 - mrSges(2,2) * t472 + t484 * t372 + t488 * t379 + pkin(2) * t493 + pkin(6) * t497 + mrSges(3,1) * t443 - mrSges(3,2) * t444 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);];
tauB = t1;
