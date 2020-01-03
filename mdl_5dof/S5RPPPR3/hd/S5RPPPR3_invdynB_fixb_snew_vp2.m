% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPPPR3
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3]';
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
% Datum: 2019-12-31 17:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPPPR3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR3_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR3_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR3_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR3_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR3_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR3_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR3_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:43:57
% EndTime: 2019-12-31 17:43:58
% DurationCPUTime: 1.56s
% Computational Cost: add. (12257->214), mult. (25705->262), div. (0->0), fcn. (14420->8), ass. (0->94)
t470 = cos(pkin(8));
t512 = (Ifges(4,6) - Ifges(5,6)) * t470;
t473 = sin(qJ(1));
t475 = cos(qJ(1));
t450 = t473 * g(1) - t475 * g(2);
t448 = qJDD(1) * pkin(1) + t450;
t451 = -t475 * g(1) - t473 * g(2);
t476 = qJD(1) ^ 2;
t449 = -t476 * pkin(1) + t451;
t469 = sin(pkin(7));
t471 = cos(pkin(7));
t431 = t469 * t448 + t471 * t449;
t511 = -t476 * pkin(2) + qJDD(1) * qJ(3) + (2 * qJD(1) * qJD(3)) + t431;
t510 = Ifges(5,4) + Ifges(4,5);
t468 = sin(pkin(8));
t465 = t468 ^ 2;
t466 = t470 ^ 2;
t500 = t466 * t476;
t508 = t465 * t476 + t500;
t467 = -g(3) + qJDD(2);
t414 = t470 * t467 - t511 * t468;
t507 = Ifges(4,1) + Ifges(5,1);
t506 = Ifges(4,4) - Ifges(5,5);
t505 = Ifges(4,2) + Ifges(5,3);
t504 = mrSges(4,2) * t468;
t503 = mrSges(5,3) * t468;
t502 = qJ(4) * t468;
t499 = t476 * qJ(3);
t444 = (-mrSges(5,1) * t470 - t503) * qJD(1);
t445 = (-mrSges(4,1) * t470 + t504) * qJD(1);
t484 = -pkin(3) * t470 - t502;
t443 = t484 * qJD(1);
t496 = t468 * qJD(1);
t410 = t443 * t496 + qJDD(4) - t414;
t406 = (-pkin(4) * t470 * t476 - pkin(6) * qJDD(1)) * t468 + t410;
t415 = t468 * t467 + t511 * t470;
t495 = t470 * qJD(1);
t412 = t443 * t495 + t415;
t493 = qJDD(1) * t470;
t407 = -pkin(4) * t500 - pkin(6) * t493 + t412;
t472 = sin(qJ(5));
t474 = cos(qJ(5));
t404 = t474 * t406 - t472 * t407;
t481 = -t468 * t472 - t470 * t474;
t437 = t481 * qJD(1);
t482 = t468 * t474 - t470 * t472;
t438 = t482 * qJD(1);
t425 = -t437 * mrSges(6,1) + t438 * mrSges(6,2);
t429 = t437 * qJD(5) + t482 * qJDD(1);
t432 = -qJD(5) * mrSges(6,2) + t437 * mrSges(6,3);
t402 = m(6) * t404 + qJDD(5) * mrSges(6,1) - t429 * mrSges(6,3) + qJD(5) * t432 - t438 * t425;
t405 = t472 * t406 + t474 * t407;
t428 = -t438 * qJD(5) + t481 * qJDD(1);
t433 = qJD(5) * mrSges(6,1) - t438 * mrSges(6,3);
t403 = m(6) * t405 - qJDD(5) * mrSges(6,2) + t428 * mrSges(6,3) - qJD(5) * t433 + t437 * t425;
t395 = t474 * t402 + t472 * t403;
t478 = -m(5) * t410 - t395;
t393 = m(4) * t414 + ((-mrSges(5,2) - mrSges(4,3)) * qJDD(1) + (-t444 - t445) * qJD(1)) * t468 + t478;
t486 = -t472 * t402 + t474 * t403;
t480 = m(5) * t412 + mrSges(5,2) * t493 + t444 * t495 + t486;
t394 = m(4) * t415 + (qJDD(1) * mrSges(4,3) + qJD(1) * t445) * t470 + t480;
t487 = -t468 * t393 + t470 * t394;
t388 = m(3) * t431 - t476 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t487;
t430 = t471 * t448 - t469 * t449;
t490 = -qJDD(3) + t430;
t479 = -0.2e1 * qJD(4) * t496 - t490;
t413 = -t499 + (-pkin(2) + t484) * qJDD(1) + t479;
t409 = (qJ(3) + (-t465 - t466) * pkin(6)) * t476 + (t502 + pkin(2) + (pkin(3) + pkin(4)) * t470) * qJDD(1) - t479;
t483 = -m(6) * t409 + t428 * mrSges(6,1) - t429 * mrSges(6,2) + t437 * t432 - t438 * t433;
t400 = m(5) * t413 - mrSges(5,1) * t493 - t508 * mrSges(5,2) - qJDD(1) * t503 + t483;
t418 = -qJDD(1) * pkin(2) - t490 - t499;
t477 = -m(4) * t418 + mrSges(4,1) * t493 + t508 * mrSges(4,3) - t400;
t399 = t477 + (mrSges(3,1) - t504) * qJDD(1) - t476 * mrSges(3,2) + m(3) * t430;
t385 = t469 * t388 + t471 * t399;
t383 = m(2) * t450 + qJDD(1) * mrSges(2,1) - t476 * mrSges(2,2) + t385;
t488 = t471 * t388 - t469 * t399;
t384 = m(2) * t451 - t476 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t488;
t498 = t475 * t383 + t473 * t384;
t389 = t470 * t393 + t468 * t394;
t497 = (t510 * t468 + t512) * qJD(1);
t492 = m(3) * t467 + t389;
t489 = -t473 * t383 + t475 * t384;
t421 = Ifges(6,1) * t438 + Ifges(6,4) * t437 + Ifges(6,5) * qJD(5);
t420 = Ifges(6,4) * t438 + Ifges(6,2) * t437 + Ifges(6,6) * qJD(5);
t419 = Ifges(6,5) * t438 + Ifges(6,6) * t437 + Ifges(6,3) * qJD(5);
t397 = mrSges(6,2) * t409 - mrSges(6,3) * t404 + Ifges(6,1) * t429 + Ifges(6,4) * t428 + Ifges(6,5) * qJDD(5) - qJD(5) * t420 + t437 * t419;
t396 = -mrSges(6,1) * t409 + mrSges(6,3) * t405 + Ifges(6,4) * t429 + Ifges(6,2) * t428 + Ifges(6,6) * qJDD(5) + qJD(5) * t421 - t438 * t419;
t379 = mrSges(4,2) * t418 + mrSges(5,2) * t410 - mrSges(4,3) * t414 - mrSges(5,3) * t413 - pkin(6) * t395 - qJ(4) * t400 - t472 * t396 + t474 * t397 + t497 * t495 + (t507 * t468 + t506 * t470) * qJDD(1);
t378 = -mrSges(4,1) * t418 + mrSges(4,3) * t415 - mrSges(5,1) * t413 + mrSges(5,2) * t412 - t472 * t397 - t474 * t396 - pkin(4) * t483 - pkin(6) * t486 - pkin(3) * t400 - t497 * t496 + (t506 * t468 + t505 * t470) * qJDD(1);
t377 = Ifges(6,3) * qJDD(5) + t476 * Ifges(3,5) - qJ(4) * t480 - pkin(3) * t478 - mrSges(3,1) * t467 + Ifges(6,5) * t429 + mrSges(3,3) * t431 - t437 * t421 + t438 * t420 + mrSges(4,2) * t415 + Ifges(6,6) * t428 + mrSges(5,1) * t410 - mrSges(5,3) * t412 - mrSges(4,1) * t414 + mrSges(6,1) * t404 - mrSges(6,2) * t405 + pkin(4) * t395 - pkin(2) * t389 + (Ifges(3,6) - t512 + (mrSges(5,2) * pkin(3) - t510) * t468) * qJDD(1) + (t506 * t466 * qJD(1) + (pkin(3) * t444 - t506 * t496 + (-t505 + t507) * t495) * t468) * qJD(1);
t376 = mrSges(3,2) * t467 - mrSges(3,3) * t430 + Ifges(3,5) * qJDD(1) - t476 * Ifges(3,6) - qJ(3) * t389 - t468 * t378 + t470 * t379;
t375 = -mrSges(2,2) * g(3) - mrSges(2,3) * t450 + Ifges(2,5) * qJDD(1) - t476 * Ifges(2,6) - qJ(2) * t385 + t471 * t376 - t469 * t377;
t374 = mrSges(2,1) * g(3) + mrSges(2,3) * t451 + t476 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t492 + qJ(2) * t488 + t469 * t376 + t471 * t377;
t1 = [-m(1) * g(1) + t489; -m(1) * g(2) + t498; (-m(1) - m(2)) * g(3) + t492; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t498 - t473 * t374 + t475 * t375; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t489 + t475 * t374 + t473 * t375; pkin(1) * t385 + mrSges(2,1) * t450 - mrSges(2,2) * t451 + t470 * t378 + pkin(2) * t477 + qJ(3) * t487 + mrSges(3,1) * t430 - mrSges(3,2) * t431 + t468 * t379 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-pkin(2) * t504 + Ifges(2,3) + Ifges(3,3)) * qJDD(1);];
tauB = t1;
