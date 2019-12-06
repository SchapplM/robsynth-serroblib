% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PRPRP5
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
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PRPRP5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP5_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP5_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP5_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP5_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP5_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP5_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP5_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:37:47
% EndTime: 2019-12-05 15:37:51
% DurationCPUTime: 1.90s
% Computational Cost: add. (15535->215), mult. (34261->263), div. (0->0), fcn. (21974->8), ass. (0->98)
t524 = Ifges(5,1) + Ifges(6,1);
t518 = Ifges(5,4) - Ifges(6,5);
t517 = Ifges(5,5) + Ifges(6,4);
t523 = Ifges(5,2) + Ifges(6,3);
t516 = Ifges(5,6) - Ifges(6,6);
t522 = -Ifges(5,3) - Ifges(6,2);
t485 = qJD(2) ^ 2;
t521 = cos(qJ(4));
t480 = cos(pkin(8));
t520 = pkin(3) * t480;
t519 = -mrSges(5,3) - mrSges(6,2);
t478 = sin(pkin(8));
t515 = mrSges(4,2) * t478;
t514 = cos(pkin(7));
t474 = t480 ^ 2;
t513 = t474 * t485;
t479 = sin(pkin(7));
t463 = -t514 * g(1) - t479 * g(2);
t477 = -g(3) + qJDD(1);
t482 = sin(qJ(2));
t483 = cos(qJ(2));
t453 = t483 * t463 + t482 * t477;
t445 = -t485 * pkin(2) + qJDD(2) * qJ(3) + t453;
t462 = t479 * g(1) - t514 * g(2);
t503 = qJD(2) * qJD(3);
t507 = -t480 * t462 - 0.2e1 * t478 * t503;
t421 = (-pkin(6) * qJDD(2) + t485 * t520 - t445) * t478 + t507;
t424 = -t478 * t462 + (t445 + 0.2e1 * t503) * t480;
t501 = qJDD(2) * t480;
t422 = -pkin(3) * t513 + pkin(6) * t501 + t424;
t481 = sin(qJ(4));
t418 = t481 * t421 + t521 * t422;
t499 = t480 * t521;
t502 = qJDD(2) * t478;
t488 = t521 * t478 + t480 * t481;
t455 = t488 * qJD(2);
t505 = t455 * qJD(4);
t441 = -qJDD(2) * t499 + t481 * t502 + t505;
t449 = qJD(4) * mrSges(5,1) - t455 * mrSges(5,3);
t504 = t478 * qJD(2);
t454 = -qJD(2) * t499 + t481 * t504;
t435 = t454 * pkin(4) - t455 * qJ(5);
t484 = qJD(4) ^ 2;
t413 = -t484 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) - t454 * t435 + t418;
t450 = -qJD(4) * mrSges(6,1) + t455 * mrSges(6,2);
t500 = m(6) * t413 + qJDD(4) * mrSges(6,3) + qJD(4) * t450;
t436 = t454 * mrSges(6,1) - t455 * mrSges(6,3);
t508 = -t454 * mrSges(5,1) - t455 * mrSges(5,2) - t436;
t409 = m(5) * t418 - qJDD(4) * mrSges(5,2) - qJD(4) * t449 + t519 * t441 + t508 * t454 + t500;
t417 = t521 * t421 - t481 * t422;
t506 = t454 * qJD(4);
t442 = t488 * qJDD(2) - t506;
t448 = -qJD(4) * mrSges(5,2) - t454 * mrSges(5,3);
t414 = -qJDD(4) * pkin(4) - t484 * qJ(5) + t455 * t435 + qJDD(5) - t417;
t451 = -t454 * mrSges(6,2) + qJD(4) * mrSges(6,3);
t494 = -m(6) * t414 + qJDD(4) * mrSges(6,1) + qJD(4) * t451;
t410 = m(5) * t417 + qJDD(4) * mrSges(5,1) + qJD(4) * t448 + t519 * t442 + t508 * t455 + t494;
t403 = t481 * t409 + t521 * t410;
t423 = -t478 * t445 + t507;
t489 = mrSges(4,3) * qJDD(2) + t485 * (-mrSges(4,1) * t480 + t515);
t401 = m(4) * t423 - t489 * t478 + t403;
t495 = t521 * t409 - t481 * t410;
t402 = m(4) * t424 + t489 * t480 + t495;
t496 = -t478 * t401 + t480 * t402;
t394 = m(3) * t453 - t485 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t496;
t452 = -t482 * t463 + t483 * t477;
t490 = qJDD(3) - t452;
t444 = -qJDD(2) * pkin(2) - t485 * qJ(3) + t490;
t473 = t478 ^ 2;
t425 = (-pkin(2) - t520) * qJDD(2) + (-qJ(3) + (-t473 - t474) * pkin(6)) * t485 + t490;
t416 = -0.2e1 * qJD(5) * t455 + (-t442 + t506) * qJ(5) + (t441 + t505) * pkin(4) + t425;
t411 = m(6) * t416 + t441 * mrSges(6,1) - t442 * mrSges(6,3) - t455 * t450 + t454 * t451;
t487 = m(5) * t425 + t441 * mrSges(5,1) + t442 * mrSges(5,2) + t454 * t448 + t455 * t449 + t411;
t486 = -m(4) * t444 + mrSges(4,1) * t501 - t487 + (t473 * t485 + t513) * mrSges(4,3);
t405 = t486 + m(3) * t452 - t485 * mrSges(3,2) + (mrSges(3,1) - t515) * qJDD(2);
t497 = t483 * t394 - t482 * t405;
t390 = m(2) * t463 + t497;
t397 = t480 * t401 + t478 * t402;
t396 = (m(2) + m(3)) * t462 - t397;
t512 = t479 * t390 + t514 * t396;
t391 = t482 * t394 + t483 * t405;
t511 = -t516 * qJD(4) + t523 * t454 - t518 * t455;
t510 = t522 * qJD(4) + t516 * t454 - t517 * t455;
t509 = t517 * qJD(4) - t518 * t454 + t524 * t455;
t498 = t514 * t390 - t479 * t396;
t493 = Ifges(4,1) * t478 + Ifges(4,4) * t480;
t492 = Ifges(4,4) * t478 + Ifges(4,2) * t480;
t491 = Ifges(4,5) * t478 + Ifges(4,6) * t480;
t461 = t491 * qJD(2);
t399 = mrSges(5,2) * t425 + mrSges(6,2) * t414 - mrSges(5,3) * t417 - mrSges(6,3) * t416 - qJ(5) * t411 + t511 * qJD(4) + t517 * qJDD(4) - t518 * t441 + t524 * t442 + t510 * t454;
t398 = -mrSges(5,1) * t425 - mrSges(6,1) * t416 + mrSges(6,2) * t413 + mrSges(5,3) * t418 - pkin(4) * t411 + t509 * qJD(4) + t516 * qJDD(4) - t523 * t441 + t518 * t442 + t510 * t455;
t387 = t480 * qJD(2) * t461 + mrSges(4,2) * t444 - mrSges(4,3) * t423 - pkin(6) * t403 + t493 * qJDD(2) - t481 * t398 + t521 * t399;
t386 = -mrSges(4,1) * t444 + mrSges(4,3) * t424 - pkin(3) * t487 + pkin(6) * t495 + t492 * qJDD(2) + t521 * t398 + t481 * t399 - t461 * t504;
t385 = mrSges(3,3) * t453 - mrSges(6,3) * t413 + mrSges(6,1) * t414 - mrSges(5,1) * t417 + mrSges(5,2) * t418 - mrSges(4,1) * t423 + mrSges(4,2) * t424 + mrSges(3,1) * t462 - qJ(5) * t500 - pkin(4) * t494 - pkin(3) * t403 - pkin(2) * t397 + (pkin(4) * t436 + t511) * t455 + (qJ(5) * t436 - t509) * t454 + (pkin(4) * mrSges(6,2) - t517) * t442 + (qJ(5) * mrSges(6,2) + t516) * t441 + t522 * qJDD(4) + (Ifges(3,6) - t491) * qJDD(2) + (-t478 * t492 + t480 * t493 + Ifges(3,5)) * t485;
t384 = -mrSges(3,2) * t462 - mrSges(3,3) * t452 + Ifges(3,5) * qJDD(2) - t485 * Ifges(3,6) - qJ(3) * t397 - t478 * t386 + t480 * t387;
t383 = -mrSges(2,1) * t477 + mrSges(2,3) * t463 - Ifges(3,3) * qJDD(2) - mrSges(3,1) * t452 + mrSges(3,2) * t453 - t478 * t387 - t480 * t386 - pkin(2) * (-mrSges(4,2) * t502 + t486) - qJ(3) * t496 - pkin(1) * t391;
t382 = mrSges(2,2) * t477 - mrSges(2,3) * t462 - pkin(5) * t391 + t483 * t384 - t482 * t385;
t1 = [-m(1) * g(1) + t498; -m(1) * g(2) + t512; -m(1) * g(3) + m(2) * t477 + t391; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t512 + t514 * t382 - t479 * t383; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t498 + t479 * t382 + t514 * t383; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + mrSges(2,1) * t462 - mrSges(2,2) * t463 + t482 * t384 + t483 * t385 + pkin(1) * (m(3) * t462 - t397) + pkin(5) * t497;];
tauB = t1;
