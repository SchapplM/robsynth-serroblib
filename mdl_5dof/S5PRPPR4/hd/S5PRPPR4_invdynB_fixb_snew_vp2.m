% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PRPPR4
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
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
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
% Datum: 2019-12-31 17:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PRPPR4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR4_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR4_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPPR4_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR4_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR4_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR4_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPPR4_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:36:51
% EndTime: 2019-12-31 17:36:53
% DurationCPUTime: 1.46s
% Computational Cost: add. (10975->203), mult. (23775->254), div. (0->0), fcn. (14420->8), ass. (0->92)
t446 = cos(pkin(8));
t488 = (Ifges(4,6) - Ifges(5,6)) * t446;
t445 = sin(pkin(7));
t447 = cos(pkin(7));
t427 = g(1) * t445 - g(2) * t447;
t428 = -g(1) * t447 - g(2) * t445;
t449 = sin(qJ(2));
t451 = cos(qJ(2));
t410 = t449 * t427 + t451 * t428;
t452 = qJD(2) ^ 2;
t487 = -pkin(2) * t452 + qJDD(2) * qJ(3) + (2 * qJD(2) * qJD(3)) + t410;
t486 = Ifges(5,4) + Ifges(4,5);
t444 = sin(pkin(8));
t441 = t444 ^ 2;
t442 = t446 ^ 2;
t476 = t442 * t452;
t484 = t441 * t452 + t476;
t443 = -g(3) + qJDD(1);
t393 = t443 * t446 - t487 * t444;
t483 = Ifges(4,1) + Ifges(5,1);
t482 = Ifges(4,4) - Ifges(5,5);
t481 = Ifges(4,2) + Ifges(5,3);
t480 = mrSges(4,2) * t444;
t479 = qJ(3) * t452;
t478 = qJ(4) * t444;
t423 = (-mrSges(5,1) * t446 - mrSges(5,3) * t444) * qJD(2);
t424 = (-mrSges(4,1) * t446 + t480) * qJD(2);
t460 = -pkin(3) * t446 - t478;
t422 = t460 * qJD(2);
t473 = t444 * qJD(2);
t389 = t422 * t473 + qJDD(4) - t393;
t385 = (-pkin(4) * t446 * t452 - pkin(6) * qJDD(2)) * t444 + t389;
t394 = t444 * t443 + t487 * t446;
t472 = t446 * qJD(2);
t391 = t422 * t472 + t394;
t469 = qJDD(2) * t446;
t386 = -pkin(4) * t476 - pkin(6) * t469 + t391;
t448 = sin(qJ(5));
t450 = cos(qJ(5));
t383 = t385 * t450 - t386 * t448;
t457 = -t444 * t448 - t446 * t450;
t416 = t457 * qJD(2);
t458 = t444 * t450 - t446 * t448;
t417 = t458 * qJD(2);
t402 = -mrSges(6,1) * t416 + mrSges(6,2) * t417;
t407 = qJD(5) * t416 + t458 * qJDD(2);
t411 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t416;
t381 = m(6) * t383 + qJDD(5) * mrSges(6,1) - mrSges(6,3) * t407 + qJD(5) * t411 - t402 * t417;
t384 = t385 * t448 + t386 * t450;
t406 = -qJD(5) * t417 + t457 * qJDD(2);
t412 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t417;
t382 = m(6) * t384 - qJDD(5) * mrSges(6,2) + mrSges(6,3) * t406 - qJD(5) * t412 + t402 * t416;
t374 = t381 * t450 + t382 * t448;
t454 = -m(5) * t389 - t374;
t372 = m(4) * t393 + ((-mrSges(5,2) - mrSges(4,3)) * qJDD(2) + (-t423 - t424) * qJD(2)) * t444 + t454;
t462 = -t381 * t448 + t450 * t382;
t456 = m(5) * t391 + mrSges(5,2) * t469 + t423 * t472 + t462;
t373 = m(4) * t394 + (qJDD(2) * mrSges(4,3) + qJD(2) * t424) * t446 + t456;
t463 = -t372 * t444 + t446 * t373;
t367 = m(3) * t410 - mrSges(3,1) * t452 - qJDD(2) * mrSges(3,2) + t463;
t409 = t451 * t427 - t449 * t428;
t466 = -qJDD(3) + t409;
t455 = -0.2e1 * qJD(4) * t473 - t466;
t392 = -t479 + (-pkin(2) + t460) * qJDD(2) + t455;
t388 = (qJ(3) + (-t441 - t442) * pkin(6)) * t452 + (t478 + pkin(2) + (pkin(3) + pkin(4)) * t446) * qJDD(2) - t455;
t459 = -m(6) * t388 + t406 * mrSges(6,1) - t407 * mrSges(6,2) + t416 * t411 - t417 * t412;
t470 = qJDD(2) * t444;
t379 = m(5) * t392 - mrSges(5,1) * t469 - t484 * mrSges(5,2) - mrSges(5,3) * t470 + t459;
t405 = -qJDD(2) * pkin(2) - t466 - t479;
t453 = -m(4) * t405 + mrSges(4,1) * t469 + t484 * mrSges(4,3) - t379;
t378 = t453 + m(3) * t409 - mrSges(3,2) * t452 + (mrSges(3,1) - t480) * qJDD(2);
t364 = t449 * t367 + t451 * t378;
t362 = m(2) * t427 + t364;
t464 = t451 * t367 - t449 * t378;
t363 = m(2) * t428 + t464;
t475 = t447 * t362 + t445 * t363;
t368 = t446 * t372 + t444 * t373;
t474 = (t486 * t444 + t488) * qJD(2);
t468 = m(3) * t443 + t368;
t465 = -t362 * t445 + t447 * t363;
t397 = Ifges(6,1) * t417 + Ifges(6,4) * t416 + Ifges(6,5) * qJD(5);
t396 = Ifges(6,4) * t417 + Ifges(6,2) * t416 + Ifges(6,6) * qJD(5);
t395 = Ifges(6,5) * t417 + Ifges(6,6) * t416 + Ifges(6,3) * qJD(5);
t376 = mrSges(6,2) * t388 - mrSges(6,3) * t383 + Ifges(6,1) * t407 + Ifges(6,4) * t406 + Ifges(6,5) * qJDD(5) - qJD(5) * t396 + t395 * t416;
t375 = -mrSges(6,1) * t388 + mrSges(6,3) * t384 + Ifges(6,4) * t407 + Ifges(6,2) * t406 + Ifges(6,6) * qJDD(5) + qJD(5) * t397 - t395 * t417;
t358 = mrSges(4,2) * t405 + mrSges(5,2) * t389 - mrSges(4,3) * t393 - mrSges(5,3) * t392 - pkin(6) * t374 - qJ(4) * t379 - t375 * t448 + t376 * t450 + t474 * t472 + (t483 * t444 + t482 * t446) * qJDD(2);
t357 = -mrSges(4,1) * t405 + mrSges(4,3) * t394 - mrSges(5,1) * t392 + mrSges(5,2) * t391 - t448 * t376 - t450 * t375 - pkin(4) * t459 - pkin(6) * t462 - pkin(3) * t379 - t474 * t473 + (t482 * t444 + t481 * t446) * qJDD(2);
t356 = Ifges(6,3) * qJDD(5) - pkin(3) * t454 - qJ(4) * t456 + t452 * Ifges(3,5) - mrSges(3,1) * t443 - t416 * t397 + t417 * t396 + Ifges(6,6) * t406 + Ifges(6,5) * t407 + mrSges(3,3) * t410 - mrSges(5,3) * t391 - mrSges(4,1) * t393 + mrSges(4,2) * t394 + mrSges(6,1) * t383 - mrSges(6,2) * t384 + mrSges(5,1) * t389 + pkin(4) * t374 - pkin(2) * t368 + (Ifges(3,6) - t488 + (mrSges(5,2) * pkin(3) - t486) * t444) * qJDD(2) + (t482 * t442 * qJD(2) + (pkin(3) * t423 - t482 * t473 + (-t481 + t483) * t472) * t444) * qJD(2);
t355 = mrSges(3,2) * t443 - mrSges(3,3) * t409 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t452 - qJ(3) * t368 - t357 * t444 + t358 * t446;
t354 = mrSges(2,2) * t443 - mrSges(2,3) * t427 - pkin(5) * t364 + t355 * t451 - t356 * t449;
t353 = -mrSges(2,1) * t443 + mrSges(2,3) * t428 - pkin(1) * t468 + pkin(5) * t464 + t449 * t355 + t451 * t356;
t1 = [-m(1) * g(1) + t465; -m(1) * g(2) + t475; -m(1) * g(3) + m(2) * t443 + t468; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t475 - t445 * t353 + t447 * t354; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t465 + t447 * t353 + t445 * t354; pkin(1) * t364 + mrSges(2,1) * t427 - mrSges(2,2) * t428 + qJ(3) * t463 + t444 * t358 + t446 * t357 + pkin(2) * (-mrSges(4,2) * t470 + t453) + mrSges(3,1) * t409 - mrSges(3,2) * t410 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(3,3) * qJDD(2);];
tauB = t1;
