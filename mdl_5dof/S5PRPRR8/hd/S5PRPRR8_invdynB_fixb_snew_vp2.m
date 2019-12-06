% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PRPRR8
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1]';
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
% Datum: 2019-12-05 16:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PRPRR8_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR8_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR8_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR8_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR8_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR8_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR8_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR8_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR8_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:02:40
% EndTime: 2019-12-05 16:02:44
% DurationCPUTime: 1.99s
% Computational Cost: add. (21056->219), mult. (38022->277), div. (0->0), fcn. (23471->10), ass. (0->99)
t440 = sin(pkin(9));
t442 = cos(pkin(9));
t428 = g(1) * t440 - g(2) * t442;
t429 = -g(1) * t442 - g(2) * t440;
t437 = -g(3) + qJDD(1);
t449 = cos(qJ(2));
t443 = cos(pkin(5));
t446 = sin(qJ(2));
t470 = t443 * t446;
t441 = sin(pkin(5));
t471 = t441 * t446;
t398 = t428 * t470 + t449 * t429 + t437 * t471;
t478 = -qJDD(2) * qJ(3) - (2 * qJD(3) * qJD(2)) - t398;
t397 = -t446 * t429 + (t428 * t443 + t437 * t441) * t449;
t477 = -pkin(2) - pkin(7);
t476 = mrSges(3,1) - mrSges(4,2);
t475 = (-Ifges(4,4) + Ifges(3,5));
t474 = Ifges(4,5) - Ifges(3,6);
t451 = qJD(2) ^ 2;
t453 = -qJ(3) * t451 + qJDD(3) - t397;
t394 = t477 * qJDD(2) + t453;
t411 = -t428 * t441 + t437 * t443;
t445 = sin(qJ(4));
t448 = cos(qJ(4));
t390 = t445 * t394 + t448 * t411;
t424 = (mrSges(5,1) * t445 + mrSges(5,2) * t448) * qJD(2);
t466 = qJD(2) * qJD(4);
t463 = t448 * t466;
t426 = -qJDD(2) * t445 - t463;
t467 = qJD(2) * t448;
t431 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t467;
t425 = (pkin(4) * t445 - pkin(8) * t448) * qJD(2);
t450 = qJD(4) ^ 2;
t468 = qJD(2) * t445;
t387 = -pkin(4) * t450 + qJDD(4) * pkin(8) - t425 * t468 + t390;
t393 = t477 * t451 - t478;
t464 = t445 * t466;
t427 = qJDD(2) * t448 - t464;
t388 = (-t427 + t464) * pkin(8) + (-t426 + t463) * pkin(4) + t393;
t444 = sin(qJ(5));
t447 = cos(qJ(5));
t384 = -t387 * t444 + t388 * t447;
t422 = qJD(4) * t447 - t444 * t467;
t405 = qJD(5) * t422 + qJDD(4) * t444 + t427 * t447;
t423 = qJD(4) * t444 + t447 * t467;
t406 = -mrSges(6,1) * t422 + mrSges(6,2) * t423;
t433 = qJD(5) + t468;
t408 = -mrSges(6,2) * t433 + mrSges(6,3) * t422;
t419 = qJDD(5) - t426;
t382 = m(6) * t384 + mrSges(6,1) * t419 - mrSges(6,3) * t405 - t406 * t423 + t408 * t433;
t385 = t387 * t447 + t388 * t444;
t404 = -qJD(5) * t423 + qJDD(4) * t447 - t427 * t444;
t409 = mrSges(6,1) * t433 - mrSges(6,3) * t423;
t383 = m(6) * t385 - mrSges(6,2) * t419 + mrSges(6,3) * t404 + t406 * t422 - t409 * t433;
t460 = -t382 * t444 + t447 * t383;
t374 = m(5) * t390 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t426 - qJD(4) * t431 - t424 * t468 + t460;
t472 = t411 * t445;
t389 = t394 * t448 - t472;
t430 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t468;
t386 = -qJDD(4) * pkin(4) - pkin(8) * t450 + t472 + (qJD(2) * t425 - t394) * t448;
t454 = -m(6) * t386 + t404 * mrSges(6,1) - mrSges(6,2) * t405 + t422 * t408 - t409 * t423;
t378 = m(5) * t389 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t427 + qJD(4) * t430 - t424 * t467 + t454;
t368 = t374 * t445 + t378 * t448;
t396 = -qJDD(2) * pkin(2) + t453;
t456 = -m(4) * t396 + (t451 * mrSges(4,3)) - t368;
t364 = m(3) * t397 - (mrSges(3,2) * t451) + t476 * qJDD(2) + t456;
t473 = t364 * t449;
t461 = t448 * t374 - t378 * t445;
t367 = m(4) * t411 + t461;
t366 = m(3) * t411 + t367;
t395 = pkin(2) * t451 + t478;
t375 = t447 * t382 + t444 * t383;
t455 = -m(5) * t393 + t426 * mrSges(5,1) - t427 * mrSges(5,2) - t430 * t468 - t431 * t467 - t375;
t452 = -m(4) * t395 + (t451 * mrSges(4,2)) + qJDD(2) * mrSges(4,3) - t455;
t372 = m(3) * t398 - (mrSges(3,1) * t451) - qJDD(2) * mrSges(3,2) + t452;
t355 = -t366 * t441 + t372 * t470 + t443 * t473;
t353 = m(2) * t428 + t355;
t359 = -t364 * t446 + t449 * t372;
t358 = m(2) * t429 + t359;
t469 = t442 * t353 + t440 * t358;
t354 = t443 * t366 + t372 * t471 + t441 * t473;
t462 = -t353 * t440 + t442 * t358;
t399 = Ifges(6,5) * t423 + Ifges(6,6) * t422 + Ifges(6,3) * t433;
t401 = Ifges(6,1) * t423 + Ifges(6,4) * t422 + Ifges(6,5) * t433;
t376 = -mrSges(6,1) * t386 + mrSges(6,3) * t385 + Ifges(6,4) * t405 + Ifges(6,2) * t404 + Ifges(6,6) * t419 - t399 * t423 + t401 * t433;
t400 = Ifges(6,4) * t423 + Ifges(6,2) * t422 + Ifges(6,6) * t433;
t377 = mrSges(6,2) * t386 - mrSges(6,3) * t384 + Ifges(6,1) * t405 + Ifges(6,4) * t404 + Ifges(6,5) * t419 + t399 * t422 - t400 * t433;
t413 = (Ifges(5,3) * qJD(4)) + (Ifges(5,5) * t448 - Ifges(5,6) * t445) * qJD(2);
t414 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t448 - Ifges(5,2) * t445) * qJD(2);
t360 = mrSges(5,2) * t393 - mrSges(5,3) * t389 + Ifges(5,1) * t427 + Ifges(5,4) * t426 + Ifges(5,5) * qJDD(4) - pkin(8) * t375 - qJD(4) * t414 - t376 * t444 + t377 * t447 - t413 * t468;
t415 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t448 - Ifges(5,4) * t445) * qJD(2);
t361 = -mrSges(5,1) * t393 - mrSges(6,1) * t384 + mrSges(6,2) * t385 + mrSges(5,3) * t390 + Ifges(5,4) * t427 - Ifges(6,5) * t405 + Ifges(5,2) * t426 + Ifges(5,6) * qJDD(4) - Ifges(6,6) * t404 - Ifges(6,3) * t419 - pkin(4) * t375 + qJD(4) * t415 - t400 * t423 + t401 * t422 - t413 * t467;
t350 = -mrSges(4,1) * t395 + mrSges(3,3) * t398 - pkin(2) * t367 - pkin(3) * t455 - pkin(7) * t461 - t474 * qJDD(2) - t445 * t360 - t448 * t361 - t476 * t411 + (t475 * t451);
t351 = -qJ(3) * t367 - mrSges(3,3) * t397 + pkin(3) * t368 + mrSges(4,1) * t396 + pkin(4) * t454 + pkin(8) * t460 + Ifges(5,6) * t426 + Ifges(5,3) * qJDD(4) + mrSges(5,1) * t389 - mrSges(5,2) * t390 + t444 * t377 + t447 * t376 + Ifges(5,5) * t427 + t474 * t451 + (mrSges(3,2) - mrSges(4,3)) * t411 + t475 * qJDD(2) + (t414 * t448 + t415 * t445) * qJD(2);
t457 = pkin(6) * t359 + t350 * t449 + t351 * t446;
t349 = mrSges(3,1) * t397 - mrSges(3,2) * t398 + mrSges(4,2) * t396 - mrSges(4,3) * t395 + t448 * t360 - t445 * t361 - pkin(7) * t368 + pkin(2) * t456 + qJ(3) * t452 + (-mrSges(4,2) * pkin(2) + Ifges(4,1) + Ifges(3,3)) * qJDD(2);
t348 = mrSges(2,2) * t437 - mrSges(2,3) * t428 - t350 * t446 + t351 * t449 + (-t354 * t441 - t355 * t443) * pkin(6);
t347 = -mrSges(2,1) * t437 + mrSges(2,3) * t429 - pkin(1) * t354 - t349 * t441 + t457 * t443;
t1 = [-m(1) * g(1) + t462; -m(1) * g(2) + t469; -m(1) * g(3) + m(2) * t437 + t354; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t469 - t440 * t347 + t442 * t348; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t462 + t442 * t347 + t440 * t348; -mrSges(1,1) * g(2) + mrSges(2,1) * t428 + mrSges(1,2) * g(1) - mrSges(2,2) * t429 + pkin(1) * t355 + t349 * t443 + t457 * t441;];
tauB = t1;
