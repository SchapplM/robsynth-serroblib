% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPRPR3
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
% Datum: 2019-12-05 17:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPRPR3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR3_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR3_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR3_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR3_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR3_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:51:26
% EndTime: 2019-12-05 17:51:28
% DurationCPUTime: 2.54s
% Computational Cost: add. (32310->206), mult. (45201->270), div. (0->0), fcn. (25113->10), ass. (0->99)
t461 = 2 * qJD(4);
t428 = sin(qJ(1));
t431 = cos(qJ(1));
t408 = t431 * g(2) + t428 * g(3);
t402 = qJDD(1) * pkin(1) + t408;
t407 = t428 * g(2) - t431 * g(3);
t432 = qJD(1) ^ 2;
t403 = -t432 * pkin(1) + t407;
t423 = sin(pkin(8));
t425 = cos(pkin(8));
t388 = t425 * t402 - t423 * t403;
t383 = qJDD(1) * pkin(2) + t388;
t389 = t423 * t402 + t425 * t403;
t384 = -t432 * pkin(2) + t389;
t427 = sin(qJ(3));
t430 = cos(qJ(3));
t379 = t427 * t383 + t430 * t384;
t420 = (qJD(1) + qJD(3));
t418 = t420 ^ 2;
t419 = qJDD(1) + qJDD(3);
t377 = -t418 * pkin(3) + t419 * qJ(4) + t379;
t460 = (t420 * t461) + t377;
t422 = sin(pkin(9));
t459 = mrSges(5,2) * t422;
t457 = mrSges(5,3) * t419;
t456 = t422 * t420;
t426 = sin(qJ(5));
t455 = t422 * t426;
t429 = cos(qJ(5));
t454 = t422 * t429;
t424 = cos(pkin(9));
t453 = t424 * t419;
t452 = t424 * t420;
t421 = -g(1) + qJDD(2);
t451 = t424 * t421;
t373 = t422 * t421 + t460 * t424;
t396 = (-mrSges(5,1) * t424 + t459) * t420;
t439 = -pkin(4) * t424 - pkin(7) * t422;
t398 = t439 * t420;
t371 = t398 * t452 + t373;
t378 = t430 * t383 - t427 * t384;
t434 = -t418 * qJ(4) + qJDD(4) - t378;
t374 = (-pkin(3) + t439) * t419 + t434;
t368 = -t426 * t371 + t429 * t374;
t405 = qJD(5) - t452;
t448 = t420 * t455;
t391 = -t405 * mrSges(6,2) - mrSges(6,3) * t448;
t393 = (mrSges(6,1) * t426 + mrSges(6,2) * t429) * t456;
t449 = qJD(5) * t420;
t395 = (t419 * t429 - t426 * t449) * t422;
t404 = qJDD(5) - t453;
t447 = t420 * t454;
t366 = m(6) * t368 + t404 * mrSges(6,1) - t395 * mrSges(6,3) + t405 * t391 - t393 * t447;
t369 = t429 * t371 + t426 * t374;
t392 = t405 * mrSges(6,1) - mrSges(6,3) * t447;
t394 = (-t419 * t426 - t429 * t449) * t422;
t367 = m(6) * t369 - t404 * mrSges(6,2) + t394 * mrSges(6,3) - t405 * t392 - t393 * t448;
t441 = -t426 * t366 + t429 * t367;
t359 = m(5) * t373 + (t396 * t420 + t457) * t424 + t441;
t372 = -t460 * t422 + t451;
t370 = -t451 + (t377 + (t461 + t398) * t420) * t422;
t435 = -m(6) * t370 + t394 * mrSges(6,1) - t395 * mrSges(6,2);
t364 = m(5) * t372 + (-t457 + (-t391 * t426 - t392 * t429 - t396) * t420) * t422 + t435;
t442 = t424 * t359 - t422 * t364;
t353 = m(4) * t379 - t418 * mrSges(4,1) - t419 * mrSges(4,2) + t442;
t360 = t429 * t366 + t426 * t367;
t376 = -t419 * pkin(3) + t434;
t433 = -m(5) * t376 + mrSges(5,1) * t453 - t360 + (t422 ^ 2 + t424 ^ 2) * mrSges(5,3) * t418;
t356 = m(4) * t378 - t418 * mrSges(4,2) + (mrSges(4,1) - t459) * t419 + t433;
t348 = t427 * t353 + t430 * t356;
t346 = m(3) * t388 + qJDD(1) * mrSges(3,1) - t432 * mrSges(3,2) + t348;
t443 = t430 * t353 - t427 * t356;
t347 = m(3) * t389 - t432 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t443;
t341 = t425 * t346 + t423 * t347;
t354 = t422 * t359 + t424 * t364;
t446 = m(4) * t421 + t354;
t444 = -t423 * t346 + t425 * t347;
t339 = m(2) * t407 - t432 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t444;
t340 = m(2) * t408 + qJDD(1) * mrSges(2,1) - t432 * mrSges(2,2) + t341;
t445 = t431 * t339 - t428 * t340;
t440 = m(3) * t421 + t446;
t438 = Ifges(5,1) * t422 + Ifges(5,4) * t424;
t437 = Ifges(5,5) * t422 + Ifges(5,6) * t424;
t436 = -t428 * t339 - t431 * t340;
t397 = t437 * t420;
t387 = Ifges(6,5) * t405 + (Ifges(6,1) * t429 - Ifges(6,4) * t426) * t456;
t386 = Ifges(6,6) * t405 + (Ifges(6,4) * t429 - Ifges(6,2) * t426) * t456;
t385 = Ifges(6,3) * t405 + (Ifges(6,5) * t429 - Ifges(6,6) * t426) * t456;
t362 = mrSges(6,2) * t370 - mrSges(6,3) * t368 + Ifges(6,1) * t395 + Ifges(6,4) * t394 + Ifges(6,5) * t404 - t385 * t448 - t405 * t386;
t361 = -mrSges(6,1) * t370 + mrSges(6,3) * t369 + Ifges(6,4) * t395 + Ifges(6,2) * t394 + Ifges(6,6) * t404 - t385 * t447 + t405 * t387;
t350 = Ifges(5,2) * t453 - mrSges(5,1) * t376 - mrSges(6,1) * t368 + mrSges(6,2) * t369 + mrSges(5,3) * t373 - Ifges(6,5) * t395 - Ifges(6,6) * t394 - Ifges(6,3) * t404 - pkin(4) * t360 + (Ifges(5,4) * t419 + (-t386 * t429 - t387 * t426 - t397) * t420) * t422;
t349 = mrSges(5,2) * t376 - mrSges(5,3) * t372 - pkin(7) * t360 - t426 * t361 + t429 * t362 + t397 * t452 + t438 * t419;
t342 = t418 * Ifges(4,5) - mrSges(4,1) * t421 + mrSges(4,3) * t379 - mrSges(5,1) * t372 + mrSges(5,2) * t373 - t426 * t362 - t429 * t361 - pkin(4) * t435 - pkin(7) * t441 - pkin(3) * t354 + (Ifges(4,6) - t437) * t419 + (-pkin(4) * (-t391 * t455 - t392 * t454) + (-t422 * (Ifges(5,4) * t422 + Ifges(5,2) * t424) + t424 * t438) * t420) * t420;
t337 = mrSges(4,2) * t421 - mrSges(4,3) * t378 + Ifges(4,5) * t419 - t418 * Ifges(4,6) - qJ(4) * t354 + t424 * t349 - t422 * t350;
t336 = mrSges(3,2) * t421 - mrSges(3,3) * t388 + Ifges(3,5) * qJDD(1) - t432 * Ifges(3,6) - pkin(6) * t348 + t430 * t337 - t427 * t342;
t335 = -mrSges(3,1) * t421 + mrSges(3,3) * t389 + t432 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t446 + pkin(6) * t443 + t427 * t337 + t430 * t342;
t334 = -mrSges(2,2) * g(1) - mrSges(2,3) * t408 + Ifges(2,5) * qJDD(1) - t432 * Ifges(2,6) - qJ(2) * t341 - t423 * t335 + t425 * t336;
t333 = mrSges(2,1) * g(1) + mrSges(2,3) * t407 + t432 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t440 + qJ(2) * t444 + t425 * t335 + t423 * t336;
t1 = [(-m(1) - m(2)) * g(1) + t440; -m(1) * g(2) + t436; -m(1) * g(3) + t445; pkin(1) * t341 + pkin(2) * t348 + mrSges(3,1) * t388 - mrSges(3,2) * t389 + t422 * t349 + t424 * t350 + pkin(3) * (-t419 * t459 + t433) + qJ(4) * t442 - mrSges(4,2) * t379 + mrSges(4,1) * t378 + mrSges(2,1) * t408 - mrSges(2,2) * t407 + Ifges(4,3) * t419 - mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + (Ifges(2,3) + Ifges(3,3)) * qJDD(1); mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t445 - t431 * t333 - t428 * t334; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t436 - t428 * t333 + t431 * t334;];
tauB = t1;
