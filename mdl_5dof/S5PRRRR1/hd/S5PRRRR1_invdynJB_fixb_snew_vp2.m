% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PRRRR1
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% Datum: 2019-12-05 17:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PRRRR1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(2,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR1_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR1_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_invdynJB_fixb_snew_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR1_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR1_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR1_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:03:03
% EndTime: 2019-12-05 17:03:05
% DurationCPUTime: 1.21s
% Computational Cost: add. (10096->223), mult. (19323->281), div. (0->0), fcn. (13744->8), ass. (0->86)
t432 = sin(qJ(4));
t433 = sin(qJ(3));
t436 = cos(qJ(4));
t437 = cos(qJ(3));
t415 = (t432 * t433 - t436 * t437) * qJD(2);
t430 = -g(3) + qJDD(1);
t434 = sin(qJ(2));
t438 = cos(qJ(2));
t422 = -t438 * g(1) + t434 * t430;
t410 = t433 * g(2) + t437 * t422;
t439 = qJD(2) ^ 2;
t404 = (-t437 ^ 2 * t439 - qJD(3) ^ 2) * pkin(2) + t410;
t409 = t437 * g(2) - t433 * t422;
t445 = (t433 * t437 * t439 + qJDD(3)) * pkin(2) + t409;
t388 = t436 * t404 + t432 * t445;
t452 = qJD(2) * qJD(3);
t450 = t433 * t452;
t420 = t437 * qJDD(2) - t450;
t421 = t434 * g(1) + t438 * t430;
t403 = (-t420 + t450) * pkin(2) - t421;
t431 = sin(qJ(5));
t435 = cos(qJ(5));
t379 = -t431 * t388 + t435 * t403;
t419 = t433 * qJDD(2) + t437 * t452;
t395 = -t415 * qJD(4) + t436 * t419 + t432 * t420;
t416 = (t432 * t437 + t433 * t436) * qJD(2);
t429 = qJD(3) + qJD(4);
t405 = -t431 * t416 + t435 * t429;
t428 = qJDD(3) + qJDD(4);
t382 = t405 * qJD(5) + t435 * t395 + t431 * t428;
t406 = t435 * t416 + t431 * t429;
t389 = -t405 * mrSges(6,1) + t406 * mrSges(6,2);
t394 = -t416 * qJD(4) - t432 * t419 + t436 * t420;
t393 = qJDD(5) - t394;
t411 = qJD(5) + t415;
t396 = -t411 * mrSges(6,2) + t405 * mrSges(6,3);
t377 = m(6) * t379 + t393 * mrSges(6,1) - t382 * mrSges(6,3) - t406 * t389 + t411 * t396;
t380 = t435 * t388 + t431 * t403;
t381 = -t406 * qJD(5) - t431 * t395 + t435 * t428;
t397 = t411 * mrSges(6,1) - t406 * mrSges(6,3);
t378 = m(6) * t380 - t393 * mrSges(6,2) + t381 * mrSges(6,3) + t405 * t389 - t411 * t397;
t401 = t415 * mrSges(5,1) + t416 * mrSges(5,2);
t408 = t429 * mrSges(5,1) - t416 * mrSges(5,3);
t369 = m(5) * t388 - t428 * mrSges(5,2) + t394 * mrSges(5,3) - t431 * t377 + t435 * t378 - t415 * t401 - t429 * t408;
t387 = t432 * t404 - t436 * t445;
t407 = -t429 * mrSges(5,2) - t415 * mrSges(5,3);
t376 = t428 * mrSges(5,1) + t381 * mrSges(6,1) - t382 * mrSges(6,2) - t395 * mrSges(5,3) + t405 * t396 - t406 * t397 - t416 * t401 + t429 * t407 + (-m(5) - m(6)) * t387;
t364 = t432 * t369 + t436 * t376;
t413 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t433 + Ifges(4,2) * t437) * qJD(2);
t414 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t433 + Ifges(4,4) * t437) * qJD(2);
t383 = Ifges(6,5) * t406 + Ifges(6,6) * t405 + Ifges(6,3) * t411;
t385 = Ifges(6,1) * t406 + Ifges(6,4) * t405 + Ifges(6,5) * t411;
t373 = -mrSges(6,1) * t387 + mrSges(6,3) * t380 + Ifges(6,4) * t382 + Ifges(6,2) * t381 + Ifges(6,6) * t393 - t406 * t383 + t411 * t385;
t384 = Ifges(6,4) * t406 + Ifges(6,2) * t405 + Ifges(6,6) * t411;
t374 = mrSges(6,2) * t387 - mrSges(6,3) * t379 + Ifges(6,1) * t382 + Ifges(6,4) * t381 + Ifges(6,5) * t393 + t405 * t383 - t411 * t384;
t399 = Ifges(5,4) * t416 - Ifges(5,2) * t415 + Ifges(5,6) * t429;
t400 = Ifges(5,1) * t416 - Ifges(5,4) * t415 + Ifges(5,5) * t429;
t444 = mrSges(5,1) * t387 + mrSges(5,2) * t388 - Ifges(5,5) * t395 - Ifges(5,6) * t394 - Ifges(5,3) * t428 - t435 * t373 - t431 * t374 - t416 * t399 - t415 * t400;
t458 = mrSges(4,1) * t409 - mrSges(4,2) * t410 + Ifges(4,5) * t419 + Ifges(4,6) * t420 + Ifges(4,3) * qJDD(3) + pkin(2) * t364 + (t433 * t413 - t437 * t414) * qJD(2) - t444;
t457 = -m(2) - m(3);
t456 = mrSges(1,3) + mrSges(2,3);
t418 = (-mrSges(4,1) * t437 + mrSges(4,2) * t433) * qJD(2);
t453 = qJD(2) * t437;
t424 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t453;
t454 = qJD(2) * t433;
t362 = m(4) * t409 + qJDD(3) * mrSges(4,1) - t419 * mrSges(4,3) + qJD(3) * t424 - t418 * t454 + t364;
t423 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t454;
t363 = m(4) * t410 - qJDD(3) * mrSges(4,2) + t420 * mrSges(4,3) - qJD(3) * t423 + t436 * t369 - t432 * t376 + t418 * t453;
t359 = m(3) * t422 - t439 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t433 * t362 + t437 * t363;
t441 = m(5) * t403 - t394 * mrSges(5,1) + t395 * mrSges(5,2) + t435 * t377 + t431 * t378 + t415 * t407 + t416 * t408;
t367 = qJDD(2) * mrSges(3,1) + t420 * mrSges(4,1) - t439 * mrSges(3,2) - t419 * mrSges(4,2) + (m(3) + m(4)) * t421 + (-t423 * t433 + t424 * t437) * qJD(2) - t441;
t455 = t434 * t359 + t438 * t367;
t451 = m(2) * t430 + t455;
t449 = t438 * t359 - t434 * t367;
t448 = -t437 * t362 - t433 * t363;
t398 = Ifges(5,5) * t416 - Ifges(5,6) * t415 + Ifges(5,3) * t429;
t365 = mrSges(5,2) * t403 + mrSges(5,3) * t387 + Ifges(5,1) * t395 + Ifges(5,4) * t394 + Ifges(5,5) * t428 - t431 * t373 + t435 * t374 - t415 * t398 - t429 * t399;
t442 = mrSges(6,1) * t379 - mrSges(6,2) * t380 + Ifges(6,5) * t382 + Ifges(6,6) * t381 + Ifges(6,3) * t393 + t406 * t384 - t405 * t385;
t370 = -mrSges(5,1) * t403 + mrSges(5,3) * t388 + Ifges(5,4) * t395 + Ifges(5,2) * t394 + Ifges(5,6) * t428 - t416 * t398 + t429 * t400 - t442;
t412 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t433 + Ifges(4,6) * t437) * qJD(2);
t356 = mrSges(4,1) * t421 + mrSges(4,3) * t410 + Ifges(4,4) * t419 + Ifges(4,2) * t420 + Ifges(4,6) * qJDD(3) - pkin(2) * t441 + qJD(3) * t414 + t432 * t365 + t436 * t370 - t412 * t454;
t361 = -mrSges(4,2) * t421 - mrSges(4,3) * t409 + Ifges(4,1) * t419 + Ifges(4,4) * t420 + Ifges(4,5) * qJDD(3) - qJD(3) * t413 + t436 * t365 - t432 * t370 + t412 * t453;
t443 = mrSges(3,1) * t421 - mrSges(3,2) * t422 + Ifges(3,3) * qJDD(2) + t437 * t356 + t433 * t361;
t360 = -mrSges(3,1) * g(2) + mrSges(3,3) * t422 + t439 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - t458;
t355 = mrSges(3,2) * g(2) - mrSges(3,3) * t421 + Ifges(3,5) * qJDD(2) - t439 * Ifges(3,6) - t433 * t356 + t437 * t361;
t1 = [(-m(1) - m(2)) * g(1) + t449; (-m(1) + t457) * g(2) + t448; -m(1) * g(3) + t451; -mrSges(1,2) * g(3) + mrSges(2,2) * t430 + t438 * t355 - t434 * t360 - qJ(1) * t448 + (-qJ(1) * t457 + t456) * g(2); qJ(1) * t449 + mrSges(1,1) * g(3) - pkin(1) * t455 - mrSges(2,1) * t430 + (-qJ(1) * m(2) - t456) * g(1) - t443; t434 * t355 + t438 * t360 + pkin(1) * t448 + (-pkin(1) * m(3) - mrSges(1,1) - mrSges(2,1)) * g(2) + (mrSges(1,2) + mrSges(2,2)) * g(1); t451; t443; t458; -t444; t442;];
tauJB = t1;
