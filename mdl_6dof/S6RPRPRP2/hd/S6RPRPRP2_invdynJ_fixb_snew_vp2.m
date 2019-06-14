% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJ [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 17:36
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRPRP2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP2_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP2_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP2_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP2_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP2_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP2_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP2_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:32:55
% EndTime: 2019-05-05 17:32:58
% DurationCPUTime: 2.06s
% Computational Cost: add. (15233->268), mult. (32625->332), div. (0->0), fcn. (21136->10), ass. (0->110)
t501 = -2 * qJD(4);
t500 = Ifges(6,1) + Ifges(7,1);
t494 = Ifges(6,4) - Ifges(7,5);
t493 = -Ifges(6,5) - Ifges(7,4);
t499 = Ifges(6,2) + Ifges(7,3);
t492 = Ifges(6,6) - Ifges(7,6);
t498 = -Ifges(6,3) - Ifges(7,2);
t464 = sin(qJ(1));
t466 = cos(qJ(1));
t480 = t464 * g(1) - t466 * g(2);
t443 = qJDD(1) * pkin(1) + t480;
t468 = qJD(1) ^ 2;
t474 = -t466 * g(1) - t464 * g(2);
t445 = -t468 * pkin(1) + t474;
t459 = sin(pkin(9));
t461 = cos(pkin(9));
t421 = t459 * t443 + t461 * t445;
t412 = -t468 * pkin(2) + qJDD(1) * pkin(7) + t421;
t457 = -g(3) + qJDD(2);
t463 = sin(qJ(3));
t465 = cos(qJ(3));
t401 = -t463 * t412 + t465 * t457;
t483 = qJD(1) * qJD(3);
t481 = t465 * t483;
t446 = t463 * qJDD(1) + t481;
t380 = (-t446 + t481) * qJ(4) + (t463 * t465 * t468 + qJDD(3)) * pkin(3) + t401;
t402 = t465 * t412 + t463 * t457;
t447 = t465 * qJDD(1) - t463 * t483;
t486 = qJD(1) * t463;
t448 = qJD(3) * pkin(3) - qJ(4) * t486;
t456 = t465 ^ 2;
t381 = -t456 * t468 * pkin(3) + t447 * qJ(4) - qJD(3) * t448 + t402;
t458 = sin(pkin(10));
t460 = cos(pkin(10));
t433 = (t465 * t458 + t463 * t460) * qJD(1);
t373 = t460 * t380 - t458 * t381 + t433 * t501;
t432 = (t463 * t458 - t465 * t460) * qJD(1);
t423 = t460 * t446 + t458 * t447;
t462 = sin(qJ(5));
t496 = cos(qJ(5));
t424 = -qJD(3) * t496 + t462 * t433;
t394 = -t424 * qJD(5) + t462 * qJDD(3) + t423 * t496;
t425 = t462 * qJD(3) + t433 * t496;
t398 = t424 * mrSges(7,1) - t425 * mrSges(7,3);
t374 = t458 * t380 + t460 * t381 + t432 * t501;
t415 = t432 * pkin(4) - t433 * pkin(8);
t467 = qJD(3) ^ 2;
t372 = -t467 * pkin(4) + qJDD(3) * pkin(8) - t432 * t415 + t374;
t420 = t461 * t443 - t459 * t445;
t472 = -qJDD(1) * pkin(2) - t420;
t382 = -t447 * pkin(3) + qJDD(4) + t448 * t486 + (-qJ(4) * t456 - pkin(7)) * t468 + t472;
t422 = -t458 * t446 + t460 * t447;
t376 = (qJD(3) * t432 - t423) * pkin(8) + (qJD(3) * t433 - t422) * pkin(4) + t382;
t368 = -t462 * t372 + t376 * t496;
t397 = t424 * pkin(5) - t425 * qJ(6);
t419 = qJDD(5) - t422;
t431 = qJD(5) + t432;
t430 = t431 ^ 2;
t366 = -t419 * pkin(5) - t430 * qJ(6) + t425 * t397 + qJDD(6) - t368;
t403 = -t424 * mrSges(7,2) + t431 * mrSges(7,3);
t475 = -m(7) * t366 + t419 * mrSges(7,1) + t431 * t403;
t362 = t394 * mrSges(7,2) + t425 * t398 - t475;
t369 = t372 * t496 + t462 * t376;
t365 = -t430 * pkin(5) + t419 * qJ(6) + 0.2e1 * qJD(6) * t431 - t424 * t397 + t369;
t393 = t425 * qJD(5) - qJDD(3) * t496 + t462 * t423;
t406 = -t431 * mrSges(7,1) + t425 * mrSges(7,2);
t482 = m(7) * t365 + t419 * mrSges(7,3) + t431 * t406;
t488 = t494 * t424 - t500 * t425 + t493 * t431;
t489 = t499 * t424 - t494 * t425 - t492 * t431;
t497 = -t393 * t492 - t394 * t493 - t498 * t419 - t424 * t488 - t425 * t489 + mrSges(6,1) * t368 - mrSges(7,1) * t366 - mrSges(6,2) * t369 + mrSges(7,3) * t365 - pkin(5) * t362 + qJ(6) * (-t393 * mrSges(7,2) - t424 * t398 + t482);
t495 = -mrSges(6,3) - mrSges(7,2);
t414 = t432 * mrSges(5,1) + t433 * mrSges(5,2);
t427 = qJD(3) * mrSges(5,1) - t433 * mrSges(5,3);
t405 = t431 * mrSges(6,1) - t425 * mrSges(6,3);
t487 = -t424 * mrSges(6,1) - t425 * mrSges(6,2) - t398;
t358 = m(6) * t369 - t419 * mrSges(6,2) + t393 * t495 - t431 * t405 + t424 * t487 + t482;
t404 = -t431 * mrSges(6,2) - t424 * mrSges(6,3);
t360 = m(6) * t368 + t419 * mrSges(6,1) + t394 * t495 + t431 * t404 + t425 * t487 + t475;
t477 = t358 * t496 - t360 * t462;
t349 = m(5) * t374 - qJDD(3) * mrSges(5,2) + mrSges(5,3) * t422 - qJD(3) * t427 - t414 * t432 + t477;
t426 = -qJD(3) * mrSges(5,2) - t432 * mrSges(5,3);
t371 = -qJDD(3) * pkin(4) - t467 * pkin(8) + t433 * t415 - t373;
t367 = -0.2e1 * qJD(6) * t425 + (t424 * t431 - t394) * qJ(6) + (t425 * t431 + t393) * pkin(5) + t371;
t363 = m(7) * t367 + t393 * mrSges(7,1) - t394 * mrSges(7,3) + t424 * t403 - t425 * t406;
t470 = -m(6) * t371 - t393 * mrSges(6,1) - t394 * mrSges(6,2) - t424 * t404 - t425 * t405 - t363;
t355 = m(5) * t373 + qJDD(3) * mrSges(5,1) - t423 * mrSges(5,3) + qJD(3) * t426 - t433 * t414 + t470;
t346 = t458 * t349 + t460 * t355;
t353 = t462 * t358 + t360 * t496;
t490 = t492 * t424 + t493 * t425 + t498 * t431;
t485 = qJD(1) * t465;
t444 = (-t465 * mrSges(4,1) + t463 * mrSges(4,2)) * qJD(1);
t450 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t485;
t344 = m(4) * t401 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t446 + qJD(3) * t450 - t444 * t486 + t346;
t449 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t486;
t478 = t460 * t349 - t355 * t458;
t345 = m(4) * t402 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t447 - qJD(3) * t449 + t444 * t485 + t478;
t479 = -t344 * t463 + t465 * t345;
t352 = m(5) * t382 - t422 * mrSges(5,1) + mrSges(5,2) * t423 + t432 * t426 + t427 * t433 + t353;
t411 = -t468 * pkin(7) + t472;
t469 = -m(4) * t411 + t447 * mrSges(4,1) - mrSges(4,2) * t446 - t449 * t486 + t450 * t485 - t352;
t439 = Ifges(4,5) * qJD(3) + (t463 * Ifges(4,1) + t465 * Ifges(4,4)) * qJD(1);
t438 = Ifges(4,6) * qJD(3) + (t463 * Ifges(4,4) + t465 * Ifges(4,2)) * qJD(1);
t410 = Ifges(5,1) * t433 - Ifges(5,4) * t432 + Ifges(5,5) * qJD(3);
t409 = Ifges(5,4) * t433 - Ifges(5,2) * t432 + Ifges(5,6) * qJD(3);
t408 = Ifges(5,5) * t433 - Ifges(5,6) * t432 + Ifges(5,3) * qJD(3);
t351 = mrSges(6,2) * t371 + mrSges(7,2) * t366 - mrSges(6,3) * t368 - mrSges(7,3) * t367 - qJ(6) * t363 - t494 * t393 + t500 * t394 - t493 * t419 + t490 * t424 + t489 * t431;
t350 = -mrSges(6,1) * t371 - mrSges(7,1) * t367 + mrSges(7,2) * t365 + mrSges(6,3) * t369 - pkin(5) * t363 - t499 * t393 + t494 * t394 + t492 * t419 + t490 * t425 - t488 * t431;
t342 = -mrSges(5,1) * t382 + mrSges(5,3) * t374 + Ifges(5,4) * t423 + Ifges(5,2) * t422 + Ifges(5,6) * qJDD(3) - pkin(4) * t353 + qJD(3) * t410 - t433 * t408 - t497;
t341 = mrSges(5,2) * t382 - mrSges(5,3) * t373 + Ifges(5,1) * t423 + Ifges(5,4) * t422 + Ifges(5,5) * qJDD(3) - pkin(8) * t353 - qJD(3) * t409 - t462 * t350 + t351 * t496 - t432 * t408;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t480 - mrSges(2,2) * t474 + Ifges(3,3) * qJDD(1) + mrSges(3,1) * t420 - mrSges(3,2) * t421 + t463 * (mrSges(4,2) * t411 - mrSges(4,3) * t401 + Ifges(4,1) * t446 + Ifges(4,4) * t447 + Ifges(4,5) * qJDD(3) - qJ(4) * t346 - qJD(3) * t438 + t460 * t341 - t458 * t342) + t465 * (-mrSges(4,1) * t411 + mrSges(4,3) * t402 + Ifges(4,4) * t446 + Ifges(4,2) * t447 + Ifges(4,6) * qJDD(3) - pkin(3) * t352 + qJ(4) * t478 + qJD(3) * t439 + t458 * t341 + t460 * t342) + pkin(2) * t469 + pkin(7) * t479 + pkin(1) * (t459 * (m(3) * t421 - mrSges(3,1) * t468 - qJDD(1) * mrSges(3,2) + t479) + t461 * (m(3) * t420 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t468 + t469)); m(3) * t457 + t344 * t465 + t345 * t463; Ifges(4,5) * t446 + Ifges(4,6) * t447 + mrSges(4,1) * t401 - mrSges(4,2) * t402 + Ifges(5,5) * t423 + Ifges(5,6) * t422 + t433 * t409 + t432 * t410 + mrSges(5,1) * t373 - mrSges(5,2) * t374 + t462 * t351 + t496 * t350 + pkin(4) * t470 + pkin(8) * t477 + pkin(3) * t346 + (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + (t438 * t463 - t439 * t465) * qJD(1); t352; t497; t362;];
tauJ  = t1;
