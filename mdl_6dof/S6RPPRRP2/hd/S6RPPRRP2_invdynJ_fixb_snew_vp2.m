% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPPRRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-05-05 14:49
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPPRRP2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP2_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP2_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP2_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP2_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP2_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP2_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP2_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:47:15
% EndTime: 2019-05-05 14:47:17
% DurationCPUTime: 1.87s
% Computational Cost: add. (13352->242), mult. (29470->295), div. (0->0), fcn. (20160->10), ass. (0->110)
t499 = Ifges(6,1) + Ifges(7,1);
t492 = Ifges(6,4) - Ifges(7,5);
t491 = -Ifges(6,5) - Ifges(7,4);
t498 = Ifges(6,2) + Ifges(7,3);
t490 = Ifges(6,6) - Ifges(7,6);
t497 = -Ifges(6,3) - Ifges(7,2);
t460 = qJD(1) ^ 2;
t450 = sin(pkin(10));
t452 = cos(pkin(10));
t455 = sin(qJ(4));
t457 = cos(qJ(4));
t466 = t450 * t455 - t452 * t457;
t428 = t466 * qJD(1);
t467 = t450 * t457 + t452 * t455;
t429 = t467 * qJD(1);
t478 = qJD(4) * t429;
t418 = -qJDD(1) * t466 - t478;
t479 = qJD(4) * t428;
t419 = qJDD(1) * t467 - t479;
t454 = sin(qJ(5));
t495 = cos(qJ(5));
t422 = -qJD(4) * t495 + t429 * t454;
t390 = -qJD(5) * t422 + qJDD(4) * t454 + t419 * t495;
t423 = qJD(4) * t454 + t429 * t495;
t398 = mrSges(7,1) * t422 - mrSges(7,3) * t423;
t456 = sin(qJ(1));
t458 = cos(qJ(1));
t474 = g(1) * t456 - g(2) * t458;
t435 = qJDD(1) * pkin(1) + t474;
t469 = -g(1) * t458 - g(2) * t456;
t436 = -pkin(1) * t460 + t469;
t451 = sin(pkin(9));
t453 = cos(pkin(9));
t421 = t435 * t451 + t436 * t453;
t410 = -pkin(2) * t460 + qJDD(1) * qJ(3) + t421;
t449 = -g(3) + qJDD(2);
t477 = qJD(1) * qJD(3);
t481 = t449 * t452 - 0.2e1 * t450 * t477;
t494 = pkin(3) * t452;
t388 = (-pkin(7) * qJDD(1) + t460 * t494 - t410) * t450 + t481;
t396 = t450 * t449 + (t410 + 0.2e1 * t477) * t452;
t476 = qJDD(1) * t452;
t448 = t452 ^ 2;
t487 = t448 * t460;
t391 = -pkin(3) * t487 + pkin(7) * t476 + t396;
t374 = t388 * t455 + t391 * t457;
t417 = pkin(4) * t428 - pkin(8) * t429;
t459 = qJD(4) ^ 2;
t370 = -pkin(4) * t459 + qJDD(4) * pkin(8) - t417 * t428 + t374;
t447 = t450 ^ 2;
t420 = t435 * t453 - t436 * t451;
t468 = qJDD(3) - t420;
t394 = (-pkin(2) - t494) * qJDD(1) + (-qJ(3) + (-t447 - t448) * pkin(7)) * t460 + t468;
t372 = (-t419 + t479) * pkin(8) + (-t418 + t478) * pkin(4) + t394;
t366 = -t370 * t454 + t372 * t495;
t397 = pkin(5) * t422 - qJ(6) * t423;
t416 = qJDD(5) - t418;
t427 = qJD(5) + t428;
t426 = t427 ^ 2;
t364 = -pkin(5) * t416 - qJ(6) * t426 + t397 * t423 + qJDD(6) - t366;
t401 = -mrSges(7,2) * t422 + mrSges(7,3) * t427;
t470 = -m(7) * t364 + mrSges(7,1) * t416 + t401 * t427;
t360 = mrSges(7,2) * t390 + t398 * t423 - t470;
t367 = t370 * t495 + t372 * t454;
t363 = -pkin(5) * t426 + qJ(6) * t416 + 0.2e1 * qJD(6) * t427 - t397 * t422 + t367;
t389 = qJD(5) * t423 - qJDD(4) * t495 + t419 * t454;
t404 = -mrSges(7,1) * t427 + mrSges(7,2) * t423;
t475 = m(7) * t363 + mrSges(7,3) * t416 + t404 * t427;
t483 = t422 * t492 - t423 * t499 + t427 * t491;
t484 = t422 * t498 - t423 * t492 - t427 * t490;
t496 = -t490 * t389 - t491 * t390 - t497 * t416 - t483 * t422 - t484 * t423 + mrSges(6,1) * t366 - mrSges(7,1) * t364 - mrSges(6,2) * t367 + mrSges(7,3) * t363 - pkin(5) * t360 + qJ(6) * (-mrSges(7,2) * t389 - t398 * t422 + t475);
t493 = -mrSges(6,3) - mrSges(7,2);
t488 = mrSges(4,2) * t450;
t414 = mrSges(5,1) * t428 + mrSges(5,2) * t429;
t425 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t429;
t403 = mrSges(6,1) * t427 - mrSges(6,3) * t423;
t482 = -mrSges(6,1) * t422 - mrSges(6,2) * t423 - t398;
t356 = m(6) * t367 - mrSges(6,2) * t416 + t389 * t493 - t403 * t427 + t422 * t482 + t475;
t402 = -mrSges(6,2) * t427 - mrSges(6,3) * t422;
t358 = m(6) * t366 + mrSges(6,1) * t416 + t390 * t493 + t402 * t427 + t423 * t482 + t470;
t471 = t356 * t495 - t358 * t454;
t348 = m(5) * t374 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t418 - qJD(4) * t425 - t414 * t428 + t471;
t373 = t388 * t457 - t391 * t455;
t424 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t428;
t369 = -qJDD(4) * pkin(4) - pkin(8) * t459 + t417 * t429 - t373;
t365 = -0.2e1 * qJD(6) * t423 + (t422 * t427 - t390) * qJ(6) + (t423 * t427 + t389) * pkin(5) + t369;
t361 = m(7) * t365 + mrSges(7,1) * t389 - mrSges(7,3) * t390 + t401 * t422 - t404 * t423;
t461 = -m(6) * t369 - mrSges(6,1) * t389 - mrSges(6,2) * t390 - t402 * t422 - t403 * t423 - t361;
t353 = m(5) * t373 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t419 + qJD(4) * t424 - t414 * t429 + t461;
t486 = t348 * t455 + t353 * t457;
t351 = t356 * t454 + t358 * t495;
t485 = t422 * t490 + t423 * t491 + t427 * t497;
t395 = -t410 * t450 + t481;
t465 = mrSges(4,3) * qJDD(1) + t460 * (-mrSges(4,1) * t452 + t488);
t343 = m(4) * t395 - t450 * t465 + t486;
t472 = t348 * t457 - t353 * t455;
t344 = m(4) * t396 + t452 * t465 + t472;
t473 = -t343 * t450 + t344 * t452;
t464 = m(5) * t394 - mrSges(5,1) * t418 + mrSges(5,2) * t419 + t424 * t428 + t425 * t429 + t351;
t406 = -qJDD(1) * pkin(2) - qJ(3) * t460 + t468;
t462 = -m(4) * t406 + mrSges(4,1) * t476 - t464 + (t447 * t460 + t487) * mrSges(4,3);
t409 = Ifges(5,1) * t429 - Ifges(5,4) * t428 + Ifges(5,5) * qJD(4);
t408 = Ifges(5,4) * t429 - Ifges(5,2) * t428 + Ifges(5,6) * qJD(4);
t407 = Ifges(5,5) * t429 - Ifges(5,6) * t428 + Ifges(5,3) * qJD(4);
t350 = mrSges(6,2) * t369 + mrSges(7,2) * t364 - mrSges(6,3) * t366 - mrSges(7,3) * t365 - qJ(6) * t361 - t389 * t492 + t390 * t499 - t416 * t491 + t422 * t485 + t427 * t484;
t349 = -mrSges(6,1) * t369 - mrSges(7,1) * t365 + mrSges(7,2) * t363 + mrSges(6,3) * t367 - pkin(5) * t361 - t389 * t498 + t390 * t492 + t416 * t490 + t423 * t485 - t427 * t483;
t345 = qJDD(1) * t488 - t462;
t341 = -mrSges(5,1) * t394 + mrSges(5,3) * t374 + Ifges(5,4) * t419 + Ifges(5,2) * t418 + Ifges(5,6) * qJDD(4) - pkin(4) * t351 + qJD(4) * t409 - t429 * t407 - t496;
t340 = mrSges(5,2) * t394 - mrSges(5,3) * t373 + Ifges(5,1) * t419 + Ifges(5,4) * t418 + Ifges(5,5) * qJDD(4) - pkin(8) * t351 - qJD(4) * t408 - t349 * t454 + t350 * t495 - t407 * t428;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t474 - mrSges(2,2) * t469 + Ifges(3,3) * qJDD(1) + mrSges(3,1) * t420 - mrSges(3,2) * t421 + t450 * (mrSges(4,2) * t406 - mrSges(4,3) * t395 + t457 * t340 - t455 * t341 - pkin(7) * t486 + (Ifges(4,1) * t450 + Ifges(4,4) * t452) * qJDD(1)) + t452 * (-mrSges(4,1) * t406 + mrSges(4,3) * t396 + t455 * t340 + t457 * t341 - pkin(3) * t464 + pkin(7) * t472 + (Ifges(4,4) * t450 + Ifges(4,2) * t452) * qJDD(1)) - pkin(2) * t345 + qJ(3) * t473 + pkin(1) * (t451 * (m(3) * t421 - mrSges(3,1) * t460 - qJDD(1) * mrSges(3,2) + t473) + t453 * (t462 + (mrSges(3,1) - t488) * qJDD(1) + m(3) * t420 - mrSges(3,2) * t460)); m(3) * t449 + t343 * t452 + t344 * t450; t345; mrSges(5,1) * t373 - mrSges(5,2) * t374 + Ifges(5,5) * t419 + Ifges(5,6) * t418 + Ifges(5,3) * qJDD(4) + pkin(4) * t461 + pkin(8) * t471 + t349 * t495 + t454 * t350 + t429 * t408 + t428 * t409; t496; t360;];
tauJ  = t1;
