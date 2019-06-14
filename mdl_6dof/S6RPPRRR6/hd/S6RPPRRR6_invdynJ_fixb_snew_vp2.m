% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPPRRR6
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
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
% Datum: 2019-05-05 15:58
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPPRRR6_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR6_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR6_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR6_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR6_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR6_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR6_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR6_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR6_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:57:17
% EndTime: 2019-05-05 15:57:19
% DurationCPUTime: 1.21s
% Computational Cost: add. (9766->232), mult. (18403->285), div. (0->0), fcn. (10576->8), ass. (0->91)
t425 = qJD(1) ^ 2;
t419 = sin(qJ(1));
t423 = cos(qJ(1));
t443 = t419 * g(1) - t423 * g(2);
t385 = -qJDD(1) * pkin(1) - t425 * qJ(2) + qJDD(2) - t443;
t379 = -qJDD(1) * qJ(3) - (2 * qJD(3) * qJD(1)) + t385;
t376 = -t425 * pkin(7) - t379;
t422 = cos(qJ(4));
t441 = qJD(1) * qJD(4);
t404 = t422 * t441;
t418 = sin(qJ(4));
t398 = -t418 * qJDD(1) - t404;
t438 = t418 * t441;
t399 = t422 * qJDD(1) - t438;
t357 = (-t399 + t438) * pkin(8) + (-t398 + t404) * pkin(4) + t376;
t435 = -t423 * g(1) - t419 * g(2);
t445 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t435;
t380 = qJDD(3) + (-pkin(1) - qJ(3)) * t425 + t445;
t377 = -qJDD(1) * pkin(7) + t380;
t371 = -t422 * g(3) + t418 * t377;
t397 = (t418 * pkin(4) - t422 * pkin(8)) * qJD(1);
t406 = t418 * qJD(1);
t424 = qJD(4) ^ 2;
t360 = -t424 * pkin(4) + qJDD(4) * pkin(8) - t397 * t406 + t371;
t417 = sin(qJ(5));
t421 = cos(qJ(5));
t342 = t421 * t357 - t417 * t360;
t442 = t422 * qJD(1);
t394 = t421 * qJD(4) - t417 * t442;
t369 = t394 * qJD(5) + t417 * qJDD(4) + t421 * t399;
t393 = qJDD(5) - t398;
t395 = t417 * qJD(4) + t421 * t442;
t403 = t406 + qJD(5);
t339 = (t394 * t403 - t369) * pkin(9) + (t394 * t395 + t393) * pkin(5) + t342;
t343 = t417 * t357 + t421 * t360;
t368 = -t395 * qJD(5) + t421 * qJDD(4) - t417 * t399;
t383 = t403 * pkin(5) - t395 * pkin(9);
t392 = t394 ^ 2;
t340 = -t392 * pkin(5) + t368 * pkin(9) - t403 * t383 + t343;
t416 = sin(qJ(6));
t420 = cos(qJ(6));
t337 = t420 * t339 - t416 * t340;
t372 = t420 * t394 - t416 * t395;
t349 = t372 * qJD(6) + t416 * t368 + t420 * t369;
t373 = t416 * t394 + t420 * t395;
t356 = -t372 * mrSges(7,1) + t373 * mrSges(7,2);
t402 = qJD(6) + t403;
t361 = -t402 * mrSges(7,2) + t372 * mrSges(7,3);
t391 = qJDD(6) + t393;
t334 = m(7) * t337 + t391 * mrSges(7,1) - t349 * mrSges(7,3) - t373 * t356 + t402 * t361;
t338 = t416 * t339 + t420 * t340;
t348 = -t373 * qJD(6) + t420 * t368 - t416 * t369;
t362 = t402 * mrSges(7,1) - t373 * mrSges(7,3);
t335 = m(7) * t338 - t391 * mrSges(7,2) + t348 * mrSges(7,3) + t372 * t356 - t402 * t362;
t327 = t420 * t334 + t416 * t335;
t375 = -t394 * mrSges(6,1) + t395 * mrSges(6,2);
t381 = -t403 * mrSges(6,2) + t394 * mrSges(6,3);
t325 = m(6) * t342 + t393 * mrSges(6,1) - t369 * mrSges(6,3) - t395 * t375 + t403 * t381 + t327;
t382 = t403 * mrSges(6,1) - t395 * mrSges(6,3);
t436 = -t416 * t334 + t420 * t335;
t326 = m(6) * t343 - t393 * mrSges(6,2) + t368 * mrSges(6,3) + t394 * t375 - t403 * t382 + t436;
t321 = t421 * t325 + t417 * t326;
t400 = -(qJD(4) * mrSges(5,2)) - mrSges(5,3) * t406;
t401 = (qJD(4) * mrSges(5,1)) - mrSges(5,3) * t442;
t446 = m(4) * t379 - m(5) * t376 + t398 * mrSges(5,1) - t399 * mrSges(5,2) - t321 + (-t400 * t418 - t401 * t422) * qJD(1);
t370 = t418 * g(3) + t422 * t377;
t396 = (t418 * mrSges(5,1) + t422 * mrSges(5,2)) * qJD(1);
t359 = -qJDD(4) * pkin(4) - t424 * pkin(8) + t397 * t442 - t370;
t341 = -t368 * pkin(5) - t392 * pkin(9) + t395 * t383 + t359;
t430 = m(7) * t341 - t348 * mrSges(7,1) + t349 * mrSges(7,2) - t372 * t361 + t373 * t362;
t427 = -m(6) * t359 + t368 * mrSges(6,1) - t369 * mrSges(6,2) + t394 * t381 - t395 * t382 - t430;
t437 = -t417 * t325 + t421 * t326;
t444 = t418 * (m(5) * t371 - qJDD(4) * mrSges(5,2) + t398 * mrSges(5,3) - qJD(4) * t401 - t396 * t406 + t437) + t422 * (m(5) * t370 + qJDD(4) * mrSges(5,1) - t399 * mrSges(5,3) + qJD(4) * t400 - t396 * t442 + t427);
t432 = m(4) * t380 + qJDD(1) * mrSges(4,2) - t425 * mrSges(4,3) + t444;
t351 = Ifges(7,4) * t373 + Ifges(7,2) * t372 + Ifges(7,6) * t402;
t352 = Ifges(7,1) * t373 + Ifges(7,4) * t372 + Ifges(7,5) * t402;
t429 = -mrSges(7,1) * t337 + mrSges(7,2) * t338 - Ifges(7,5) * t349 - Ifges(7,6) * t348 - Ifges(7,3) * t391 - t373 * t351 + t372 * t352;
t364 = Ifges(6,4) * t395 + Ifges(6,2) * t394 + Ifges(6,6) * t403;
t365 = Ifges(6,1) * t395 + Ifges(6,4) * t394 + Ifges(6,5) * t403;
t426 = mrSges(6,1) * t342 - mrSges(6,2) * t343 + Ifges(6,5) * t369 + Ifges(6,6) * t368 + Ifges(6,3) * t393 + pkin(5) * t327 + t395 * t364 - t394 * t365 - t429;
t390 = (Ifges(5,5) * qJD(4)) + (t422 * Ifges(5,1) - t418 * Ifges(5,4)) * qJD(1);
t389 = (Ifges(5,6) * qJD(4)) + (t422 * Ifges(5,4) - t418 * Ifges(5,2)) * qJD(1);
t384 = t425 * pkin(1) - t445;
t363 = Ifges(6,5) * t395 + Ifges(6,6) * t394 + Ifges(6,3) * t403;
t350 = Ifges(7,5) * t373 + Ifges(7,6) * t372 + Ifges(7,3) * t402;
t329 = mrSges(7,2) * t341 - mrSges(7,3) * t337 + Ifges(7,1) * t349 + Ifges(7,4) * t348 + Ifges(7,5) * t391 + t372 * t350 - t402 * t351;
t328 = -mrSges(7,1) * t341 + mrSges(7,3) * t338 + Ifges(7,4) * t349 + Ifges(7,2) * t348 + Ifges(7,6) * t391 - t373 * t350 + t402 * t352;
t319 = m(3) * t385 + (-mrSges(4,2) - mrSges(3,3)) * t425 + (mrSges(3,2) - mrSges(4,3)) * qJDD(1) + t446;
t318 = mrSges(6,2) * t359 - mrSges(6,3) * t342 + Ifges(6,1) * t369 + Ifges(6,4) * t368 + Ifges(6,5) * t393 - pkin(9) * t327 - t416 * t328 + t420 * t329 + t394 * t363 - t403 * t364;
t317 = -mrSges(6,1) * t359 + mrSges(6,3) * t343 + Ifges(6,4) * t369 + Ifges(6,2) * t368 + Ifges(6,6) * t393 - pkin(5) * t430 + pkin(9) * t436 + t420 * t328 + t416 * t329 - t395 * t363 + t403 * t365;
t1 = [qJ(2) * (-m(3) * t384 + t425 * mrSges(3,2) + t432) - pkin(1) * t319 + mrSges(2,1) * t443 - mrSges(2,2) * t435 + mrSges(3,2) * t385 - mrSges(3,3) * t384 + t422 * (mrSges(5,2) * t376 - mrSges(5,3) * t370 + Ifges(5,1) * t399 + Ifges(5,4) * t398 + Ifges(5,5) * qJDD(4) - pkin(8) * t321 - qJD(4) * t389 - t417 * t317 + t421 * t318) - t418 * (-mrSges(5,1) * t376 + mrSges(5,3) * t371 + Ifges(5,4) * t399 + Ifges(5,2) * t398 + Ifges(5,6) * qJDD(4) - pkin(4) * t321 + qJD(4) * t390 - t426) - pkin(7) * t444 + mrSges(4,2) * t380 - mrSges(4,3) * t379 + (qJ(2) * mrSges(3,3) + Ifges(3,1) + Ifges(4,1) + Ifges(2,3)) * qJDD(1) + (t425 * mrSges(4,2) + qJDD(1) * mrSges(4,3) - t446) * qJ(3); t319; t432; Ifges(5,5) * t399 + Ifges(5,6) * t398 + Ifges(5,3) * qJDD(4) + mrSges(5,1) * t370 - mrSges(5,2) * t371 + t417 * t318 + t421 * t317 + pkin(4) * t427 + pkin(8) * t437 + (t422 * t389 + t418 * t390) * qJD(1); t426; -t429;];
tauJ  = t1;
