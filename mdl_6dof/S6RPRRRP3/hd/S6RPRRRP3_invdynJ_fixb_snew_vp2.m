% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRRRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-05-06 01:22
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRRRP3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP3_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP3_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP3_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP3_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP3_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP3_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP3_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:18:05
% EndTime: 2019-05-06 01:18:10
% DurationCPUTime: 2.58s
% Computational Cost: add. (20730->269), mult. (39893->327), div. (0->0), fcn. (26099->10), ass. (0->112)
t491 = Ifges(6,1) + Ifges(7,1);
t482 = Ifges(6,4) - Ifges(7,5);
t489 = Ifges(7,4) + Ifges(6,5);
t490 = Ifges(6,2) + Ifges(7,3);
t488 = Ifges(6,6) - Ifges(7,6);
t487 = -Ifges(6,3) - Ifges(7,2);
t456 = sin(qJ(4));
t459 = cos(qJ(4));
t457 = sin(qJ(3));
t479 = qJD(1) * t457;
t436 = qJD(3) * t459 - t456 * t479;
t437 = qJD(3) * t456 + t459 * t479;
t455 = sin(qJ(5));
t484 = cos(qJ(5));
t411 = -t436 * t484 + t455 * t437;
t412 = t455 * t436 + t437 * t484;
t460 = cos(qJ(3));
t478 = qJD(1) * t460;
t448 = qJD(4) - t478;
t447 = qJD(5) + t448;
t486 = t490 * t411 - t482 * t412 - t488 * t447;
t485 = -t482 * t411 + t491 * t412 + t489 * t447;
t483 = -mrSges(6,3) - mrSges(7,2);
t458 = sin(qJ(1));
t461 = cos(qJ(1));
t474 = t458 * g(1) - g(2) * t461;
t438 = qJDD(1) * pkin(1) + t474;
t463 = qJD(1) ^ 2;
t469 = -g(1) * t461 - g(2) * t458;
t440 = -pkin(1) * t463 + t469;
t453 = sin(pkin(10));
t454 = cos(pkin(10));
t413 = t454 * t438 - t453 * t440;
t401 = -qJDD(1) * pkin(2) - t463 * pkin(7) - t413;
t477 = qJD(1) * qJD(3);
t475 = t460 * t477;
t442 = qJDD(1) * t457 + t475;
t450 = t457 * t477;
t443 = qJDD(1) * t460 - t450;
t380 = (-t442 - t475) * pkin(8) + (-t443 + t450) * pkin(3) + t401;
t414 = t453 * t438 + t454 * t440;
t402 = -pkin(2) * t463 + qJDD(1) * pkin(7) + t414;
t452 = -g(3) + qJDD(2);
t393 = t460 * t402 + t457 * t452;
t441 = (-pkin(3) * t460 - pkin(8) * t457) * qJD(1);
t462 = qJD(3) ^ 2;
t386 = -pkin(3) * t462 + qJDD(3) * pkin(8) + t441 * t478 + t393;
t359 = t459 * t380 - t456 * t386;
t410 = qJD(4) * t436 + qJDD(3) * t456 + t442 * t459;
t435 = qJDD(4) - t443;
t355 = (t436 * t448 - t410) * pkin(9) + (t436 * t437 + t435) * pkin(4) + t359;
t360 = t456 * t380 + t459 * t386;
t409 = -qJD(4) * t437 + qJDD(3) * t459 - t442 * t456;
t418 = pkin(4) * t448 - pkin(9) * t437;
t434 = t436 ^ 2;
t357 = -pkin(4) * t434 + pkin(9) * t409 - t418 * t448 + t360;
t353 = t455 * t355 + t357 * t484;
t370 = t412 * qJD(5) - t409 * t484 + t455 * t410;
t396 = mrSges(6,1) * t447 - mrSges(6,3) * t412;
t431 = qJDD(5) + t435;
t387 = pkin(5) * t411 - qJ(6) * t412;
t446 = t447 ^ 2;
t347 = -pkin(5) * t446 + qJ(6) * t431 + 0.2e1 * qJD(6) * t447 - t387 * t411 + t353;
t397 = -mrSges(7,1) * t447 + mrSges(7,2) * t412;
t476 = m(7) * t347 + t431 * mrSges(7,3) + t447 * t397;
t388 = mrSges(7,1) * t411 - mrSges(7,3) * t412;
t480 = -mrSges(6,1) * t411 - mrSges(6,2) * t412 - t388;
t337 = m(6) * t353 - t431 * mrSges(6,2) + t370 * t483 - t447 * t396 + t411 * t480 + t476;
t352 = t355 * t484 - t455 * t357;
t371 = -t411 * qJD(5) + t455 * t409 + t410 * t484;
t395 = -mrSges(6,2) * t447 - mrSges(6,3) * t411;
t348 = -t431 * pkin(5) - t446 * qJ(6) + t412 * t387 + qJDD(6) - t352;
t394 = -mrSges(7,2) * t411 + mrSges(7,3) * t447;
t470 = -m(7) * t348 + t431 * mrSges(7,1) + t447 * t394;
t339 = m(6) * t352 + t431 * mrSges(6,1) + t371 * t483 + t447 * t395 + t412 * t480 + t470;
t333 = t455 * t337 + t339 * t484;
t481 = t488 * t411 - t489 * t412 + t487 * t447;
t439 = (-mrSges(4,1) * t460 + mrSges(4,2) * t457) * qJD(1);
t444 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t479;
t415 = -mrSges(5,1) * t436 + mrSges(5,2) * t437;
t416 = -mrSges(5,2) * t448 + mrSges(5,3) * t436;
t329 = m(5) * t359 + mrSges(5,1) * t435 - mrSges(5,3) * t410 - t415 * t437 + t416 * t448 + t333;
t417 = mrSges(5,1) * t448 - mrSges(5,3) * t437;
t471 = t337 * t484 - t339 * t455;
t330 = m(5) * t360 - mrSges(5,2) * t435 + mrSges(5,3) * t409 + t415 * t436 - t417 * t448 + t471;
t472 = -t329 * t456 + t459 * t330;
t326 = m(4) * t393 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t443 - qJD(3) * t444 + t439 * t478 + t472;
t392 = -t457 * t402 + t460 * t452;
t445 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t478;
t385 = -qJDD(3) * pkin(3) - t462 * pkin(8) + t441 * t479 - t392;
t358 = -t409 * pkin(4) - t434 * pkin(9) + t437 * t418 + t385;
t350 = -0.2e1 * qJD(6) * t412 + (t411 * t447 - t371) * qJ(6) + (t412 * t447 + t370) * pkin(5) + t358;
t344 = m(7) * t350 + t370 * mrSges(7,1) - t371 * mrSges(7,3) + t411 * t394 - t412 * t397;
t468 = m(6) * t358 + t370 * mrSges(6,1) + t371 * mrSges(6,2) + t411 * t395 + t412 * t396 + t344;
t465 = -m(5) * t385 + t409 * mrSges(5,1) - t410 * mrSges(5,2) + t436 * t416 - t437 * t417 - t468;
t334 = m(4) * t392 + qJDD(3) * mrSges(4,1) - t442 * mrSges(4,3) + qJD(3) * t445 - t439 * t479 + t465;
t473 = t460 * t326 - t334 * t457;
t327 = t329 * t459 + t330 * t456;
t467 = -m(4) * t401 + t443 * mrSges(4,1) - mrSges(4,2) * t442 - t444 * t479 + t445 * t478 - t327;
t343 = t371 * mrSges(7,2) + t412 * t388 - t470;
t466 = -mrSges(6,1) * t352 + mrSges(7,1) * t348 + mrSges(6,2) * t353 - mrSges(7,3) * t347 + pkin(5) * t343 - qJ(6) * t476 + t487 * t431 + t486 * t412 + (qJ(6) * t388 - t485) * t411 - t489 * t371 + (mrSges(7,2) * qJ(6) + t488) * t370;
t404 = Ifges(5,4) * t437 + Ifges(5,2) * t436 + Ifges(5,6) * t448;
t405 = Ifges(5,1) * t437 + Ifges(5,4) * t436 + Ifges(5,5) * t448;
t464 = mrSges(5,1) * t359 - mrSges(5,2) * t360 + Ifges(5,5) * t410 + Ifges(5,6) * t409 + Ifges(5,3) * t435 + pkin(4) * t333 + t437 * t404 - t436 * t405 - t466;
t430 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t457 + Ifges(4,4) * t460) * qJD(1);
t429 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t457 + Ifges(4,2) * t460) * qJD(1);
t403 = Ifges(5,5) * t437 + Ifges(5,6) * t436 + Ifges(5,3) * t448;
t332 = mrSges(6,2) * t358 + mrSges(7,2) * t348 - mrSges(6,3) * t352 - mrSges(7,3) * t350 - qJ(6) * t344 - t482 * t370 + t491 * t371 + t481 * t411 + t489 * t431 + t486 * t447;
t331 = -mrSges(6,1) * t358 - mrSges(7,1) * t350 + mrSges(7,2) * t347 + mrSges(6,3) * t353 - pkin(5) * t344 - t490 * t370 + t482 * t371 + t481 * t412 + t488 * t431 + t485 * t447;
t324 = mrSges(5,2) * t385 - mrSges(5,3) * t359 + Ifges(5,1) * t410 + Ifges(5,4) * t409 + Ifges(5,5) * t435 - pkin(9) * t333 - t455 * t331 + t332 * t484 + t436 * t403 - t448 * t404;
t323 = -mrSges(5,1) * t385 + mrSges(5,3) * t360 + Ifges(5,4) * t410 + Ifges(5,2) * t409 + Ifges(5,6) * t435 - pkin(4) * t468 + pkin(9) * t471 + t331 * t484 + t455 * t332 - t437 * t403 + t448 * t405;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t474 - mrSges(2,2) * t469 + Ifges(3,3) * qJDD(1) + mrSges(3,1) * t413 - mrSges(3,2) * t414 + t457 * (mrSges(4,2) * t401 - mrSges(4,3) * t392 + Ifges(4,1) * t442 + Ifges(4,4) * t443 + Ifges(4,5) * qJDD(3) - pkin(8) * t327 - qJD(3) * t429 - t456 * t323 + t459 * t324) + t460 * (-mrSges(4,1) * t401 + mrSges(4,3) * t393 + Ifges(4,4) * t442 + Ifges(4,2) * t443 + Ifges(4,6) * qJDD(3) - pkin(3) * t327 + qJD(3) * t430 - t464) + pkin(2) * t467 + pkin(7) * t473 + pkin(1) * (t453 * (m(3) * t414 - mrSges(3,1) * t463 - qJDD(1) * mrSges(3,2) + t473) + t454 * (m(3) * t413 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t463 + t467)); m(3) * t452 + t326 * t457 + t334 * t460; Ifges(4,5) * t442 + Ifges(4,6) * t443 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t392 - mrSges(4,2) * t393 + t456 * t324 + t459 * t323 + pkin(3) * t465 + pkin(8) * t472 + (t429 * t457 - t430 * t460) * qJD(1); t464; -t466; t343;];
tauJ  = t1;
