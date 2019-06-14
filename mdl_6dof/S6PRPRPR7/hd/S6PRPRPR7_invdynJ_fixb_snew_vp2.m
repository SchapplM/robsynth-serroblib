% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRPRPR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
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
% Datum: 2019-05-04 23:22
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRPRPR7_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR7_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR7_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR7_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR7_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR7_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR7_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 23:19:07
% EndTime: 2019-05-04 23:19:09
% DurationCPUTime: 1.20s
% Computational Cost: add. (3769->213), mult. (7092->256), div. (0->0), fcn. (3992->10), ass. (0->96)
t430 = sin(qJ(4));
t433 = cos(qJ(4));
t466 = Ifges(5,4) + Ifges(6,6);
t477 = t433 * (Ifges(5,1) + Ifges(6,2)) - t430 * t466;
t476 = t433 * t466 + t430 * (-Ifges(5,2) - Ifges(6,3));
t465 = Ifges(5,5) - Ifges(6,4);
t464 = Ifges(5,6) - Ifges(6,5);
t457 = qJD(2) * qJD(4);
t454 = t433 * t457;
t401 = qJDD(2) * t430 + t454;
t416 = t430 * t457;
t402 = qJDD(2) * t433 - t416;
t425 = sin(pkin(10));
t427 = cos(pkin(10));
t405 = -g(1) * t427 - g(2) * t425;
t431 = sin(qJ(2));
t434 = cos(qJ(2));
t404 = g(1) * t425 - g(2) * t427;
t422 = -g(3) + qJDD(1);
t426 = sin(pkin(6));
t428 = cos(pkin(6));
t474 = t404 * t428 + t422 * t426;
t363 = t434 * t405 + t474 * t431;
t450 = qJDD(2) * qJ(3) + (2 * qJD(3) * qJD(2)) + t363;
t458 = t433 * qJD(2);
t469 = -2 * qJD(5);
t441 = pkin(4) * t454 + t458 * t469 + t450 + (-t402 + t416) * qJ(5);
t436 = qJD(2) ^ 2;
t468 = (-pkin(2) - pkin(8));
t455 = t468 * t436;
t352 = pkin(4) * t401 + t441 + t455;
t459 = t430 * qJD(2);
t408 = mrSges(6,1) * t459 - qJD(4) * mrSges(6,3);
t409 = mrSges(6,1) * t458 + qJD(4) * mrSges(6,2);
t362 = -t431 * t405 + t434 * t474;
t439 = -qJ(3) * t436 + qJDD(3) - t362;
t359 = qJDD(2) * t468 + t439;
t378 = -t404 * t426 + t422 * t428;
t353 = t359 * t433 - t430 * t378;
t398 = (t430 * pkin(4) - t433 * qJ(5)) * qJD(2);
t435 = qJD(4) ^ 2;
t350 = -qJDD(4) * pkin(4) - qJ(5) * t435 + t398 * t458 + qJDD(5) - t353;
t347 = (t430 * t433 * t436 - qJDD(4)) * pkin(9) + (t402 + t416) * pkin(5) + t350;
t411 = pkin(5) * t458 - qJD(4) * pkin(9);
t421 = t430 ^ 2;
t348 = -t411 * t458 + (pkin(4) + pkin(9)) * t401 + (-pkin(5) * t421 + t468) * t436 + t441;
t429 = sin(qJ(6));
t432 = cos(qJ(6));
t343 = t347 * t432 - t348 * t429;
t396 = -qJD(4) * t429 + t432 * t459;
t372 = qJD(6) * t396 + qJDD(4) * t432 + t401 * t429;
t397 = qJD(4) * t432 + t429 * t459;
t373 = -mrSges(7,1) * t396 + mrSges(7,2) * t397;
t414 = qJD(6) + t458;
t376 = -mrSges(7,2) * t414 + mrSges(7,3) * t396;
t393 = qJDD(6) + t402;
t340 = m(7) * t343 + mrSges(7,1) * t393 - t372 * mrSges(7,3) - t373 * t397 + t376 * t414;
t344 = t347 * t429 + t348 * t432;
t371 = -qJD(6) * t397 - qJDD(4) * t429 + t401 * t432;
t377 = mrSges(7,1) * t414 - mrSges(7,3) * t397;
t341 = m(7) * t344 - mrSges(7,2) * t393 + t371 * mrSges(7,3) + t373 * t396 - t377 * t414;
t453 = -t429 * t340 + t432 * t341;
t475 = (-t408 * t430 - t409 * t433) * qJD(2) + m(6) * t352 - t402 * mrSges(6,3) + t453;
t471 = (t476 * qJD(2) + t464 * qJD(4)) * t433;
t357 = t455 + t450;
t360 = (pkin(2) * t436) - t450;
t406 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t459;
t407 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t458;
t467 = mrSges(5,1) - mrSges(6,2);
t470 = -m(4) * t360 + m(5) * t357 + (t436 * mrSges(4,2)) + t402 * mrSges(5,2) + qJDD(2) * mrSges(4,3) + t401 * t467 + t406 * t459 + t407 * t458 + t475;
t354 = t430 * t359 + t433 * t378;
t460 = t477 * qJD(2) + t465 * qJD(4);
t399 = (-t430 * mrSges(6,2) - t433 * mrSges(6,3)) * qJD(2);
t452 = qJD(2) * (-t399 - (t430 * mrSges(5,1) + t433 * mrSges(5,2)) * qJD(2));
t335 = t340 * t432 + t341 * t429;
t443 = -m(6) * t350 - t402 * mrSges(6,1) - t335;
t332 = m(5) * t353 - mrSges(5,3) * t402 + t467 * qJDD(4) + (t406 - t408) * qJD(4) + t433 * t452 + t443;
t440 = -pkin(4) * t435 + qJDD(4) * qJ(5) - t398 * t459 + t354;
t349 = qJD(4) * t469 - t440;
t346 = -pkin(9) * t421 * t436 - pkin(5) * t401 + ((2 * qJD(5)) + t411) * qJD(4) + t440;
t445 = -m(7) * t346 + t371 * mrSges(7,1) - t372 * mrSges(7,2) + t396 * t376 - t397 * t377;
t438 = -m(6) * t349 + qJDD(4) * mrSges(6,3) + qJD(4) * t409 - t445;
t338 = m(5) * t354 - qJDD(4) * mrSges(5,2) - qJD(4) * t407 + (-mrSges(5,3) - mrSges(6,1)) * t401 + t430 * t452 + t438;
t449 = t332 * t433 + t338 * t430;
t361 = -qJDD(2) * pkin(2) + t439;
t444 = -m(4) * t361 + t436 * mrSges(4,3) - t449;
t365 = Ifges(7,4) * t397 + Ifges(7,2) * t396 + Ifges(7,6) * t414;
t366 = Ifges(7,1) * t397 + Ifges(7,4) * t396 + Ifges(7,5) * t414;
t442 = mrSges(7,1) * t343 - mrSges(7,2) * t344 + Ifges(7,5) * t372 + Ifges(7,6) * t371 + Ifges(7,3) * t393 + t397 * t365 - t396 * t366;
t364 = Ifges(7,5) * t397 + Ifges(7,6) * t396 + Ifges(7,3) * t414;
t337 = mrSges(7,2) * t346 - mrSges(7,3) * t343 + Ifges(7,1) * t372 + Ifges(7,4) * t371 + Ifges(7,5) * t393 + t364 * t396 - t365 * t414;
t336 = -mrSges(7,1) * t346 + mrSges(7,3) * t344 + Ifges(7,4) * t372 + Ifges(7,2) * t371 + Ifges(7,6) * t393 - t364 * t397 + t366 * t414;
t334 = qJDD(4) * mrSges(6,2) + qJD(4) * t408 + t399 * t458 - t443;
t333 = -mrSges(6,2) * t401 + t475;
t331 = qJDD(2) * mrSges(4,2) - t444;
t1 = [m(2) * t422 + t428 * (-t332 * t430 + t338 * t433 + (m(3) + m(4)) * t378) + (t431 * (m(3) * t363 - (t436 * mrSges(3,1)) - qJDD(2) * mrSges(3,2) + t470) + t434 * (m(3) * t362 - mrSges(3,2) * t436 + (mrSges(3,1) - mrSges(4,2)) * qJDD(2) + t444)) * t426; mrSges(3,1) * t362 - mrSges(3,2) * t363 + mrSges(4,2) * t361 - mrSges(4,3) * t360 + t433 * (mrSges(6,1) * t350 + mrSges(5,2) * t357 - mrSges(5,3) * t353 - mrSges(6,3) * t352 + pkin(5) * t335 - qJ(5) * t333 + t442) - t430 * (-mrSges(5,1) * t357 - mrSges(6,1) * t349 + mrSges(6,2) * t352 + mrSges(5,3) * t354 - pkin(4) * t333 - pkin(5) * t445 - pkin(9) * t453 - t432 * t336 - t429 * t337) - pkin(8) * t449 - pkin(2) * t331 + t477 * t402 + (-t430 * t464 + t433 * t465) * qJDD(4) + (Ifges(3,3) + Ifges(4,1)) * qJDD(2) + (-t430 * t460 - t471) * qJD(4) - t476 * t401 + t470 * qJ(3); t331; mrSges(5,1) * t353 - mrSges(5,2) * t354 + mrSges(6,2) * t350 - mrSges(6,3) * t349 + t432 * t337 - t429 * t336 - pkin(9) * t335 - pkin(4) * t334 + qJ(5) * t438 + t465 * t402 + (-qJ(5) * mrSges(6,1) - t464) * t401 + (Ifges(5,3) + Ifges(6,1)) * qJDD(4) + (t471 + (-qJ(5) * t399 + t460) * t430) * qJD(2); t334; t442;];
tauJ  = t1;
