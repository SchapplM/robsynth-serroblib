% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRRPR11
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
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
% tauJ [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:36
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRRPR11_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR11_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR11_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR11_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR11_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR11_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR11_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR11_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR11_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:32:41
% EndTime: 2019-12-31 21:32:44
% DurationCPUTime: 1.65s
% Computational Cost: add. (10526->258), mult. (20648->313), div. (0->0), fcn. (13035->8), ass. (0->105)
t455 = Ifges(4,1) + Ifges(5,1);
t447 = Ifges(4,4) - Ifges(5,5);
t446 = Ifges(4,5) + Ifges(5,4);
t454 = -Ifges(4,2) - Ifges(5,3);
t445 = Ifges(4,6) - Ifges(5,6);
t453 = Ifges(4,3) + Ifges(5,2);
t414 = sin(qJ(3));
t415 = sin(qJ(2));
t438 = qJD(1) * t415;
t449 = cos(qJ(3));
t396 = -t449 * qJD(2) + t414 * t438;
t418 = cos(qJ(2));
t436 = qJD(1) * qJD(2);
t434 = t418 * t436;
t400 = qJDD(1) * t415 + t434;
t366 = -t396 * qJD(3) + t414 * qJDD(2) + t449 * t400;
t421 = qJD(1) ^ 2;
t416 = sin(qJ(1));
t419 = cos(qJ(1));
t429 = -g(1) * t419 - g(2) * t416;
t389 = -pkin(1) * t421 + qJDD(1) * pkin(6) + t429;
t378 = -t418 * g(3) - t415 * t389;
t399 = (-pkin(2) * t418 - pkin(7) * t415) * qJD(1);
t420 = qJD(2) ^ 2;
t426 = qJDD(2) * pkin(2) + pkin(7) * t420 - t399 * t438 + t378;
t437 = qJD(1) * t418;
t407 = qJD(3) - t437;
t443 = t396 * t407;
t452 = (-t366 + t443) * qJ(4) - t426;
t397 = t414 * qJD(2) + t449 * t438;
t372 = mrSges(5,1) * t396 - mrSges(5,3) * t397;
t433 = g(1) * t416 - t419 * g(2);
t388 = -qJDD(1) * pkin(1) - pkin(6) * t421 - t433;
t435 = t415 * t436;
t401 = t418 * qJDD(1) - t435;
t346 = (-t400 - t434) * pkin(7) + (-t401 + t435) * pkin(2) + t388;
t379 = -g(3) * t415 + t418 * t389;
t350 = -pkin(2) * t420 + qJDD(2) * pkin(7) + t399 * t437 + t379;
t337 = t449 * t346 - t414 * t350;
t371 = pkin(3) * t396 - qJ(4) * t397;
t395 = qJDD(3) - t401;
t406 = t407 ^ 2;
t329 = -t395 * pkin(3) - t406 * qJ(4) + t397 * t371 + qJDD(4) - t337;
t322 = (-t366 - t443) * pkin(8) + (t396 * t397 - t395) * pkin(4) + t329;
t338 = t414 * t346 + t449 * t350;
t450 = 2 * qJD(4);
t327 = -pkin(3) * t406 + t395 * qJ(4) - t396 * t371 + t407 * t450 + t338;
t365 = qJD(3) * t397 - t449 * qJDD(2) + t400 * t414;
t380 = -pkin(4) * t407 - pkin(8) * t397;
t394 = t396 ^ 2;
t324 = -pkin(4) * t394 + pkin(8) * t365 + t380 * t407 + t327;
t413 = sin(qJ(5));
t417 = cos(qJ(5));
t320 = t322 * t417 - t324 * t413;
t367 = t396 * t417 - t397 * t413;
t336 = qJD(5) * t367 + t365 * t413 + t366 * t417;
t368 = t396 * t413 + t397 * t417;
t344 = -mrSges(6,1) * t367 + mrSges(6,2) * t368;
t405 = qJD(5) - t407;
t351 = -mrSges(6,2) * t405 + mrSges(6,3) * t367;
t393 = qJDD(5) - t395;
t317 = m(6) * t320 + mrSges(6,1) * t393 - mrSges(6,3) * t336 - t344 * t368 + t351 * t405;
t321 = t322 * t413 + t324 * t417;
t335 = -qJD(5) * t368 + t365 * t417 - t366 * t413;
t352 = mrSges(6,1) * t405 - mrSges(6,3) * t368;
t318 = m(6) * t321 - mrSges(6,2) * t393 + mrSges(6,3) * t335 + t344 * t367 - t352 * t405;
t312 = t317 * t417 + t318 * t413;
t377 = -mrSges(5,2) * t396 + mrSges(5,3) * t407;
t425 = -m(5) * t329 + t395 * mrSges(5,1) + t407 * t377 - t312;
t311 = mrSges(5,2) * t366 + t372 * t397 - t425;
t340 = Ifges(6,4) * t368 + Ifges(6,2) * t367 + Ifges(6,6) * t405;
t341 = Ifges(6,1) * t368 + Ifges(6,4) * t367 + Ifges(6,5) * t405;
t424 = -mrSges(6,1) * t320 + mrSges(6,2) * t321 - Ifges(6,5) * t336 - Ifges(6,6) * t335 - Ifges(6,3) * t393 - t368 * t340 + t367 * t341;
t376 = -mrSges(5,1) * t407 + mrSges(5,2) * t397;
t431 = -t317 * t413 + t417 * t318;
t427 = m(5) * t327 + t395 * mrSges(5,3) + t407 * t376 + t431;
t440 = -t447 * t396 + t455 * t397 + t446 * t407;
t441 = t454 * t396 + t447 * t397 + t445 * t407;
t451 = -t445 * t365 + t446 * t366 + t453 * t395 + t440 * t396 + t441 * t397 + mrSges(4,1) * t337 - mrSges(5,1) * t329 - mrSges(4,2) * t338 + mrSges(5,3) * t327 - pkin(3) * t311 - pkin(4) * t312 + qJ(4) * (-mrSges(5,2) * t365 - t372 * t396 + t427) + t424;
t448 = -mrSges(4,3) - mrSges(5,2);
t442 = t445 * t396 - t446 * t397 - t453 * t407;
t439 = -mrSges(4,1) * t396 - mrSges(4,2) * t397 - t372;
t375 = mrSges(4,1) * t407 - mrSges(4,3) * t397;
t308 = m(4) * t338 - mrSges(4,2) * t395 + t448 * t365 - t375 * t407 + t439 * t396 + t427;
t374 = -mrSges(4,2) * t407 - mrSges(4,3) * t396;
t309 = m(4) * t337 + mrSges(4,1) * t395 + t448 * t366 + t374 * t407 + t439 * t397 + t425;
t432 = t449 * t308 - t309 * t414;
t325 = -pkin(8) * t394 + (-pkin(3) - pkin(4)) * t365 + (-pkin(3) * t407 + t380 + t450) * t397 - t452;
t428 = -m(6) * t325 + t335 * mrSges(6,1) - t336 * mrSges(6,2) + t367 * t351 - t368 * t352;
t306 = t414 * t308 + t449 * t309;
t328 = -0.2e1 * qJD(4) * t397 + (t397 * t407 + t365) * pkin(3) + t452;
t315 = m(5) * t328 + t365 * mrSges(5,1) - t366 * mrSges(5,3) - t397 * t376 + t396 * t377 + t428;
t423 = m(4) * t426 - t365 * mrSges(4,1) - t366 * mrSges(4,2) - t396 * t374 - t397 * t375 - t315;
t403 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t437;
t402 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t438;
t398 = (-mrSges(3,1) * t418 + mrSges(3,2) * t415) * qJD(1);
t387 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t415 + Ifges(3,4) * t418) * qJD(1);
t386 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t415 + Ifges(3,2) * t418) * qJD(1);
t385 = Ifges(3,3) * qJD(2) + (t415 * Ifges(3,5) + t418 * Ifges(3,6)) * qJD(1);
t339 = Ifges(6,5) * t368 + Ifges(6,6) * t367 + Ifges(6,3) * t405;
t314 = mrSges(6,2) * t325 - mrSges(6,3) * t320 + Ifges(6,1) * t336 + Ifges(6,4) * t335 + Ifges(6,5) * t393 + t339 * t367 - t340 * t405;
t313 = -mrSges(6,1) * t325 + mrSges(6,3) * t321 + Ifges(6,4) * t336 + Ifges(6,2) * t335 + Ifges(6,6) * t393 - t339 * t368 + t341 * t405;
t305 = -mrSges(4,2) * t426 + mrSges(5,2) * t329 - mrSges(4,3) * t337 - mrSges(5,3) * t328 - pkin(8) * t312 - qJ(4) * t315 - t313 * t413 + t314 * t417 - t447 * t365 + t455 * t366 + t446 * t395 + t442 * t396 - t441 * t407;
t304 = mrSges(4,1) * t426 - mrSges(5,1) * t328 + mrSges(5,2) * t327 + mrSges(4,3) * t338 - pkin(3) * t315 - pkin(4) * t428 - pkin(8) * t431 - t417 * t313 - t413 * t314 + t454 * t365 + t447 * t366 + t445 * t395 + t442 * t397 + t440 * t407;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t433 - mrSges(2,2) * t429 + t415 * (mrSges(3,2) * t388 - mrSges(3,3) * t378 + Ifges(3,1) * t400 + Ifges(3,4) * t401 + Ifges(3,5) * qJDD(2) - pkin(7) * t306 - qJD(2) * t386 - t414 * t304 + t449 * t305 + t385 * t437) + t418 * (-mrSges(3,1) * t388 + mrSges(3,3) * t379 + Ifges(3,4) * t400 + Ifges(3,2) * t401 + Ifges(3,6) * qJDD(2) - pkin(2) * t306 + qJD(2) * t387 - t385 * t438 - t451) + pkin(1) * (-m(3) * t388 + t401 * mrSges(3,1) - t400 * mrSges(3,2) + (-t402 * t415 + t403 * t418) * qJD(1) - t306) + pkin(6) * (t418 * (m(3) * t379 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t401 - qJD(2) * t402 + t398 * t437 + t432) - t415 * (m(3) * t378 + qJDD(2) * mrSges(3,1) - t400 * mrSges(3,3) + qJD(2) * t403 - t398 * t438 + t423)); Ifges(3,5) * t400 + Ifges(3,6) * t401 + Ifges(3,3) * qJDD(2) + mrSges(3,1) * t378 - mrSges(3,2) * t379 + t414 * t305 + t449 * t304 + pkin(2) * t423 + pkin(7) * t432 + (t415 * t386 - t418 * t387) * qJD(1); t451; t311; -t424;];
tauJ = t1;
