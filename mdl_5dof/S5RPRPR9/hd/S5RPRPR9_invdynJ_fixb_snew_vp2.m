% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPRPR9
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% Datum: 2019-12-31 18:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPRPR9_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR9_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR9_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR9_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR9_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR9_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR9_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR9_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR9_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:23:50
% EndTime: 2019-12-31 18:23:52
% DurationCPUTime: 1.11s
% Computational Cost: add. (2821->202), mult. (5499->244), div. (0->0), fcn. (2772->8), ass. (0->92)
t357 = sin(qJ(3));
t360 = cos(qJ(3));
t389 = Ifges(4,4) + Ifges(5,6);
t401 = t357 * t389 + t360 * (Ifges(4,2) + Ifges(5,3));
t400 = t357 * (Ifges(4,1) + Ifges(5,2)) + t360 * t389;
t388 = Ifges(4,5) - Ifges(5,4);
t387 = Ifges(4,6) - Ifges(5,5);
t355 = cos(pkin(8));
t397 = pkin(1) * t355 + pkin(2);
t396 = (t401 * qJD(1) + t387 * qJD(3)) * t357 - t360 * (t400 * qJD(1) + t388 * qJD(3));
t395 = -2 * qJD(4);
t394 = -pkin(3) - pkin(7);
t363 = qJD(1) ^ 2;
t392 = pkin(6) * t363;
t391 = pkin(7) * t363;
t353 = -g(3) + qJDD(2);
t386 = t353 * t360;
t380 = qJD(1) * qJD(3);
t378 = t357 * t380;
t338 = qJDD(1) * t360 - t378;
t382 = qJD(1) * t357;
t344 = pkin(4) * t382 - qJD(3) * pkin(7);
t352 = t360 ^ 2;
t379 = t360 * t380;
t337 = qJDD(1) * t357 + t379;
t358 = sin(qJ(1));
t361 = cos(qJ(1));
t377 = t358 * g(1) - g(2) * t361;
t332 = qJDD(1) * pkin(1) + t377;
t375 = -g(1) * t361 - g(2) * t358;
t336 = -pkin(1) * t363 + t375;
t354 = sin(pkin(8));
t307 = t332 * t355 - t354 * t336;
t372 = -qJDD(1) * pkin(2) - t307;
t364 = pkin(3) * t378 + t382 * t395 + (-t337 - t379) * qJ(4) + t372;
t284 = -t344 * t382 + (-pkin(4) * t352 - pkin(6)) * t363 + t394 * t338 + t364;
t308 = t354 * t332 + t355 * t336;
t297 = -pkin(2) * t363 + qJDD(1) * pkin(6) + t308;
t294 = t357 * t297;
t333 = (-pkin(3) * t360 - qJ(4) * t357) * qJD(1);
t362 = qJD(3) ^ 2;
t373 = -qJ(4) * t362 + t333 * t382 + qJDD(4) + t294;
t287 = pkin(4) * t337 + t394 * qJDD(3) + (-pkin(4) * t380 - t357 * t391 - t353) * t360 + t373;
t356 = sin(qJ(5));
t359 = cos(qJ(5));
t282 = -t284 * t356 + t287 * t359;
t381 = qJD(1) * t360;
t330 = -qJD(3) * t356 - t359 * t381;
t306 = qJD(5) * t330 + qJDD(3) * t359 - t338 * t356;
t331 = qJD(3) * t359 - t356 * t381;
t309 = -mrSges(6,1) * t330 + mrSges(6,2) * t331;
t346 = qJD(5) + t382;
t310 = -mrSges(6,2) * t346 + mrSges(6,3) * t330;
t329 = qJDD(5) + t337;
t279 = m(6) * t282 + mrSges(6,1) * t329 - mrSges(6,3) * t306 - t309 * t331 + t310 * t346;
t283 = t284 * t359 + t287 * t356;
t305 = -qJD(5) * t331 - qJDD(3) * t356 - t338 * t359;
t311 = mrSges(6,1) * t346 - mrSges(6,3) * t331;
t280 = m(6) * t283 - mrSges(6,2) * t329 + mrSges(6,3) * t305 + t309 * t330 - t311 * t346;
t385 = -t356 * t279 + t359 * t280;
t293 = t360 * t297 + t357 * t353;
t292 = -t294 + t386;
t334 = (mrSges(5,2) * t360 - mrSges(5,3) * t357) * qJD(1);
t335 = (-mrSges(4,1) * t360 + mrSges(4,2) * t357) * qJD(1);
t341 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t381;
t342 = -mrSges(5,1) * t381 - qJD(3) * mrSges(5,3);
t272 = t279 * t359 + t280 * t356;
t290 = -qJDD(3) * pkin(3) + t373 - t386;
t369 = -m(5) * t290 - t337 * mrSges(5,1) - t272;
t269 = m(4) * t292 - mrSges(4,3) * t337 + (mrSges(4,1) - mrSges(5,2)) * qJDD(3) + (t341 - t342) * qJD(3) + (-t334 - t335) * t382 + t369;
t340 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t382;
t366 = -pkin(3) * t362 + qJDD(3) * qJ(4) + t333 * t381 + t293;
t289 = qJD(3) * t395 - t366;
t343 = mrSges(5,1) * t382 + qJD(3) * mrSges(5,2);
t286 = -t352 * t391 + pkin(4) * t338 + ((2 * qJD(4)) + t344) * qJD(3) + t366;
t370 = -m(6) * t286 + mrSges(6,1) * t305 - t306 * mrSges(6,2) + t310 * t330 - t331 * t311;
t365 = -m(5) * t289 + qJDD(3) * mrSges(5,3) + qJD(3) * t343 + t334 * t381 - t370;
t276 = t335 * t381 + m(4) * t293 - qJDD(3) * mrSges(4,2) - qJD(3) * t340 + (mrSges(4,3) + mrSges(5,1)) * t338 + t365;
t376 = -t269 * t357 + t360 * t276;
t288 = -pkin(3) * t338 + t364 - t392;
t371 = m(5) * t288 + t338 * mrSges(5,2) - t343 * t382 + t385;
t299 = Ifges(6,4) * t331 + Ifges(6,2) * t330 + Ifges(6,6) * t346;
t300 = Ifges(6,1) * t331 + Ifges(6,4) * t330 + Ifges(6,5) * t346;
t368 = mrSges(6,1) * t282 - mrSges(6,2) * t283 + Ifges(6,5) * t306 + Ifges(6,6) * t305 + Ifges(6,3) * t329 + t331 * t299 - t330 * t300;
t296 = t372 - t392;
t367 = -m(4) * t296 + t338 * mrSges(4,1) + t341 * t381 - t371;
t298 = Ifges(6,5) * t331 + Ifges(6,6) * t330 + Ifges(6,3) * t346;
t274 = mrSges(6,2) * t286 - mrSges(6,3) * t282 + Ifges(6,1) * t306 + Ifges(6,4) * t305 + Ifges(6,5) * t329 + t298 * t330 - t299 * t346;
t273 = -mrSges(6,1) * t286 + mrSges(6,3) * t283 + Ifges(6,4) * t306 + Ifges(6,2) * t305 + Ifges(6,6) * t329 - t298 * t331 + t300 * t346;
t271 = qJDD(3) * mrSges(5,2) + qJD(3) * t342 + t334 * t382 - t369;
t270 = -mrSges(5,3) * t337 + t342 * t381 + t371;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t377 - mrSges(2,2) * t375 + Ifges(3,3) * qJDD(1) + mrSges(3,1) * t307 - mrSges(3,2) * t308 + t357 * (mrSges(5,1) * t290 + mrSges(4,2) * t296 - mrSges(4,3) * t292 - mrSges(5,3) * t288 + pkin(4) * t272 - qJ(4) * t270 + t368) + t360 * (-mrSges(4,1) * t296 - mrSges(5,1) * t289 + mrSges(5,2) * t288 + mrSges(4,3) * t293 - pkin(3) * t270 - pkin(4) * t370 - pkin(7) * t385 - t359 * t273 - t356 * t274) + pkin(2) * t367 + pkin(6) * t376 + pkin(1) * (t354 * (m(3) * t308 - mrSges(3,1) * t363 - qJDD(1) * mrSges(3,2) + t376) + t355 * (m(3) * t307 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t363 + t367)) + t401 * t338 + (t357 * t388 + t360 * t387) * qJDD(3) - t396 * qJD(3) + (t397 * (-mrSges(4,2) + mrSges(5,3)) + t400) * t337 + t397 * qJD(1) * (-t340 * t357 - t342 * t360); m(3) * t353 + t269 * t360 + t276 * t357; mrSges(4,1) * t292 - mrSges(4,2) * t293 + mrSges(5,2) * t290 - mrSges(5,3) * t289 + t359 * t274 - t356 * t273 - pkin(7) * t272 - pkin(3) * t271 + qJ(4) * t365 + (qJ(4) * mrSges(5,1) + t387) * t338 + t388 * t337 + (Ifges(4,3) + Ifges(5,1)) * qJDD(3) + t396 * qJD(1); t271; t368;];
tauJ = t1;
