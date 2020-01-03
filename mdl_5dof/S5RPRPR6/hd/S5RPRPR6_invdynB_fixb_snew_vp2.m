% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPRPR6
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPRPR6_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR6_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR6_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR6_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR6_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR6_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR6_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR6_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR6_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:17:46
% EndTime: 2019-12-31 18:17:47
% DurationCPUTime: 0.98s
% Computational Cost: add. (13029->178), mult. (17306->215), div. (0->0), fcn. (8388->8), ass. (0->77)
t371 = -pkin(3) - pkin(7);
t370 = mrSges(4,1) - mrSges(5,2);
t369 = -Ifges(5,4) + Ifges(4,5);
t368 = Ifges(5,5) - Ifges(4,6);
t342 = qJD(1) + qJD(3);
t346 = sin(qJ(5));
t367 = t342 * t346;
t349 = cos(qJ(5));
t366 = t342 * t349;
t348 = sin(qJ(1));
t351 = cos(qJ(1));
t331 = t348 * g(1) - g(2) * t351;
t327 = qJDD(1) * pkin(1) + t331;
t332 = -g(1) * t351 - g(2) * t348;
t352 = qJD(1) ^ 2;
t328 = -pkin(1) * t352 + t332;
t344 = sin(pkin(8));
t345 = cos(pkin(8));
t313 = t345 * t327 - t328 * t344;
t311 = qJDD(1) * pkin(2) + t313;
t314 = t344 * t327 + t345 * t328;
t312 = -pkin(2) * t352 + t314;
t347 = sin(qJ(3));
t350 = cos(qJ(3));
t306 = t350 * t311 - t347 * t312;
t340 = t342 ^ 2;
t341 = qJDD(1) + qJDD(3);
t355 = -t340 * qJ(4) + qJDD(4) - t306;
t303 = t341 * t371 + t355;
t343 = -g(3) + qJDD(2);
t299 = t303 * t349 - t343 * t346;
t321 = (mrSges(6,1) * t346 + mrSges(6,2) * t349) * t342;
t364 = qJD(5) * t342;
t323 = t341 * t349 - t346 * t364;
t329 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t367;
t297 = m(6) * t299 + qJDD(5) * mrSges(6,1) - mrSges(6,3) * t323 + qJD(5) * t329 - t321 * t366;
t300 = t303 * t346 + t343 * t349;
t322 = -t341 * t346 - t349 * t364;
t330 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t366;
t298 = m(6) * t300 - qJDD(5) * mrSges(6,2) + mrSges(6,3) * t322 - qJD(5) * t330 - t321 * t367;
t293 = t349 * t297 + t346 * t298;
t305 = -pkin(3) * t341 + t355;
t354 = -m(5) * t305 + t340 * mrSges(5,3) - t293;
t288 = m(4) * t306 - t340 * mrSges(4,2) + t341 * t370 + t354;
t307 = t347 * t311 + t350 * t312;
t356 = t341 * qJ(4) + 0.2e1 * qJD(4) * t342 + t307;
t304 = pkin(3) * t340 - t356;
t302 = t340 * t371 + t356;
t358 = -m(6) * t302 + mrSges(6,1) * t322 - t323 * mrSges(6,2) - t329 * t367 - t330 * t366;
t353 = -m(5) * t304 + t340 * mrSges(5,2) + t341 * mrSges(5,3) - t358;
t291 = m(4) * t307 - mrSges(4,1) * t340 - mrSges(4,2) * t341 + t353;
t286 = t350 * t288 + t347 * t291;
t284 = m(3) * t313 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t352 + t286;
t361 = -t288 * t347 + t350 * t291;
t285 = m(3) * t314 - mrSges(3,1) * t352 - qJDD(1) * mrSges(3,2) + t361;
t278 = t345 * t284 + t344 * t285;
t276 = m(2) * t331 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t352 + t278;
t362 = -t284 * t344 + t345 * t285;
t277 = m(2) * t332 - mrSges(2,1) * t352 - qJDD(1) * mrSges(2,2) + t362;
t365 = t351 * t276 + t348 * t277;
t363 = -t276 * t348 + t351 * t277;
t360 = -t346 * t297 + t349 * t298;
t292 = m(5) * t343 + t360;
t359 = m(4) * t343 + t292;
t357 = m(3) * t343 + t359;
t317 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t349 - Ifges(6,4) * t346) * t342;
t316 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t349 - Ifges(6,2) * t346) * t342;
t315 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t349 - Ifges(6,6) * t346) * t342;
t295 = mrSges(6,2) * t302 - mrSges(6,3) * t299 + Ifges(6,1) * t323 + Ifges(6,4) * t322 + Ifges(6,5) * qJDD(5) - qJD(5) * t316 - t315 * t367;
t294 = -mrSges(6,1) * t302 + mrSges(6,3) * t300 + Ifges(6,4) * t323 + Ifges(6,2) * t322 + Ifges(6,6) * qJDD(5) + qJD(5) * t317 - t315 * t366;
t280 = mrSges(5,1) * t305 + mrSges(6,1) * t299 - mrSges(6,2) * t300 - mrSges(4,3) * t306 + Ifges(6,5) * t323 + Ifges(6,6) * t322 + Ifges(6,3) * qJDD(5) + pkin(4) * t293 - qJ(4) * t292 + (mrSges(4,2) - mrSges(5,3)) * t343 + (t316 * t349 + t317 * t346) * t342 + t369 * t341 + t368 * t340;
t279 = -mrSges(5,1) * t304 + mrSges(4,3) * t307 - pkin(3) * t292 - pkin(4) * t358 - pkin(7) * t360 - t349 * t294 - t346 * t295 + t340 * t369 - t341 * t368 - t343 * t370;
t272 = mrSges(3,2) * t343 - mrSges(3,3) * t313 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t352 - pkin(6) * t286 - t279 * t347 + t280 * t350;
t271 = -mrSges(3,1) * t343 + mrSges(3,3) * t314 + t352 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t359 + pkin(6) * t361 + t350 * t279 + t347 * t280;
t270 = -mrSges(2,2) * g(3) - mrSges(2,3) * t331 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t352 - qJ(2) * t278 - t271 * t344 + t272 * t345;
t269 = mrSges(2,1) * g(3) + mrSges(2,3) * t332 + t352 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t357 + qJ(2) * t362 + t345 * t271 + t344 * t272;
t1 = [-m(1) * g(1) + t363; -m(1) * g(2) + t365; (-m(1) - m(2)) * g(3) + t357; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t365 - t348 * t269 + t351 * t270; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t363 + t351 * t269 + t348 * t270; pkin(1) * t278 + mrSges(2,1) * t331 - mrSges(2,2) * t332 + pkin(2) * t286 + mrSges(3,1) * t313 - mrSges(3,2) * t314 + pkin(3) * t354 + qJ(4) * t353 + mrSges(4,1) * t306 - mrSges(4,2) * t307 + t349 * t295 - t346 * t294 - pkin(7) * t293 + mrSges(5,2) * t305 - mrSges(5,3) * t304 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-mrSges(5,2) * pkin(3) + Ifges(5,1) + Ifges(4,3)) * t341 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);];
tauB = t1;
