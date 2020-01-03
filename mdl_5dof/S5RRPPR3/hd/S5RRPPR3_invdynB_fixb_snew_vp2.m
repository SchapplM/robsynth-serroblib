% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRPPR3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRPPR3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR3_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR3_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR3_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR3_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR3_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR3_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR3_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:26:32
% EndTime: 2019-12-31 19:26:33
% DurationCPUTime: 0.98s
% Computational Cost: add. (13939->179), mult. (17306->215), div. (0->0), fcn. (8388->8), ass. (0->76)
t370 = -pkin(3) - pkin(7);
t369 = mrSges(4,1) - mrSges(5,2);
t368 = -Ifges(5,4) + Ifges(4,5);
t367 = Ifges(5,5) - Ifges(4,6);
t342 = qJD(1) + qJD(2);
t346 = sin(qJ(5));
t366 = t342 * t346;
t349 = cos(qJ(5));
t365 = t342 * t349;
t348 = sin(qJ(1));
t351 = cos(qJ(1));
t332 = t348 * g(1) - t351 * g(2);
t328 = qJDD(1) * pkin(1) + t332;
t333 = -t351 * g(1) - t348 * g(2);
t352 = qJD(1) ^ 2;
t329 = -t352 * pkin(1) + t333;
t347 = sin(qJ(2));
t350 = cos(qJ(2));
t314 = t350 * t328 - t347 * t329;
t341 = qJDD(1) + qJDD(2);
t312 = t341 * pkin(2) + t314;
t315 = t347 * t328 + t350 * t329;
t340 = t342 ^ 2;
t313 = -t340 * pkin(2) + t315;
t344 = sin(pkin(8));
t345 = cos(pkin(8));
t307 = t345 * t312 - t344 * t313;
t355 = -t340 * qJ(4) + qJDD(4) - t307;
t304 = t370 * t341 + t355;
t343 = -g(3) + qJDD(3);
t300 = t349 * t304 - t346 * t343;
t322 = (mrSges(6,1) * t346 + mrSges(6,2) * t349) * t342;
t363 = qJD(5) * t342;
t324 = t349 * t341 - t346 * t363;
t330 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t366;
t298 = m(6) * t300 + qJDD(5) * mrSges(6,1) - t324 * mrSges(6,3) + qJD(5) * t330 - t322 * t365;
t301 = t346 * t304 + t349 * t343;
t323 = -t346 * t341 - t349 * t363;
t331 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t365;
t299 = m(6) * t301 - qJDD(5) * mrSges(6,2) + t323 * mrSges(6,3) - qJD(5) * t331 - t322 * t366;
t294 = t349 * t298 + t346 * t299;
t306 = -t341 * pkin(3) + t355;
t354 = -m(5) * t306 + t340 * mrSges(5,3) - t294;
t289 = m(4) * t307 - t340 * mrSges(4,2) + t369 * t341 + t354;
t308 = t344 * t312 + t345 * t313;
t356 = t341 * qJ(4) + 0.2e1 * qJD(4) * t342 + t308;
t305 = t340 * pkin(3) - t356;
t303 = t370 * t340 + t356;
t357 = -m(6) * t303 + t323 * mrSges(6,1) - t324 * mrSges(6,2) - t330 * t366 - t331 * t365;
t353 = -m(5) * t305 + t340 * mrSges(5,2) + t341 * mrSges(5,3) - t357;
t292 = m(4) * t308 - t340 * mrSges(4,1) - t341 * mrSges(4,2) + t353;
t287 = t345 * t289 + t344 * t292;
t285 = m(3) * t314 + t341 * mrSges(3,1) - t340 * mrSges(3,2) + t287;
t360 = -t344 * t289 + t345 * t292;
t286 = m(3) * t315 - t340 * mrSges(3,1) - t341 * mrSges(3,2) + t360;
t279 = t350 * t285 + t347 * t286;
t277 = m(2) * t332 + qJDD(1) * mrSges(2,1) - t352 * mrSges(2,2) + t279;
t361 = -t347 * t285 + t350 * t286;
t278 = m(2) * t333 - t352 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t361;
t364 = t351 * t277 + t348 * t278;
t362 = -t348 * t277 + t351 * t278;
t359 = -t346 * t298 + t349 * t299;
t293 = m(5) * t343 + t359;
t358 = m(4) * t343 + t293;
t318 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t349 - Ifges(6,4) * t346) * t342;
t317 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t349 - Ifges(6,2) * t346) * t342;
t316 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t349 - Ifges(6,6) * t346) * t342;
t296 = mrSges(6,2) * t303 - mrSges(6,3) * t300 + Ifges(6,1) * t324 + Ifges(6,4) * t323 + Ifges(6,5) * qJDD(5) - qJD(5) * t317 - t316 * t366;
t295 = -mrSges(6,1) * t303 + mrSges(6,3) * t301 + Ifges(6,4) * t324 + Ifges(6,2) * t323 + Ifges(6,6) * qJDD(5) + qJD(5) * t318 - t316 * t365;
t281 = mrSges(5,1) * t306 + mrSges(6,1) * t300 - mrSges(6,2) * t301 - mrSges(4,3) * t307 + Ifges(6,5) * t324 + Ifges(6,6) * t323 + Ifges(6,3) * qJDD(5) + pkin(4) * t294 - qJ(4) * t293 + (mrSges(4,2) - mrSges(5,3)) * t343 + (t349 * t317 + t346 * t318) * t342 + t368 * t341 + t367 * t340;
t280 = -mrSges(5,1) * t305 + mrSges(4,3) * t308 - pkin(3) * t293 - pkin(4) * t357 - pkin(7) * t359 - t349 * t295 - t346 * t296 + t368 * t340 - t367 * t341 - t369 * t343;
t273 = -mrSges(3,2) * g(3) - mrSges(3,3) * t314 + Ifges(3,5) * t341 - t340 * Ifges(3,6) - qJ(3) * t287 - t344 * t280 + t345 * t281;
t272 = mrSges(3,1) * g(3) + mrSges(3,3) * t315 + t340 * Ifges(3,5) + Ifges(3,6) * t341 - pkin(2) * t358 + qJ(3) * t360 + t345 * t280 + t344 * t281;
t271 = -mrSges(2,2) * g(3) - mrSges(2,3) * t332 + Ifges(2,5) * qJDD(1) - t352 * Ifges(2,6) - pkin(6) * t279 - t347 * t272 + t350 * t273;
t270 = Ifges(2,6) * qJDD(1) + t352 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t333 + t347 * t273 + t350 * t272 - pkin(1) * (-m(3) * g(3) + t358) + pkin(6) * t361;
t1 = [-m(1) * g(1) + t362; -m(1) * g(2) + t364; (-m(1) - m(2) - m(3)) * g(3) + t358; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t364 - t348 * t270 + t351 * t271; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t362 + t351 * t270 + t348 * t271; pkin(1) * t279 + mrSges(2,1) * t332 - mrSges(2,2) * t333 + pkin(2) * t287 - mrSges(3,2) * t315 + mrSges(3,1) * t314 + pkin(3) * t354 + qJ(4) * t353 + mrSges(4,1) * t307 - mrSges(4,2) * t308 + t349 * t296 - t346 * t295 - pkin(7) * t294 + mrSges(5,2) * t306 - mrSges(5,3) * t305 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(2,3) * qJDD(1) + (-pkin(3) * mrSges(5,2) + Ifges(5,1) + Ifges(3,3) + Ifges(4,3)) * t341;];
tauB = t1;
