% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PRRRR7
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PRRRR7_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR7_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR7_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR7_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR7_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR7_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR7_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR7_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR7_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:12:09
% EndTime: 2019-12-05 17:12:10
% DurationCPUTime: 1.03s
% Computational Cost: add. (8401->212), mult. (17473->276), div. (0->0), fcn. (11886->10), ass. (0->89)
t355 = sin(pkin(9));
t356 = cos(pkin(9));
t342 = -t356 * g(1) - t355 * g(2);
t354 = -g(3) + qJDD(1);
t360 = sin(qJ(2));
t364 = cos(qJ(2));
t324 = t364 * t342 + t360 * t354;
t365 = qJD(2) ^ 2;
t319 = -t365 * pkin(2) + qJDD(2) * pkin(6) + t324;
t341 = -t355 * g(1) + t356 * g(2);
t359 = sin(qJ(3));
t363 = cos(qJ(3));
t308 = -t359 * t319 + t363 * t341;
t376 = qJD(2) * qJD(3);
t375 = t363 * t376;
t339 = t359 * qJDD(2) + t375;
t297 = (-t339 + t375) * pkin(7) + (t359 * t363 * t365 + qJDD(3)) * pkin(3) + t308;
t309 = t363 * t319 + t359 * t341;
t340 = t363 * qJDD(2) - t359 * t376;
t377 = t359 * qJD(2);
t345 = qJD(3) * pkin(3) - pkin(7) * t377;
t353 = t363 ^ 2;
t298 = -t353 * t365 * pkin(3) + t340 * pkin(7) - qJD(3) * t345 + t309;
t358 = sin(qJ(4));
t362 = cos(qJ(4));
t279 = t362 * t297 - t358 * t298;
t329 = (-t359 * t358 + t363 * t362) * qJD(2);
t305 = t329 * qJD(4) + t362 * t339 + t358 * t340;
t330 = (t363 * t358 + t359 * t362) * qJD(2);
t351 = qJDD(3) + qJDD(4);
t352 = qJD(3) + qJD(4);
t274 = (t329 * t352 - t305) * pkin(8) + (t329 * t330 + t351) * pkin(4) + t279;
t280 = t358 * t297 + t362 * t298;
t304 = -t330 * qJD(4) - t358 * t339 + t362 * t340;
t322 = t352 * pkin(4) - t330 * pkin(8);
t325 = t329 ^ 2;
t275 = -t325 * pkin(4) + t304 * pkin(8) - t352 * t322 + t280;
t357 = sin(qJ(5));
t361 = cos(qJ(5));
t272 = t361 * t274 - t357 * t275;
t314 = t361 * t329 - t357 * t330;
t286 = t314 * qJD(5) + t357 * t304 + t361 * t305;
t315 = t357 * t329 + t361 * t330;
t293 = -t314 * mrSges(6,1) + t315 * mrSges(6,2);
t350 = qJD(5) + t352;
t306 = -t350 * mrSges(6,2) + t314 * mrSges(6,3);
t349 = qJDD(5) + t351;
t269 = m(6) * t272 + t349 * mrSges(6,1) - t286 * mrSges(6,3) - t315 * t293 + t350 * t306;
t273 = t357 * t274 + t361 * t275;
t285 = -t315 * qJD(5) + t361 * t304 - t357 * t305;
t307 = t350 * mrSges(6,1) - t315 * mrSges(6,3);
t270 = m(6) * t273 - t349 * mrSges(6,2) + t285 * mrSges(6,3) + t314 * t293 - t350 * t307;
t263 = t361 * t269 + t357 * t270;
t316 = -t329 * mrSges(5,1) + t330 * mrSges(5,2);
t320 = -t352 * mrSges(5,2) + t329 * mrSges(5,3);
t260 = m(5) * t279 + t351 * mrSges(5,1) - t305 * mrSges(5,3) - t330 * t316 + t352 * t320 + t263;
t321 = t352 * mrSges(5,1) - t330 * mrSges(5,3);
t372 = -t357 * t269 + t361 * t270;
t261 = m(5) * t280 - t351 * mrSges(5,2) + t304 * mrSges(5,3) + t329 * t316 - t352 * t321 + t372;
t256 = t362 * t260 + t358 * t261;
t378 = qJD(2) * t363;
t338 = (-t363 * mrSges(4,1) + t359 * mrSges(4,2)) * qJD(2);
t343 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t377;
t344 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t378;
t373 = -t358 * t260 + t362 * t261;
t374 = -t359 * (m(4) * t308 + qJDD(3) * mrSges(4,1) - t339 * mrSges(4,3) + qJD(3) * t344 - t338 * t377 + t256) + t363 * (m(4) * t309 - qJDD(3) * mrSges(4,2) + t340 * mrSges(4,3) - qJD(3) * t343 + t338 * t378 + t373);
t323 = -t360 * t342 + t364 * t354;
t370 = -qJDD(2) * pkin(2) - t323;
t299 = -t340 * pkin(3) + t345 * t377 + (-pkin(7) * t353 - pkin(6)) * t365 + t370;
t277 = -t304 * pkin(4) - t325 * pkin(8) + t330 * t322 + t299;
t371 = m(6) * t277 - t285 * mrSges(6,1) + t286 * mrSges(6,2) - t314 * t306 + t315 * t307;
t289 = Ifges(6,4) * t315 + Ifges(6,2) * t314 + Ifges(6,6) * t350;
t290 = Ifges(6,1) * t315 + Ifges(6,4) * t314 + Ifges(6,5) * t350;
t369 = mrSges(6,1) * t272 - mrSges(6,2) * t273 + Ifges(6,5) * t286 + Ifges(6,6) * t285 + Ifges(6,3) * t349 + t315 * t289 - t314 * t290;
t368 = m(5) * t299 - t304 * mrSges(5,1) + t305 * mrSges(5,2) - t329 * t320 + t330 * t321 + t371;
t311 = Ifges(5,4) * t330 + Ifges(5,2) * t329 + Ifges(5,6) * t352;
t312 = Ifges(5,1) * t330 + Ifges(5,4) * t329 + Ifges(5,5) * t352;
t367 = mrSges(5,1) * t279 - mrSges(5,2) * t280 + Ifges(5,5) * t305 + Ifges(5,6) * t304 + Ifges(5,3) * t351 + pkin(4) * t263 + t330 * t311 - t329 * t312 + t369;
t318 = -t365 * pkin(6) + t370;
t366 = -m(4) * t318 + t340 * mrSges(4,1) - t339 * mrSges(4,2) - t343 * t377 + t344 * t378 - t368;
t328 = Ifges(4,5) * qJD(3) + (t359 * Ifges(4,1) + t363 * Ifges(4,4)) * qJD(2);
t327 = Ifges(4,6) * qJD(3) + (t359 * Ifges(4,4) + t363 * Ifges(4,2)) * qJD(2);
t310 = Ifges(5,5) * t330 + Ifges(5,6) * t329 + Ifges(5,3) * t352;
t288 = Ifges(6,5) * t315 + Ifges(6,6) * t314 + Ifges(6,3) * t350;
t265 = mrSges(6,2) * t277 - mrSges(6,3) * t272 + Ifges(6,1) * t286 + Ifges(6,4) * t285 + Ifges(6,5) * t349 + t314 * t288 - t350 * t289;
t264 = -mrSges(6,1) * t277 + mrSges(6,3) * t273 + Ifges(6,4) * t286 + Ifges(6,2) * t285 + Ifges(6,6) * t349 - t315 * t288 + t350 * t290;
t253 = mrSges(5,2) * t299 - mrSges(5,3) * t279 + Ifges(5,1) * t305 + Ifges(5,4) * t304 + Ifges(5,5) * t351 - pkin(8) * t263 - t357 * t264 + t361 * t265 + t329 * t310 - t352 * t311;
t252 = -mrSges(5,1) * t299 + mrSges(5,3) * t280 + Ifges(5,4) * t305 + Ifges(5,2) * t304 + Ifges(5,6) * t351 - pkin(4) * t371 + pkin(8) * t372 + t361 * t264 + t357 * t265 - t330 * t310 + t352 * t312;
t1 = [m(2) * t354 + t360 * (m(3) * t324 - t365 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t374) + t364 * (m(3) * t323 + qJDD(2) * mrSges(3,1) - t365 * mrSges(3,2) + t366); Ifges(3,3) * qJDD(2) + mrSges(3,1) * t323 - mrSges(3,2) * t324 + t359 * (mrSges(4,2) * t318 - mrSges(4,3) * t308 + Ifges(4,1) * t339 + Ifges(4,4) * t340 + Ifges(4,5) * qJDD(3) - pkin(7) * t256 - qJD(3) * t327 - t358 * t252 + t362 * t253) + t363 * (-mrSges(4,1) * t318 + mrSges(4,3) * t309 + Ifges(4,4) * t339 + Ifges(4,2) * t340 + Ifges(4,6) * qJDD(3) - pkin(3) * t368 + pkin(7) * t373 + qJD(3) * t328 + t362 * t252 + t358 * t253) + pkin(2) * t366 + pkin(6) * t374; t367 + Ifges(4,3) * qJDD(3) + (t359 * t327 - t363 * t328) * qJD(2) + pkin(3) * t256 + mrSges(4,1) * t308 - mrSges(4,2) * t309 + Ifges(4,5) * t339 + Ifges(4,6) * t340; t367; t369;];
tauJ = t1;
