% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPPRR6
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-12-31 17:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPPRR6_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR6_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR6_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR6_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR6_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR6_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR6_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR6_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR6_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:57:53
% EndTime: 2019-12-31 17:57:54
% DurationCPUTime: 0.91s
% Computational Cost: add. (6021->192), mult. (13170->247), div. (0->0), fcn. (8655->10), ass. (0->91)
t358 = qJD(1) ^ 2;
t347 = sin(pkin(9));
t349 = cos(pkin(9));
t352 = sin(qJ(4));
t355 = cos(qJ(4));
t364 = t347 * t352 - t349 * t355;
t327 = t364 * qJD(1);
t365 = t347 * t355 + t349 * t352;
t328 = t365 * qJD(1);
t374 = qJD(4) * t328;
t318 = -t364 * qJDD(1) - t374;
t381 = t349 * pkin(3);
t380 = mrSges(4,2) * t347;
t345 = t349 ^ 2;
t379 = t345 * t358;
t353 = sin(qJ(1));
t356 = cos(qJ(1));
t371 = t353 * g(1) - g(2) * t356;
t334 = qJDD(1) * pkin(1) + t371;
t367 = -g(1) * t356 - g(2) * t353;
t335 = -pkin(1) * t358 + t367;
t348 = sin(pkin(8));
t350 = cos(pkin(8));
t321 = t348 * t334 + t350 * t335;
t312 = -pkin(2) * t358 + qJDD(1) * qJ(3) + t321;
t346 = -g(3) + qJDD(2);
t373 = qJD(1) * qJD(3);
t377 = t349 * t346 - 0.2e1 * t347 * t373;
t297 = (-pkin(6) * qJDD(1) + t358 * t381 - t312) * t347 + t377;
t303 = t347 * t346 + (t312 + 0.2e1 * t373) * t349;
t372 = t349 * qJDD(1);
t300 = -pkin(3) * t379 + pkin(6) * t372 + t303;
t289 = t352 * t297 + t355 * t300;
t314 = mrSges(5,1) * t327 + mrSges(5,2) * t328;
t325 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t328;
t317 = pkin(4) * t327 - pkin(7) * t328;
t357 = qJD(4) ^ 2;
t286 = -pkin(4) * t357 + qJDD(4) * pkin(7) - t317 * t327 + t289;
t344 = t347 ^ 2;
t320 = t334 * t350 - t348 * t335;
t366 = qJDD(3) - t320;
t301 = (-pkin(2) - t381) * qJDD(1) + (-qJ(3) + (-t344 - t345) * pkin(6)) * t358 + t366;
t375 = qJD(4) * t327;
t319 = t365 * qJDD(1) - t375;
t287 = (-t319 + t375) * pkin(7) + (-t318 + t374) * pkin(4) + t301;
t351 = sin(qJ(5));
t354 = cos(qJ(5));
t283 = -t286 * t351 + t287 * t354;
t322 = qJD(4) * t354 - t328 * t351;
t299 = qJD(5) * t322 + qJDD(4) * t351 + t319 * t354;
t323 = qJD(4) * t351 + t328 * t354;
t304 = -mrSges(6,1) * t322 + mrSges(6,2) * t323;
t326 = qJD(5) + t327;
t305 = -mrSges(6,2) * t326 + mrSges(6,3) * t322;
t316 = qJDD(5) - t318;
t281 = m(6) * t283 + mrSges(6,1) * t316 - mrSges(6,3) * t299 - t304 * t323 + t305 * t326;
t284 = t286 * t354 + t287 * t351;
t298 = -qJD(5) * t323 + qJDD(4) * t354 - t319 * t351;
t306 = mrSges(6,1) * t326 - mrSges(6,3) * t323;
t282 = m(6) * t284 - mrSges(6,2) * t316 + mrSges(6,3) * t298 + t304 * t322 - t306 * t326;
t368 = -t281 * t351 + t354 * t282;
t272 = m(5) * t289 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t318 - qJD(4) * t325 - t314 * t327 + t368;
t288 = t297 * t355 - t300 * t352;
t324 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t327;
t285 = -qJDD(4) * pkin(4) - pkin(7) * t357 + t317 * t328 - t288;
t362 = -m(6) * t285 + t298 * mrSges(6,1) - mrSges(6,2) * t299 + t322 * t305 - t306 * t323;
t277 = m(5) * t288 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t319 + qJD(4) * t324 - t314 * t328 + t362;
t378 = t352 * t272 + t355 * t277;
t273 = t354 * t281 + t351 * t282;
t302 = -t312 * t347 + t377;
t363 = mrSges(4,3) * qJDD(1) + t358 * (-t349 * mrSges(4,1) + t380);
t267 = m(4) * t302 - t363 * t347 + t378;
t369 = t355 * t272 - t277 * t352;
t268 = m(4) * t303 + t363 * t349 + t369;
t370 = -t347 * t267 + t349 * t268;
t361 = m(5) * t301 - t318 * mrSges(5,1) + mrSges(5,2) * t319 + t327 * t324 + t325 * t328 + t273;
t308 = -qJDD(1) * pkin(2) - qJ(3) * t358 + t366;
t360 = -m(4) * t308 + mrSges(4,1) * t372 - t361 + (t344 * t358 + t379) * mrSges(4,3);
t292 = Ifges(6,4) * t323 + Ifges(6,2) * t322 + Ifges(6,6) * t326;
t293 = Ifges(6,1) * t323 + Ifges(6,4) * t322 + Ifges(6,5) * t326;
t359 = mrSges(6,1) * t283 - mrSges(6,2) * t284 + Ifges(6,5) * t299 + Ifges(6,6) * t298 + Ifges(6,3) * t316 + t292 * t323 - t293 * t322;
t311 = Ifges(5,1) * t328 - Ifges(5,4) * t327 + Ifges(5,5) * qJD(4);
t310 = Ifges(5,4) * t328 - Ifges(5,2) * t327 + Ifges(5,6) * qJD(4);
t309 = Ifges(5,5) * t328 - Ifges(5,6) * t327 + Ifges(5,3) * qJD(4);
t291 = Ifges(6,5) * t323 + Ifges(6,6) * t322 + Ifges(6,3) * t326;
t275 = mrSges(6,2) * t285 - mrSges(6,3) * t283 + Ifges(6,1) * t299 + Ifges(6,4) * t298 + Ifges(6,5) * t316 + t291 * t322 - t292 * t326;
t274 = -mrSges(6,1) * t285 + mrSges(6,3) * t284 + Ifges(6,4) * t299 + Ifges(6,2) * t298 + Ifges(6,6) * t316 - t291 * t323 + t293 * t326;
t269 = qJDD(1) * t380 - t360;
t265 = -mrSges(5,1) * t301 + mrSges(5,3) * t289 + Ifges(5,4) * t319 + Ifges(5,2) * t318 + Ifges(5,6) * qJDD(4) - pkin(4) * t273 + qJD(4) * t311 - t309 * t328 - t359;
t264 = mrSges(5,2) * t301 - mrSges(5,3) * t288 + Ifges(5,1) * t319 + Ifges(5,4) * t318 + Ifges(5,5) * qJDD(4) - pkin(7) * t273 - qJD(4) * t310 - t274 * t351 + t275 * t354 - t309 * t327;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t371 - mrSges(2,2) * t367 + Ifges(3,3) * qJDD(1) + mrSges(3,1) * t320 - mrSges(3,2) * t321 + t347 * (mrSges(4,2) * t308 - mrSges(4,3) * t302 + t355 * t264 - t352 * t265 - pkin(6) * t378 + (Ifges(4,1) * t347 + Ifges(4,4) * t349) * qJDD(1)) + t349 * (-mrSges(4,1) * t308 + mrSges(4,3) * t303 + t352 * t264 + t355 * t265 - pkin(3) * t361 + pkin(6) * t369 + (Ifges(4,4) * t347 + Ifges(4,2) * t349) * qJDD(1)) - pkin(2) * t269 + qJ(3) * t370 + pkin(1) * (t348 * (m(3) * t321 - mrSges(3,1) * t358 - qJDD(1) * mrSges(3,2) + t370) + t350 * (t360 + m(3) * t320 - mrSges(3,2) * t358 + (mrSges(3,1) - t380) * qJDD(1))); m(3) * t346 + t267 * t349 + t268 * t347; t269; mrSges(5,1) * t288 - mrSges(5,2) * t289 + Ifges(5,5) * t319 + Ifges(5,6) * t318 + Ifges(5,3) * qJDD(4) + pkin(4) * t362 + pkin(7) * t368 + t354 * t274 + t351 * t275 + t328 * t310 + t327 * t311; t359;];
tauJ = t1;
