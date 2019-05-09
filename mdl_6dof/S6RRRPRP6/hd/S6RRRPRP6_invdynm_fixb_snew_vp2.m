% Calculate vector of cutting torques with Newton-Euler for
% S6RRRPRP6
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
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
% m [3x7]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-07 08:11
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRPRP6_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP6_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP6_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP6_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP6_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP6_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP6_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP6_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP6_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 08:00:04
% EndTime: 2019-05-07 08:01:18
% DurationCPUTime: 38.11s
% Computational Cost: add. (686622->393), mult. (1514332->500), div. (0->0), fcn. (1210941->12), ass. (0->155)
t341 = sin(pkin(6));
t346 = sin(qJ(2));
t350 = cos(qJ(2));
t373 = qJD(1) * qJD(2);
t328 = (-qJDD(1) * t350 + t346 * t373) * t341;
t376 = qJD(1) * t341;
t326 = (-pkin(2) * t350 - pkin(9) * t346) * t376;
t343 = cos(pkin(6));
t336 = qJD(1) * t343 + qJD(2);
t334 = t336 ^ 2;
t335 = qJDD(1) * t343 + qJDD(2);
t375 = qJD(1) * t350;
t347 = sin(qJ(1));
t351 = cos(qJ(1));
t332 = t347 * g(1) - g(2) * t351;
t352 = qJD(1) ^ 2;
t386 = pkin(8) * t341;
t323 = qJDD(1) * pkin(1) + t352 * t386 + t332;
t333 = -g(1) * t351 - g(2) * t347;
t324 = -pkin(1) * t352 + qJDD(1) * t386 + t333;
t381 = t343 * t346;
t377 = t323 * t381 + t350 * t324;
t282 = -t334 * pkin(2) + t335 * pkin(9) + (-g(3) * t346 + t326 * t375) * t341 + t377;
t327 = (qJDD(1) * t346 + t350 * t373) * t341;
t385 = t343 * g(3);
t283 = t328 * pkin(2) - t327 * pkin(9) - t385 + (-t323 + (pkin(2) * t346 - pkin(9) * t350) * t336 * qJD(1)) * t341;
t345 = sin(qJ(3));
t349 = cos(qJ(3));
t242 = -t345 * t282 + t349 * t283;
t370 = t346 * t376;
t315 = t336 * t349 - t345 * t370;
t296 = qJD(3) * t315 + t327 * t349 + t335 * t345;
t316 = t336 * t345 + t349 * t370;
t320 = qJDD(3) + t328;
t369 = t341 * t375;
t331 = qJD(3) - t369;
t222 = (t315 * t331 - t296) * qJ(4) + (t315 * t316 + t320) * pkin(3) + t242;
t243 = t349 * t282 + t345 * t283;
t295 = -qJD(3) * t316 - t327 * t345 + t335 * t349;
t306 = pkin(3) * t331 - qJ(4) * t316;
t314 = t315 ^ 2;
t225 = -pkin(3) * t314 + qJ(4) * t295 - t306 * t331 + t243;
t340 = sin(pkin(11));
t342 = cos(pkin(11));
t303 = t315 * t340 + t316 * t342;
t216 = -0.2e1 * qJD(4) * t303 + t222 * t342 - t340 * t225;
t271 = t295 * t340 + t296 * t342;
t344 = sin(qJ(5));
t348 = cos(qJ(5));
t285 = -t303 * t344 + t331 * t348;
t240 = qJD(5) * t285 + t271 * t348 + t320 * t344;
t286 = t303 * t348 + t331 * t344;
t255 = -mrSges(7,1) * t285 + mrSges(7,2) * t286;
t302 = t315 * t342 - t316 * t340;
t217 = 0.2e1 * qJD(4) * t302 + t340 * t222 + t342 * t225;
t277 = -pkin(4) * t302 - pkin(10) * t303;
t330 = t331 ^ 2;
t214 = -pkin(4) * t330 + pkin(10) * t320 + t277 * t302 + t217;
t380 = t343 * t350;
t382 = t341 * t350;
t297 = -g(3) * t382 + t323 * t380 - t346 * t324;
t281 = -t335 * pkin(2) - t334 * pkin(9) + t326 * t370 - t297;
t229 = -t295 * pkin(3) - t314 * qJ(4) + t316 * t306 + qJDD(4) + t281;
t270 = t295 * t342 - t296 * t340;
t220 = (-t302 * t331 - t271) * pkin(10) + (t303 * t331 - t270) * pkin(4) + t229;
t208 = -t344 * t214 + t348 * t220;
t269 = qJDD(5) - t270;
t301 = qJD(5) - t302;
t204 = -0.2e1 * qJD(6) * t286 + (t285 * t301 - t240) * qJ(6) + (t285 * t286 + t269) * pkin(5) + t208;
t257 = -mrSges(7,2) * t301 + mrSges(7,3) * t285;
t372 = m(7) * t204 + t269 * mrSges(7,1) + t301 * t257;
t201 = -t240 * mrSges(7,3) - t286 * t255 + t372;
t209 = t348 * t214 + t344 * t220;
t239 = -qJD(5) * t286 - t271 * t344 + t320 * t348;
t247 = Ifges(6,4) * t286 + Ifges(6,2) * t285 + Ifges(6,6) * t301;
t248 = Ifges(7,1) * t286 + Ifges(7,4) * t285 + Ifges(7,5) * t301;
t249 = Ifges(6,1) * t286 + Ifges(6,4) * t285 + Ifges(6,5) * t301;
t259 = pkin(5) * t301 - qJ(6) * t286;
t284 = t285 ^ 2;
t207 = -pkin(5) * t284 + qJ(6) * t239 + 0.2e1 * qJD(6) * t285 - t259 * t301 + t209;
t246 = Ifges(7,4) * t286 + Ifges(7,2) * t285 + Ifges(7,6) * t301;
t361 = -mrSges(7,1) * t204 + mrSges(7,2) * t207 - Ifges(7,5) * t240 - Ifges(7,6) * t239 - Ifges(7,3) * t269 - t286 * t246;
t387 = mrSges(6,1) * t208 - mrSges(6,2) * t209 + Ifges(6,5) * t240 + Ifges(6,6) * t239 + Ifges(6,3) * t269 + pkin(5) * t201 + t286 * t247 - (t249 + t248) * t285 - t361;
t384 = -mrSges(6,2) - mrSges(7,2);
t383 = t341 * t346;
t276 = -mrSges(5,1) * t302 + mrSges(5,2) * t303;
t288 = mrSges(5,1) * t331 - mrSges(5,3) * t303;
t256 = -mrSges(6,1) * t285 + mrSges(6,2) * t286;
t258 = -mrSges(6,2) * t301 + mrSges(6,3) * t285;
t193 = m(6) * t208 + t269 * mrSges(6,1) + t301 * t258 + (-t255 - t256) * t286 + (-mrSges(6,3) - mrSges(7,3)) * t240 + t372;
t371 = m(7) * t207 + t239 * mrSges(7,3) + t285 * t255;
t260 = mrSges(7,1) * t301 - mrSges(7,3) * t286;
t378 = -mrSges(6,1) * t301 + mrSges(6,3) * t286 - t260;
t198 = m(6) * t209 + t239 * mrSges(6,3) + t285 * t256 + t384 * t269 + t378 * t301 + t371;
t366 = -t193 * t344 + t348 * t198;
t186 = m(5) * t217 - mrSges(5,2) * t320 + mrSges(5,3) * t270 + t276 * t302 - t288 * t331 + t366;
t287 = -mrSges(5,2) * t331 + mrSges(5,3) * t302;
t213 = -pkin(4) * t320 - pkin(10) * t330 + t303 * t277 - t216;
t211 = -pkin(5) * t239 - qJ(6) * t284 + t259 * t286 + qJDD(6) + t213;
t364 = -m(7) * t211 + t239 * mrSges(7,1) + t285 * t257;
t356 = -m(6) * t213 + t239 * mrSges(6,1) + t384 * t240 + t285 * t258 + t378 * t286 + t364;
t195 = m(5) * t216 + t320 * mrSges(5,1) - t271 * mrSges(5,3) - t303 * t276 + t331 * t287 + t356;
t179 = t340 * t186 + t342 * t195;
t304 = -mrSges(4,1) * t315 + mrSges(4,2) * t316;
t305 = -mrSges(4,2) * t331 + mrSges(4,3) * t315;
t177 = m(4) * t242 + mrSges(4,1) * t320 - mrSges(4,3) * t296 - t304 * t316 + t305 * t331 + t179;
t307 = mrSges(4,1) * t331 - mrSges(4,3) * t316;
t367 = t342 * t186 - t195 * t340;
t178 = m(4) * t243 - mrSges(4,2) * t320 + mrSges(4,3) * t295 + t304 * t315 - t307 * t331 + t367;
t172 = t349 * t177 + t345 * t178;
t190 = t348 * t193 + t344 * t198;
t298 = -g(3) * t383 + t377;
t321 = mrSges(3,1) * t336 - mrSges(3,3) * t370;
t325 = (-mrSges(3,1) * t350 + mrSges(3,2) * t346) * t376;
t368 = -t177 * t345 + t349 * t178;
t170 = m(3) * t298 - mrSges(3,2) * t335 - mrSges(3,3) * t328 - t321 * t336 + t325 * t369 + t368;
t322 = -mrSges(3,2) * t336 + mrSges(3,3) * t369;
t358 = m(5) * t229 - t270 * mrSges(5,1) + mrSges(5,2) * t271 - t302 * t287 + t288 * t303 + t190;
t355 = -m(4) * t281 + t295 * mrSges(4,1) - mrSges(4,2) * t296 + t315 * t305 - t307 * t316 - t358;
t183 = m(3) * t297 + mrSges(3,1) * t335 - mrSges(3,3) * t327 + t322 * t336 - t325 * t370 + t355;
t166 = t350 * t170 - t183 * t346;
t311 = -t341 * t323 - t385;
t171 = m(3) * t311 + t328 * mrSges(3,1) + t327 * mrSges(3,2) + (t321 * t346 - t322 * t350) * t376 + t172;
t162 = t170 * t381 - t171 * t341 + t183 * t380;
t362 = -mrSges(7,1) * t211 + mrSges(7,3) * t207 + Ifges(7,4) * t240 + Ifges(7,2) * t239 + Ifges(7,6) * t269 + t301 * t248;
t244 = Ifges(7,5) * t286 + Ifges(7,6) * t285 + Ifges(7,3) * t301;
t360 = mrSges(7,2) * t211 - mrSges(7,3) * t204 + Ifges(7,1) * t240 + Ifges(7,4) * t239 + Ifges(7,5) * t269 + t285 * t244;
t245 = Ifges(6,5) * t286 + Ifges(6,6) * t285 + Ifges(6,3) * t301;
t181 = Ifges(6,4) * t240 + Ifges(6,2) * t239 + Ifges(6,6) * t269 + t301 * t249 - mrSges(6,1) * t213 + mrSges(6,3) * t209 - pkin(5) * (t240 * mrSges(7,2) - t364) + qJ(6) * (-t269 * mrSges(7,2) - t301 * t260 + t371) + (-pkin(5) * t260 - t244 - t245) * t286 + t362;
t188 = mrSges(6,2) * t213 - mrSges(6,3) * t208 + Ifges(6,1) * t240 + Ifges(6,4) * t239 + Ifges(6,5) * t269 - qJ(6) * t201 + t285 * t245 + (-t246 - t247) * t301 + t360;
t272 = Ifges(5,5) * t303 + Ifges(5,6) * t302 + Ifges(5,3) * t331;
t273 = Ifges(5,4) * t303 + Ifges(5,2) * t302 + Ifges(5,6) * t331;
t167 = mrSges(5,2) * t229 - mrSges(5,3) * t216 + Ifges(5,1) * t271 + Ifges(5,4) * t270 + Ifges(5,5) * t320 - pkin(10) * t190 - t181 * t344 + t188 * t348 + t272 * t302 - t273 * t331;
t274 = Ifges(5,1) * t303 + Ifges(5,4) * t302 + Ifges(5,5) * t331;
t173 = -mrSges(5,1) * t229 + mrSges(5,3) * t217 + Ifges(5,4) * t271 + Ifges(5,2) * t270 + Ifges(5,6) * t320 - pkin(4) * t190 - t303 * t272 + t331 * t274 - t387;
t289 = Ifges(4,5) * t316 + Ifges(4,6) * t315 + Ifges(4,3) * t331;
t291 = Ifges(4,1) * t316 + Ifges(4,4) * t315 + Ifges(4,5) * t331;
t158 = -mrSges(4,1) * t281 + mrSges(4,3) * t243 + Ifges(4,4) * t296 + Ifges(4,2) * t295 + Ifges(4,6) * t320 - pkin(3) * t358 + qJ(4) * t367 + t340 * t167 + t342 * t173 - t316 * t289 + t331 * t291;
t290 = Ifges(4,4) * t316 + Ifges(4,2) * t315 + Ifges(4,6) * t331;
t163 = mrSges(4,2) * t281 - mrSges(4,3) * t242 + Ifges(4,1) * t296 + Ifges(4,4) * t295 + Ifges(4,5) * t320 - qJ(4) * t179 + t167 * t342 - t173 * t340 + t289 * t315 - t290 * t331;
t309 = Ifges(3,6) * t336 + (Ifges(3,4) * t346 + Ifges(3,2) * t350) * t376;
t310 = Ifges(3,5) * t336 + (Ifges(3,1) * t346 + Ifges(3,4) * t350) * t376;
t153 = Ifges(3,5) * t327 - Ifges(3,6) * t328 + Ifges(3,3) * t335 + mrSges(3,1) * t297 - mrSges(3,2) * t298 + t345 * t163 + t349 * t158 + pkin(2) * t355 + pkin(9) * t368 + (t309 * t346 - t310 * t350) * t376;
t308 = Ifges(3,3) * t336 + (Ifges(3,5) * t346 + Ifges(3,6) * t350) * t376;
t155 = mrSges(3,2) * t311 - mrSges(3,3) * t297 + Ifges(3,1) * t327 - Ifges(3,4) * t328 + Ifges(3,5) * t335 - pkin(9) * t172 - t158 * t345 + t163 * t349 + t308 * t369 - t309 * t336;
t357 = -mrSges(5,1) * t216 + mrSges(5,2) * t217 - Ifges(5,5) * t271 - Ifges(5,6) * t270 - Ifges(5,3) * t320 - pkin(4) * t356 - pkin(10) * t366 - t348 * t181 - t344 * t188 - t303 * t273 + t302 * t274;
t353 = mrSges(4,1) * t242 - mrSges(4,2) * t243 + Ifges(4,5) * t296 + Ifges(4,6) * t295 + Ifges(4,3) * t320 + pkin(3) * t179 + t316 * t290 - t315 * t291 - t357;
t157 = -mrSges(3,1) * t311 + mrSges(3,3) * t298 + Ifges(3,4) * t327 - Ifges(3,2) * t328 + Ifges(3,6) * t335 - pkin(2) * t172 - t308 * t370 + t336 * t310 - t353;
t359 = mrSges(2,1) * t332 - mrSges(2,2) * t333 + Ifges(2,3) * qJDD(1) + pkin(1) * t162 + t343 * t153 + t155 * t383 + t157 * t382 + t166 * t386;
t164 = m(2) * t333 - mrSges(2,1) * t352 - qJDD(1) * mrSges(2,2) + t166;
t161 = t343 * t171 + (t170 * t346 + t183 * t350) * t341;
t159 = m(2) * t332 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t352 + t162;
t151 = -mrSges(2,2) * g(3) - mrSges(2,3) * t332 + Ifges(2,5) * qJDD(1) - t352 * Ifges(2,6) + t350 * t155 - t346 * t157 + (-t161 * t341 - t162 * t343) * pkin(8);
t150 = mrSges(2,1) * g(3) + mrSges(2,3) * t333 + t352 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t161 - t341 * t153 + (pkin(8) * t166 + t155 * t346 + t157 * t350) * t343;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t351 * t151 - t347 * t150 - pkin(7) * (t159 * t351 + t164 * t347), t151, t155, t163, t167, t188, -t246 * t301 + t360; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t347 * t151 + t351 * t150 + pkin(7) * (-t159 * t347 + t164 * t351), t150, t157, t158, t173, t181, -t286 * t244 + t362; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t359, t359, t153, t353, -t357, t387, -t285 * t248 - t361;];
m_new  = t1;
