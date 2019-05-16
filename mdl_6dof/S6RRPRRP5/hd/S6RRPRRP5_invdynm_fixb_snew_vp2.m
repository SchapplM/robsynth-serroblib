% Calculate vector of cutting torques with Newton-Euler for
% S6RRPRRP5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-05-06 17:56
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPRRP5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP5_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP5_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP5_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP5_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP5_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP5_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP5_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 17:45:44
% EndTime: 2019-05-06 17:46:50
% DurationCPUTime: 34.00s
% Computational Cost: add. (557510->391), mult. (1458449->499), div. (0->0), fcn. (1148323->12), ass. (0->155)
t386 = -2 * qJD(3);
t340 = sin(pkin(11));
t342 = cos(pkin(11));
t346 = sin(qJ(2));
t350 = cos(qJ(2));
t341 = sin(pkin(6));
t375 = qJD(1) * t341;
t318 = (t340 * t346 - t342 * t350) * t375;
t373 = qJD(1) * qJD(2);
t327 = (qJDD(1) * t346 + t350 * t373) * t341;
t343 = cos(pkin(6));
t334 = t343 * qJDD(1) + qJDD(2);
t335 = t343 * qJD(1) + qJD(2);
t347 = sin(qJ(1));
t351 = cos(qJ(1));
t331 = t347 * g(1) - t351 * g(2);
t352 = qJD(1) ^ 2;
t384 = pkin(8) * t341;
t324 = qJDD(1) * pkin(1) + t352 * t384 + t331;
t332 = -t351 * g(1) - t347 * g(2);
t325 = -t352 * pkin(1) + qJDD(1) * t384 + t332;
t378 = t343 * t350;
t364 = t324 * t378 - t346 * t325;
t382 = t341 ^ 2 * t352;
t258 = t334 * pkin(2) - t327 * qJ(3) + (pkin(2) * t346 * t382 + (qJ(3) * qJD(1) * t335 - g(3)) * t341) * t350 + t364;
t379 = t343 * t346;
t381 = t341 * t346;
t293 = -g(3) * t381 + t324 * t379 + t350 * t325;
t369 = t346 * t375;
t321 = t335 * pkin(2) - qJ(3) * t369;
t328 = (qJDD(1) * t350 - t346 * t373) * t341;
t372 = t350 ^ 2 * t382;
t264 = -pkin(2) * t372 + t328 * qJ(3) - t335 * t321 + t293;
t319 = (t340 * t350 + t342 * t346) * t375;
t224 = t342 * t258 - t340 * t264 + t319 * t386;
t301 = t342 * t327 + t340 * t328;
t345 = sin(qJ(4));
t349 = cos(qJ(4));
t303 = -t345 * t319 + t349 * t335;
t276 = t303 * qJD(4) + t349 * t301 + t345 * t334;
t304 = t349 * t319 + t345 * t335;
t317 = qJD(4) + t318;
t344 = sin(qJ(5));
t348 = cos(qJ(5));
t284 = -t344 * t304 + t348 * t317;
t300 = -t340 * t327 + t342 * t328;
t299 = qJDD(4) - t300;
t240 = t284 * qJD(5) + t348 * t276 + t344 * t299;
t285 = t348 * t304 + t344 * t317;
t252 = -t284 * mrSges(7,1) + t285 * mrSges(7,2);
t225 = t340 * t258 + t342 * t264 + t318 * t386;
t295 = t318 * pkin(3) - t319 * pkin(9);
t333 = t335 ^ 2;
t222 = -t333 * pkin(3) + t334 * pkin(9) - t318 * t295 + t225;
t310 = -t343 * g(3) - t341 * t324;
t278 = -t328 * pkin(2) - qJ(3) * t372 + t321 * t369 + qJDD(3) + t310;
t229 = (t318 * t335 - t301) * pkin(9) + (t319 * t335 - t300) * pkin(3) + t278;
t218 = t349 * t222 + t345 * t229;
t281 = -t303 * pkin(4) - t304 * pkin(10);
t316 = t317 ^ 2;
t213 = -t316 * pkin(4) + t299 * pkin(10) + t303 * t281 + t218;
t221 = -t334 * pkin(3) - t333 * pkin(9) + t319 * t295 - t224;
t275 = -t304 * qJD(4) - t345 * t301 + t349 * t334;
t216 = (-t303 * t317 - t276) * pkin(10) + (t304 * t317 - t275) * pkin(4) + t221;
t207 = -t344 * t213 + t348 * t216;
t274 = qJDD(5) - t275;
t302 = qJD(5) - t303;
t203 = -0.2e1 * qJD(6) * t285 + (t284 * t302 - t240) * qJ(6) + (t284 * t285 + t274) * pkin(5) + t207;
t259 = -t302 * mrSges(7,2) + t284 * mrSges(7,3);
t371 = m(7) * t203 + t274 * mrSges(7,1) + t302 * t259;
t200 = -t240 * mrSges(7,3) - t285 * t252 + t371;
t208 = t348 * t213 + t344 * t216;
t239 = -t285 * qJD(5) - t344 * t276 + t348 * t299;
t246 = Ifges(6,4) * t285 + Ifges(6,2) * t284 + Ifges(6,6) * t302;
t247 = Ifges(7,1) * t285 + Ifges(7,4) * t284 + Ifges(7,5) * t302;
t248 = Ifges(6,1) * t285 + Ifges(6,4) * t284 + Ifges(6,5) * t302;
t261 = t302 * pkin(5) - t285 * qJ(6);
t283 = t284 ^ 2;
t206 = -t283 * pkin(5) + t239 * qJ(6) + 0.2e1 * qJD(6) * t284 - t302 * t261 + t208;
t245 = Ifges(7,4) * t285 + Ifges(7,2) * t284 + Ifges(7,6) * t302;
t360 = -mrSges(7,1) * t203 + mrSges(7,2) * t206 - Ifges(7,5) * t240 - Ifges(7,6) * t239 - Ifges(7,3) * t274 - t285 * t245;
t385 = mrSges(6,1) * t207 - mrSges(6,2) * t208 + Ifges(6,5) * t240 + Ifges(6,6) * t239 + Ifges(6,3) * t274 + pkin(5) * t200 + t285 * t246 - (t248 + t247) * t284 - t360;
t383 = -mrSges(6,2) - mrSges(7,2);
t380 = t341 * t350;
t294 = t318 * mrSges(4,1) + t319 * mrSges(4,2);
t306 = t335 * mrSges(4,1) - t319 * mrSges(4,3);
t253 = -t284 * mrSges(6,1) + t285 * mrSges(6,2);
t260 = -t302 * mrSges(6,2) + t284 * mrSges(6,3);
t194 = m(6) * t207 + t274 * mrSges(6,1) + t302 * t260 + (-t252 - t253) * t285 + (-mrSges(6,3) - mrSges(7,3)) * t240 + t371;
t370 = m(7) * t206 + t239 * mrSges(7,3) + t284 * t252;
t262 = t302 * mrSges(7,1) - t285 * mrSges(7,3);
t376 = -t302 * mrSges(6,1) + t285 * mrSges(6,3) - t262;
t196 = m(6) * t208 + t239 * mrSges(6,3) + t284 * t253 + t274 * t383 + t302 * t376 + t370;
t193 = -t344 * t194 + t348 * t196;
t280 = -t303 * mrSges(5,1) + t304 * mrSges(5,2);
t287 = t317 * mrSges(5,1) - t304 * mrSges(5,3);
t189 = m(5) * t218 - t299 * mrSges(5,2) + t275 * mrSges(5,3) + t303 * t280 - t317 * t287 + t193;
t217 = -t345 * t222 + t349 * t229;
t212 = -t299 * pkin(4) - t316 * pkin(10) + t304 * t281 - t217;
t210 = -t239 * pkin(5) - t283 * qJ(6) + t285 * t261 + qJDD(6) + t212;
t363 = -m(7) * t210 + t239 * mrSges(7,1) + t284 * t259;
t199 = -m(6) * t212 + t239 * mrSges(6,1) + t240 * t383 + t284 * t260 + t285 * t376 + t363;
t286 = -t317 * mrSges(5,2) + t303 * mrSges(5,3);
t198 = m(5) * t217 + t299 * mrSges(5,1) - t276 * mrSges(5,3) - t304 * t280 + t317 * t286 + t199;
t366 = t349 * t189 - t345 * t198;
t180 = m(4) * t225 - t334 * mrSges(4,2) + t300 * mrSges(4,3) - t318 * t294 - t335 * t306 + t366;
t305 = -t335 * mrSges(4,2) - t318 * mrSges(4,3);
t192 = t348 * t194 + t344 * t196;
t355 = -m(5) * t221 + t275 * mrSges(5,1) - t276 * mrSges(5,2) + t303 * t286 - t304 * t287 - t192;
t186 = m(4) * t224 + t334 * mrSges(4,1) - t301 * mrSges(4,3) - t319 * t294 + t335 * t305 + t355;
t175 = t340 * t180 + t342 * t186;
t183 = t345 * t189 + t349 * t198;
t368 = t350 * t375;
t292 = -g(3) * t380 + t364;
t323 = -t335 * mrSges(3,2) + mrSges(3,3) * t368;
t326 = (-mrSges(3,1) * t350 + mrSges(3,2) * t346) * t375;
t173 = m(3) * t292 + t334 * mrSges(3,1) - t327 * mrSges(3,3) + t335 * t323 - t326 * t369 + t175;
t322 = t335 * mrSges(3,1) - mrSges(3,3) * t369;
t367 = t342 * t180 - t340 * t186;
t174 = m(3) * t293 - t334 * mrSges(3,2) + t328 * mrSges(3,3) - t335 * t322 + t326 * t368 + t367;
t167 = -t346 * t173 + t350 * t174;
t357 = m(4) * t278 - t300 * mrSges(4,1) + t301 * mrSges(4,2) + t318 * t305 + t319 * t306 + t183;
t181 = m(3) * t310 - t328 * mrSges(3,1) + t327 * mrSges(3,2) + (t322 * t346 - t323 * t350) * t375 + t357;
t163 = t173 * t378 + t174 * t379 - t341 * t181;
t361 = -mrSges(7,1) * t210 + mrSges(7,3) * t206 + Ifges(7,4) * t240 + Ifges(7,2) * t239 + Ifges(7,6) * t274 + t302 * t247;
t243 = Ifges(7,5) * t285 + Ifges(7,6) * t284 + Ifges(7,3) * t302;
t359 = mrSges(7,2) * t210 - mrSges(7,3) * t203 + Ifges(7,1) * t240 + Ifges(7,4) * t239 + Ifges(7,5) * t274 + t284 * t243;
t244 = Ifges(6,5) * t285 + Ifges(6,6) * t284 + Ifges(6,3) * t302;
t184 = Ifges(6,4) * t240 + Ifges(6,2) * t239 + Ifges(6,6) * t274 + t302 * t248 - mrSges(6,1) * t212 + mrSges(6,3) * t208 - pkin(5) * (t240 * mrSges(7,2) - t363) + qJ(6) * (-t274 * mrSges(7,2) - t302 * t262 + t370) + (-pkin(5) * t262 - t243 - t244) * t285 + t361;
t191 = mrSges(6,2) * t212 - mrSges(6,3) * t207 + Ifges(6,1) * t240 + Ifges(6,4) * t239 + Ifges(6,5) * t274 - qJ(6) * t200 + t284 * t244 + (-t245 - t246) * t302 + t359;
t265 = Ifges(5,5) * t304 + Ifges(5,6) * t303 + Ifges(5,3) * t317;
t266 = Ifges(5,4) * t304 + Ifges(5,2) * t303 + Ifges(5,6) * t317;
t169 = mrSges(5,2) * t221 - mrSges(5,3) * t217 + Ifges(5,1) * t276 + Ifges(5,4) * t275 + Ifges(5,5) * t299 - pkin(10) * t192 - t344 * t184 + t348 * t191 + t303 * t265 - t317 * t266;
t267 = Ifges(5,1) * t304 + Ifges(5,4) * t303 + Ifges(5,5) * t317;
t177 = -mrSges(5,1) * t221 + mrSges(5,3) * t218 + Ifges(5,4) * t276 + Ifges(5,2) * t275 + Ifges(5,6) * t299 - pkin(4) * t192 - t304 * t265 + t317 * t267 - t385;
t288 = Ifges(4,5) * t319 - Ifges(4,6) * t318 + Ifges(4,3) * t335;
t289 = Ifges(4,4) * t319 - Ifges(4,2) * t318 + Ifges(4,6) * t335;
t159 = mrSges(4,2) * t278 - mrSges(4,3) * t224 + Ifges(4,1) * t301 + Ifges(4,4) * t300 + Ifges(4,5) * t334 - pkin(9) * t183 + t349 * t169 - t345 * t177 - t318 * t288 - t335 * t289;
t290 = Ifges(4,1) * t319 - Ifges(4,4) * t318 + Ifges(4,5) * t335;
t353 = mrSges(5,1) * t217 - mrSges(5,2) * t218 + Ifges(5,5) * t276 + Ifges(5,6) * t275 + Ifges(5,3) * t299 + pkin(4) * t199 + pkin(10) * t193 + t348 * t184 + t344 * t191 + t304 * t266 - t303 * t267;
t164 = -mrSges(4,1) * t278 + mrSges(4,3) * t225 + Ifges(4,4) * t301 + Ifges(4,2) * t300 + Ifges(4,6) * t334 - pkin(3) * t183 - t319 * t288 + t335 * t290 - t353;
t307 = Ifges(3,3) * t335 + (Ifges(3,5) * t346 + Ifges(3,6) * t350) * t375;
t309 = Ifges(3,5) * t335 + (Ifges(3,1) * t346 + Ifges(3,4) * t350) * t375;
t154 = -mrSges(3,1) * t310 + mrSges(3,3) * t293 + Ifges(3,4) * t327 + Ifges(3,2) * t328 + Ifges(3,6) * t334 - pkin(2) * t357 + qJ(3) * t367 + t340 * t159 + t342 * t164 - t307 * t369 + t335 * t309;
t308 = Ifges(3,6) * t335 + (Ifges(3,4) * t346 + Ifges(3,2) * t350) * t375;
t156 = mrSges(3,2) * t310 - mrSges(3,3) * t292 + Ifges(3,1) * t327 + Ifges(3,4) * t328 + Ifges(3,5) * t334 - qJ(3) * t175 + t342 * t159 - t340 * t164 + t307 * t368 - t335 * t308;
t356 = mrSges(4,1) * t224 - mrSges(4,2) * t225 + Ifges(4,5) * t301 + Ifges(4,6) * t300 + Ifges(4,3) * t334 + pkin(3) * t355 + pkin(9) * t366 + t345 * t169 + t349 * t177 + t319 * t289 + t318 * t290;
t158 = t356 + (t308 * t346 - t309 * t350) * t375 + pkin(2) * t175 + mrSges(3,1) * t292 - mrSges(3,2) * t293 + Ifges(3,5) * t327 + Ifges(3,6) * t328 + Ifges(3,3) * t334;
t358 = mrSges(2,1) * t331 - mrSges(2,2) * t332 + Ifges(2,3) * qJDD(1) + pkin(1) * t163 + t154 * t380 + t156 * t381 + t343 * t158 + t167 * t384;
t165 = m(2) * t332 - t352 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t167;
t162 = t343 * t181 + (t173 * t350 + t174 * t346) * t341;
t160 = m(2) * t331 + qJDD(1) * mrSges(2,1) - t352 * mrSges(2,2) + t163;
t152 = -mrSges(2,2) * g(3) - mrSges(2,3) * t331 + Ifges(2,5) * qJDD(1) - t352 * Ifges(2,6) - t346 * t154 + t350 * t156 + (-t162 * t341 - t163 * t343) * pkin(8);
t151 = mrSges(2,1) * g(3) + mrSges(2,3) * t332 + t352 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t162 - t341 * t158 + (pkin(8) * t167 + t154 * t350 + t156 * t346) * t343;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t351 * t152 - t347 * t151 - pkin(7) * (t351 * t160 + t347 * t165), t152, t156, t159, t169, t191, -t302 * t245 + t359; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t347 * t152 + t351 * t151 + pkin(7) * (-t347 * t160 + t351 * t165), t151, t154, t164, t177, t184, -t285 * t243 + t361; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t358, t358, t158, t356, t353, t385, -t284 * t247 - t360;];
m_new  = t1;
