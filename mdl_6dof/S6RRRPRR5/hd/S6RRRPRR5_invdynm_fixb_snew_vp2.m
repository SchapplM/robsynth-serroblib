% Calculate vector of cutting torques with Newton-Euler for
% S6RRRPRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
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
% Datum: 2019-05-07 10:45
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRPRR5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR5_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR5_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR5_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR5_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR5_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR5_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR5_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 10:39:31
% EndTime: 2019-05-07 10:40:13
% DurationCPUTime: 17.24s
% Computational Cost: add. (288197->388), mult. (592407->465), div. (0->0), fcn. (414725->10), ass. (0->147)
t341 = sin(qJ(2));
t345 = cos(qJ(2));
t369 = qJD(1) * qJD(2);
t318 = t341 * qJDD(1) + t345 * t369;
t342 = sin(qJ(1));
t346 = cos(qJ(1));
t325 = -t346 * g(1) - t342 * g(2);
t347 = qJD(1) ^ 2;
t312 = -t347 * pkin(1) + qJDD(1) * pkin(7) + t325;
t375 = t341 * t312;
t378 = pkin(2) * t347;
t253 = qJDD(2) * pkin(2) - t318 * pkin(8) - t375 + (pkin(8) * t369 + t341 * t378 - g(3)) * t345;
t294 = -t341 * g(3) + t345 * t312;
t319 = t345 * qJDD(1) - t341 * t369;
t372 = qJD(1) * t341;
t323 = qJD(2) * pkin(2) - pkin(8) * t372;
t337 = t345 ^ 2;
t254 = t319 * pkin(8) - qJD(2) * t323 - t337 * t378 + t294;
t340 = sin(qJ(3));
t379 = cos(qJ(3));
t231 = t253 * t379 - t340 * t254;
t232 = t340 * t253 + t254 * t379;
t310 = (t340 * t345 + t341 * t379) * qJD(1);
t270 = t310 * qJD(3) + t340 * t318 - t319 * t379;
t371 = qJD(1) * t345;
t309 = t340 * t372 - t371 * t379;
t271 = -t309 * qJD(3) + t318 * t379 + t340 * t319;
t335 = qJD(2) + qJD(3);
t279 = Ifges(4,4) * t310 - Ifges(4,2) * t309 + Ifges(4,6) * t335;
t287 = -mrSges(5,2) * t309 - mrSges(5,3) * t310;
t297 = mrSges(5,1) * t309 - mrSges(5,3) * t335;
t334 = qJDD(2) + qJDD(3);
t285 = pkin(3) * t309 - qJ(4) * t310;
t333 = t335 ^ 2;
t382 = -2 * qJD(4);
t222 = t333 * pkin(3) - t334 * qJ(4) + t309 * t285 + t335 * t382 - t232;
t299 = pkin(4) * t310 - pkin(9) * t335;
t305 = t309 ^ 2;
t213 = -pkin(4) * t270 - pkin(9) * t305 + t335 * t299 - t222;
t339 = sin(qJ(5));
t344 = cos(qJ(5));
t292 = t339 * t309 + t344 * t335;
t237 = -t292 * qJD(5) + t344 * t270 - t339 * t334;
t291 = t344 * t309 - t339 * t335;
t238 = t291 * qJD(5) + t339 * t270 + t344 * t334;
t304 = qJD(5) + t310;
t273 = -mrSges(6,2) * t304 + mrSges(6,3) * t291;
t274 = mrSges(6,1) * t304 - mrSges(6,3) * t292;
t275 = pkin(5) * t304 - pkin(10) * t292;
t290 = t291 ^ 2;
t206 = -pkin(5) * t237 - pkin(10) * t290 + t275 * t292 + t213;
t338 = sin(qJ(6));
t343 = cos(qJ(6));
t246 = t338 * t291 + t343 * t292;
t217 = -t246 * qJD(6) + t343 * t237 - t338 * t238;
t245 = t343 * t291 - t338 * t292;
t218 = t245 * qJD(6) + t338 * t237 + t343 * t238;
t302 = qJD(6) + t304;
t239 = -mrSges(7,2) * t302 + mrSges(7,3) * t245;
t240 = mrSges(7,1) * t302 - mrSges(7,3) * t246;
t361 = m(7) * t206 - t217 * mrSges(7,1) + t218 * mrSges(7,2) - t239 * t245 + t246 * t240;
t196 = -m(6) * t213 + mrSges(6,1) * t237 - t238 * mrSges(6,2) + t273 * t291 - t292 * t274 - t361;
t298 = mrSges(5,1) * t310 + mrSges(5,2) * t335;
t354 = -m(5) * t222 + t334 * mrSges(5,3) + t335 * t298 - t196;
t324 = t342 * g(1) - t346 * g(2);
t363 = -qJDD(1) * pkin(1) - t324;
t272 = -t319 * pkin(2) + t323 * t372 + (-pkin(8) * t337 - pkin(7)) * t347 + t363;
t376 = t309 * t335;
t352 = (-t271 + t376) * qJ(4) + t272 + (t335 * pkin(3) + t382) * t310;
t208 = -t305 * pkin(4) - t310 * t299 + (pkin(3) + pkin(9)) * t270 + t352;
t224 = -t334 * pkin(3) - t333 * qJ(4) + t310 * t285 + qJDD(4) - t231;
t211 = (t309 * t310 - t334) * pkin(9) + (t271 + t376) * pkin(4) + t224;
t203 = -t208 * t339 + t344 * t211;
t269 = qJDD(5) + t271;
t200 = (t291 * t304 - t238) * pkin(10) + (t291 * t292 + t269) * pkin(5) + t203;
t204 = t344 * t208 + t339 * t211;
t201 = -pkin(5) * t290 + pkin(10) * t237 - t275 * t304 + t204;
t199 = t338 * t200 + t343 * t201;
t225 = Ifges(7,5) * t246 + Ifges(7,6) * t245 + Ifges(7,3) * t302;
t227 = Ifges(7,1) * t246 + Ifges(7,4) * t245 + Ifges(7,5) * t302;
t257 = qJDD(6) + t269;
t184 = -mrSges(7,1) * t206 + mrSges(7,3) * t199 + Ifges(7,4) * t218 + Ifges(7,2) * t217 + Ifges(7,6) * t257 - t225 * t246 + t227 * t302;
t198 = t343 * t200 - t338 * t201;
t226 = Ifges(7,4) * t246 + Ifges(7,2) * t245 + Ifges(7,6) * t302;
t185 = mrSges(7,2) * t206 - mrSges(7,3) * t198 + Ifges(7,1) * t218 + Ifges(7,4) * t217 + Ifges(7,5) * t257 + t225 * t245 - t226 * t302;
t241 = Ifges(6,5) * t292 + Ifges(6,6) * t291 + Ifges(6,3) * t304;
t243 = Ifges(6,1) * t292 + Ifges(6,4) * t291 + Ifges(6,5) * t304;
t229 = -mrSges(7,1) * t245 + mrSges(7,2) * t246;
t192 = m(7) * t198 + mrSges(7,1) * t257 - t218 * mrSges(7,3) - t229 * t246 + t239 * t302;
t193 = m(7) * t199 - mrSges(7,2) * t257 + t217 * mrSges(7,3) + t229 * t245 - t240 * t302;
t365 = -t192 * t338 + t343 * t193;
t166 = -mrSges(6,1) * t213 + mrSges(6,3) * t204 + Ifges(6,4) * t238 + Ifges(6,2) * t237 + Ifges(6,6) * t269 - pkin(5) * t361 + pkin(10) * t365 + t343 * t184 + t338 * t185 - t292 * t241 + t304 * t243;
t183 = t343 * t192 + t338 * t193;
t242 = Ifges(6,4) * t292 + Ifges(6,2) * t291 + Ifges(6,6) * t304;
t168 = mrSges(6,2) * t213 - mrSges(6,3) * t203 + Ifges(6,1) * t238 + Ifges(6,4) * t237 + Ifges(6,5) * t269 - pkin(10) * t183 - t338 * t184 + t343 * t185 + t291 * t241 - t304 * t242;
t250 = -mrSges(6,1) * t291 + mrSges(6,2) * t292;
t180 = m(6) * t203 + mrSges(6,1) * t269 - mrSges(6,3) * t238 - t250 * t292 + t273 * t304 + t183;
t181 = m(6) * t204 - mrSges(6,2) * t269 + mrSges(6,3) * t237 + t250 * t291 - t274 * t304 + t365;
t176 = t344 * t180 + t339 * t181;
t276 = Ifges(5,5) * t335 - Ifges(5,6) * t310 + Ifges(5,3) * t309;
t357 = -mrSges(5,2) * t224 + mrSges(5,3) * t222 - Ifges(5,1) * t334 + Ifges(5,4) * t271 - Ifges(5,5) * t270 + pkin(9) * t176 + t339 * t166 - t344 * t168 + t310 * t276;
t359 = -m(5) * t224 - t271 * mrSges(5,1) - t310 * t287 - t176;
t278 = Ifges(5,4) * t335 - Ifges(5,2) * t310 + Ifges(5,6) * t309;
t373 = Ifges(4,1) * t310 - Ifges(4,4) * t309 + Ifges(4,5) * t335 - t278;
t383 = t373 * t309 - mrSges(4,2) * t232 + pkin(3) * (-t334 * mrSges(5,2) - t335 * t297 + t359) + qJ(4) * (-mrSges(5,1) * t270 - t287 * t309 + t354) + mrSges(4,1) * t231 + t310 * t279 - Ifges(4,6) * t270 + Ifges(4,5) * t271 + Ifges(4,3) * t334 - t357;
t286 = mrSges(4,1) * t309 + mrSges(4,2) * t310;
t295 = -mrSges(4,2) * t335 - mrSges(4,3) * t309;
t172 = m(4) * t231 - t271 * mrSges(4,3) - t310 * t286 + (t295 - t297) * t335 + (mrSges(4,1) - mrSges(5,2)) * t334 + t359;
t296 = mrSges(4,1) * t335 - mrSges(4,3) * t310;
t188 = t354 + m(4) * t232 - mrSges(4,2) * t334 - t296 * t335 + (-mrSges(4,3) - mrSges(5,1)) * t270 + (-t286 - t287) * t309;
t165 = t172 * t379 + t340 * t188;
t293 = -t345 * g(3) - t375;
t307 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t341 + Ifges(3,2) * t345) * qJD(1);
t308 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t341 + Ifges(3,4) * t345) * qJD(1);
t381 = mrSges(3,1) * t293 - mrSges(3,2) * t294 + Ifges(3,5) * t318 + Ifges(3,6) * t319 + Ifges(3,3) * qJDD(2) + pkin(2) * t165 + (t341 * t307 - t345 * t308) * qJD(1) + t383;
t377 = Ifges(4,4) + Ifges(5,6);
t177 = -t339 * t180 + t344 * t181;
t280 = Ifges(5,1) * t335 - Ifges(5,4) * t310 + Ifges(5,5) * t309;
t374 = -Ifges(4,5) * t310 + Ifges(4,6) * t309 - Ifges(4,3) * t335 - t280;
t317 = (-mrSges(3,1) * t345 + mrSges(3,2) * t341) * qJD(1);
t322 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t371;
t163 = m(3) * t293 + qJDD(2) * mrSges(3,1) - t318 * mrSges(3,3) + qJD(2) * t322 - t317 * t372 + t165;
t321 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t372;
t366 = -t340 * t172 + t188 * t379;
t164 = m(3) * t294 - qJDD(2) * mrSges(3,2) + t319 * mrSges(3,3) - qJD(2) * t321 + t317 * t371 + t366;
t367 = -t341 * t163 + t345 * t164;
t220 = t270 * pkin(3) + t352;
t173 = m(5) * t220 - t270 * mrSges(5,2) - t271 * mrSges(5,3) - t309 * t297 - t310 * t298 + t177;
t356 = -mrSges(5,1) * t222 + mrSges(5,2) * t220 - pkin(4) * t196 - pkin(9) * t177 - t344 * t166 - t339 * t168;
t157 = -mrSges(4,1) * t272 + mrSges(4,3) * t232 - pkin(3) * t173 + t373 * t335 + (Ifges(4,6) - Ifges(5,5)) * t334 + t374 * t310 + t377 * t271 + (-Ifges(4,2) - Ifges(5,3)) * t270 + t356;
t358 = -mrSges(7,1) * t198 + mrSges(7,2) * t199 - Ifges(7,5) * t218 - Ifges(7,6) * t217 - Ifges(7,3) * t257 - t246 * t226 + t245 * t227;
t353 = -mrSges(6,1) * t203 + mrSges(6,2) * t204 - Ifges(6,5) * t238 - Ifges(6,6) * t237 - Ifges(6,3) * t269 - pkin(5) * t183 - t292 * t242 + t291 * t243 + t358;
t349 = -mrSges(5,1) * t224 + mrSges(5,3) * t220 - pkin(4) * t176 + t353;
t158 = (t276 - t279) * t335 + (Ifges(4,5) - Ifges(5,4)) * t334 + t374 * t309 + (Ifges(4,1) + Ifges(5,2)) * t271 - t377 * t270 + mrSges(4,2) * t272 - mrSges(4,3) * t231 - qJ(4) * t173 - t349;
t306 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t341 + Ifges(3,6) * t345) * qJD(1);
t311 = -t347 * pkin(7) + t363;
t355 = m(4) * t272 + t270 * mrSges(4,1) + t271 * mrSges(4,2) + t309 * t295 + t310 * t296 + t173;
t153 = -mrSges(3,1) * t311 + mrSges(3,3) * t294 + Ifges(3,4) * t318 + Ifges(3,2) * t319 + Ifges(3,6) * qJDD(2) - pkin(2) * t355 + pkin(8) * t366 + qJD(2) * t308 + t157 * t379 + t340 * t158 - t306 * t372;
t155 = mrSges(3,2) * t311 - mrSges(3,3) * t293 + Ifges(3,1) * t318 + Ifges(3,4) * t319 + Ifges(3,5) * qJDD(2) - pkin(8) * t165 - qJD(2) * t307 - t340 * t157 + t158 * t379 + t306 * t371;
t350 = -m(3) * t311 + t319 * mrSges(3,1) - t318 * mrSges(3,2) - t321 * t372 + t322 * t371 - t355;
t360 = mrSges(2,1) * t324 - mrSges(2,2) * t325 + Ifges(2,3) * qJDD(1) + pkin(1) * t350 + pkin(7) * t367 + t345 * t153 + t341 * t155;
t169 = m(2) * t324 + qJDD(1) * mrSges(2,1) - t347 * mrSges(2,2) + t350;
t161 = t345 * t163 + t341 * t164;
t159 = m(2) * t325 - t347 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t367;
t156 = mrSges(2,1) * g(3) + mrSges(2,3) * t325 + t347 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t161 - t381;
t151 = -mrSges(2,2) * g(3) - mrSges(2,3) * t324 + Ifges(2,5) * qJDD(1) - t347 * Ifges(2,6) - pkin(7) * t161 - t341 * t153 + t345 * t155;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t346 * t151 - t342 * t156 - pkin(6) * (t342 * t159 + t346 * t169), t151, t155, t158, -t309 * t278 - t357, t168, t185; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t342 * t151 + t346 * t156 + pkin(6) * (t346 * t159 - t342 * t169), t156, t153, t157, Ifges(5,4) * t334 - Ifges(5,2) * t271 + Ifges(5,6) * t270 - t335 * t276 + t309 * t280 + t349, t166, t184; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t360, t360, t381, t383, Ifges(5,5) * t334 - Ifges(5,6) * t271 + Ifges(5,3) * t270 + t335 * t278 + t310 * t280 - t356, -t353, -t358;];
m_new  = t1;
