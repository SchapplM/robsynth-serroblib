% Calculate vector of cutting torques with Newton-Euler for
% S6RRPRRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-05-06 20:23
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPRRR3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR3_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR3_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR3_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR3_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR3_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR3_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR3_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 20:16:34
% EndTime: 2019-05-06 20:18:11
% DurationCPUTime: 45.44s
% Computational Cost: add. (791004->386), mult. (1816945->487), div. (0->0), fcn. (1342684->12), ass. (0->153)
t338 = sin(qJ(2));
t343 = cos(qJ(2));
t366 = qJD(1) * qJD(2);
t320 = t338 * qJDD(1) + t343 * t366;
t339 = sin(qJ(1));
t344 = cos(qJ(1));
t327 = -t344 * g(1) - t339 * g(2);
t346 = qJD(1) ^ 2;
t315 = -t346 * pkin(1) + qJDD(1) * pkin(7) + t327;
t370 = t338 * t315;
t371 = pkin(2) * t346;
t271 = qJDD(2) * pkin(2) - t320 * qJ(3) - t370 + (qJ(3) * t366 + t338 * t371 - g(3)) * t343;
t299 = -t338 * g(3) + t343 * t315;
t321 = t343 * qJDD(1) - t338 * t366;
t369 = qJD(1) * t338;
t323 = qJD(2) * pkin(2) - qJ(3) * t369;
t332 = t343 ^ 2;
t272 = t321 * qJ(3) - qJD(2) * t323 - t332 * t371 + t299;
t333 = sin(pkin(11));
t334 = cos(pkin(11));
t309 = (t333 * t343 + t334 * t338) * qJD(1);
t248 = -0.2e1 * qJD(3) * t309 + t334 * t271 - t333 * t272;
t368 = qJD(1) * t343;
t308 = -t333 * t369 + t334 * t368;
t249 = 0.2e1 * qJD(3) * t308 + t333 * t271 + t334 * t272;
t284 = -t308 * mrSges(4,1) + t309 * mrSges(4,2);
t293 = -t333 * t320 + t334 * t321;
t301 = qJD(2) * mrSges(4,1) - t309 * mrSges(4,3);
t287 = -t308 * pkin(3) - t309 * pkin(8);
t345 = qJD(2) ^ 2;
t237 = -t345 * pkin(3) + qJDD(2) * pkin(8) + t287 * t308 + t249;
t326 = t339 * g(1) - t344 * g(2);
t357 = -qJDD(1) * pkin(1) - t326;
t275 = -t321 * pkin(2) + qJDD(3) + t323 * t369 + (-qJ(3) * t332 - pkin(7)) * t346 + t357;
t294 = t334 * t320 + t333 * t321;
t240 = (-qJD(2) * t308 - t294) * pkin(8) + (qJD(2) * t309 - t293) * pkin(3) + t275;
t337 = sin(qJ(4));
t342 = cos(qJ(4));
t221 = -t337 * t237 + t342 * t240;
t296 = t342 * qJD(2) - t337 * t309;
t263 = t296 * qJD(4) + t337 * qJDD(2) + t342 * t294;
t292 = qJDD(4) - t293;
t297 = t337 * qJD(2) + t342 * t309;
t307 = qJD(4) - t308;
t210 = (t296 * t307 - t263) * pkin(9) + (t296 * t297 + t292) * pkin(4) + t221;
t222 = t342 * t237 + t337 * t240;
t262 = -t297 * qJD(4) + t342 * qJDD(2) - t337 * t294;
t278 = pkin(4) * t307 - pkin(9) * t297;
t295 = t296 ^ 2;
t218 = -pkin(4) * t295 + t262 * pkin(9) - t278 * t307 + t222;
t336 = sin(qJ(5));
t341 = cos(qJ(5));
t204 = t341 * t210 - t336 * t218;
t268 = t341 * t296 - t336 * t297;
t233 = t268 * qJD(5) + t336 * t262 + t341 * t263;
t269 = t336 * t296 + t341 * t297;
t288 = qJDD(5) + t292;
t303 = qJD(5) + t307;
t201 = (t268 * t303 - t233) * pkin(10) + (t268 * t269 + t288) * pkin(5) + t204;
t205 = t336 * t210 + t341 * t218;
t232 = -t269 * qJD(5) + t341 * t262 - t336 * t263;
t255 = pkin(5) * t303 - t269 * pkin(10);
t267 = t268 ^ 2;
t202 = -t267 * pkin(5) + t232 * pkin(10) - t255 * t303 + t205;
t335 = sin(qJ(6));
t340 = cos(qJ(6));
t199 = t340 * t201 - t335 * t202;
t250 = t340 * t268 - t335 * t269;
t216 = t250 * qJD(6) + t335 * t232 + t340 * t233;
t251 = t335 * t268 + t340 * t269;
t228 = -mrSges(7,1) * t250 + mrSges(7,2) * t251;
t302 = qJD(6) + t303;
t241 = -mrSges(7,2) * t302 + t250 * mrSges(7,3);
t286 = qJDD(6) + t288;
t194 = m(7) * t199 + t286 * mrSges(7,1) - t216 * mrSges(7,3) - t251 * t228 + t241 * t302;
t200 = t335 * t201 + t340 * t202;
t215 = -t251 * qJD(6) + t340 * t232 - t335 * t233;
t242 = mrSges(7,1) * t302 - t251 * mrSges(7,3);
t195 = m(7) * t200 - t286 * mrSges(7,2) + t215 * mrSges(7,3) + t250 * t228 - t242 * t302;
t186 = t340 * t194 + t335 * t195;
t252 = -mrSges(6,1) * t268 + mrSges(6,2) * t269;
t253 = -mrSges(6,2) * t303 + t268 * mrSges(6,3);
t183 = m(6) * t204 + t288 * mrSges(6,1) - t233 * mrSges(6,3) - t269 * t252 + t253 * t303 + t186;
t254 = mrSges(6,1) * t303 - t269 * mrSges(6,3);
t361 = -t335 * t194 + t340 * t195;
t184 = m(6) * t205 - t288 * mrSges(6,2) + t232 * mrSges(6,3) + t268 * t252 - t303 * t254 + t361;
t179 = t341 * t183 + t336 * t184;
t273 = -mrSges(5,1) * t296 + mrSges(5,2) * t297;
t276 = -mrSges(5,2) * t307 + mrSges(5,3) * t296;
t177 = m(5) * t221 + mrSges(5,1) * t292 - t263 * mrSges(5,3) - t273 * t297 + t276 * t307 + t179;
t277 = mrSges(5,1) * t307 - mrSges(5,3) * t297;
t362 = -t336 * t183 + t341 * t184;
t178 = m(5) * t222 - t292 * mrSges(5,2) + t262 * mrSges(5,3) + t296 * t273 - t307 * t277 + t362;
t363 = -t337 * t177 + t342 * t178;
t168 = m(4) * t249 - qJDD(2) * mrSges(4,2) + t293 * mrSges(4,3) - qJD(2) * t301 + t308 * t284 + t363;
t300 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t308;
t236 = -qJDD(2) * pkin(3) - t345 * pkin(8) + t309 * t287 - t248;
t220 = -t262 * pkin(4) - t295 * pkin(9) + t297 * t278 + t236;
t207 = -t232 * pkin(5) - t267 * pkin(10) + t269 * t255 + t220;
t359 = m(7) * t207 - t215 * mrSges(7,1) + t216 * mrSges(7,2) - t250 * t241 + t251 * t242;
t352 = m(6) * t220 - t232 * mrSges(6,1) + t233 * mrSges(6,2) - t268 * t253 + t269 * t254 + t359;
t349 = -m(5) * t236 + t262 * mrSges(5,1) - t263 * mrSges(5,2) + t296 * t276 - t297 * t277 - t352;
t190 = m(4) * t248 + qJDD(2) * mrSges(4,1) - t294 * mrSges(4,3) + qJD(2) * t300 - t309 * t284 + t349;
t163 = t333 * t168 + t334 * t190;
t298 = -t343 * g(3) - t370;
t311 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t338 + Ifges(3,2) * t343) * qJD(1);
t312 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t338 + Ifges(3,4) * t343) * qJD(1);
t223 = Ifges(7,5) * t251 + Ifges(7,6) * t250 + Ifges(7,3) * t302;
t225 = Ifges(7,1) * t251 + Ifges(7,4) * t250 + Ifges(7,5) * t302;
t187 = -mrSges(7,1) * t207 + mrSges(7,3) * t200 + Ifges(7,4) * t216 + Ifges(7,2) * t215 + Ifges(7,6) * t286 - t251 * t223 + t225 * t302;
t224 = Ifges(7,4) * t251 + Ifges(7,2) * t250 + Ifges(7,6) * t302;
t188 = mrSges(7,2) * t207 - mrSges(7,3) * t199 + Ifges(7,1) * t216 + Ifges(7,4) * t215 + Ifges(7,5) * t286 + t250 * t223 - t224 * t302;
t243 = Ifges(6,5) * t269 + Ifges(6,6) * t268 + Ifges(6,3) * t303;
t245 = Ifges(6,1) * t269 + Ifges(6,4) * t268 + Ifges(6,5) * t303;
t172 = -mrSges(6,1) * t220 + mrSges(6,3) * t205 + Ifges(6,4) * t233 + Ifges(6,2) * t232 + Ifges(6,6) * t288 - pkin(5) * t359 + pkin(10) * t361 + t340 * t187 + t335 * t188 - t269 * t243 + t303 * t245;
t244 = Ifges(6,4) * t269 + Ifges(6,2) * t268 + Ifges(6,6) * t303;
t173 = mrSges(6,2) * t220 - mrSges(6,3) * t204 + Ifges(6,1) * t233 + Ifges(6,4) * t232 + Ifges(6,5) * t288 - pkin(10) * t186 - t335 * t187 + t340 * t188 + t268 * t243 - t303 * t244;
t256 = Ifges(5,5) * t297 + Ifges(5,6) * t296 + Ifges(5,3) * t307;
t258 = Ifges(5,1) * t297 + Ifges(5,4) * t296 + Ifges(5,5) * t307;
t157 = -mrSges(5,1) * t236 + mrSges(5,3) * t222 + Ifges(5,4) * t263 + Ifges(5,2) * t262 + Ifges(5,6) * t292 - pkin(4) * t352 + pkin(9) * t362 + t341 * t172 + t336 * t173 - t297 * t256 + t307 * t258;
t257 = Ifges(5,4) * t297 + Ifges(5,2) * t296 + Ifges(5,6) * t307;
t159 = mrSges(5,2) * t236 - mrSges(5,3) * t221 + Ifges(5,1) * t263 + Ifges(5,4) * t262 + Ifges(5,5) * t292 - pkin(9) * t179 - t336 * t172 + t341 * t173 + t296 * t256 - t307 * t257;
t280 = Ifges(4,4) * t309 + Ifges(4,2) * t308 + Ifges(4,6) * qJD(2);
t281 = Ifges(4,1) * t309 + Ifges(4,4) * t308 + Ifges(4,5) * qJD(2);
t353 = -mrSges(4,1) * t248 + mrSges(4,2) * t249 - Ifges(4,5) * t294 - Ifges(4,6) * t293 - Ifges(4,3) * qJDD(2) - pkin(3) * t349 - pkin(8) * t363 - t342 * t157 - t337 * t159 - t309 * t280 + t308 * t281;
t372 = mrSges(3,1) * t298 - mrSges(3,2) * t299 + Ifges(3,5) * t320 + Ifges(3,6) * t321 + Ifges(3,3) * qJDD(2) + pkin(2) * t163 + (t338 * t311 - t343 * t312) * qJD(1) - t353;
t170 = t342 * t177 + t337 * t178;
t319 = (-mrSges(3,1) * t343 + mrSges(3,2) * t338) * qJD(1);
t325 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t368;
t161 = m(3) * t298 + qJDD(2) * mrSges(3,1) - t320 * mrSges(3,3) + qJD(2) * t325 - t319 * t369 + t163;
t324 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t369;
t364 = t334 * t168 - t333 * t190;
t162 = m(3) * t299 - qJDD(2) * mrSges(3,2) + t321 * mrSges(3,3) - qJD(2) * t324 + t319 * t368 + t364;
t365 = -t338 * t161 + t343 * t162;
t279 = Ifges(4,5) * t309 + Ifges(4,6) * t308 + Ifges(4,3) * qJD(2);
t151 = mrSges(4,2) * t275 - mrSges(4,3) * t248 + Ifges(4,1) * t294 + Ifges(4,4) * t293 + Ifges(4,5) * qJDD(2) - pkin(8) * t170 - qJD(2) * t280 - t337 * t157 + t342 * t159 + t308 * t279;
t355 = -mrSges(7,1) * t199 + mrSges(7,2) * t200 - Ifges(7,5) * t216 - Ifges(7,6) * t215 - Ifges(7,3) * t286 - t251 * t224 + t250 * t225;
t351 = -mrSges(6,1) * t204 + mrSges(6,2) * t205 - Ifges(6,5) * t233 - Ifges(6,6) * t232 - Ifges(6,3) * t288 - pkin(5) * t186 - t269 * t244 + t268 * t245 + t355;
t347 = mrSges(5,1) * t221 - mrSges(5,2) * t222 + Ifges(5,5) * t263 + Ifges(5,6) * t262 + Ifges(5,3) * t292 + pkin(4) * t179 + t297 * t257 - t296 * t258 - t351;
t155 = -mrSges(4,1) * t275 + mrSges(4,3) * t249 + Ifges(4,4) * t294 + Ifges(4,2) * t293 + Ifges(4,6) * qJDD(2) - pkin(3) * t170 + qJD(2) * t281 - t309 * t279 - t347;
t310 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t338 + Ifges(3,6) * t343) * qJD(1);
t314 = -t346 * pkin(7) + t357;
t354 = m(4) * t275 - t293 * mrSges(4,1) + t294 * mrSges(4,2) - t308 * t300 + t309 * t301 + t170;
t147 = -mrSges(3,1) * t314 + mrSges(3,3) * t299 + Ifges(3,4) * t320 + Ifges(3,2) * t321 + Ifges(3,6) * qJDD(2) - pkin(2) * t354 + qJ(3) * t364 + qJD(2) * t312 + t333 * t151 + t334 * t155 - t310 * t369;
t150 = mrSges(3,2) * t314 - mrSges(3,3) * t298 + Ifges(3,1) * t320 + Ifges(3,4) * t321 + Ifges(3,5) * qJDD(2) - qJ(3) * t163 - qJD(2) * t311 + t334 * t151 - t333 * t155 + t310 * t368;
t350 = -m(3) * t314 + t321 * mrSges(3,1) - t320 * mrSges(3,2) - t324 * t369 + t325 * t368 - t354;
t356 = mrSges(2,1) * t326 - mrSges(2,2) * t327 + Ifges(2,3) * qJDD(1) + pkin(1) * t350 + pkin(7) * t365 + t343 * t147 + t338 * t150;
t164 = m(2) * t326 + qJDD(1) * mrSges(2,1) - t346 * mrSges(2,2) + t350;
t154 = t343 * t161 + t338 * t162;
t152 = m(2) * t327 - t346 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t365;
t148 = mrSges(2,1) * g(3) + mrSges(2,3) * t327 + t346 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t154 - t372;
t145 = -mrSges(2,2) * g(3) - mrSges(2,3) * t326 + Ifges(2,5) * qJDD(1) - t346 * Ifges(2,6) - pkin(7) * t154 - t338 * t147 + t343 * t150;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t344 * t145 - t339 * t148 - pkin(6) * (t339 * t152 + t344 * t164), t145, t150, t151, t159, t173, t188; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t339 * t145 + t344 * t148 + pkin(6) * (t344 * t152 - t339 * t164), t148, t147, t155, t157, t172, t187; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t356, t356, t372, -t353, t347, -t351, -t355;];
m_new  = t1;
