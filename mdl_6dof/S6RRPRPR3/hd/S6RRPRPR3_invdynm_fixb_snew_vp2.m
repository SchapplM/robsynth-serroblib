% Calculate vector of cutting torques with Newton-Euler for
% S6RRPRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-05-06 13:20
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPRPR3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR3_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR3_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR3_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR3_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR3_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR3_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR3_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 13:14:34
% EndTime: 2019-05-06 13:15:38
% DurationCPUTime: 44.09s
% Computational Cost: add. (759319->385), mult. (1766360->489), div. (0->0), fcn. (1292746->12), ass. (0->151)
t339 = sin(qJ(2));
t343 = cos(qJ(2));
t366 = qJD(1) * qJD(2);
t320 = t339 * qJDD(1) + t343 * t366;
t340 = sin(qJ(1));
t344 = cos(qJ(1));
t327 = -t344 * g(1) - t340 * g(2);
t346 = qJD(1) ^ 2;
t315 = -t346 * pkin(1) + qJDD(1) * pkin(7) + t327;
t370 = t339 * t315;
t371 = pkin(2) * t346;
t273 = qJDD(2) * pkin(2) - t320 * qJ(3) - t370 + (qJ(3) * t366 + t339 * t371 - g(3)) * t343;
t300 = -t339 * g(3) + t343 * t315;
t321 = t343 * qJDD(1) - t339 * t366;
t369 = qJD(1) * t339;
t323 = qJD(2) * pkin(2) - qJ(3) * t369;
t332 = t343 ^ 2;
t274 = t321 * qJ(3) - qJD(2) * t323 - t332 * t371 + t300;
t334 = sin(pkin(10));
t336 = cos(pkin(10));
t309 = (t334 * t343 + t336 * t339) * qJD(1);
t248 = -0.2e1 * qJD(3) * t309 + t336 * t273 - t334 * t274;
t368 = qJD(1) * t343;
t308 = -t334 * t369 + t336 * t368;
t249 = 0.2e1 * qJD(3) * t308 + t334 * t273 + t336 * t274;
t285 = -t308 * mrSges(4,1) + t309 * mrSges(4,2);
t294 = -t334 * t320 + t336 * t321;
t302 = qJD(2) * mrSges(4,1) - t309 * mrSges(4,3);
t287 = -t308 * pkin(3) - t309 * pkin(8);
t345 = qJD(2) ^ 2;
t232 = -t345 * pkin(3) + qJDD(2) * pkin(8) + t308 * t287 + t249;
t326 = t340 * g(1) - t344 * g(2);
t357 = -qJDD(1) * pkin(1) - t326;
t277 = -t321 * pkin(2) + qJDD(3) + t323 * t369 + (-qJ(3) * t332 - pkin(7)) * t346 + t357;
t295 = t336 * t320 + t334 * t321;
t235 = (-qJD(2) * t308 - t295) * pkin(8) + (qJD(2) * t309 - t294) * pkin(3) + t277;
t338 = sin(qJ(4));
t342 = cos(qJ(4));
t221 = -t338 * t232 + t342 * t235;
t297 = t342 * qJD(2) - t338 * t309;
t263 = t297 * qJD(4) + t338 * qJDD(2) + t342 * t295;
t293 = qJDD(4) - t294;
t298 = t338 * qJD(2) + t342 * t309;
t307 = qJD(4) - t308;
t210 = (t297 * t307 - t263) * qJ(5) + (t297 * t298 + t293) * pkin(4) + t221;
t222 = t342 * t232 + t338 * t235;
t262 = -t298 * qJD(4) + t342 * qJDD(2) - t338 * t295;
t279 = t307 * pkin(4) - t298 * qJ(5);
t296 = t297 ^ 2;
t212 = -t296 * pkin(4) + t262 * qJ(5) - t307 * t279 + t222;
t333 = sin(pkin(11));
t335 = cos(pkin(11));
t271 = t333 * t297 + t335 * t298;
t204 = -0.2e1 * qJD(5) * t271 + t335 * t210 - t333 * t212;
t242 = t333 * t262 + t335 * t263;
t270 = t335 * t297 - t333 * t298;
t201 = (t270 * t307 - t242) * pkin(9) + (t270 * t271 + t293) * pkin(5) + t204;
t205 = 0.2e1 * qJD(5) * t270 + t333 * t210 + t335 * t212;
t241 = t335 * t262 - t333 * t263;
t255 = t307 * pkin(5) - t271 * pkin(9);
t267 = t270 ^ 2;
t202 = -t267 * pkin(5) + t241 * pkin(9) - t307 * t255 + t205;
t337 = sin(qJ(6));
t341 = cos(qJ(6));
t199 = t341 * t201 - t337 * t202;
t250 = t341 * t270 - t337 * t271;
t218 = t250 * qJD(6) + t337 * t241 + t341 * t242;
t251 = t337 * t270 + t341 * t271;
t228 = -t250 * mrSges(7,1) + t251 * mrSges(7,2);
t303 = qJD(6) + t307;
t236 = -t303 * mrSges(7,2) + t250 * mrSges(7,3);
t289 = qJDD(6) + t293;
t192 = m(7) * t199 + t289 * mrSges(7,1) - t218 * mrSges(7,3) - t251 * t228 + t303 * t236;
t200 = t337 * t201 + t341 * t202;
t217 = -t251 * qJD(6) + t341 * t241 - t337 * t242;
t237 = t303 * mrSges(7,1) - t251 * mrSges(7,3);
t193 = m(7) * t200 - t289 * mrSges(7,2) + t217 * mrSges(7,3) + t250 * t228 - t303 * t237;
t186 = t341 * t192 + t337 * t193;
t252 = -t270 * mrSges(6,1) + t271 * mrSges(6,2);
t253 = -t307 * mrSges(6,2) + t270 * mrSges(6,3);
t183 = m(6) * t204 + t293 * mrSges(6,1) - t242 * mrSges(6,3) - t271 * t252 + t307 * t253 + t186;
t254 = t307 * mrSges(6,1) - t271 * mrSges(6,3);
t361 = -t337 * t192 + t341 * t193;
t184 = m(6) * t205 - t293 * mrSges(6,2) + t241 * mrSges(6,3) + t270 * t252 - t307 * t254 + t361;
t179 = t335 * t183 + t333 * t184;
t275 = -t297 * mrSges(5,1) + t298 * mrSges(5,2);
t278 = -t307 * mrSges(5,2) + t297 * mrSges(5,3);
t177 = m(5) * t221 + t293 * mrSges(5,1) - t263 * mrSges(5,3) - t298 * t275 + t307 * t278 + t179;
t280 = t307 * mrSges(5,1) - t298 * mrSges(5,3);
t362 = -t333 * t183 + t335 * t184;
t178 = m(5) * t222 - t293 * mrSges(5,2) + t262 * mrSges(5,3) + t297 * t275 - t307 * t280 + t362;
t363 = -t338 * t177 + t342 * t178;
t168 = m(4) * t249 - qJDD(2) * mrSges(4,2) + t294 * mrSges(4,3) - qJD(2) * t302 + t308 * t285 + t363;
t301 = -qJD(2) * mrSges(4,2) + t308 * mrSges(4,3);
t231 = -qJDD(2) * pkin(3) - t345 * pkin(8) + t309 * t287 - t248;
t220 = -t262 * pkin(4) - t296 * qJ(5) + t298 * t279 + qJDD(5) + t231;
t207 = -t241 * pkin(5) - t267 * pkin(9) + t271 * t255 + t220;
t359 = m(7) * t207 - t217 * mrSges(7,1) + t218 * mrSges(7,2) - t250 * t236 + t251 * t237;
t352 = m(6) * t220 - t241 * mrSges(6,1) + t242 * mrSges(6,2) - t270 * t253 + t271 * t254 + t359;
t349 = -m(5) * t231 + t262 * mrSges(5,1) - t263 * mrSges(5,2) + t297 * t278 - t298 * t280 - t352;
t195 = m(4) * t248 + qJDD(2) * mrSges(4,1) - t295 * mrSges(4,3) + qJD(2) * t301 - t309 * t285 + t349;
t163 = t334 * t168 + t336 * t195;
t299 = -t343 * g(3) - t370;
t311 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t339 + Ifges(3,2) * t343) * qJD(1);
t312 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t339 + Ifges(3,4) * t343) * qJD(1);
t223 = Ifges(7,5) * t251 + Ifges(7,6) * t250 + Ifges(7,3) * t303;
t225 = Ifges(7,1) * t251 + Ifges(7,4) * t250 + Ifges(7,5) * t303;
t187 = -mrSges(7,1) * t207 + mrSges(7,3) * t200 + Ifges(7,4) * t218 + Ifges(7,2) * t217 + Ifges(7,6) * t289 - t251 * t223 + t303 * t225;
t224 = Ifges(7,4) * t251 + Ifges(7,2) * t250 + Ifges(7,6) * t303;
t188 = mrSges(7,2) * t207 - mrSges(7,3) * t199 + Ifges(7,1) * t218 + Ifges(7,4) * t217 + Ifges(7,5) * t289 + t250 * t223 - t303 * t224;
t243 = Ifges(6,5) * t271 + Ifges(6,6) * t270 + Ifges(6,3) * t307;
t245 = Ifges(6,1) * t271 + Ifges(6,4) * t270 + Ifges(6,5) * t307;
t172 = -mrSges(6,1) * t220 + mrSges(6,3) * t205 + Ifges(6,4) * t242 + Ifges(6,2) * t241 + Ifges(6,6) * t293 - pkin(5) * t359 + pkin(9) * t361 + t341 * t187 + t337 * t188 - t271 * t243 + t307 * t245;
t244 = Ifges(6,4) * t271 + Ifges(6,2) * t270 + Ifges(6,6) * t307;
t173 = mrSges(6,2) * t220 - mrSges(6,3) * t204 + Ifges(6,1) * t242 + Ifges(6,4) * t241 + Ifges(6,5) * t293 - pkin(9) * t186 - t337 * t187 + t341 * t188 + t270 * t243 - t307 * t244;
t256 = Ifges(5,5) * t298 + Ifges(5,6) * t297 + Ifges(5,3) * t307;
t258 = Ifges(5,1) * t298 + Ifges(5,4) * t297 + Ifges(5,5) * t307;
t157 = -mrSges(5,1) * t231 + mrSges(5,3) * t222 + Ifges(5,4) * t263 + Ifges(5,2) * t262 + Ifges(5,6) * t293 - pkin(4) * t352 + qJ(5) * t362 + t335 * t172 + t333 * t173 - t298 * t256 + t307 * t258;
t257 = Ifges(5,4) * t298 + Ifges(5,2) * t297 + Ifges(5,6) * t307;
t159 = mrSges(5,2) * t231 - mrSges(5,3) * t221 + Ifges(5,1) * t263 + Ifges(5,4) * t262 + Ifges(5,5) * t293 - qJ(5) * t179 - t333 * t172 + t335 * t173 + t297 * t256 - t307 * t257;
t282 = Ifges(4,4) * t309 + Ifges(4,2) * t308 + Ifges(4,6) * qJD(2);
t283 = Ifges(4,1) * t309 + Ifges(4,4) * t308 + Ifges(4,5) * qJD(2);
t353 = -mrSges(4,1) * t248 + mrSges(4,2) * t249 - Ifges(4,5) * t295 - Ifges(4,6) * t294 - Ifges(4,3) * qJDD(2) - pkin(3) * t349 - pkin(8) * t363 - t342 * t157 - t338 * t159 - t309 * t282 + t308 * t283;
t372 = mrSges(3,1) * t299 - mrSges(3,2) * t300 + Ifges(3,5) * t320 + Ifges(3,6) * t321 + Ifges(3,3) * qJDD(2) + pkin(2) * t163 + (t339 * t311 - t343 * t312) * qJD(1) - t353;
t170 = t342 * t177 + t338 * t178;
t319 = (-mrSges(3,1) * t343 + mrSges(3,2) * t339) * qJD(1);
t325 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t368;
t161 = m(3) * t299 + qJDD(2) * mrSges(3,1) - t320 * mrSges(3,3) + qJD(2) * t325 - t319 * t369 + t163;
t324 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t369;
t364 = t336 * t168 - t334 * t195;
t162 = m(3) * t300 - qJDD(2) * mrSges(3,2) + t321 * mrSges(3,3) - qJD(2) * t324 + t319 * t368 + t364;
t365 = -t339 * t161 + t343 * t162;
t281 = Ifges(4,5) * t309 + Ifges(4,6) * t308 + Ifges(4,3) * qJD(2);
t151 = mrSges(4,2) * t277 - mrSges(4,3) * t248 + Ifges(4,1) * t295 + Ifges(4,4) * t294 + Ifges(4,5) * qJDD(2) - pkin(8) * t170 - qJD(2) * t282 - t338 * t157 + t342 * t159 + t308 * t281;
t355 = -mrSges(7,1) * t199 + mrSges(7,2) * t200 - Ifges(7,5) * t218 - Ifges(7,6) * t217 - Ifges(7,3) * t289 - t251 * t224 + t250 * t225;
t351 = -mrSges(6,1) * t204 + mrSges(6,2) * t205 - Ifges(6,5) * t242 - Ifges(6,6) * t241 - Ifges(6,3) * t293 - pkin(5) * t186 - t271 * t244 + t270 * t245 + t355;
t347 = mrSges(5,1) * t221 - mrSges(5,2) * t222 + Ifges(5,5) * t263 + Ifges(5,6) * t262 + Ifges(5,3) * t293 + pkin(4) * t179 + t298 * t257 - t297 * t258 - t351;
t155 = -mrSges(4,1) * t277 + mrSges(4,3) * t249 + Ifges(4,4) * t295 + Ifges(4,2) * t294 + Ifges(4,6) * qJDD(2) - pkin(3) * t170 + qJD(2) * t283 - t309 * t281 - t347;
t310 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t339 + Ifges(3,6) * t343) * qJD(1);
t314 = -t346 * pkin(7) + t357;
t354 = m(4) * t277 - t294 * mrSges(4,1) + t295 * mrSges(4,2) - t308 * t301 + t309 * t302 + t170;
t147 = -mrSges(3,1) * t314 + mrSges(3,3) * t300 + Ifges(3,4) * t320 + Ifges(3,2) * t321 + Ifges(3,6) * qJDD(2) - pkin(2) * t354 + qJ(3) * t364 + qJD(2) * t312 + t334 * t151 + t336 * t155 - t310 * t369;
t150 = mrSges(3,2) * t314 - mrSges(3,3) * t299 + Ifges(3,1) * t320 + Ifges(3,4) * t321 + Ifges(3,5) * qJDD(2) - qJ(3) * t163 - qJD(2) * t311 + t336 * t151 - t334 * t155 + t310 * t368;
t350 = -m(3) * t314 + t321 * mrSges(3,1) - t320 * mrSges(3,2) - t324 * t369 + t325 * t368 - t354;
t356 = mrSges(2,1) * t326 - mrSges(2,2) * t327 + Ifges(2,3) * qJDD(1) + pkin(1) * t350 + pkin(7) * t365 + t343 * t147 + t339 * t150;
t164 = m(2) * t326 + qJDD(1) * mrSges(2,1) - t346 * mrSges(2,2) + t350;
t154 = t343 * t161 + t339 * t162;
t152 = m(2) * t327 - t346 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t365;
t148 = mrSges(2,1) * g(3) + mrSges(2,3) * t327 + t346 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t154 - t372;
t145 = -mrSges(2,2) * g(3) - mrSges(2,3) * t326 + Ifges(2,5) * qJDD(1) - t346 * Ifges(2,6) - pkin(7) * t154 - t339 * t147 + t343 * t150;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t344 * t145 - t340 * t148 - pkin(6) * (t340 * t152 + t344 * t164), t145, t150, t151, t159, t173, t188; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t340 * t145 + t344 * t148 + pkin(6) * (t344 * t152 - t340 * t164), t148, t147, t155, t157, t172, t187; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t356, t356, t372, -t353, t347, -t351, -t355;];
m_new  = t1;
