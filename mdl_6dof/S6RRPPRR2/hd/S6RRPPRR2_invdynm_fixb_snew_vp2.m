% Calculate vector of cutting torques with Newton-Euler for
% S6RRPPRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3,theta4]';
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
% Datum: 2019-05-06 09:49
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPPRR2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR2_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR2_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR2_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR2_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR2_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR2_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR2_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 09:44:05
% EndTime: 2019-05-06 09:45:13
% DurationCPUTime: 42.93s
% Computational Cost: add. (727520->386), mult. (1728582->488), div. (0->0), fcn. (1261750->12), ass. (0->152)
t374 = -2 * qJD(3);
t339 = sin(qJ(2));
t343 = cos(qJ(2));
t366 = qJD(1) * qJD(2);
t321 = t339 * qJDD(1) + t343 * t366;
t340 = sin(qJ(1));
t344 = cos(qJ(1));
t328 = -t344 * g(1) - t340 * g(2);
t346 = qJD(1) ^ 2;
t316 = -t346 * pkin(1) + qJDD(1) * pkin(7) + t328;
t370 = t339 * t316;
t372 = pkin(2) * t346;
t269 = qJDD(2) * pkin(2) - t321 * qJ(3) - t370 + (qJ(3) * t366 + t339 * t372 - g(3)) * t343;
t301 = -t339 * g(3) + t343 * t316;
t322 = t343 * qJDD(1) - t339 * t366;
t369 = qJD(1) * t339;
t324 = qJD(2) * pkin(2) - qJ(3) * t369;
t333 = t343 ^ 2;
t271 = t322 * qJ(3) - qJD(2) * t324 - t333 * t372 + t301;
t335 = sin(pkin(10));
t371 = cos(pkin(10));
t310 = (t335 * t343 + t371 * t339) * qJD(1);
t249 = t371 * t269 - t335 * t271 + t310 * t374;
t368 = qJD(1) * t343;
t309 = t335 * t369 - t371 * t368;
t250 = t335 * t269 + t371 * t271 + t309 * t374;
t285 = t309 * mrSges(4,1) + t310 * mrSges(4,2);
t293 = t335 * t321 - t371 * t322;
t303 = qJD(2) * mrSges(4,1) - t310 * mrSges(4,3);
t284 = t309 * pkin(3) - t310 * qJ(4);
t345 = qJD(2) ^ 2;
t233 = -t345 * pkin(3) + qJDD(2) * qJ(4) - t309 * t284 + t250;
t327 = t340 * g(1) - t344 * g(2);
t357 = -qJDD(1) * pkin(1) - t327;
t273 = -t322 * pkin(2) + qJDD(3) + t324 * t369 + (-qJ(3) * t333 - pkin(7)) * t346 + t357;
t294 = t371 * t321 + t335 * t322;
t236 = (qJD(2) * t309 - t294) * qJ(4) + (qJD(2) * t310 + t293) * pkin(3) + t273;
t334 = sin(pkin(11));
t336 = cos(pkin(11));
t299 = t334 * qJD(2) + t336 * t310;
t221 = -0.2e1 * qJD(4) * t299 - t334 * t233 + t336 * t236;
t279 = t334 * qJDD(2) + t336 * t294;
t298 = t336 * qJD(2) - t334 * t310;
t211 = (t298 * t309 - t279) * pkin(8) + (t298 * t299 + t293) * pkin(4) + t221;
t222 = 0.2e1 * qJD(4) * t298 + t336 * t233 + t334 * t236;
t276 = t309 * pkin(4) - t299 * pkin(8);
t278 = t336 * qJDD(2) - t334 * t294;
t297 = t298 ^ 2;
t213 = -t297 * pkin(4) + t278 * pkin(8) - t309 * t276 + t222;
t338 = sin(qJ(5));
t342 = cos(qJ(5));
t205 = t342 * t211 - t338 * t213;
t264 = t342 * t298 - t338 * t299;
t241 = t264 * qJD(5) + t338 * t278 + t342 * t279;
t265 = t338 * t298 + t342 * t299;
t292 = qJDD(5) + t293;
t308 = qJD(5) + t309;
t202 = (t264 * t308 - t241) * pkin(9) + (t264 * t265 + t292) * pkin(5) + t205;
t206 = t338 * t211 + t342 * t213;
t240 = -t265 * qJD(5) + t342 * t278 - t338 * t279;
t256 = t308 * pkin(5) - t265 * pkin(9);
t263 = t264 ^ 2;
t203 = -t263 * pkin(5) + t240 * pkin(9) - t308 * t256 + t206;
t337 = sin(qJ(6));
t341 = cos(qJ(6));
t200 = t341 * t202 - t337 * t203;
t251 = t341 * t264 - t337 * t265;
t219 = t251 * qJD(6) + t337 * t240 + t341 * t241;
t252 = t337 * t264 + t341 * t265;
t229 = -t251 * mrSges(7,1) + t252 * mrSges(7,2);
t304 = qJD(6) + t308;
t242 = -t304 * mrSges(7,2) + t251 * mrSges(7,3);
t288 = qJDD(6) + t292;
t193 = m(7) * t200 + t288 * mrSges(7,1) - t219 * mrSges(7,3) - t252 * t229 + t304 * t242;
t201 = t337 * t202 + t341 * t203;
t218 = -t252 * qJD(6) + t341 * t240 - t337 * t241;
t243 = t304 * mrSges(7,1) - t252 * mrSges(7,3);
t194 = m(7) * t201 - t288 * mrSges(7,2) + t218 * mrSges(7,3) + t251 * t229 - t304 * t243;
t187 = t341 * t193 + t337 * t194;
t253 = -t264 * mrSges(6,1) + t265 * mrSges(6,2);
t254 = -t308 * mrSges(6,2) + t264 * mrSges(6,3);
t184 = m(6) * t205 + t292 * mrSges(6,1) - t241 * mrSges(6,3) - t265 * t253 + t308 * t254 + t187;
t255 = t308 * mrSges(6,1) - t265 * mrSges(6,3);
t361 = -t337 * t193 + t341 * t194;
t185 = m(6) * t206 - t292 * mrSges(6,2) + t240 * mrSges(6,3) + t264 * t253 - t308 * t255 + t361;
t180 = t342 * t184 + t338 * t185;
t270 = -t298 * mrSges(5,1) + t299 * mrSges(5,2);
t274 = -t309 * mrSges(5,2) + t298 * mrSges(5,3);
t178 = m(5) * t221 + t293 * mrSges(5,1) - t279 * mrSges(5,3) - t299 * t270 + t309 * t274 + t180;
t275 = t309 * mrSges(5,1) - t299 * mrSges(5,3);
t362 = -t338 * t184 + t342 * t185;
t179 = m(5) * t222 - t293 * mrSges(5,2) + t278 * mrSges(5,3) + t298 * t270 - t309 * t275 + t362;
t363 = -t334 * t178 + t336 * t179;
t169 = m(4) * t250 - qJDD(2) * mrSges(4,2) - t293 * mrSges(4,3) - qJD(2) * t303 - t309 * t285 + t363;
t302 = -qJD(2) * mrSges(4,2) - t309 * mrSges(4,3);
t232 = -qJDD(2) * pkin(3) - t345 * qJ(4) + t310 * t284 + qJDD(4) - t249;
t223 = -t278 * pkin(4) - t297 * pkin(8) + t299 * t276 + t232;
t208 = -t240 * pkin(5) - t263 * pkin(9) + t265 * t256 + t223;
t359 = m(7) * t208 - t218 * mrSges(7,1) + t219 * mrSges(7,2) - t251 * t242 + t252 * t243;
t352 = m(6) * t223 - t240 * mrSges(6,1) + t241 * mrSges(6,2) - t264 * t254 + t265 * t255 + t359;
t349 = -m(5) * t232 + t278 * mrSges(5,1) - t279 * mrSges(5,2) + t298 * t274 - t299 * t275 - t352;
t196 = m(4) * t249 + qJDD(2) * mrSges(4,1) - t294 * mrSges(4,3) + qJD(2) * t302 - t310 * t285 + t349;
t164 = t335 * t169 + t371 * t196;
t300 = -t343 * g(3) - t370;
t312 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t339 + Ifges(3,2) * t343) * qJD(1);
t313 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t339 + Ifges(3,4) * t343) * qJD(1);
t224 = Ifges(7,5) * t252 + Ifges(7,6) * t251 + Ifges(7,3) * t304;
t226 = Ifges(7,1) * t252 + Ifges(7,4) * t251 + Ifges(7,5) * t304;
t188 = -mrSges(7,1) * t208 + mrSges(7,3) * t201 + Ifges(7,4) * t219 + Ifges(7,2) * t218 + Ifges(7,6) * t288 - t252 * t224 + t304 * t226;
t225 = Ifges(7,4) * t252 + Ifges(7,2) * t251 + Ifges(7,6) * t304;
t189 = mrSges(7,2) * t208 - mrSges(7,3) * t200 + Ifges(7,1) * t219 + Ifges(7,4) * t218 + Ifges(7,5) * t288 + t251 * t224 - t304 * t225;
t244 = Ifges(6,5) * t265 + Ifges(6,6) * t264 + Ifges(6,3) * t308;
t246 = Ifges(6,1) * t265 + Ifges(6,4) * t264 + Ifges(6,5) * t308;
t173 = -mrSges(6,1) * t223 + mrSges(6,3) * t206 + Ifges(6,4) * t241 + Ifges(6,2) * t240 + Ifges(6,6) * t292 - pkin(5) * t359 + pkin(9) * t361 + t341 * t188 + t337 * t189 - t265 * t244 + t308 * t246;
t245 = Ifges(6,4) * t265 + Ifges(6,2) * t264 + Ifges(6,6) * t308;
t174 = mrSges(6,2) * t223 - mrSges(6,3) * t205 + Ifges(6,1) * t241 + Ifges(6,4) * t240 + Ifges(6,5) * t292 - pkin(9) * t187 - t337 * t188 + t341 * t189 + t264 * t244 - t308 * t245;
t257 = Ifges(5,5) * t299 + Ifges(5,6) * t298 + Ifges(5,3) * t309;
t259 = Ifges(5,1) * t299 + Ifges(5,4) * t298 + Ifges(5,5) * t309;
t158 = -mrSges(5,1) * t232 + mrSges(5,3) * t222 + Ifges(5,4) * t279 + Ifges(5,2) * t278 + Ifges(5,6) * t293 - pkin(4) * t352 + pkin(8) * t362 + t342 * t173 + t338 * t174 - t299 * t257 + t309 * t259;
t258 = Ifges(5,4) * t299 + Ifges(5,2) * t298 + Ifges(5,6) * t309;
t160 = mrSges(5,2) * t232 - mrSges(5,3) * t221 + Ifges(5,1) * t279 + Ifges(5,4) * t278 + Ifges(5,5) * t293 - pkin(8) * t180 - t338 * t173 + t342 * t174 + t298 * t257 - t309 * t258;
t281 = Ifges(4,4) * t310 - Ifges(4,2) * t309 + Ifges(4,6) * qJD(2);
t282 = Ifges(4,1) * t310 - Ifges(4,4) * t309 + Ifges(4,5) * qJD(2);
t353 = -mrSges(4,1) * t249 + mrSges(4,2) * t250 - Ifges(4,5) * t294 + Ifges(4,6) * t293 - Ifges(4,3) * qJDD(2) - pkin(3) * t349 - qJ(4) * t363 - t336 * t158 - t334 * t160 - t310 * t281 - t309 * t282;
t373 = mrSges(3,1) * t300 - mrSges(3,2) * t301 + Ifges(3,5) * t321 + Ifges(3,6) * t322 + Ifges(3,3) * qJDD(2) + pkin(2) * t164 + (t339 * t312 - t343 * t313) * qJD(1) - t353;
t171 = t336 * t178 + t334 * t179;
t320 = (-mrSges(3,1) * t343 + mrSges(3,2) * t339) * qJD(1);
t326 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t368;
t162 = m(3) * t300 + qJDD(2) * mrSges(3,1) - t321 * mrSges(3,3) + qJD(2) * t326 - t320 * t369 + t164;
t325 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t369;
t364 = t371 * t169 - t335 * t196;
t163 = m(3) * t301 - qJDD(2) * mrSges(3,2) + t322 * mrSges(3,3) - qJD(2) * t325 + t320 * t368 + t364;
t365 = -t339 * t162 + t343 * t163;
t280 = Ifges(4,5) * t310 - Ifges(4,6) * t309 + Ifges(4,3) * qJD(2);
t152 = mrSges(4,2) * t273 - mrSges(4,3) * t249 + Ifges(4,1) * t294 - Ifges(4,4) * t293 + Ifges(4,5) * qJDD(2) - qJ(4) * t171 - qJD(2) * t281 - t334 * t158 + t336 * t160 - t309 * t280;
t355 = -mrSges(7,1) * t200 + mrSges(7,2) * t201 - Ifges(7,5) * t219 - Ifges(7,6) * t218 - Ifges(7,3) * t288 - t252 * t225 + t251 * t226;
t351 = -mrSges(6,1) * t205 + mrSges(6,2) * t206 - Ifges(6,5) * t241 - Ifges(6,6) * t240 - Ifges(6,3) * t292 - pkin(5) * t187 - t265 * t245 + t264 * t246 + t355;
t347 = mrSges(5,1) * t221 - mrSges(5,2) * t222 + Ifges(5,5) * t279 + Ifges(5,6) * t278 + pkin(4) * t180 + t299 * t258 - t298 * t259 - t351;
t156 = -t347 + (-Ifges(5,3) - Ifges(4,2)) * t293 - t310 * t280 + Ifges(4,4) * t294 + qJD(2) * t282 - mrSges(4,1) * t273 + mrSges(4,3) * t250 - pkin(3) * t171 + Ifges(4,6) * qJDD(2);
t311 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t339 + Ifges(3,6) * t343) * qJD(1);
t315 = -t346 * pkin(7) + t357;
t354 = m(4) * t273 + t293 * mrSges(4,1) + t294 * mrSges(4,2) + t309 * t302 + t310 * t303 + t171;
t148 = -mrSges(3,1) * t315 + mrSges(3,3) * t301 + Ifges(3,4) * t321 + Ifges(3,2) * t322 + Ifges(3,6) * qJDD(2) - pkin(2) * t354 + qJ(3) * t364 + qJD(2) * t313 + t335 * t152 + t371 * t156 - t311 * t369;
t151 = mrSges(3,2) * t315 - mrSges(3,3) * t300 + Ifges(3,1) * t321 + Ifges(3,4) * t322 + Ifges(3,5) * qJDD(2) - qJ(3) * t164 - qJD(2) * t312 + t371 * t152 - t335 * t156 + t311 * t368;
t350 = -m(3) * t315 + t322 * mrSges(3,1) - t321 * mrSges(3,2) - t325 * t369 + t326 * t368 - t354;
t356 = mrSges(2,1) * t327 - mrSges(2,2) * t328 + Ifges(2,3) * qJDD(1) + pkin(1) * t350 + pkin(7) * t365 + t343 * t148 + t339 * t151;
t165 = m(2) * t327 + qJDD(1) * mrSges(2,1) - t346 * mrSges(2,2) + t350;
t155 = t343 * t162 + t339 * t163;
t153 = m(2) * t328 - t346 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t365;
t149 = mrSges(2,1) * g(3) + mrSges(2,3) * t328 + t346 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t155 - t373;
t146 = -mrSges(2,2) * g(3) - mrSges(2,3) * t327 + Ifges(2,5) * qJDD(1) - t346 * Ifges(2,6) - pkin(7) * t155 - t339 * t148 + t343 * t151;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t344 * t146 - t340 * t149 - pkin(6) * (t340 * t153 + t344 * t165), t146, t151, t152, t160, t174, t189; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t340 * t146 + t344 * t149 + pkin(6) * (t344 * t153 - t340 * t165), t149, t148, t156, t158, t173, t188; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t356, t356, t373, -t353, Ifges(5,3) * t293 + t347, -t351, -t355;];
m_new  = t1;
