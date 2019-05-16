% Calculate vector of cutting torques with Newton-Euler for
% S6PRRPPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
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
% Datum: 2019-05-05 03:15
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRRPPR4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR4_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR4_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPPR4_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR4_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR4_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR4_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPPR4_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 03:08:30
% EndTime: 2019-05-05 03:08:50
% DurationCPUTime: 10.75s
% Computational Cost: add. (181342->340), mult. (378111->417), div. (0->0), fcn. (253103->12), ass. (0->137)
t357 = -2 * qJD(4);
t307 = sin(pkin(10));
t309 = cos(pkin(10));
t294 = g(1) * t307 - g(2) * t309;
t295 = -g(1) * t309 - g(2) * t307;
t305 = -g(3) + qJDD(1);
t316 = cos(qJ(2));
t310 = cos(pkin(6));
t313 = sin(qJ(2));
t344 = t310 * t313;
t308 = sin(pkin(6));
t345 = t308 * t313;
t238 = t294 * t344 + t316 * t295 + t305 * t345;
t318 = qJD(2) ^ 2;
t232 = -pkin(2) * t318 + qJDD(2) * pkin(8) + t238;
t267 = -t294 * t308 + t305 * t310;
t312 = sin(qJ(3));
t315 = cos(qJ(3));
t226 = t315 * t232 + t312 * t267;
t289 = (-pkin(3) * t315 - qJ(4) * t312) * qJD(2);
t317 = qJD(3) ^ 2;
t339 = qJD(2) * t315;
t211 = -pkin(3) * t317 + qJDD(3) * qJ(4) + t289 * t339 + t226;
t237 = -t313 * t295 + (t294 * t310 + t305 * t308) * t316;
t231 = -qJDD(2) * pkin(2) - t318 * pkin(8) - t237;
t337 = qJD(2) * qJD(3);
t335 = t315 * t337;
t291 = qJDD(2) * t312 + t335;
t301 = t312 * t337;
t292 = qJDD(2) * t315 - t301;
t217 = (-t291 - t335) * qJ(4) + (-t292 + t301) * pkin(3) + t231;
t306 = sin(pkin(11));
t340 = qJD(2) * t312;
t349 = cos(pkin(11));
t279 = -t349 * qJD(3) + t306 * t340;
t204 = t349 * t211 + t306 * t217 + t279 * t357;
t280 = t306 * qJD(3) + t349 * t340;
t262 = -mrSges(5,1) * t339 - mrSges(5,3) * t280;
t263 = mrSges(6,1) * t339 + mrSges(6,2) * t280;
t264 = -t349 * qJDD(3) + t291 * t306;
t203 = -t306 * t211 + t349 * t217 + t280 * t357;
t249 = pkin(4) * t279 - qJ(5) * t280;
t346 = t315 ^ 2 * t318;
t201 = t292 * pkin(4) - qJ(5) * t346 + t280 * t249 + qJDD(5) - t203;
t265 = t306 * qJDD(3) + t349 * t291;
t336 = t279 * t339;
t194 = (-t265 + t336) * pkin(9) + (t279 * t280 + t292) * pkin(5) + t201;
t352 = -2 * qJD(5);
t198 = -pkin(4) * t346 - t292 * qJ(5) - t279 * t249 + t339 * t352 + t204;
t266 = pkin(5) * t339 - pkin(9) * t280;
t277 = t279 ^ 2;
t195 = -pkin(5) * t277 + pkin(9) * t264 - t266 * t339 + t198;
t311 = sin(qJ(6));
t314 = cos(qJ(6));
t191 = t194 * t314 - t195 * t311;
t247 = t279 * t314 - t280 * t311;
t219 = qJD(6) * t247 + t264 * t311 + t265 * t314;
t248 = t279 * t311 + t280 * t314;
t227 = -mrSges(7,1) * t247 + mrSges(7,2) * t248;
t299 = qJD(6) + t339;
t233 = -mrSges(7,2) * t299 + mrSges(7,3) * t247;
t286 = qJDD(6) + t292;
t186 = m(7) * t191 + mrSges(7,1) * t286 - mrSges(7,3) * t219 - t227 * t248 + t233 * t299;
t192 = t194 * t311 + t195 * t314;
t218 = -qJD(6) * t248 + t264 * t314 - t265 * t311;
t234 = mrSges(7,1) * t299 - mrSges(7,3) * t248;
t187 = m(7) * t192 - mrSges(7,2) * t286 + mrSges(7,3) * t218 + t227 * t247 - t234 * t299;
t178 = -t311 * t186 + t314 * t187;
t329 = m(6) * t198 - t292 * mrSges(6,3) + t178;
t250 = mrSges(6,1) * t279 - mrSges(6,3) * t280;
t341 = -mrSges(5,1) * t279 - mrSges(5,2) * t280 - t250;
t350 = -mrSges(5,3) - mrSges(6,2);
t173 = m(5) * t204 + t292 * mrSges(5,2) + t341 * t279 + t350 * t264 + (t262 - t263) * t339 + t329;
t260 = -mrSges(6,2) * t279 - mrSges(6,3) * t339;
t261 = mrSges(5,2) * t339 - mrSges(5,3) * t279;
t177 = t314 * t186 + t311 * t187;
t327 = -m(6) * t201 - t292 * mrSges(6,1) - t177;
t174 = m(5) * t203 - t292 * mrSges(5,1) + t341 * t280 + t350 * t265 + (-t260 - t261) * t339 + t327;
t171 = t349 * t173 - t174 * t306;
t290 = (-mrSges(4,1) * t315 + mrSges(4,2) * t312) * qJD(2);
t296 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t340;
t169 = m(4) * t226 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t292 - qJD(3) * t296 + t290 * t339 + t171;
t225 = -t312 * t232 + t315 * t267;
t325 = qJDD(3) * pkin(3) + t317 * qJ(4) - t289 * t340 - qJDD(4) + t225;
t355 = -(t265 + t336) * qJ(5) - t325;
t202 = -t277 * pkin(9) + (-pkin(4) - pkin(5)) * t264 + (pkin(4) * t339 + (2 * qJD(5)) + t266) * t280 - t355;
t193 = -m(7) * t202 + t218 * mrSges(7,1) - t219 * mrSges(7,2) + t247 * t233 - t248 * t234;
t205 = t280 * t352 + (-t280 * t339 + t264) * pkin(4) + t355;
t188 = m(6) * t205 + t264 * mrSges(6,1) - t265 * mrSges(6,3) + t279 * t260 - t280 * t263 + t193;
t184 = m(5) * t325 - t264 * mrSges(5,1) - t265 * mrSges(5,2) - t279 * t261 - t280 * t262 - t188;
t297 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t339;
t183 = m(4) * t225 + qJDD(3) * mrSges(4,1) - t291 * mrSges(4,3) + qJD(3) * t297 - t290 * t340 + t184;
t164 = t312 * t169 + t315 * t183;
t243 = Ifges(6,1) * t280 - Ifges(6,4) * t339 + Ifges(6,5) * t279;
t244 = Ifges(5,1) * t280 - Ifges(5,4) * t279 - Ifges(5,5) * t339;
t220 = Ifges(7,5) * t248 + Ifges(7,6) * t247 + Ifges(7,3) * t299;
t222 = Ifges(7,1) * t248 + Ifges(7,4) * t247 + Ifges(7,5) * t299;
t180 = -mrSges(7,1) * t202 + mrSges(7,3) * t192 + Ifges(7,4) * t219 + Ifges(7,2) * t218 + Ifges(7,6) * t286 - t220 * t248 + t222 * t299;
t221 = Ifges(7,4) * t248 + Ifges(7,2) * t247 + Ifges(7,6) * t299;
t181 = mrSges(7,2) * t202 - mrSges(7,3) * t191 + Ifges(7,1) * t219 + Ifges(7,4) * t218 + Ifges(7,5) * t286 + t220 * t247 - t221 * t299;
t323 = -mrSges(6,1) * t205 + mrSges(6,2) * t198 - pkin(5) * t193 - pkin(9) * t178 - t314 * t180 - t311 * t181;
t241 = Ifges(6,4) * t280 - Ifges(6,2) * t339 + Ifges(6,6) * t279;
t343 = -Ifges(5,5) * t280 + Ifges(5,6) * t279 + Ifges(5,3) * t339 - t241;
t158 = mrSges(5,1) * t325 + mrSges(5,3) * t204 - pkin(4) * t188 + (-Ifges(5,6) + Ifges(6,6)) * t292 + t343 * t280 + (Ifges(5,4) - Ifges(6,5)) * t265 + (-Ifges(5,2) - Ifges(6,3)) * t264 + (-t243 - t244) * t339 + t323;
t324 = mrSges(6,2) * t201 - mrSges(6,3) * t205 + Ifges(6,1) * t265 - Ifges(6,4) * t292 + Ifges(6,5) * t264 - pkin(9) * t177 - t180 * t311 + t314 * t181;
t239 = Ifges(6,5) * t280 - Ifges(6,6) * t339 + Ifges(6,3) * t279;
t342 = Ifges(5,4) * t280 - Ifges(5,2) * t279 - Ifges(5,6) * t339 - t239;
t162 = -mrSges(5,2) * t325 - mrSges(5,3) * t203 + Ifges(5,1) * t265 - Ifges(5,4) * t264 - Ifges(5,5) * t292 - qJ(5) * t188 + t343 * t279 + t342 * t339 + t324;
t273 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t312 + Ifges(4,2) * t315) * qJD(2);
t274 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t312 + Ifges(4,4) * t315) * qJD(2);
t353 = mrSges(4,1) * t225 - mrSges(4,2) * t226 + Ifges(4,5) * t291 + Ifges(4,6) * t292 + Ifges(4,3) * qJDD(3) + pkin(3) * t184 + qJ(4) * t171 + (t273 * t312 - t274 * t315) * qJD(2) + t349 * t158 + t306 * t162;
t148 = -mrSges(3,1) * t267 + mrSges(3,3) * t238 + t318 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t164 - t353;
t334 = t315 * t169 - t183 * t312;
t161 = m(3) * t238 - mrSges(3,1) * t318 - qJDD(2) * mrSges(3,2) + t334;
t170 = t306 * t173 + t349 * t174;
t321 = -m(4) * t231 + t292 * mrSges(4,1) - t291 * mrSges(4,2) - t296 * t340 + t297 * t339 - t170;
t166 = m(3) * t237 + qJDD(2) * mrSges(3,1) - t318 * mrSges(3,2) + t321;
t157 = t316 * t161 - t166 * t313;
t356 = pkin(7) * t157 + t148 * t316;
t328 = mrSges(7,1) * t191 - mrSges(7,2) * t192 + Ifges(7,5) * t219 + Ifges(7,6) * t218 + Ifges(7,3) * t286 + t248 * t221 - t247 * t222;
t322 = mrSges(6,1) * t201 - mrSges(6,3) * t198 - Ifges(6,4) * t265 + Ifges(6,2) * t292 - Ifges(6,6) * t264 + pkin(5) * t177 - t279 * t243 + t328;
t354 = -t342 * t280 - mrSges(5,1) * t203 + mrSges(5,2) * t204 - Ifges(5,5) * t265 + Ifges(5,6) * t264 - pkin(4) * (-t265 * mrSges(6,2) - t280 * t250 - t260 * t339 + t327) - qJ(5) * (-t264 * mrSges(6,2) - t279 * t250 - t263 * t339 + t329) - t279 * t244 + t322;
t347 = t166 * t316;
t163 = m(3) * t267 + t164;
t153 = t161 * t344 - t163 * t308 + t310 * t347;
t272 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t312 + Ifges(4,6) * t315) * qJD(2);
t149 = mrSges(4,2) * t231 - mrSges(4,3) * t225 + Ifges(4,1) * t291 + Ifges(4,4) * t292 + Ifges(4,5) * qJDD(3) - qJ(4) * t170 - qJD(3) * t273 - t306 * t158 + t349 * t162 + t272 * t339;
t154 = -t272 * t340 + Ifges(4,6) * qJDD(3) - pkin(3) * t170 + (Ifges(4,2) + Ifges(5,3)) * t292 + mrSges(4,3) * t226 - mrSges(4,1) * t231 + qJD(3) * t274 + Ifges(4,4) * t291 + t354;
t144 = mrSges(3,1) * t237 - mrSges(3,2) * t238 + Ifges(3,3) * qJDD(2) + pkin(2) * t321 + pkin(8) * t334 + t312 * t149 + t315 * t154;
t146 = mrSges(3,2) * t267 - mrSges(3,3) * t237 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t318 - pkin(8) * t164 + t149 * t315 - t154 * t312;
t326 = mrSges(2,1) * t294 - mrSges(2,2) * t295 + pkin(1) * t153 + t310 * t144 + t146 * t345 + t308 * t356;
t155 = m(2) * t295 + t157;
t152 = t310 * t163 + (t161 * t313 + t347) * t308;
t150 = m(2) * t294 + t153;
t142 = mrSges(2,2) * t305 - mrSges(2,3) * t294 + t316 * t146 - t313 * t148 + (-t152 * t308 - t153 * t310) * pkin(7);
t141 = -mrSges(2,1) * t305 + mrSges(2,3) * t295 - pkin(1) * t152 - t308 * t144 + (t146 * t313 + t356) * t310;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t309 * t142 - t307 * t141 - qJ(1) * (t150 * t309 + t155 * t307), t142, t146, t149, t162, -t239 * t339 - t241 * t279 + t324, t181; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t307 * t142 + t309 * t141 + qJ(1) * (-t150 * t307 + t155 * t309), t141, t148, t154, t158, -t280 * t239 - t322, t180; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t326, t326, t144, t353, -Ifges(5,3) * t292 - t354, Ifges(6,5) * t265 - Ifges(6,6) * t292 + Ifges(6,3) * t264 + t280 * t241 + t243 * t339 - t323, t328;];
m_new  = t1;
