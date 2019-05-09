% Calculate vector of cutting torques with Newton-Euler for
% S6PRRRPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
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
% Datum: 2019-05-05 07:33
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRRRPR3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR3_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR3_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR3_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR3_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR3_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR3_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR3_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 07:27:39
% EndTime: 2019-05-05 07:28:02
% DurationCPUTime: 12.22s
% Computational Cost: add. (215765->342), mult. (431551->419), div. (0->0), fcn. (299044->12), ass. (0->139)
t306 = sin(pkin(11));
t308 = cos(pkin(11));
t291 = g(1) * t306 - g(2) * t308;
t292 = -g(1) * t308 - g(2) * t306;
t305 = -g(3) + qJDD(1);
t316 = cos(qJ(2));
t309 = cos(pkin(6));
t313 = sin(qJ(2));
t343 = t309 * t313;
t307 = sin(pkin(6));
t344 = t307 * t313;
t253 = t291 * t343 + t292 * t316 + t305 * t344;
t317 = qJD(2) ^ 2;
t245 = -pkin(2) * t317 + qJDD(2) * pkin(8) + t253;
t270 = -t291 * t307 + t305 * t309;
t312 = sin(qJ(3));
t315 = cos(qJ(3));
t222 = -t312 * t245 + t270 * t315;
t338 = qJD(2) * qJD(3);
t337 = t315 * t338;
t288 = qJDD(2) * t312 + t337;
t209 = (-t288 + t337) * pkin(9) + (t312 * t315 * t317 + qJDD(3)) * pkin(3) + t222;
t223 = t245 * t315 + t270 * t312;
t289 = qJDD(2) * t315 - t312 * t338;
t340 = qJD(2) * t312;
t296 = qJD(3) * pkin(3) - pkin(9) * t340;
t304 = t315 ^ 2;
t210 = -pkin(3) * t304 * t317 + pkin(9) * t289 - qJD(3) * t296 + t223;
t311 = sin(qJ(4));
t350 = cos(qJ(4));
t204 = t209 * t350 - t210 * t311;
t205 = t209 * t311 + t210 * t350;
t280 = (t311 * t315 + t312 * t350) * qJD(2);
t240 = qJD(4) * t280 + t288 * t311 - t289 * t350;
t339 = qJD(2) * t315;
t279 = t311 * t340 - t339 * t350;
t241 = -qJD(4) * t279 + t288 * t350 + t289 * t311;
t303 = qJD(3) + qJD(4);
t249 = Ifges(5,4) * t280 - Ifges(5,2) * t279 + Ifges(5,6) * t303;
t259 = -mrSges(6,2) * t279 - mrSges(6,3) * t280;
t267 = mrSges(6,1) * t279 - mrSges(6,3) * t303;
t302 = qJDD(3) + qJDD(4);
t257 = pkin(4) * t279 - qJ(5) * t280;
t301 = t303 ^ 2;
t200 = -t302 * pkin(4) - t301 * qJ(5) + t257 * t280 + qJDD(5) - t204;
t345 = t279 * t303;
t194 = (t279 * t280 - t302) * pkin(10) + (t241 + t345) * pkin(5) + t200;
t269 = pkin(5) * t280 - pkin(10) * t303;
t275 = t279 ^ 2;
t252 = -t313 * t292 + (t291 * t309 + t305 * t307) * t316;
t327 = -qJDD(2) * pkin(2) - t252;
t217 = -t289 * pkin(3) + t296 * t340 + (-pkin(9) * t304 - pkin(8)) * t317 + t327;
t351 = -2 * qJD(5);
t319 = (-t241 + t345) * qJ(5) + t217 + (pkin(4) * t303 + t351) * t280;
t197 = -t275 * pkin(5) - t280 * t269 + (pkin(4) + pkin(10)) * t240 + t319;
t310 = sin(qJ(6));
t314 = cos(qJ(6));
t191 = t194 * t314 - t197 * t310;
t261 = t279 * t314 - t303 * t310;
t216 = qJD(6) * t261 + t240 * t310 + t302 * t314;
t262 = t279 * t310 + t303 * t314;
t224 = -mrSges(7,1) * t261 + mrSges(7,2) * t262;
t238 = qJDD(6) + t241;
t273 = qJD(6) + t280;
t242 = -mrSges(7,2) * t273 + mrSges(7,3) * t261;
t188 = m(7) * t191 + mrSges(7,1) * t238 - mrSges(7,3) * t216 - t224 * t262 + t242 * t273;
t192 = t194 * t310 + t197 * t314;
t215 = -qJD(6) * t262 + t240 * t314 - t302 * t310;
t243 = mrSges(7,1) * t273 - mrSges(7,3) * t262;
t189 = m(7) * t192 - mrSges(7,2) * t238 + mrSges(7,3) * t215 + t224 * t261 - t243 * t273;
t176 = t314 * t188 + t310 * t189;
t329 = -t301 * pkin(4) + t302 * qJ(5) - t279 * t257 + t205;
t196 = -t240 * pkin(5) - t275 * pkin(10) + ((2 * qJD(5)) + t269) * t303 + t329;
t218 = Ifges(7,5) * t262 + Ifges(7,6) * t261 + Ifges(7,3) * t273;
t220 = Ifges(7,1) * t262 + Ifges(7,4) * t261 + Ifges(7,5) * t273;
t179 = -mrSges(7,1) * t196 + mrSges(7,3) * t192 + Ifges(7,4) * t216 + Ifges(7,2) * t215 + Ifges(7,6) * t238 - t218 * t262 + t220 * t273;
t219 = Ifges(7,4) * t262 + Ifges(7,2) * t261 + Ifges(7,6) * t273;
t180 = mrSges(7,2) * t196 - mrSges(7,3) * t191 + Ifges(7,1) * t216 + Ifges(7,4) * t215 + Ifges(7,5) * t238 + t218 * t261 - t219 * t273;
t198 = t303 * t351 - t329;
t246 = Ifges(6,5) * t303 - Ifges(6,6) * t280 + Ifges(6,3) * t279;
t325 = -mrSges(6,2) * t200 + mrSges(6,3) * t198 - Ifges(6,1) * t302 + Ifges(6,4) * t241 - Ifges(6,5) * t240 + pkin(10) * t176 + t310 * t179 - t180 * t314 + t246 * t280;
t193 = -m(7) * t196 + t215 * mrSges(7,1) - mrSges(7,2) * t216 + t261 * t242 - t243 * t262;
t268 = mrSges(6,1) * t280 + mrSges(6,2) * t303;
t326 = -m(6) * t198 + mrSges(6,3) * t302 + t268 * t303 - t193;
t330 = -m(6) * t200 - mrSges(6,1) * t241 - t259 * t280 - t176;
t248 = Ifges(6,4) * t303 - Ifges(6,2) * t280 + Ifges(6,6) * t279;
t341 = Ifges(5,1) * t280 - Ifges(5,4) * t279 + Ifges(5,5) * t303 - t248;
t355 = -mrSges(5,2) * t205 + pkin(4) * (-t302 * mrSges(6,2) - t303 * t267 + t330) + qJ(5) * (-t240 * mrSges(6,1) - t279 * t259 + t326) + mrSges(5,1) * t204 + t280 * t249 - Ifges(5,6) * t240 + Ifges(5,5) * t241 + Ifges(5,3) * t302 - t325 + t341 * t279;
t258 = mrSges(5,1) * t279 + mrSges(5,2) * t280;
t265 = -mrSges(5,2) * t303 - mrSges(5,3) * t279;
t172 = m(5) * t204 - t241 * mrSges(5,3) - t280 * t258 + (t265 - t267) * t303 + (mrSges(5,1) - mrSges(6,2)) * t302 + t330;
t266 = mrSges(5,1) * t303 - mrSges(5,3) * t280;
t183 = m(5) * t205 - t302 * mrSges(5,2) - t303 * t266 + (-t258 - t259) * t279 + (-mrSges(5,3) - mrSges(6,1)) * t240 + t326;
t168 = t172 * t350 + t183 * t311;
t277 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t312 + Ifges(4,2) * t315) * qJD(2);
t278 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t312 + Ifges(4,4) * t315) * qJD(2);
t354 = mrSges(4,1) * t222 - mrSges(4,2) * t223 + Ifges(4,5) * t288 + Ifges(4,6) * t289 + Ifges(4,3) * qJDD(3) + pkin(3) * t168 + (t277 * t312 - t278 * t315) * qJD(2) + t355;
t287 = (-mrSges(4,1) * t315 + mrSges(4,2) * t312) * qJD(2);
t294 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t339;
t166 = m(4) * t222 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t288 + qJD(3) * t294 - t287 * t340 + t168;
t293 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t340;
t334 = -t172 * t311 + t183 * t350;
t167 = m(4) * t223 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t289 - qJD(3) * t293 + t287 * t339 + t334;
t161 = t166 * t315 + t167 * t312;
t147 = -mrSges(3,1) * t270 + mrSges(3,3) * t253 + t317 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t161 - t354;
t335 = -t166 * t312 + t167 * t315;
t159 = m(3) * t253 - mrSges(3,1) * t317 - qJDD(2) * mrSges(3,2) + t335;
t244 = -t317 * pkin(8) + t327;
t177 = -t188 * t310 + t189 * t314;
t202 = t240 * pkin(4) + t319;
t174 = m(6) * t202 - mrSges(6,2) * t240 - mrSges(6,3) * t241 - t267 * t279 - t268 * t280 + t177;
t323 = m(5) * t217 + mrSges(5,1) * t240 + mrSges(5,2) * t241 + t265 * t279 + t266 * t280 + t174;
t320 = -m(4) * t244 + mrSges(4,1) * t289 - mrSges(4,2) * t288 - t293 * t340 + t294 * t339 - t323;
t170 = m(3) * t252 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t317 + t320;
t155 = t159 * t316 - t170 * t313;
t353 = pkin(7) * t155 + t147 * t316;
t348 = Ifges(5,4) + Ifges(6,6);
t346 = t170 * t316;
t250 = Ifges(6,1) * t303 - Ifges(6,4) * t280 + Ifges(6,5) * t279;
t342 = -Ifges(5,5) * t280 + Ifges(5,6) * t279 - Ifges(5,3) * t303 - t250;
t160 = m(3) * t270 + t161;
t152 = t159 * t343 - t160 * t307 + t309 * t346;
t324 = -mrSges(6,1) * t198 + mrSges(6,2) * t202 - pkin(5) * t193 - pkin(10) * t177 - t314 * t179 - t310 * t180;
t156 = -mrSges(5,1) * t217 + mrSges(5,3) * t205 - pkin(4) * t174 + t341 * t303 + (Ifges(5,6) - Ifges(6,5)) * t302 + t342 * t280 + t348 * t241 + (-Ifges(5,2) - Ifges(6,3)) * t240 + t324;
t328 = mrSges(7,1) * t191 - mrSges(7,2) * t192 + Ifges(7,5) * t216 + Ifges(7,6) * t215 + Ifges(7,3) * t238 + t219 * t262 - t261 * t220;
t322 = mrSges(6,1) * t200 - mrSges(6,3) * t202 + pkin(5) * t176 + t328;
t162 = (-t249 + t246) * t303 + (Ifges(5,5) - Ifges(6,4)) * t302 + t342 * t279 + (Ifges(5,1) + Ifges(6,2)) * t241 - t348 * t240 + mrSges(5,2) * t217 - mrSges(5,3) * t204 - qJ(5) * t174 + t322;
t276 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t312 + Ifges(4,6) * t315) * qJD(2);
t145 = -mrSges(4,1) * t244 + mrSges(4,3) * t223 + Ifges(4,4) * t288 + Ifges(4,2) * t289 + Ifges(4,6) * qJDD(3) - pkin(3) * t323 + pkin(9) * t334 + qJD(3) * t278 + t156 * t350 + t311 * t162 - t276 * t340;
t148 = mrSges(4,2) * t244 - mrSges(4,3) * t222 + Ifges(4,1) * t288 + Ifges(4,4) * t289 + Ifges(4,5) * qJDD(3) - pkin(9) * t168 - qJD(3) * t277 - t156 * t311 + t162 * t350 + t276 * t339;
t142 = mrSges(3,1) * t252 - mrSges(3,2) * t253 + Ifges(3,3) * qJDD(2) + pkin(2) * t320 + pkin(8) * t335 + t315 * t145 + t312 * t148;
t144 = mrSges(3,2) * t270 - mrSges(3,3) * t252 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t317 - pkin(8) * t161 - t145 * t312 + t148 * t315;
t331 = mrSges(2,1) * t291 - mrSges(2,2) * t292 + pkin(1) * t152 + t309 * t142 + t144 * t344 + t307 * t353;
t153 = m(2) * t292 + t155;
t151 = t309 * t160 + (t159 * t313 + t346) * t307;
t149 = m(2) * t291 + t152;
t140 = mrSges(2,2) * t305 - mrSges(2,3) * t291 + t316 * t144 - t313 * t147 + (-t151 * t307 - t152 * t309) * pkin(7);
t139 = -mrSges(2,1) * t305 + mrSges(2,3) * t292 - pkin(1) * t151 - t307 * t142 + (t144 * t313 + t353) * t309;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t308 * t140 - t306 * t139 - qJ(1) * (t149 * t308 + t153 * t306), t140, t144, t148, t162, -t279 * t248 - t325, t180; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t306 * t140 + t308 * t139 + qJ(1) * (-t149 * t306 + t153 * t308), t139, t147, t145, t156, Ifges(6,4) * t302 - Ifges(6,2) * t241 + Ifges(6,6) * t240 - t303 * t246 + t279 * t250 - t322, t179; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t331, t331, t142, t354, t355, Ifges(6,5) * t302 - Ifges(6,6) * t241 + Ifges(6,3) * t240 + t303 * t248 + t280 * t250 - t324, t328;];
m_new  = t1;
