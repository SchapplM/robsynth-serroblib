% Calculate vector of cutting torques with Newton-Euler for
% S6PRRRRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-05-05 09:33
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRRRRP1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP1_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP1_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRP1_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP1_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP1_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP1_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP1_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 09:27:52
% EndTime: 2019-05-05 09:28:19
% DurationCPUTime: 15.00s
% Computational Cost: add. (281230->336), mult. (549722->421), div. (0->0), fcn. (392124->12), ass. (0->136)
t303 = sin(pkin(11));
t305 = cos(pkin(11));
t291 = t303 * g(1) - t305 * g(2);
t292 = -t305 * g(1) - t303 * g(2);
t302 = -g(3) + qJDD(1);
t314 = cos(qJ(2));
t306 = cos(pkin(6));
t310 = sin(qJ(2));
t342 = t306 * t310;
t304 = sin(pkin(6));
t343 = t304 * t310;
t263 = t291 * t342 + t314 * t292 + t302 * t343;
t315 = qJD(2) ^ 2;
t258 = -t315 * pkin(2) + qJDD(2) * pkin(8) + t263;
t274 = -t304 * t291 + t306 * t302;
t309 = sin(qJ(3));
t313 = cos(qJ(3));
t234 = -t309 * t258 + t313 * t274;
t337 = qJD(2) * qJD(3);
t334 = t313 * t337;
t288 = t309 * qJDD(2) + t334;
t208 = (-t288 + t334) * pkin(9) + (t309 * t313 * t315 + qJDD(3)) * pkin(3) + t234;
t235 = t313 * t258 + t309 * t274;
t289 = t313 * qJDD(2) - t309 * t337;
t339 = qJD(2) * t309;
t296 = qJD(3) * pkin(3) - pkin(9) * t339;
t301 = t313 ^ 2;
t209 = -t301 * t315 * pkin(3) + t289 * pkin(9) - qJD(3) * t296 + t235;
t308 = sin(qJ(4));
t312 = cos(qJ(4));
t204 = t308 * t208 + t312 * t209;
t281 = (t308 * t313 + t309 * t312) * qJD(2);
t250 = -t281 * qJD(4) - t308 * t288 + t312 * t289;
t280 = (t308 * t309 - t312 * t313) * qJD(2);
t265 = t280 * mrSges(5,1) + t281 * mrSges(5,2);
t300 = qJD(3) + qJD(4);
t273 = t300 * mrSges(5,1) - t281 * mrSges(5,3);
t299 = qJDD(3) + qJDD(4);
t266 = t280 * pkin(4) - t281 * pkin(10);
t298 = t300 ^ 2;
t198 = -t298 * pkin(4) + t299 * pkin(10) - t280 * t266 + t204;
t262 = -t310 * t292 + (t291 * t306 + t302 * t304) * t314;
t322 = -qJDD(2) * pkin(2) - t262;
t225 = -t289 * pkin(3) + t296 * t339 + (-pkin(9) * t301 - pkin(8)) * t315 + t322;
t251 = -t280 * qJD(4) + t312 * t288 + t308 * t289;
t201 = (t280 * t300 - t251) * pkin(10) + (t281 * t300 - t250) * pkin(4) + t225;
t307 = sin(qJ(5));
t311 = cos(qJ(5));
t192 = -t307 * t198 + t311 * t201;
t268 = -t307 * t281 + t311 * t300;
t223 = t268 * qJD(5) + t311 * t251 + t307 * t299;
t269 = t311 * t281 + t307 * t300;
t237 = -t268 * mrSges(7,1) + t269 * mrSges(7,2);
t238 = -t268 * mrSges(6,1) + t269 * mrSges(6,2);
t248 = qJDD(5) - t250;
t275 = qJD(5) + t280;
t253 = -t275 * mrSges(6,2) + t268 * mrSges(6,3);
t188 = -0.2e1 * qJD(6) * t269 + (t268 * t275 - t223) * qJ(6) + (t268 * t269 + t248) * pkin(5) + t192;
t252 = -t275 * mrSges(7,2) + t268 * mrSges(7,3);
t336 = m(7) * t188 + t248 * mrSges(7,1) + t275 * t252;
t177 = m(6) * t192 + t248 * mrSges(6,1) + t275 * t253 + (-t237 - t238) * t269 + (-mrSges(6,3) - mrSges(7,3)) * t223 + t336;
t193 = t311 * t198 + t307 * t201;
t222 = -t269 * qJD(5) - t307 * t251 + t311 * t299;
t254 = t275 * pkin(5) - t269 * qJ(6);
t267 = t268 ^ 2;
t191 = -t267 * pkin(5) + t222 * qJ(6) + 0.2e1 * qJD(6) * t268 - t275 * t254 + t193;
t335 = m(7) * t191 + t222 * mrSges(7,3) + t268 * t237;
t255 = t275 * mrSges(7,1) - t269 * mrSges(7,3);
t340 = -t275 * mrSges(6,1) + t269 * mrSges(6,3) - t255;
t346 = -mrSges(6,2) - mrSges(7,2);
t180 = m(6) * t193 + t222 * mrSges(6,3) + t268 * t238 + t248 * t346 + t275 * t340 + t335;
t331 = -t307 * t177 + t311 * t180;
t170 = m(5) * t204 - t299 * mrSges(5,2) + t250 * mrSges(5,3) - t280 * t265 - t300 * t273 + t331;
t203 = t312 * t208 - t308 * t209;
t272 = -t300 * mrSges(5,2) - t280 * mrSges(5,3);
t197 = -t299 * pkin(4) - t298 * pkin(10) + t281 * t266 - t203;
t195 = -t222 * pkin(5) - t267 * qJ(6) + t269 * t254 + qJDD(6) + t197;
t330 = -m(7) * t195 + t222 * mrSges(7,1) + t268 * t252;
t319 = -m(6) * t197 + t222 * mrSges(6,1) + t223 * t346 + t268 * t253 + t269 * t340 + t330;
t182 = m(5) * t203 + t299 * mrSges(5,1) - t251 * mrSges(5,3) - t281 * t265 + t300 * t272 + t319;
t163 = t308 * t170 + t312 * t182;
t287 = (-mrSges(4,1) * t313 + mrSges(4,2) * t309) * qJD(2);
t338 = qJD(2) * t313;
t294 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t338;
t161 = m(4) * t234 + qJDD(3) * mrSges(4,1) - t288 * mrSges(4,3) + qJD(3) * t294 - t287 * t339 + t163;
t293 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t339;
t332 = t312 * t170 - t308 * t182;
t162 = m(4) * t235 - qJDD(3) * mrSges(4,2) + t289 * mrSges(4,3) - qJD(3) * t293 + t287 * t338 + t332;
t156 = t313 * t161 + t309 * t162;
t278 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t309 + Ifges(4,2) * t313) * qJD(2);
t279 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t309 + Ifges(4,4) * t313) * qJD(2);
t226 = Ifges(7,5) * t269 + Ifges(7,6) * t268 + Ifges(7,3) * t275;
t227 = Ifges(6,5) * t269 + Ifges(6,6) * t268 + Ifges(6,3) * t275;
t231 = Ifges(6,1) * t269 + Ifges(6,4) * t268 + Ifges(6,5) * t275;
t230 = Ifges(7,1) * t269 + Ifges(7,4) * t268 + Ifges(7,5) * t275;
t326 = -mrSges(7,1) * t195 + mrSges(7,3) * t191 + Ifges(7,4) * t223 + Ifges(7,2) * t222 + Ifges(7,6) * t248 + t275 * t230;
t165 = Ifges(6,4) * t223 + Ifges(6,2) * t222 + Ifges(6,6) * t248 + t275 * t231 - mrSges(6,1) * t197 + mrSges(6,3) * t193 - pkin(5) * (t223 * mrSges(7,2) - t330) + qJ(6) * (-t248 * mrSges(7,2) - t275 * t255 + t335) + (-pkin(5) * t255 - t226 - t227) * t269 + t326;
t185 = -t223 * mrSges(7,3) - t269 * t237 + t336;
t228 = Ifges(7,4) * t269 + Ifges(7,2) * t268 + Ifges(7,6) * t275;
t229 = Ifges(6,4) * t269 + Ifges(6,2) * t268 + Ifges(6,6) * t275;
t324 = mrSges(7,2) * t195 - mrSges(7,3) * t188 + Ifges(7,1) * t223 + Ifges(7,4) * t222 + Ifges(7,5) * t248 + t268 * t226;
t172 = mrSges(6,2) * t197 - mrSges(6,3) * t192 + Ifges(6,1) * t223 + Ifges(6,4) * t222 + Ifges(6,5) * t248 - qJ(6) * t185 + t268 * t227 + (-t228 - t229) * t275 + t324;
t260 = Ifges(5,4) * t281 - Ifges(5,2) * t280 + Ifges(5,6) * t300;
t261 = Ifges(5,1) * t281 - Ifges(5,4) * t280 + Ifges(5,5) * t300;
t320 = -mrSges(5,1) * t203 + mrSges(5,2) * t204 - Ifges(5,5) * t251 - Ifges(5,6) * t250 - Ifges(5,3) * t299 - pkin(4) * t319 - pkin(10) * t331 - t311 * t165 - t307 * t172 - t281 * t260 - t280 * t261;
t348 = mrSges(4,1) * t234 - mrSges(4,2) * t235 + Ifges(4,5) * t288 + Ifges(4,6) * t289 + Ifges(4,3) * qJDD(3) + pkin(3) * t163 + (t309 * t278 - t313 * t279) * qJD(2) - t320;
t141 = -mrSges(3,1) * t274 + mrSges(3,3) * t263 + t315 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t156 - t348;
t333 = -t309 * t161 + t313 * t162;
t154 = m(3) * t263 - t315 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t333;
t257 = -t315 * pkin(8) + t322;
t174 = t311 * t177 + t307 * t180;
t321 = m(5) * t225 - t250 * mrSges(5,1) + t251 * mrSges(5,2) + t280 * t272 + t281 * t273 + t174;
t318 = -m(4) * t257 + t289 * mrSges(4,1) - t288 * mrSges(4,2) - t293 * t339 + t294 * t338 - t321;
t167 = m(3) * t262 + qJDD(2) * mrSges(3,1) - t315 * mrSges(3,2) + t318;
t150 = t314 * t154 - t310 * t167;
t350 = pkin(7) * t150 + t141 * t314;
t325 = -mrSges(7,1) * t188 + mrSges(7,2) * t191 - Ifges(7,5) * t223 - Ifges(7,6) * t222 - Ifges(7,3) * t248 - t269 * t228;
t349 = mrSges(6,1) * t192 - mrSges(6,2) * t193 + Ifges(6,5) * t223 + Ifges(6,6) * t222 + Ifges(6,3) * t248 + pkin(5) * t185 + t269 * t229 - (t231 + t230) * t268 - t325;
t344 = t167 * t314;
t155 = m(3) * t274 + t156;
t146 = t154 * t342 - t304 * t155 + t306 * t344;
t259 = Ifges(5,5) * t281 - Ifges(5,6) * t280 + Ifges(5,3) * t300;
t151 = mrSges(5,2) * t225 - mrSges(5,3) * t203 + Ifges(5,1) * t251 + Ifges(5,4) * t250 + Ifges(5,5) * t299 - pkin(10) * t174 - t307 * t165 + t311 * t172 - t280 * t259 - t300 * t260;
t157 = -mrSges(5,1) * t225 + mrSges(5,3) * t204 + Ifges(5,4) * t251 + Ifges(5,2) * t250 + Ifges(5,6) * t299 - pkin(4) * t174 - t281 * t259 + t300 * t261 - t349;
t277 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t309 + Ifges(4,6) * t313) * qJD(2);
t142 = -mrSges(4,1) * t257 + mrSges(4,3) * t235 + Ifges(4,4) * t288 + Ifges(4,2) * t289 + Ifges(4,6) * qJDD(3) - pkin(3) * t321 + pkin(9) * t332 + qJD(3) * t279 + t308 * t151 + t312 * t157 - t277 * t339;
t147 = mrSges(4,2) * t257 - mrSges(4,3) * t234 + Ifges(4,1) * t288 + Ifges(4,4) * t289 + Ifges(4,5) * qJDD(3) - pkin(9) * t163 - qJD(3) * t278 + t312 * t151 - t308 * t157 + t277 * t338;
t137 = mrSges(3,1) * t262 - mrSges(3,2) * t263 + Ifges(3,3) * qJDD(2) + pkin(2) * t318 + pkin(8) * t333 + t313 * t142 + t309 * t147;
t139 = mrSges(3,2) * t274 - mrSges(3,3) * t262 + Ifges(3,5) * qJDD(2) - t315 * Ifges(3,6) - pkin(8) * t156 - t309 * t142 + t313 * t147;
t323 = mrSges(2,1) * t291 - mrSges(2,2) * t292 + pkin(1) * t146 + t306 * t137 + t139 * t343 + t304 * t350;
t148 = m(2) * t292 + t150;
t145 = t306 * t155 + (t154 * t310 + t344) * t304;
t143 = m(2) * t291 + t146;
t135 = mrSges(2,2) * t302 - mrSges(2,3) * t291 + t314 * t139 - t310 * t141 + (-t145 * t304 - t146 * t306) * pkin(7);
t134 = -mrSges(2,1) * t302 + mrSges(2,3) * t292 - pkin(1) * t145 - t304 * t137 + (t139 * t310 + t350) * t306;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t305 * t135 - t303 * t134 - qJ(1) * (t305 * t143 + t303 * t148), t135, t139, t147, t151, t172, -t275 * t228 + t324; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t303 * t135 + t305 * t134 + qJ(1) * (-t303 * t143 + t305 * t148), t134, t141, t142, t157, t165, -t269 * t226 + t326; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t323, t323, t137, t348, -t320, t349, -t268 * t230 - t325;];
m_new  = t1;
