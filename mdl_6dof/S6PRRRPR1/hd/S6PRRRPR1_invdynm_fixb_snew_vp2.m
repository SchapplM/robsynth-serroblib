% Calculate vector of cutting torques with Newton-Euler for
% S6PRRRPR1
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-05-05 07:03
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRRRPR1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR1_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR1_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR1_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR1_invdynm_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR1_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR1_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR1_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 06:58:10
% EndTime: 2019-05-05 06:58:54
% DurationCPUTime: 33.61s
% Computational Cost: add. (621374->341), mult. (1327579->441), div. (0->0), fcn. (987207->14), ass. (0->146)
t305 = sin(pkin(11));
t308 = cos(pkin(11));
t291 = t305 * g(1) - t308 * g(2);
t292 = -t308 * g(1) - t305 * g(2);
t303 = -g(3) + qJDD(1);
t317 = cos(qJ(2));
t309 = cos(pkin(6));
t313 = sin(qJ(2));
t340 = t309 * t313;
t306 = sin(pkin(6));
t341 = t306 * t313;
t263 = t291 * t340 + t317 * t292 + t303 * t341;
t318 = qJD(2) ^ 2;
t256 = -t318 * pkin(2) + qJDD(2) * pkin(8) + t263;
t274 = -t306 * t291 + t309 * t303;
t312 = sin(qJ(3));
t316 = cos(qJ(3));
t242 = -t312 * t256 + t316 * t274;
t337 = qJD(2) * qJD(3);
t336 = t316 * t337;
t288 = t312 * qJDD(2) + t336;
t227 = (-t288 + t336) * pkin(9) + (t312 * t316 * t318 + qJDD(3)) * pkin(3) + t242;
t243 = t316 * t256 + t312 * t274;
t289 = t316 * qJDD(2) - t312 * t337;
t339 = qJD(2) * t312;
t296 = qJD(3) * pkin(3) - pkin(9) * t339;
t302 = t316 ^ 2;
t230 = -t302 * t318 * pkin(3) + t289 * pkin(9) - qJD(3) * t296 + t243;
t311 = sin(qJ(4));
t315 = cos(qJ(4));
t207 = t315 * t227 - t311 * t230;
t280 = (-t311 * t312 + t315 * t316) * qJD(2);
t252 = t280 * qJD(4) + t315 * t288 + t311 * t289;
t281 = (t311 * t316 + t312 * t315) * qJD(2);
t300 = qJDD(3) + qJDD(4);
t301 = qJD(3) + qJD(4);
t203 = (t280 * t301 - t252) * qJ(5) + (t280 * t281 + t300) * pkin(4) + t207;
t208 = t311 * t227 + t315 * t230;
t251 = -t281 * qJD(4) - t311 * t288 + t315 * t289;
t272 = t301 * pkin(4) - t281 * qJ(5);
t276 = t280 ^ 2;
t205 = -t276 * pkin(4) + t251 * qJ(5) - t301 * t272 + t208;
t304 = sin(pkin(12));
t307 = cos(pkin(12));
t266 = t307 * t280 - t304 * t281;
t345 = 2 * qJD(5);
t200 = t304 * t203 + t307 * t205 + t266 * t345;
t228 = t307 * t251 - t304 * t252;
t267 = t304 * t280 + t307 * t281;
t239 = -t266 * mrSges(6,1) + t267 * mrSges(6,2);
t254 = t301 * mrSges(6,1) - t267 * mrSges(6,3);
t240 = -t266 * pkin(5) - t267 * pkin(10);
t299 = t301 ^ 2;
t197 = -t299 * pkin(5) + t300 * pkin(10) + t266 * t240 + t200;
t262 = -t313 * t292 + (t291 * t309 + t303 * t306) * t317;
t325 = -qJDD(2) * pkin(2) - t262;
t241 = -t289 * pkin(3) + t296 * t339 + (-pkin(9) * t302 - pkin(8)) * t318 + t325;
t210 = -t251 * pkin(4) - t276 * qJ(5) + t281 * t272 + qJDD(5) + t241;
t229 = t304 * t251 + t307 * t252;
t201 = (-t266 * t301 - t229) * pkin(10) + (t267 * t301 - t228) * pkin(5) + t210;
t310 = sin(qJ(6));
t314 = cos(qJ(6));
t194 = -t310 * t197 + t314 * t201;
t248 = -t310 * t267 + t314 * t301;
t213 = t248 * qJD(6) + t314 * t229 + t310 * t300;
t226 = qJDD(6) - t228;
t249 = t314 * t267 + t310 * t301;
t231 = -t248 * mrSges(7,1) + t249 * mrSges(7,2);
t258 = qJD(6) - t266;
t232 = -t258 * mrSges(7,2) + t248 * mrSges(7,3);
t190 = m(7) * t194 + t226 * mrSges(7,1) - t213 * mrSges(7,3) - t249 * t231 + t258 * t232;
t195 = t314 * t197 + t310 * t201;
t212 = -t249 * qJD(6) - t310 * t229 + t314 * t300;
t233 = t258 * mrSges(7,1) - t249 * mrSges(7,3);
t191 = m(7) * t195 - t226 * mrSges(7,2) + t212 * mrSges(7,3) + t248 * t231 - t258 * t233;
t332 = -t310 * t190 + t314 * t191;
t177 = m(6) * t200 - t300 * mrSges(6,2) + t228 * mrSges(6,3) + t266 * t239 - t301 * t254 + t332;
t331 = -t307 * t203 + t304 * t205;
t199 = -0.2e1 * qJD(5) * t267 - t331;
t253 = -t301 * mrSges(6,2) + t266 * mrSges(6,3);
t196 = -t300 * pkin(5) - t299 * pkin(10) + (t345 + t240) * t267 + t331;
t326 = -m(7) * t196 + t212 * mrSges(7,1) - t213 * mrSges(7,2) + t248 * t232 - t249 * t233;
t186 = m(6) * t199 + t300 * mrSges(6,1) - t229 * mrSges(6,3) - t267 * t239 + t301 * t253 + t326;
t172 = t304 * t177 + t307 * t186;
t268 = -t280 * mrSges(5,1) + t281 * mrSges(5,2);
t271 = -t301 * mrSges(5,2) + t280 * mrSges(5,3);
t169 = m(5) * t207 + t300 * mrSges(5,1) - t252 * mrSges(5,3) - t281 * t268 + t301 * t271 + t172;
t273 = t301 * mrSges(5,1) - t281 * mrSges(5,3);
t333 = t307 * t177 - t304 * t186;
t170 = m(5) * t208 - t300 * mrSges(5,2) + t251 * mrSges(5,3) + t280 * t268 - t301 * t273 + t333;
t163 = t315 * t169 + t311 * t170;
t287 = (-mrSges(4,1) * t316 + mrSges(4,2) * t312) * qJD(2);
t338 = qJD(2) * t316;
t294 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t338;
t161 = m(4) * t242 + qJDD(3) * mrSges(4,1) - t288 * mrSges(4,3) + qJD(3) * t294 - t287 * t339 + t163;
t293 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t339;
t334 = -t311 * t169 + t315 * t170;
t162 = m(4) * t243 - qJDD(3) * mrSges(4,2) + t289 * mrSges(4,3) - qJD(3) * t293 + t287 * t338 + t334;
t156 = t316 * t161 + t312 * t162;
t278 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t312 + Ifges(4,2) * t316) * qJD(2);
t279 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t312 + Ifges(4,4) * t316) * qJD(2);
t260 = Ifges(5,4) * t281 + Ifges(5,2) * t280 + Ifges(5,6) * t301;
t261 = Ifges(5,1) * t281 + Ifges(5,4) * t280 + Ifges(5,5) * t301;
t214 = Ifges(7,5) * t249 + Ifges(7,6) * t248 + Ifges(7,3) * t258;
t216 = Ifges(7,1) * t249 + Ifges(7,4) * t248 + Ifges(7,5) * t258;
t183 = -mrSges(7,1) * t196 + mrSges(7,3) * t195 + Ifges(7,4) * t213 + Ifges(7,2) * t212 + Ifges(7,6) * t226 - t249 * t214 + t258 * t216;
t215 = Ifges(7,4) * t249 + Ifges(7,2) * t248 + Ifges(7,6) * t258;
t184 = mrSges(7,2) * t196 - mrSges(7,3) * t194 + Ifges(7,1) * t213 + Ifges(7,4) * t212 + Ifges(7,5) * t226 + t248 * t214 - t258 * t215;
t235 = Ifges(6,4) * t267 + Ifges(6,2) * t266 + Ifges(6,6) * t301;
t236 = Ifges(6,1) * t267 + Ifges(6,4) * t266 + Ifges(6,5) * t301;
t324 = -mrSges(6,1) * t199 + mrSges(6,2) * t200 - Ifges(6,5) * t229 - Ifges(6,6) * t228 - Ifges(6,3) * t300 - pkin(5) * t326 - pkin(10) * t332 - t314 * t183 - t310 * t184 - t267 * t235 + t266 * t236;
t321 = -mrSges(5,1) * t207 + mrSges(5,2) * t208 - Ifges(5,5) * t252 - Ifges(5,6) * t251 - Ifges(5,3) * t300 - pkin(4) * t172 - t281 * t260 + t280 * t261 + t324;
t346 = mrSges(4,1) * t242 - mrSges(4,2) * t243 + Ifges(4,5) * t288 + Ifges(4,6) * t289 + Ifges(4,3) * qJDD(3) + pkin(3) * t163 + (t312 * t278 - t316 * t279) * qJD(2) - t321;
t143 = -mrSges(3,1) * t274 + mrSges(3,3) * t263 + t318 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t156 - t346;
t335 = -t312 * t161 + t316 * t162;
t154 = m(3) * t263 - t318 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t335;
t255 = -t318 * pkin(8) + t325;
t179 = t314 * t190 + t310 * t191;
t328 = m(6) * t210 - t228 * mrSges(6,1) + t229 * mrSges(6,2) - t266 * t253 + t267 * t254 + t179;
t323 = m(5) * t241 - t251 * mrSges(5,1) + t252 * mrSges(5,2) - t280 * t271 + t281 * t273 + t328;
t320 = -m(4) * t255 + t289 * mrSges(4,1) - t288 * mrSges(4,2) - t293 * t339 + t294 * t338 - t323;
t174 = m(3) * t262 + qJDD(2) * mrSges(3,1) - t318 * mrSges(3,2) + t320;
t150 = t317 * t154 - t313 * t174;
t347 = pkin(7) * t150 + t143 * t317;
t342 = t174 * t317;
t155 = m(3) * t274 + t156;
t147 = t154 * t340 - t306 * t155 + t309 * t342;
t234 = Ifges(6,5) * t267 + Ifges(6,6) * t266 + Ifges(6,3) * t301;
t164 = mrSges(6,2) * t210 - mrSges(6,3) * t199 + Ifges(6,1) * t229 + Ifges(6,4) * t228 + Ifges(6,5) * t300 - pkin(10) * t179 - t310 * t183 + t314 * t184 + t266 * t234 - t301 * t235;
t322 = mrSges(7,1) * t194 - mrSges(7,2) * t195 + Ifges(7,5) * t213 + Ifges(7,6) * t212 + Ifges(7,3) * t226 + t249 * t215 - t248 * t216;
t165 = -mrSges(6,1) * t210 + mrSges(6,3) * t200 + Ifges(6,4) * t229 + Ifges(6,2) * t228 + Ifges(6,6) * t300 - pkin(5) * t179 - t267 * t234 + t301 * t236 - t322;
t259 = Ifges(5,5) * t281 + Ifges(5,6) * t280 + Ifges(5,3) * t301;
t151 = -mrSges(5,1) * t241 + mrSges(5,3) * t208 + Ifges(5,4) * t252 + Ifges(5,2) * t251 + Ifges(5,6) * t300 - pkin(4) * t328 + qJ(5) * t333 + t304 * t164 + t307 * t165 - t281 * t259 + t301 * t261;
t157 = mrSges(5,2) * t241 - mrSges(5,3) * t207 + Ifges(5,1) * t252 + Ifges(5,4) * t251 + Ifges(5,5) * t300 - qJ(5) * t172 + t307 * t164 - t304 * t165 + t280 * t259 - t301 * t260;
t277 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t312 + Ifges(4,6) * t316) * qJD(2);
t140 = -mrSges(4,1) * t255 + mrSges(4,3) * t243 + Ifges(4,4) * t288 + Ifges(4,2) * t289 + Ifges(4,6) * qJDD(3) - pkin(3) * t323 + pkin(9) * t334 + qJD(3) * t279 + t315 * t151 + t311 * t157 - t277 * t339;
t141 = mrSges(4,2) * t255 - mrSges(4,3) * t242 + Ifges(4,1) * t288 + Ifges(4,4) * t289 + Ifges(4,5) * qJDD(3) - pkin(9) * t163 - qJD(3) * t278 - t311 * t151 + t315 * t157 + t277 * t338;
t137 = mrSges(3,1) * t262 - mrSges(3,2) * t263 + Ifges(3,3) * qJDD(2) + pkin(2) * t320 + pkin(8) * t335 + t316 * t140 + t312 * t141;
t139 = mrSges(3,2) * t274 - mrSges(3,3) * t262 + Ifges(3,5) * qJDD(2) - t318 * Ifges(3,6) - pkin(8) * t156 - t312 * t140 + t316 * t141;
t327 = mrSges(2,1) * t291 - mrSges(2,2) * t292 + pkin(1) * t147 + t309 * t137 + t139 * t341 + t347 * t306;
t148 = m(2) * t292 + t150;
t146 = t309 * t155 + (t154 * t313 + t342) * t306;
t144 = m(2) * t291 + t147;
t135 = mrSges(2,2) * t303 - mrSges(2,3) * t291 + t317 * t139 - t313 * t143 + (-t146 * t306 - t147 * t309) * pkin(7);
t134 = -mrSges(2,1) * t303 + mrSges(2,3) * t292 - pkin(1) * t146 - t306 * t137 + (t139 * t313 + t347) * t309;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t308 * t135 - t305 * t134 - qJ(1) * (t308 * t144 + t305 * t148), t135, t139, t141, t157, t164, t184; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t305 * t135 + t308 * t134 + qJ(1) * (-t305 * t144 + t308 * t148), t134, t143, t140, t151, t165, t183; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t327, t327, t137, t346, -t321, -t324, t322;];
m_new  = t1;
