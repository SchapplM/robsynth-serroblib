% Calculate vector of cutting torques with Newton-Euler for
% S5RRRRR15
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d1,d2,d3,d4,d5]';
% m [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% m [3x6]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 22:28
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRRRR15_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR15_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR15_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR15_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR15_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR15_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR15_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR15_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR15_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 22:24:31
% EndTime: 2024-09-27 22:25:17
% DurationCPUTime: 24.52s
% Computational Cost: add. (928673->338), mult. (2613850->458), div. (0->0), fcn. (2118667->14), ass. (0->149)
t305 = sin(pkin(5));
t344 = pkin(8) * t305;
t304 = sin(pkin(6));
t343 = pkin(11) * t304;
t306 = cos(pkin(6));
t342 = pkin(11) * t306;
t318 = qJD(1) ^ 2;
t341 = t305 ^ 2 * t318;
t308 = sin(qJ(5));
t340 = t304 * t308;
t313 = cos(qJ(5));
t339 = t304 * t313;
t311 = sin(qJ(2));
t338 = t305 * t311;
t316 = cos(qJ(2));
t337 = t305 * t316;
t307 = cos(pkin(5));
t336 = t307 * t311;
t335 = t307 * t316;
t333 = qJD(1) * qJD(2);
t285 = (qJDD(1) * t311 + t316 * t333) * t305;
t296 = t307 * qJDD(1) + qJDD(2);
t297 = t307 * qJD(1) + qJD(2);
t312 = sin(qJ(1));
t317 = cos(qJ(1));
t293 = t312 * g(1) - g(2) * t317;
t281 = qJDD(1) * pkin(1) + t318 * t344 + t293;
t294 = -g(1) * t317 - g(2) * t312;
t282 = -pkin(1) * t318 + qJDD(1) * t344 + t294;
t327 = t281 * t335 - t311 * t282;
t238 = pkin(2) * t296 - t285 * pkin(9) + (pkin(2) * t311 * t341 + (pkin(9) * qJD(1) * t297 - g(3)) * t305) * t316 + t327;
t261 = -g(3) * t338 + t281 * t336 + t316 * t282;
t334 = qJD(1) * t305;
t331 = t311 * t334;
t284 = pkin(2) * t297 - pkin(9) * t331;
t286 = (qJDD(1) * t316 - t311 * t333) * t305;
t332 = t316 ^ 2 * t341;
t239 = -pkin(2) * t332 + pkin(9) * t286 - t284 * t297 + t261;
t310 = sin(qJ(3));
t315 = cos(qJ(3));
t218 = t315 * t238 - t310 * t239;
t277 = (-t310 * t311 + t315 * t316) * t334;
t252 = qJD(3) * t277 + t285 * t315 + t286 * t310;
t278 = (t310 * t316 + t311 * t315) * t334;
t292 = qJDD(3) + t296;
t295 = qJD(3) + t297;
t208 = (t277 * t295 - t252) * pkin(10) + (t277 * t278 + t292) * pkin(3) + t218;
t219 = t310 * t238 + t315 * t239;
t251 = -qJD(3) * t278 - t285 * t310 + t286 * t315;
t267 = pkin(3) * t295 - pkin(10) * t278;
t275 = t277 ^ 2;
t210 = -pkin(3) * t275 + pkin(10) * t251 - t267 * t295 + t219;
t309 = sin(qJ(4));
t314 = cos(qJ(4));
t203 = t208 * t309 + t314 * t210;
t262 = t277 * t314 - t278 * t309;
t263 = t277 * t309 + t278 * t314;
t240 = -pkin(4) * t262 - t263 * t343;
t291 = qJD(4) + t295;
t250 = pkin(4) * t291 - t263 * t342;
t226 = -qJD(4) * t263 + t251 * t314 - t252 * t309;
t290 = qJDD(4) + t292;
t325 = t226 * t306 + t290 * t304;
t199 = pkin(11) * t325 + t262 * t240 - t291 * t250 + t203;
t202 = t208 * t314 - t210 * t309;
t227 = qJD(4) * t262 + t251 * t309 + t252 * t314;
t324 = t262 * t306 + t291 * t304;
t245 = t324 * pkin(11);
t198 = pkin(4) * t290 - t227 * t342 - t240 * t263 + t245 * t291 + t202;
t271 = -g(3) * t307 - t281 * t305;
t243 = -pkin(2) * t286 - pkin(9) * t332 + t284 * t331 + t271;
t216 = -pkin(3) * t251 - pkin(10) * t275 + t278 * t267 + t243;
t200 = -pkin(4) * t226 - t227 * t343 - t245 * t262 + t250 * t263 + t216;
t326 = t198 * t306 + t200 * t304;
t195 = -t308 * t199 + t313 * t326;
t229 = -t308 * t263 + t313 * t324;
t205 = t229 * qJD(5) + t313 * t227 + t308 * t325;
t230 = t313 * t263 + t308 * t324;
t214 = -mrSges(6,1) * t229 + mrSges(6,2) * t230;
t220 = -t226 * t304 + t290 * t306 + qJDD(5);
t246 = -t262 * t304 + t291 * t306 + qJD(5);
t221 = -mrSges(6,2) * t246 + mrSges(6,3) * t229;
t191 = m(6) * t195 + mrSges(6,1) * t220 - mrSges(6,3) * t205 - t214 * t230 + t221 * t246;
t196 = t313 * t199 + t308 * t326;
t204 = -t230 * qJD(5) - t308 * t227 + t313 * t325;
t222 = mrSges(6,1) * t246 - mrSges(6,3) * t230;
t192 = m(6) * t196 - mrSges(6,2) * t220 + mrSges(6,3) * t204 + t214 * t229 - t222 * t246;
t197 = -t198 * t304 + t200 * t306;
t194 = m(6) * t197 - mrSges(6,1) * t204 + mrSges(6,2) * t205 - t221 * t229 + t222 * t230;
t174 = -t194 * t304 + (t191 * t313 + t192 * t308) * t306;
t241 = -mrSges(5,1) * t262 + mrSges(5,2) * t263;
t253 = -mrSges(5,2) * t291 + mrSges(5,3) * t262;
t171 = m(5) * t202 + mrSges(5,1) * t290 - mrSges(5,3) * t227 - t241 * t263 + t253 * t291 + t174;
t179 = -t191 * t308 + t192 * t313;
t254 = mrSges(5,1) * t291 - mrSges(5,3) * t263;
t177 = m(5) * t203 - mrSges(5,2) * t290 + mrSges(5,3) * t226 + t241 * t262 - t254 * t291 + t179;
t168 = t171 * t314 + t177 * t309;
t264 = -mrSges(4,1) * t277 + mrSges(4,2) * t278;
t265 = -mrSges(4,2) * t295 + mrSges(4,3) * t277;
t165 = m(4) * t218 + mrSges(4,1) * t292 - mrSges(4,3) * t252 - t264 * t278 + t265 * t295 + t168;
t266 = mrSges(4,1) * t295 - mrSges(4,3) * t278;
t328 = -t171 * t309 + t177 * t314;
t166 = m(4) * t219 - mrSges(4,2) * t292 + mrSges(4,3) * t251 + t264 * t277 - t266 * t295 + t328;
t159 = t165 * t315 + t166 * t310;
t173 = t191 * t339 + t192 * t340 + t194 * t306;
t330 = t316 * t334;
t260 = -g(3) * t337 + t327;
t280 = -mrSges(3,2) * t297 + mrSges(3,3) * t330;
t283 = (-mrSges(3,1) * t316 + mrSges(3,2) * t311) * t334;
t157 = m(3) * t260 + mrSges(3,1) * t296 - mrSges(3,3) * t285 + t280 * t297 - t283 * t331 + t159;
t279 = mrSges(3,1) * t297 - mrSges(3,3) * t331;
t329 = -t165 * t310 + t166 * t315;
t158 = m(3) * t261 - mrSges(3,2) * t296 + mrSges(3,3) * t286 - t279 * t297 + t283 * t330 + t329;
t153 = -t157 * t311 + t158 * t316;
t322 = m(5) * t216 - mrSges(5,1) * t226 + t227 * mrSges(5,2) - t253 * t262 + t263 * t254 + t173;
t320 = m(4) * t243 - mrSges(4,1) * t251 + t252 * mrSges(4,2) - t265 * t277 + t278 * t266 + t322;
t169 = (t279 * t311 - t280 * t316) * t334 + t320 + m(3) * t271 - mrSges(3,1) * t286 + mrSges(3,2) * t285;
t150 = t157 * t335 + t158 * t336 - t169 * t305;
t212 = Ifges(6,4) * t230 + Ifges(6,2) * t229 + Ifges(6,6) * t246;
t213 = Ifges(6,1) * t230 + Ifges(6,4) * t229 + Ifges(6,5) * t246;
t181 = mrSges(6,1) * t195 - mrSges(6,2) * t196 + Ifges(6,5) * t205 + Ifges(6,6) * t204 + Ifges(6,3) * t220 + t212 * t230 - t213 * t229;
t211 = Ifges(6,5) * t230 + Ifges(6,6) * t229 + Ifges(6,3) * t246;
t184 = -mrSges(6,1) * t197 + mrSges(6,3) * t196 + Ifges(6,4) * t205 + Ifges(6,2) * t204 + Ifges(6,6) * t220 - t211 * t230 + t213 * t246;
t185 = mrSges(6,2) * t197 - mrSges(6,3) * t195 + Ifges(6,1) * t205 + Ifges(6,4) * t204 + Ifges(6,5) * t220 + t211 * t229 - t212 * t246;
t231 = Ifges(5,5) * t263 + Ifges(5,6) * t262 + Ifges(5,3) * t291;
t233 = Ifges(5,1) * t263 + Ifges(5,4) * t262 + Ifges(5,5) * t291;
t160 = -mrSges(5,1) * t216 + mrSges(5,3) * t203 + Ifges(5,4) * t227 + Ifges(5,2) * t226 + Ifges(5,6) * t290 - pkin(4) * t173 - t181 * t304 - t231 * t263 + t233 * t291 + (pkin(11) * t179 + t184 * t313 + t185 * t308) * t306;
t232 = Ifges(5,4) * t263 + Ifges(5,2) * t262 + Ifges(5,6) * t291;
t161 = mrSges(5,2) * t216 - mrSges(5,3) * t202 + Ifges(5,1) * t227 + Ifges(5,4) * t226 + Ifges(5,5) * t290 - t184 * t308 + t185 * t313 + t231 * t262 - t232 * t291 + (-t173 * t304 - t174 * t306) * pkin(11);
t255 = Ifges(4,5) * t278 + Ifges(4,6) * t277 + Ifges(4,3) * t295;
t257 = Ifges(4,1) * t278 + Ifges(4,4) * t277 + Ifges(4,5) * t295;
t143 = -mrSges(4,1) * t243 + mrSges(4,3) * t219 + Ifges(4,4) * t252 + Ifges(4,2) * t251 + Ifges(4,6) * t292 - pkin(3) * t322 + pkin(10) * t328 + t314 * t160 + t309 * t161 - t278 * t255 + t295 * t257;
t256 = Ifges(4,4) * t278 + Ifges(4,2) * t277 + Ifges(4,6) * t295;
t146 = mrSges(4,2) * t243 - mrSges(4,3) * t218 + Ifges(4,1) * t252 + Ifges(4,4) * t251 + Ifges(4,5) * t292 - pkin(10) * t168 - t160 * t309 + t161 * t314 + t255 * t277 - t256 * t295;
t268 = Ifges(3,3) * t297 + (Ifges(3,5) * t311 + Ifges(3,6) * t316) * t334;
t270 = Ifges(3,5) * t297 + (Ifges(3,1) * t311 + Ifges(3,4) * t316) * t334;
t140 = -mrSges(3,1) * t271 + mrSges(3,3) * t261 + Ifges(3,4) * t285 + Ifges(3,2) * t286 + Ifges(3,6) * t296 - pkin(2) * t320 + pkin(9) * t329 + t315 * t143 + t310 * t146 - t268 * t331 + t297 * t270;
t269 = Ifges(3,6) * t297 + (Ifges(3,4) * t311 + Ifges(3,2) * t316) * t334;
t142 = mrSges(3,2) * t271 - mrSges(3,3) * t260 + Ifges(3,1) * t285 + Ifges(3,4) * t286 + Ifges(3,5) * t296 - pkin(9) * t159 - t143 * t310 + t146 * t315 + t268 * t330 - t269 * t297;
t321 = mrSges(5,1) * t202 - mrSges(5,2) * t203 + Ifges(5,5) * t227 + Ifges(5,6) * t226 + Ifges(5,3) * t290 + pkin(4) * t174 + t179 * t343 + t181 * t306 + t184 * t339 + t185 * t340 + t263 * t232 - t262 * t233;
t319 = mrSges(4,1) * t218 - mrSges(4,2) * t219 + Ifges(4,5) * t252 + Ifges(4,6) * t251 + Ifges(4,3) * t292 + pkin(3) * t168 + t278 * t256 - t277 * t257 + t321;
t145 = pkin(2) * t159 + t319 + (t269 * t311 - t270 * t316) * t334 + Ifges(3,3) * t296 + Ifges(3,5) * t285 + Ifges(3,6) * t286 + mrSges(3,1) * t260 - mrSges(3,2) * t261;
t323 = mrSges(2,1) * t293 - mrSges(2,2) * t294 + Ifges(2,3) * qJDD(1) + pkin(1) * t150 + t140 * t337 + t142 * t338 + t145 * t307 + t153 * t344;
t151 = m(2) * t294 - mrSges(2,1) * t318 - qJDD(1) * mrSges(2,2) + t153;
t149 = t307 * t169 + (t157 * t316 + t158 * t311) * t305;
t147 = m(2) * t293 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t318 + t150;
t138 = -mrSges(2,2) * g(3) - mrSges(2,3) * t293 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t318 - t140 * t311 + t142 * t316 + (-t149 * t305 - t150 * t307) * pkin(8);
t137 = mrSges(2,1) * g(3) + mrSges(2,3) * t294 + Ifges(2,5) * t318 + Ifges(2,6) * qJDD(1) - pkin(1) * t149 - t145 * t305 + (pkin(8) * t153 + t140 * t316 + t142 * t311) * t307;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t317 * t138 - t312 * t137 - pkin(7) * (t147 * t317 + t151 * t312), t138, t142, t146, t161, t185; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t312 * t138 + t317 * t137 + pkin(7) * (-t147 * t312 + t151 * t317), t137, t140, t143, t160, t184; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t323, t323, t145, t319, t321, t181;];
m_new = t1;
