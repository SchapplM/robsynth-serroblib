% Calculate vector of cutting torques with Newton-Euler for
% S6PRRRRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-05-05 10:47
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRRRRR2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR2_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR2_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRR2_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR2_invdynm_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR2_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR2_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR2_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 10:41:50
% EndTime: 2019-05-05 10:42:36
% DurationCPUTime: 32.52s
% Computational Cost: add. (636350->341), mult. (1244135->438), div. (0->0), fcn. (915570->14), ass. (0->146)
t308 = sin(pkin(12));
t310 = cos(pkin(12));
t295 = g(1) * t308 - g(2) * t310;
t296 = -g(1) * t310 - g(2) * t308;
t307 = -g(3) + qJDD(1);
t321 = cos(qJ(2));
t311 = cos(pkin(6));
t316 = sin(qJ(2));
t343 = t311 * t316;
t309 = sin(pkin(6));
t344 = t309 * t316;
t264 = t295 * t343 + t321 * t296 + t307 * t344;
t322 = qJD(2) ^ 2;
t259 = -pkin(2) * t322 + qJDD(2) * pkin(8) + t264;
t275 = -t295 * t309 + t307 * t311;
t315 = sin(qJ(3));
t320 = cos(qJ(3));
t239 = -t315 * t259 + t320 * t275;
t340 = qJD(2) * qJD(3);
t339 = t320 * t340;
t292 = qJDD(2) * t315 + t339;
t227 = (-t292 + t339) * pkin(9) + (t315 * t320 * t322 + qJDD(3)) * pkin(3) + t239;
t240 = t320 * t259 + t315 * t275;
t293 = qJDD(2) * t320 - t315 * t340;
t342 = qJD(2) * t315;
t300 = qJD(3) * pkin(3) - pkin(9) * t342;
t306 = t320 ^ 2;
t228 = -pkin(3) * t306 * t322 + pkin(9) * t293 - qJD(3) * t300 + t240;
t314 = sin(qJ(4));
t319 = cos(qJ(4));
t217 = t314 * t227 + t319 * t228;
t284 = (t314 * t320 + t315 * t319) * qJD(2);
t253 = -t284 * qJD(4) - t314 * t292 + t293 * t319;
t341 = qJD(2) * t320;
t283 = -t314 * t342 + t319 * t341;
t266 = -mrSges(5,1) * t283 + mrSges(5,2) * t284;
t305 = qJD(3) + qJD(4);
t274 = mrSges(5,1) * t305 - mrSges(5,3) * t284;
t304 = qJDD(3) + qJDD(4);
t267 = -pkin(4) * t283 - pkin(10) * t284;
t303 = t305 ^ 2;
t206 = -pkin(4) * t303 + pkin(10) * t304 + t267 * t283 + t217;
t263 = -t316 * t296 + (t295 * t311 + t307 * t309) * t321;
t329 = -qJDD(2) * pkin(2) - t263;
t234 = -t293 * pkin(3) + t300 * t342 + (-pkin(9) * t306 - pkin(8)) * t322 + t329;
t254 = qJD(4) * t283 + t292 * t319 + t293 * t314;
t214 = (-t283 * t305 - t254) * pkin(10) + (t284 * t305 - t253) * pkin(4) + t234;
t313 = sin(qJ(5));
t318 = cos(qJ(5));
t201 = -t313 * t206 + t318 * t214;
t269 = -t284 * t313 + t305 * t318;
t231 = qJD(5) * t269 + t254 * t318 + t304 * t313;
t251 = qJDD(5) - t253;
t270 = t284 * t318 + t305 * t313;
t278 = qJD(5) - t283;
t199 = (t269 * t278 - t231) * pkin(11) + (t269 * t270 + t251) * pkin(5) + t201;
t202 = t318 * t206 + t313 * t214;
t230 = -qJD(5) * t270 - t254 * t313 + t304 * t318;
t257 = pkin(5) * t278 - pkin(11) * t270;
t268 = t269 ^ 2;
t200 = -pkin(5) * t268 + pkin(11) * t230 - t257 * t278 + t202;
t312 = sin(qJ(6));
t317 = cos(qJ(6));
t197 = t199 * t317 - t200 * t312;
t241 = t269 * t317 - t270 * t312;
t211 = qJD(6) * t241 + t230 * t312 + t231 * t317;
t242 = t269 * t312 + t270 * t317;
t223 = -mrSges(7,1) * t241 + mrSges(7,2) * t242;
t276 = qJD(6) + t278;
t232 = -mrSges(7,2) * t276 + mrSges(7,3) * t241;
t246 = qJDD(6) + t251;
t192 = m(7) * t197 + mrSges(7,1) * t246 - mrSges(7,3) * t211 - t223 * t242 + t232 * t276;
t198 = t199 * t312 + t200 * t317;
t210 = -qJD(6) * t242 + t230 * t317 - t231 * t312;
t233 = mrSges(7,1) * t276 - mrSges(7,3) * t242;
t193 = m(7) * t198 - mrSges(7,2) * t246 + t210 * mrSges(7,3) + t223 * t241 - t233 * t276;
t184 = t317 * t192 + t312 * t193;
t243 = -mrSges(6,1) * t269 + mrSges(6,2) * t270;
t255 = -mrSges(6,2) * t278 + mrSges(6,3) * t269;
t182 = m(6) * t201 + mrSges(6,1) * t251 - mrSges(6,3) * t231 - t243 * t270 + t255 * t278 + t184;
t256 = mrSges(6,1) * t278 - mrSges(6,3) * t270;
t335 = -t192 * t312 + t317 * t193;
t183 = m(6) * t202 - mrSges(6,2) * t251 + mrSges(6,3) * t230 + t243 * t269 - t256 * t278 + t335;
t336 = -t182 * t313 + t318 * t183;
t175 = m(5) * t217 - mrSges(5,2) * t304 + mrSges(5,3) * t253 + t266 * t283 - t274 * t305 + t336;
t216 = t227 * t319 - t314 * t228;
t273 = -mrSges(5,2) * t305 + mrSges(5,3) * t283;
t205 = -pkin(4) * t304 - pkin(10) * t303 + t284 * t267 - t216;
t203 = -pkin(5) * t230 - pkin(11) * t268 + t257 * t270 + t205;
t331 = m(7) * t203 - t210 * mrSges(7,1) + mrSges(7,2) * t211 - t241 * t232 + t233 * t242;
t326 = -m(6) * t205 + t230 * mrSges(6,1) - mrSges(6,2) * t231 + t269 * t255 - t256 * t270 - t331;
t188 = m(5) * t216 + mrSges(5,1) * t304 - mrSges(5,3) * t254 - t266 * t284 + t273 * t305 + t326;
t166 = t314 * t175 + t319 * t188;
t291 = (-mrSges(4,1) * t320 + mrSges(4,2) * t315) * qJD(2);
t298 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t341;
t164 = m(4) * t239 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t292 + qJD(3) * t298 - t291 * t342 + t166;
t297 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t342;
t337 = t319 * t175 - t188 * t314;
t165 = m(4) * t240 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t293 - qJD(3) * t297 + t291 * t341 + t337;
t159 = t320 * t164 + t315 * t165;
t281 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t315 + Ifges(4,2) * t320) * qJD(2);
t282 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t315 + Ifges(4,4) * t320) * qJD(2);
t219 = Ifges(7,5) * t242 + Ifges(7,6) * t241 + Ifges(7,3) * t276;
t221 = Ifges(7,1) * t242 + Ifges(7,4) * t241 + Ifges(7,5) * t276;
t185 = -mrSges(7,1) * t203 + mrSges(7,3) * t198 + Ifges(7,4) * t211 + Ifges(7,2) * t210 + Ifges(7,6) * t246 - t219 * t242 + t221 * t276;
t220 = Ifges(7,4) * t242 + Ifges(7,2) * t241 + Ifges(7,6) * t276;
t186 = mrSges(7,2) * t203 - mrSges(7,3) * t197 + Ifges(7,1) * t211 + Ifges(7,4) * t210 + Ifges(7,5) * t246 + t219 * t241 - t220 * t276;
t235 = Ifges(6,5) * t270 + Ifges(6,6) * t269 + Ifges(6,3) * t278;
t237 = Ifges(6,1) * t270 + Ifges(6,4) * t269 + Ifges(6,5) * t278;
t168 = -mrSges(6,1) * t205 + mrSges(6,3) * t202 + Ifges(6,4) * t231 + Ifges(6,2) * t230 + Ifges(6,6) * t251 - pkin(5) * t331 + pkin(11) * t335 + t317 * t185 + t312 * t186 - t270 * t235 + t278 * t237;
t236 = Ifges(6,4) * t270 + Ifges(6,2) * t269 + Ifges(6,6) * t278;
t170 = mrSges(6,2) * t205 - mrSges(6,3) * t201 + Ifges(6,1) * t231 + Ifges(6,4) * t230 + Ifges(6,5) * t251 - pkin(11) * t184 - t185 * t312 + t186 * t317 + t235 * t269 - t236 * t278;
t261 = Ifges(5,4) * t284 + Ifges(5,2) * t283 + Ifges(5,6) * t305;
t262 = Ifges(5,1) * t284 + Ifges(5,4) * t283 + Ifges(5,5) * t305;
t327 = -mrSges(5,1) * t216 + mrSges(5,2) * t217 - Ifges(5,5) * t254 - Ifges(5,6) * t253 - Ifges(5,3) * t304 - pkin(4) * t326 - pkin(10) * t336 - t318 * t168 - t313 * t170 - t284 * t261 + t283 * t262;
t348 = mrSges(4,1) * t239 - mrSges(4,2) * t240 + Ifges(4,5) * t292 + Ifges(4,6) * t293 + Ifges(4,3) * qJDD(3) + pkin(3) * t166 + (t281 * t315 - t282 * t320) * qJD(2) - t327;
t145 = -mrSges(3,1) * t275 + mrSges(3,3) * t264 + t322 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t159 - t348;
t338 = -t164 * t315 + t320 * t165;
t157 = m(3) * t264 - mrSges(3,1) * t322 - qJDD(2) * mrSges(3,2) + t338;
t258 = -t322 * pkin(8) + t329;
t177 = t318 * t182 + t313 * t183;
t328 = m(5) * t234 - t253 * mrSges(5,1) + mrSges(5,2) * t254 - t283 * t273 + t274 * t284 + t177;
t325 = -m(4) * t258 + t293 * mrSges(4,1) - mrSges(4,2) * t292 - t297 * t342 + t298 * t341 - t328;
t172 = m(3) * t263 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t322 + t325;
t153 = t321 * t157 - t172 * t316;
t349 = pkin(7) * t153 + t145 * t321;
t345 = t172 * t321;
t158 = m(3) * t275 + t159;
t150 = t157 * t343 - t158 * t309 + t311 * t345;
t260 = Ifges(5,5) * t284 + Ifges(5,6) * t283 + Ifges(5,3) * t305;
t154 = mrSges(5,2) * t234 - mrSges(5,3) * t216 + Ifges(5,1) * t254 + Ifges(5,4) * t253 + Ifges(5,5) * t304 - pkin(10) * t177 - t168 * t313 + t170 * t318 + t260 * t283 - t261 * t305;
t330 = -mrSges(7,1) * t197 + mrSges(7,2) * t198 - Ifges(7,5) * t211 - Ifges(7,6) * t210 - Ifges(7,3) * t246 - t242 * t220 + t241 * t221;
t323 = mrSges(6,1) * t201 - mrSges(6,2) * t202 + Ifges(6,5) * t231 + Ifges(6,6) * t230 + Ifges(6,3) * t251 + pkin(5) * t184 + t270 * t236 - t269 * t237 - t330;
t160 = -mrSges(5,1) * t234 + mrSges(5,3) * t217 + Ifges(5,4) * t254 + Ifges(5,2) * t253 + Ifges(5,6) * t304 - pkin(4) * t177 - t284 * t260 + t305 * t262 - t323;
t280 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t315 + Ifges(4,6) * t320) * qJD(2);
t143 = -mrSges(4,1) * t258 + mrSges(4,3) * t240 + Ifges(4,4) * t292 + Ifges(4,2) * t293 + Ifges(4,6) * qJDD(3) - pkin(3) * t328 + pkin(9) * t337 + qJD(3) * t282 + t314 * t154 + t319 * t160 - t280 * t342;
t146 = mrSges(4,2) * t258 - mrSges(4,3) * t239 + Ifges(4,1) * t292 + Ifges(4,4) * t293 + Ifges(4,5) * qJDD(3) - pkin(9) * t166 - qJD(3) * t281 + t154 * t319 - t160 * t314 + t280 * t341;
t140 = mrSges(3,1) * t263 - mrSges(3,2) * t264 + Ifges(3,3) * qJDD(2) + pkin(2) * t325 + pkin(8) * t338 + t320 * t143 + t315 * t146;
t142 = mrSges(3,2) * t275 - mrSges(3,3) * t263 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t322 - pkin(8) * t159 - t143 * t315 + t146 * t320;
t332 = mrSges(2,1) * t295 - mrSges(2,2) * t296 + pkin(1) * t150 + t311 * t140 + t142 * t344 + t349 * t309;
t151 = m(2) * t296 + t153;
t149 = t311 * t158 + (t157 * t316 + t345) * t309;
t147 = m(2) * t295 + t150;
t138 = mrSges(2,2) * t307 - mrSges(2,3) * t295 + t321 * t142 - t316 * t145 + (-t149 * t309 - t150 * t311) * pkin(7);
t137 = -mrSges(2,1) * t307 + mrSges(2,3) * t296 - pkin(1) * t149 - t309 * t140 + (t142 * t316 + t349) * t311;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t310 * t138 - t308 * t137 - qJ(1) * (t147 * t310 + t151 * t308), t138, t142, t146, t154, t170, t186; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t308 * t138 + t310 * t137 + qJ(1) * (-t147 * t308 + t151 * t310), t137, t145, t143, t160, t168, t185; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t332, t332, t140, t348, -t327, t323, -t330;];
m_new  = t1;
