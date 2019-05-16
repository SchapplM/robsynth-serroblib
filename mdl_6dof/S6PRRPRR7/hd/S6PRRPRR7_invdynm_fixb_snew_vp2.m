% Calculate vector of cutting torques with Newton-Euler for
% S6PRRPRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
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
% Datum: 2019-05-05 06:06
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRRPRR7_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR7_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR7_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR7_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR7_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR7_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR7_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR7_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR7_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 06:00:42
% EndTime: 2019-05-05 06:01:01
% DurationCPUTime: 11.36s
% Computational Cost: add. (204002->348), mult. (405087->425), div. (0->0), fcn. (256680->12), ass. (0->145)
t365 = -2 * qJD(4);
t311 = sin(pkin(11));
t313 = cos(pkin(11));
t291 = g(1) * t311 - g(2) * t313;
t292 = -g(1) * t313 - g(2) * t311;
t310 = -g(3) + qJDD(1);
t322 = cos(qJ(2));
t314 = cos(pkin(6));
t318 = sin(qJ(2));
t353 = t314 * t318;
t312 = sin(pkin(6));
t354 = t312 * t318;
t238 = t291 * t353 + t322 * t292 + t310 * t354;
t324 = qJD(2) ^ 2;
t234 = -pkin(2) * t324 + qJDD(2) * pkin(8) + t238;
t317 = sin(qJ(3));
t231 = t317 * t234;
t255 = -t291 * t312 + t310 * t314;
t321 = cos(qJ(3));
t352 = t321 * t255;
t227 = -t231 + t352;
t284 = (mrSges(5,2) * t321 - mrSges(5,3) * t317) * qJD(2);
t285 = (-mrSges(4,1) * t321 + mrSges(4,2) * t317) * qJD(2);
t347 = qJD(2) * qJD(3);
t345 = t321 * t347;
t286 = qJDD(2) * t317 + t345;
t348 = qJD(2) * t321;
t294 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t348;
t295 = -mrSges(5,1) * t348 - qJD(3) * mrSges(5,3);
t304 = t317 * qJD(2);
t283 = (-pkin(3) * t321 - qJ(4) * t317) * qJD(2);
t323 = qJD(3) ^ 2;
t340 = -t323 * qJ(4) + t283 * t304 + qJDD(4) + t231;
t359 = pkin(9) * t324;
t361 = -pkin(3) - pkin(9);
t210 = t286 * pkin(4) + t361 * qJDD(3) + (-pkin(4) * t347 - t317 * t359 - t255) * t321 + t340;
t344 = t317 * t347;
t287 = qJDD(2) * t321 - t344;
t298 = pkin(4) * t304 - qJD(3) * pkin(9);
t309 = t321 ^ 2;
t237 = -t318 * t292 + (t291 * t314 + t310 * t312) * t322;
t333 = -qJDD(2) * pkin(2) - t237;
t330 = pkin(3) * t344 + t304 * t365 + (-t286 - t345) * qJ(4) + t333;
t212 = -t298 * t304 + (-pkin(4) * t309 - pkin(8)) * t324 + t361 * t287 + t330;
t316 = sin(qJ(5));
t320 = cos(qJ(5));
t202 = t320 * t210 - t316 * t212;
t281 = -qJD(3) * t316 - t320 * t348;
t247 = qJD(5) * t281 + qJDD(3) * t320 - t287 * t316;
t278 = qJDD(5) + t286;
t282 = qJD(3) * t320 - t316 * t348;
t301 = t304 + qJD(5);
t199 = (t281 * t301 - t247) * pkin(10) + (t281 * t282 + t278) * pkin(5) + t202;
t203 = t316 * t210 + t320 * t212;
t246 = -qJD(5) * t282 - qJDD(3) * t316 - t287 * t320;
t254 = pkin(5) * t301 - pkin(10) * t282;
t277 = t281 ^ 2;
t200 = -pkin(5) * t277 + pkin(10) * t246 - t254 * t301 + t203;
t315 = sin(qJ(6));
t319 = cos(qJ(6));
t197 = t199 * t319 - t200 * t315;
t248 = t281 * t319 - t282 * t315;
t218 = qJD(6) * t248 + t246 * t315 + t247 * t319;
t249 = t281 * t315 + t282 * t319;
t229 = -mrSges(7,1) * t248 + mrSges(7,2) * t249;
t299 = qJD(6) + t301;
t235 = -mrSges(7,2) * t299 + mrSges(7,3) * t248;
t269 = qJDD(6) + t278;
t193 = m(7) * t197 + mrSges(7,1) * t269 - t218 * mrSges(7,3) - t229 * t249 + t235 * t299;
t198 = t199 * t315 + t200 * t319;
t217 = -qJD(6) * t249 + t246 * t319 - t247 * t315;
t236 = mrSges(7,1) * t299 - mrSges(7,3) * t249;
t194 = m(7) * t198 - mrSges(7,2) * t269 + t217 * mrSges(7,3) + t229 * t248 - t236 * t299;
t183 = t319 * t193 + t315 * t194;
t250 = -mrSges(6,1) * t281 + mrSges(6,2) * t282;
t252 = -mrSges(6,2) * t301 + mrSges(6,3) * t281;
t180 = m(6) * t202 + mrSges(6,1) * t278 - mrSges(6,3) * t247 - t250 * t282 + t252 * t301 + t183;
t253 = mrSges(6,1) * t301 - mrSges(6,3) * t282;
t342 = -t193 * t315 + t319 * t194;
t181 = m(6) * t203 - mrSges(6,2) * t278 + mrSges(6,3) * t246 + t250 * t281 - t253 * t301 + t342;
t176 = t320 * t180 + t316 * t181;
t221 = -qJDD(3) * pkin(3) + t340 - t352;
t336 = -m(5) * t221 - t286 * mrSges(5,1) - t176;
t173 = m(4) * t227 - t286 * mrSges(4,3) + (mrSges(4,1) - mrSges(5,2)) * qJDD(3) + (t294 - t295) * qJD(3) + (-t284 - t285) * t304 + t336;
t228 = t321 * t234 + t317 * t255;
t293 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t304;
t219 = pkin(3) * t323 - qJDD(3) * qJ(4) + qJD(3) * t365 - t283 * t348 - t228;
t209 = pkin(4) * t287 + qJD(3) * t298 - t309 * t359 - t219;
t205 = -pkin(5) * t246 - pkin(10) * t277 + t254 * t282 + t209;
t337 = m(7) * t205 - t217 * mrSges(7,1) + t218 * mrSges(7,2) - t248 * t235 + t249 * t236;
t195 = -m(6) * t209 + t246 * mrSges(6,1) - t247 * mrSges(6,2) + t281 * t252 - t282 * t253 - t337;
t296 = mrSges(5,1) * t304 + qJD(3) * mrSges(5,2);
t328 = -m(5) * t219 + qJDD(3) * mrSges(5,3) + qJD(3) * t296 + t284 * t348 - t195;
t188 = -qJDD(3) * mrSges(4,2) + t285 * t348 - qJD(3) * t293 + m(4) * t228 + (mrSges(4,3) + mrSges(5,1)) * t287 + t328;
t166 = t321 * t173 + t317 * t188;
t264 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t317 + Ifges(4,4) * t321) * qJD(2);
t223 = Ifges(7,5) * t249 + Ifges(7,6) * t248 + Ifges(7,3) * t299;
t225 = Ifges(7,1) * t249 + Ifges(7,4) * t248 + Ifges(7,5) * t299;
t184 = -mrSges(7,1) * t205 + mrSges(7,3) * t198 + Ifges(7,4) * t218 + Ifges(7,2) * t217 + Ifges(7,6) * t269 - t223 * t249 + t225 * t299;
t224 = Ifges(7,4) * t249 + Ifges(7,2) * t248 + Ifges(7,6) * t299;
t185 = mrSges(7,2) * t205 - mrSges(7,3) * t197 + Ifges(7,1) * t218 + Ifges(7,4) * t217 + Ifges(7,5) * t269 + t223 * t248 - t224 * t299;
t239 = Ifges(6,5) * t282 + Ifges(6,6) * t281 + Ifges(6,3) * t301;
t241 = Ifges(6,1) * t282 + Ifges(6,4) * t281 + Ifges(6,5) * t301;
t167 = -mrSges(6,1) * t209 + mrSges(6,3) * t203 + Ifges(6,4) * t247 + Ifges(6,2) * t246 + Ifges(6,6) * t278 - pkin(5) * t337 + pkin(10) * t342 + t319 * t184 + t315 * t185 - t282 * t239 + t301 * t241;
t240 = Ifges(6,4) * t282 + Ifges(6,2) * t281 + Ifges(6,6) * t301;
t169 = mrSges(6,2) * t209 - mrSges(6,3) * t202 + Ifges(6,1) * t247 + Ifges(6,4) * t246 + Ifges(6,5) * t278 - pkin(10) * t183 - t184 * t315 + t185 * t319 + t239 * t281 - t240 * t301;
t266 = Ifges(5,4) * qJD(3) + (-Ifges(5,2) * t317 - Ifges(5,6) * t321) * qJD(2);
t332 = -mrSges(5,2) * t221 + mrSges(5,3) * t219 - Ifges(5,1) * qJDD(3) + Ifges(5,4) * t286 + Ifges(5,5) * t287 + pkin(9) * t176 + t316 * t167 - t320 * t169 - t266 * t348;
t265 = Ifges(5,5) * qJD(3) + (-Ifges(5,6) * t317 - Ifges(5,3) * t321) * qJD(2);
t349 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t317 + Ifges(4,2) * t321) * qJD(2) - t265;
t363 = qJD(2) * (-t321 * t264 + t349 * t317) + mrSges(4,1) * t227 - mrSges(4,2) * t228 + Ifges(4,5) * t286 + Ifges(4,6) * t287 + Ifges(4,3) * qJDD(3) + pkin(3) * (-qJDD(3) * mrSges(5,2) - qJD(3) * t295 - t284 * t304 + t336) + qJ(4) * (mrSges(5,1) * t287 + t328) - t332;
t152 = -mrSges(3,1) * t255 + mrSges(3,3) * t238 + t324 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t166 - t363;
t343 = -t173 * t317 + t321 * t188;
t164 = m(3) * t238 - mrSges(3,1) * t324 - qJDD(2) * mrSges(3,2) + t343;
t358 = t324 * pkin(8);
t233 = t333 - t358;
t177 = -t316 * t180 + t320 * t181;
t222 = -t287 * pkin(3) + t330 - t358;
t339 = -m(5) * t222 - t287 * mrSges(5,2) + t296 * t304 - t177;
t327 = -m(4) * t233 + t294 * t348 + t287 * mrSges(4,1) + (-mrSges(4,2) + mrSges(5,3)) * t286 + (-t293 * t317 - t295 * t321) * qJD(2) + t339;
t171 = m(3) * t237 + qJDD(2) * mrSges(3,1) - t324 * mrSges(3,2) + t327;
t161 = t322 * t164 - t171 * t318;
t364 = pkin(7) * t161 + t152 * t322;
t357 = Ifges(4,4) + Ifges(5,6);
t355 = t171 * t322;
t267 = Ifges(5,1) * qJD(3) + (-Ifges(5,4) * t317 - Ifges(5,5) * t321) * qJD(2);
t350 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t317 + Ifges(4,6) * t321) * qJD(2) + t267;
t165 = m(3) * t255 + t166;
t157 = t164 * t353 - t165 * t312 + t314 * t355;
t174 = -t286 * mrSges(5,3) + t295 * t348 - t339;
t331 = -mrSges(5,1) * t219 + mrSges(5,2) * t222 - pkin(4) * t195 - pkin(9) * t177 - t320 * t167 - t316 * t169;
t153 = -mrSges(4,1) * t233 + mrSges(4,3) * t228 - pkin(3) * t174 + (Ifges(4,2) + Ifges(5,3)) * t287 + t357 * t286 + (Ifges(4,6) - Ifges(5,5)) * qJDD(3) + (t264 - t266) * qJD(3) - t350 * t304 + t331;
t334 = -mrSges(7,1) * t197 + mrSges(7,2) * t198 - Ifges(7,5) * t218 - Ifges(7,6) * t217 - Ifges(7,3) * t269 - t249 * t224 + t248 * t225;
t329 = -mrSges(6,1) * t202 + mrSges(6,2) * t203 - Ifges(6,5) * t247 - Ifges(6,6) * t246 - Ifges(6,3) * t278 - pkin(5) * t183 - t282 * t240 + t281 * t241 + t334;
t326 = -mrSges(5,1) * t221 + mrSges(5,3) * t222 - pkin(4) * t176 + t329;
t158 = -t326 + t357 * t287 + (Ifges(5,2) + Ifges(4,1)) * t286 + (Ifges(4,5) - Ifges(5,4)) * qJDD(3) - t349 * qJD(3) + mrSges(4,2) * t233 - mrSges(4,3) * t227 - qJ(4) * t174 + t350 * t348;
t148 = mrSges(3,1) * t237 - mrSges(3,2) * t238 + Ifges(3,3) * qJDD(2) + pkin(2) * t327 + pkin(8) * t343 + t321 * t153 + t317 * t158;
t150 = mrSges(3,2) * t255 - mrSges(3,3) * t237 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t324 - pkin(8) * t166 - t153 * t317 + t158 * t321;
t335 = mrSges(2,1) * t291 - mrSges(2,2) * t292 + pkin(1) * t157 + t314 * t148 + t150 * t354 + t312 * t364;
t159 = m(2) * t292 + t161;
t156 = t314 * t165 + (t164 * t318 + t355) * t312;
t154 = m(2) * t291 + t157;
t146 = mrSges(2,2) * t310 - mrSges(2,3) * t291 + t322 * t150 - t318 * t152 + (-t156 * t312 - t157 * t314) * pkin(7);
t145 = -mrSges(2,1) * t310 + mrSges(2,3) * t292 - pkin(1) * t156 - t312 * t148 + (t150 * t318 + t364) * t314;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t313 * t146 - t311 * t145 - qJ(1) * (t154 * t313 + t159 * t311), t146, t150, t158, -t265 * t304 - t332, t169, t185; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t311 * t146 + t313 * t145 + qJ(1) * (-t154 * t311 + t159 * t313), t145, t152, t153, Ifges(5,4) * qJDD(3) - Ifges(5,2) * t286 - Ifges(5,6) * t287 - qJD(3) * t265 - t267 * t348 + t326, t167, t184; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t335, t335, t148, t363, Ifges(5,5) * qJDD(3) - Ifges(5,6) * t286 - Ifges(5,3) * t287 + qJD(3) * t266 + t267 * t304 - t331, -t329, -t334;];
m_new  = t1;
