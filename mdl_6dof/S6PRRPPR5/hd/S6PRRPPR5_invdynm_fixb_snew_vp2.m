% Calculate vector of cutting torques with Newton-Euler for
% S6PRRPPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta5]';
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
% Datum: 2019-05-05 03:29
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRRPPR5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR5_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR5_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPPR5_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR5_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR5_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR5_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPPR5_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 03:24:40
% EndTime: 2019-05-05 03:24:57
% DurationCPUTime: 10.97s
% Computational Cost: add. (189225->347), mult. (395013->427), div. (0->0), fcn. (247317->12), ass. (0->143)
t365 = -2 * qJD(4);
t311 = sin(pkin(10));
t314 = cos(pkin(10));
t292 = g(1) * t311 - g(2) * t314;
t293 = -g(1) * t314 - g(2) * t311;
t309 = -g(3) + qJDD(1);
t321 = cos(qJ(2));
t315 = cos(pkin(6));
t318 = sin(qJ(2));
t353 = t315 * t318;
t312 = sin(pkin(6));
t354 = t312 * t318;
t238 = t292 * t353 + t321 * t293 + t309 * t354;
t323 = qJD(2) ^ 2;
t234 = -pkin(2) * t323 + qJDD(2) * pkin(8) + t238;
t317 = sin(qJ(3));
t231 = t317 * t234;
t255 = -t292 * t312 + t309 * t315;
t320 = cos(qJ(3));
t352 = t320 * t255;
t227 = -t231 + t352;
t285 = (mrSges(5,2) * t320 - mrSges(5,3) * t317) * qJD(2);
t286 = (-mrSges(4,1) * t320 + mrSges(4,2) * t317) * qJD(2);
t346 = qJD(2) * qJD(3);
t344 = t320 * t346;
t287 = qJDD(2) * t317 + t344;
t347 = qJD(2) * t320;
t295 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t347;
t297 = -mrSges(5,1) * t347 - qJD(3) * mrSges(5,3);
t284 = (-pkin(3) * t320 - qJ(4) * t317) * qJD(2);
t322 = qJD(3) ^ 2;
t348 = qJD(2) * t317;
t339 = -t322 * qJ(4) + t284 * t348 + qJDD(4) + t231;
t357 = qJ(5) * t323;
t358 = -pkin(3) - qJ(5);
t210 = t287 * pkin(4) + t358 * qJDD(3) + (-pkin(4) * t346 - t317 * t357 - t255) * t320 + t339;
t343 = t317 * t346;
t288 = qJDD(2) * t320 - t343;
t296 = pkin(4) * t348 - qJD(3) * qJ(5);
t308 = t320 ^ 2;
t237 = -t318 * t293 + (t292 * t315 + t309 * t312) * t321;
t332 = -qJDD(2) * pkin(2) - t237;
t329 = pkin(3) * t343 + t348 * t365 + (-t287 - t344) * qJ(4) + t332;
t212 = -t296 * t348 + (-pkin(4) * t308 - pkin(8)) * t323 + t358 * t288 + t329;
t310 = sin(pkin(11));
t313 = cos(pkin(11));
t274 = qJD(3) * t313 - t310 * t347;
t202 = -0.2e1 * qJD(5) * t274 + t313 * t210 - t310 * t212;
t253 = qJDD(3) * t313 - t288 * t310;
t273 = -qJD(3) * t310 - t313 * t347;
t199 = (t273 * t348 - t253) * pkin(9) + (t273 * t274 + t287) * pkin(5) + t202;
t203 = 0.2e1 * qJD(5) * t273 + t310 * t210 + t313 * t212;
t252 = -qJDD(3) * t310 - t288 * t313;
t254 = pkin(5) * t348 - pkin(9) * t274;
t272 = t273 ^ 2;
t200 = -pkin(5) * t272 + pkin(9) * t252 - t254 * t348 + t203;
t316 = sin(qJ(6));
t319 = cos(qJ(6));
t197 = t199 * t319 - t200 * t316;
t243 = t273 * t319 - t274 * t316;
t222 = qJD(6) * t243 + t252 * t316 + t253 * t319;
t244 = t273 * t316 + t274 * t319;
t229 = -mrSges(7,1) * t243 + mrSges(7,2) * t244;
t301 = qJD(6) + t348;
t235 = -mrSges(7,2) * t301 + mrSges(7,3) * t243;
t281 = qJDD(6) + t287;
t192 = m(7) * t197 + mrSges(7,1) * t281 - mrSges(7,3) * t222 - t229 * t244 + t235 * t301;
t198 = t199 * t316 + t200 * t319;
t221 = -qJD(6) * t244 + t252 * t319 - t253 * t316;
t236 = mrSges(7,1) * t301 - mrSges(7,3) * t244;
t193 = m(7) * t198 - mrSges(7,2) * t281 + mrSges(7,3) * t221 + t229 * t243 - t236 * t301;
t183 = t319 * t192 + t316 * t193;
t245 = -mrSges(6,1) * t273 + mrSges(6,2) * t274;
t250 = -mrSges(6,2) * t348 + mrSges(6,3) * t273;
t180 = m(6) * t202 + mrSges(6,1) * t287 - mrSges(6,3) * t253 - t245 * t274 + t250 * t348 + t183;
t251 = mrSges(6,1) * t348 - mrSges(6,3) * t274;
t341 = -t192 * t316 + t319 * t193;
t181 = m(6) * t203 - mrSges(6,2) * t287 + mrSges(6,3) * t252 + t245 * t273 - t251 * t348 + t341;
t176 = t313 * t180 + t310 * t181;
t216 = -qJDD(3) * pkin(3) + t339 - t352;
t335 = -m(5) * t216 - t287 * mrSges(5,1) - t176;
t173 = m(4) * t227 - t287 * mrSges(4,3) + (mrSges(4,1) - mrSges(5,2)) * qJDD(3) + (t295 - t297) * qJD(3) + (-t285 - t286) * t348 + t335;
t228 = t320 * t234 + t317 * t255;
t294 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t348;
t214 = pkin(3) * t322 - qJDD(3) * qJ(4) + qJD(3) * t365 - t284 * t347 - t228;
t209 = pkin(4) * t288 + qJD(3) * t296 - t308 * t357 + qJDD(5) - t214;
t205 = -pkin(5) * t252 - pkin(9) * t272 + t254 * t274 + t209;
t336 = m(7) * t205 - t221 * mrSges(7,1) + t222 * mrSges(7,2) - t243 * t235 + t244 * t236;
t195 = -m(6) * t209 + t252 * mrSges(6,1) - t253 * mrSges(6,2) + t273 * t250 - t274 * t251 - t336;
t298 = mrSges(5,1) * t348 + qJD(3) * mrSges(5,2);
t327 = -m(5) * t214 + qJDD(3) * mrSges(5,3) + qJD(3) * t298 + t285 * t347 - t195;
t188 = -qJDD(3) * mrSges(4,2) + t286 * t347 + t327 - qJD(3) * t294 + m(4) * t228 + (mrSges(4,3) + mrSges(5,1)) * t288;
t166 = t320 * t173 + t317 * t188;
t265 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t317 + Ifges(4,4) * t320) * qJD(2);
t223 = Ifges(7,5) * t244 + Ifges(7,6) * t243 + Ifges(7,3) * t301;
t225 = Ifges(7,1) * t244 + Ifges(7,4) * t243 + Ifges(7,5) * t301;
t184 = -mrSges(7,1) * t205 + mrSges(7,3) * t198 + Ifges(7,4) * t222 + Ifges(7,2) * t221 + Ifges(7,6) * t281 - t223 * t244 + t225 * t301;
t224 = Ifges(7,4) * t244 + Ifges(7,2) * t243 + Ifges(7,6) * t301;
t185 = mrSges(7,2) * t205 - mrSges(7,3) * t197 + Ifges(7,1) * t222 + Ifges(7,4) * t221 + Ifges(7,5) * t281 + t223 * t243 - t224 * t301;
t239 = Ifges(6,5) * t274 + Ifges(6,6) * t273 + Ifges(6,3) * t348;
t241 = Ifges(6,1) * t274 + Ifges(6,4) * t273 + Ifges(6,5) * t348;
t167 = -mrSges(6,1) * t209 + mrSges(6,3) * t203 + Ifges(6,4) * t253 + Ifges(6,2) * t252 + Ifges(6,6) * t287 - pkin(5) * t336 + pkin(9) * t341 + t319 * t184 + t316 * t185 - t274 * t239 + t241 * t348;
t240 = Ifges(6,4) * t274 + Ifges(6,2) * t273 + Ifges(6,6) * t348;
t169 = mrSges(6,2) * t209 - mrSges(6,3) * t202 + Ifges(6,1) * t253 + Ifges(6,4) * t252 + Ifges(6,5) * t287 - pkin(9) * t183 - t184 * t316 + t185 * t319 + t239 * t273 - t240 * t348;
t267 = Ifges(5,4) * qJD(3) + (-Ifges(5,2) * t317 - Ifges(5,6) * t320) * qJD(2);
t331 = -mrSges(5,2) * t216 + mrSges(5,3) * t214 - Ifges(5,1) * qJDD(3) + Ifges(5,4) * t287 + Ifges(5,5) * t288 + qJ(5) * t176 + t310 * t167 - t313 * t169 - t267 * t347;
t266 = Ifges(5,5) * qJD(3) + (-Ifges(5,6) * t317 - Ifges(5,3) * t320) * qJD(2);
t349 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t317 + Ifges(4,2) * t320) * qJD(2) - t266;
t363 = (-t320 * t265 + t349 * t317) * qJD(2) + mrSges(4,1) * t227 - mrSges(4,2) * t228 + Ifges(4,5) * t287 + Ifges(4,6) * t288 + Ifges(4,3) * qJDD(3) + pkin(3) * (-qJDD(3) * mrSges(5,2) - qJD(3) * t297 - t285 * t348 + t335) + qJ(4) * (mrSges(5,1) * t288 + t327) - t331;
t152 = -mrSges(3,1) * t255 + mrSges(3,3) * t238 + t323 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t166 - t363;
t342 = -t173 * t317 + t320 * t188;
t164 = m(3) * t238 - mrSges(3,1) * t323 - qJDD(2) * mrSges(3,2) + t342;
t360 = t323 * pkin(8);
t233 = t332 - t360;
t177 = -t310 * t180 + t313 * t181;
t217 = -t288 * pkin(3) + t329 - t360;
t338 = -m(5) * t217 - t288 * mrSges(5,2) + t298 * t348 - t177;
t326 = -m(4) * t233 + t295 * t347 + t288 * mrSges(4,1) + (-mrSges(4,2) + mrSges(5,3)) * t287 + (-t294 * t317 - t297 * t320) * qJD(2) + t338;
t171 = m(3) * t237 + qJDD(2) * mrSges(3,1) - t323 * mrSges(3,2) + t326;
t161 = t321 * t164 - t171 * t318;
t364 = pkin(7) * t161 + t152 * t321;
t359 = Ifges(4,4) + Ifges(5,6);
t355 = t171 * t321;
t268 = Ifges(5,1) * qJD(3) + (-Ifges(5,4) * t317 - Ifges(5,5) * t320) * qJD(2);
t350 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t317 + Ifges(4,6) * t320) * qJD(2) + t268;
t165 = m(3) * t255 + t166;
t157 = t164 * t353 - t165 * t312 + t315 * t355;
t174 = -t287 * mrSges(5,3) + t297 * t347 - t338;
t330 = -mrSges(5,1) * t214 + mrSges(5,2) * t217 - pkin(4) * t195 - qJ(5) * t177 - t313 * t167 - t310 * t169;
t153 = -mrSges(4,1) * t233 + mrSges(4,3) * t228 - pkin(3) * t174 + (Ifges(4,2) + Ifges(5,3)) * t288 + t359 * t287 + (Ifges(4,6) - Ifges(5,5)) * qJDD(3) + (t265 - t267) * qJD(3) - t350 * t348 + t330;
t333 = -mrSges(7,1) * t197 + mrSges(7,2) * t198 - Ifges(7,5) * t222 - Ifges(7,6) * t221 - Ifges(7,3) * t281 - t244 * t224 + t243 * t225;
t328 = -mrSges(6,1) * t202 + mrSges(6,2) * t203 - Ifges(6,5) * t253 - Ifges(6,6) * t252 - Ifges(6,3) * t287 - pkin(5) * t183 - t274 * t240 + t273 * t241 + t333;
t325 = -mrSges(5,1) * t216 + mrSges(5,3) * t217 - pkin(4) * t176 + t328;
t158 = -t325 + t350 * t347 + t359 * t288 + (Ifges(4,1) + Ifges(5,2)) * t287 + (Ifges(4,5) - Ifges(5,4)) * qJDD(3) - t349 * qJD(3) + mrSges(4,2) * t233 - mrSges(4,3) * t227 - qJ(4) * t174;
t148 = mrSges(3,1) * t237 - mrSges(3,2) * t238 + Ifges(3,3) * qJDD(2) + pkin(2) * t326 + pkin(8) * t342 + t320 * t153 + t317 * t158;
t150 = mrSges(3,2) * t255 - mrSges(3,3) * t237 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t323 - pkin(8) * t166 - t153 * t317 + t158 * t320;
t334 = mrSges(2,1) * t292 - mrSges(2,2) * t293 + pkin(1) * t157 + t315 * t148 + t150 * t354 + t312 * t364;
t159 = m(2) * t293 + t161;
t156 = t315 * t165 + (t164 * t318 + t355) * t312;
t154 = m(2) * t292 + t157;
t146 = mrSges(2,2) * t309 - mrSges(2,3) * t292 + t321 * t150 - t318 * t152 + (-t156 * t312 - t157 * t315) * pkin(7);
t145 = -mrSges(2,1) * t309 + mrSges(2,3) * t293 - pkin(1) * t156 - t312 * t148 + (t150 * t318 + t364) * t315;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t314 * t146 - t311 * t145 - qJ(1) * (t154 * t314 + t159 * t311), t146, t150, t158, -t266 * t348 - t331, t169, t185; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t311 * t146 + t314 * t145 + qJ(1) * (-t154 * t311 + t159 * t314), t145, t152, t153, Ifges(5,4) * qJDD(3) - Ifges(5,2) * t287 - Ifges(5,6) * t288 - qJD(3) * t266 - t268 * t347 + t325, t167, t184; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t334, t334, t148, t363, Ifges(5,5) * qJDD(3) - Ifges(5,6) * t287 - Ifges(5,3) * t288 + qJD(3) * t267 + t268 * t348 - t330, -t328, -t333;];
m_new  = t1;
