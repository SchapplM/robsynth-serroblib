% Calculate vector of cutting torques with Newton-Euler for
% S6RPRPRR8
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
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
% Datum: 2019-05-05 19:23
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRPRR8_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR8_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR8_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR8_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR8_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR8_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR8_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR8_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR8_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 19:19:50
% EndTime: 2019-05-05 19:20:14
% DurationCPUTime: 12.40s
% Computational Cost: add. (222054->342), mult. (487016->420), div. (0->0), fcn. (328268->10), ass. (0->136)
t310 = sin(qJ(1));
t314 = cos(qJ(1));
t288 = t310 * g(1) - t314 * g(2);
t316 = qJD(1) ^ 2;
t330 = -t316 * qJ(2) + qJDD(2) - t288;
t346 = -pkin(1) - pkin(7);
t259 = t346 * qJDD(1) + t330;
t309 = sin(qJ(3));
t313 = cos(qJ(3));
t248 = t309 * g(3) + t313 * t259;
t339 = qJD(1) * qJD(3);
t337 = t309 * t339;
t283 = t313 * qJDD(1) - t337;
t221 = (-t283 - t337) * qJ(4) + (-t309 * t313 * t316 + qJDD(3)) * pkin(3) + t248;
t249 = -t313 * g(3) + t309 * t259;
t282 = -t309 * qJDD(1) - t313 * t339;
t341 = qJD(1) * t313;
t286 = qJD(3) * pkin(3) - qJ(4) * t341;
t302 = t309 ^ 2;
t222 = -t302 * t316 * pkin(3) + t282 * qJ(4) - qJD(3) * t286 + t249;
t305 = sin(pkin(10));
t306 = cos(pkin(10));
t342 = qJD(1) * t309;
t270 = -t305 * t342 + t306 * t341;
t202 = -0.2e1 * qJD(4) * t270 + t306 * t221 - t305 * t222;
t289 = -t314 * g(1) - t310 * g(2);
t331 = qJDD(1) * qJ(2) + 0.2e1 * qJD(2) * qJD(1) + t289;
t345 = mrSges(2,1) - mrSges(3,2);
t344 = Ifges(2,5) - Ifges(3,4);
t343 = -Ifges(2,6) + Ifges(3,5);
t269 = -t305 * t341 - t306 * t342;
t203 = 0.2e1 * qJD(4) * t269 + t305 * t221 + t306 * t222;
t237 = -t269 * mrSges(5,1) + t270 * mrSges(5,2);
t245 = t306 * t282 - t305 * t283;
t258 = qJD(3) * mrSges(5,1) - t270 * mrSges(5,3);
t239 = -t269 * pkin(4) - t270 * pkin(8);
t315 = qJD(3) ^ 2;
t191 = -t315 * pkin(4) + qJDD(3) * pkin(8) + t269 * t239 + t203;
t224 = -t282 * pkin(3) + qJDD(4) + t286 * t341 + (-qJ(4) * t302 + t346) * t316 + t331;
t246 = t305 * t282 + t306 * t283;
t200 = (-qJD(3) * t269 - t246) * pkin(8) + (qJD(3) * t270 - t245) * pkin(4) + t224;
t308 = sin(qJ(5));
t312 = cos(qJ(5));
t186 = -t308 * t191 + t312 * t200;
t251 = t312 * qJD(3) - t308 * t270;
t217 = t251 * qJD(5) + t308 * qJDD(3) + t312 * t246;
t244 = qJDD(5) - t245;
t252 = t308 * qJD(3) + t312 * t270;
t267 = qJD(5) - t269;
t184 = (t251 * t267 - t217) * pkin(9) + (t251 * t252 + t244) * pkin(5) + t186;
t187 = t312 * t191 + t308 * t200;
t216 = -t252 * qJD(5) + t312 * qJDD(3) - t308 * t246;
t232 = t267 * pkin(5) - t252 * pkin(9);
t250 = t251 ^ 2;
t185 = -t250 * pkin(5) + t216 * pkin(9) - t267 * t232 + t187;
t307 = sin(qJ(6));
t311 = cos(qJ(6));
t182 = t311 * t184 - t307 * t185;
t225 = t311 * t251 - t307 * t252;
t196 = t225 * qJD(6) + t307 * t216 + t311 * t217;
t226 = t307 * t251 + t311 * t252;
t208 = -t225 * mrSges(7,1) + t226 * mrSges(7,2);
t263 = qJD(6) + t267;
t209 = -t263 * mrSges(7,2) + t225 * mrSges(7,3);
t240 = qJDD(6) + t244;
t177 = m(7) * t182 + t240 * mrSges(7,1) - t196 * mrSges(7,3) - t226 * t208 + t263 * t209;
t183 = t307 * t184 + t311 * t185;
t195 = -t226 * qJD(6) + t311 * t216 - t307 * t217;
t210 = t263 * mrSges(7,1) - t226 * mrSges(7,3);
t178 = m(7) * t183 - t240 * mrSges(7,2) + t195 * mrSges(7,3) + t225 * t208 - t263 * t210;
t169 = t311 * t177 + t307 * t178;
t228 = -t251 * mrSges(6,1) + t252 * mrSges(6,2);
t230 = -t267 * mrSges(6,2) + t251 * mrSges(6,3);
t167 = m(6) * t186 + t244 * mrSges(6,1) - t217 * mrSges(6,3) - t252 * t228 + t267 * t230 + t169;
t231 = t267 * mrSges(6,1) - t252 * mrSges(6,3);
t334 = -t307 * t177 + t311 * t178;
t168 = m(6) * t187 - t244 * mrSges(6,2) + t216 * mrSges(6,3) + t251 * t228 - t267 * t231 + t334;
t335 = -t308 * t167 + t312 * t168;
t160 = m(5) * t203 - qJDD(3) * mrSges(5,2) + t245 * mrSges(5,3) - qJD(3) * t258 + t269 * t237 + t335;
t257 = -qJD(3) * mrSges(5,2) + t269 * mrSges(5,3);
t190 = -qJDD(3) * pkin(4) - t315 * pkin(8) + t270 * t239 - t202;
t188 = -t216 * pkin(5) - t250 * pkin(9) + t252 * t232 + t190;
t328 = m(7) * t188 - t195 * mrSges(7,1) + t196 * mrSges(7,2) - t225 * t209 + t226 * t210;
t321 = -m(6) * t190 + t216 * mrSges(6,1) - t217 * mrSges(6,2) + t251 * t230 - t252 * t231 - t328;
t173 = m(5) * t202 + qJDD(3) * mrSges(5,1) - t246 * mrSges(5,3) + qJD(3) * t257 - t270 * t237 + t321;
t150 = t305 * t160 + t306 * t173;
t162 = t312 * t167 + t308 * t168;
t281 = (mrSges(4,1) * t309 + mrSges(4,2) * t313) * qJD(1);
t285 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t342;
t147 = m(4) * t248 + qJDD(3) * mrSges(4,1) - t283 * mrSges(4,3) + qJD(3) * t285 - t281 * t341 + t150;
t287 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t341;
t336 = t306 * t160 - t305 * t173;
t148 = m(4) * t249 - qJDD(3) * mrSges(4,2) + t282 * mrSges(4,3) - qJD(3) * t287 - t281 * t342 + t336;
t144 = -t309 * t147 + t313 * t148;
t143 = t313 * t147 + t309 * t148;
t268 = -qJDD(1) * pkin(1) + t330;
t329 = -m(3) * t268 + t316 * mrSges(3,3) - t143;
t327 = m(5) * t224 - t245 * mrSges(5,1) + t246 * mrSges(5,2) - t269 * t257 + t270 * t258 + t162;
t205 = Ifges(7,4) * t226 + Ifges(7,2) * t225 + Ifges(7,6) * t263;
t206 = Ifges(7,1) * t226 + Ifges(7,4) * t225 + Ifges(7,5) * t263;
t326 = -mrSges(7,1) * t182 + mrSges(7,2) * t183 - Ifges(7,5) * t196 - Ifges(7,6) * t195 - Ifges(7,3) * t240 - t226 * t205 + t225 * t206;
t204 = Ifges(7,5) * t226 + Ifges(7,6) * t225 + Ifges(7,3) * t263;
t170 = -mrSges(7,1) * t188 + mrSges(7,3) * t183 + Ifges(7,4) * t196 + Ifges(7,2) * t195 + Ifges(7,6) * t240 - t226 * t204 + t263 * t206;
t171 = mrSges(7,2) * t188 - mrSges(7,3) * t182 + Ifges(7,1) * t196 + Ifges(7,4) * t195 + Ifges(7,5) * t240 + t225 * t204 - t263 * t205;
t211 = Ifges(6,5) * t252 + Ifges(6,6) * t251 + Ifges(6,3) * t267;
t213 = Ifges(6,1) * t252 + Ifges(6,4) * t251 + Ifges(6,5) * t267;
t152 = -mrSges(6,1) * t190 + mrSges(6,3) * t187 + Ifges(6,4) * t217 + Ifges(6,2) * t216 + Ifges(6,6) * t244 - pkin(5) * t328 + pkin(9) * t334 + t311 * t170 + t307 * t171 - t252 * t211 + t267 * t213;
t212 = Ifges(6,4) * t252 + Ifges(6,2) * t251 + Ifges(6,6) * t267;
t154 = mrSges(6,2) * t190 - mrSges(6,3) * t186 + Ifges(6,1) * t217 + Ifges(6,4) * t216 + Ifges(6,5) * t244 - pkin(9) * t169 - t307 * t170 + t311 * t171 + t251 * t211 - t267 * t212;
t233 = Ifges(5,5) * t270 + Ifges(5,6) * t269 + Ifges(5,3) * qJD(3);
t234 = Ifges(5,4) * t270 + Ifges(5,2) * t269 + Ifges(5,6) * qJD(3);
t139 = mrSges(5,2) * t224 - mrSges(5,3) * t202 + Ifges(5,1) * t246 + Ifges(5,4) * t245 + Ifges(5,5) * qJDD(3) - pkin(8) * t162 - qJD(3) * t234 - t308 * t152 + t312 * t154 + t269 * t233;
t235 = Ifges(5,1) * t270 + Ifges(5,4) * t269 + Ifges(5,5) * qJD(3);
t317 = mrSges(6,1) * t186 - mrSges(6,2) * t187 + Ifges(6,5) * t217 + Ifges(6,6) * t216 + Ifges(6,3) * t244 + pkin(5) * t169 + t252 * t212 - t251 * t213 - t326;
t145 = -mrSges(5,1) * t224 + mrSges(5,3) * t203 + Ifges(5,4) * t246 + Ifges(5,2) * t245 + Ifges(5,6) * qJDD(3) - pkin(4) * t162 + qJD(3) * t235 - t270 * t233 - t317;
t256 = t346 * t316 + t331;
t271 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t313 - Ifges(4,6) * t309) * qJD(1);
t273 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t313 - Ifges(4,4) * t309) * qJD(1);
t136 = -mrSges(4,1) * t256 + mrSges(4,3) * t249 + Ifges(4,4) * t283 + Ifges(4,2) * t282 + Ifges(4,6) * qJDD(3) - pkin(3) * t327 + qJ(4) * t336 + qJD(3) * t273 + t305 * t139 + t306 * t145 - t271 * t341;
t272 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t313 - Ifges(4,2) * t309) * qJD(1);
t138 = mrSges(4,2) * t256 - mrSges(4,3) * t248 + Ifges(4,1) * t283 + Ifges(4,4) * t282 + Ifges(4,5) * qJDD(3) - qJ(4) * t150 - qJD(3) * t272 + t306 * t139 - t305 * t145 - t271 * t342;
t262 = t316 * pkin(1) - t331;
t325 = mrSges(3,2) * t268 - mrSges(3,3) * t262 + Ifges(3,1) * qJDD(1) - pkin(7) * t143 - t309 * t136 + t313 * t138;
t157 = -m(4) * t256 + t282 * mrSges(4,1) - t283 * mrSges(4,2) - t285 * t342 - t287 * t341 - t327;
t324 = -mrSges(3,1) * t262 - pkin(2) * t157 - pkin(7) * t144 - t313 * t136 - t309 * t138;
t323 = -mrSges(5,1) * t202 + mrSges(5,2) * t203 - Ifges(5,5) * t246 - Ifges(5,6) * t245 - Ifges(5,3) * qJDD(3) - pkin(4) * t321 - pkin(8) * t335 - t312 * t152 - t308 * t154 - t270 * t234 + t269 * t235;
t320 = -m(3) * t262 + t316 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t157;
t322 = -mrSges(2,2) * t289 + pkin(1) * (-qJDD(1) * mrSges(3,2) + t329) + qJ(2) * t320 + mrSges(2,1) * t288 + Ifges(2,3) * qJDD(1) + t325;
t319 = -mrSges(4,1) * t248 + mrSges(4,2) * t249 - Ifges(4,5) * t283 - Ifges(4,6) * t282 - Ifges(4,3) * qJDD(3) - pkin(3) * t150 - t272 * t341 - t273 * t342 + t323;
t318 = -mrSges(3,1) * t268 - pkin(2) * t143 + t319;
t155 = m(2) * t289 - t316 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t320;
t142 = -m(3) * g(3) + t144;
t140 = m(2) * t288 - t316 * mrSges(2,2) + t345 * qJDD(1) + t329;
t135 = -t318 - qJ(2) * t142 + t343 * t316 + (-mrSges(2,2) + mrSges(3,3)) * g(3) + t344 * qJDD(1) - mrSges(2,3) * t288;
t134 = mrSges(2,3) * t289 - pkin(1) * t142 + t345 * g(3) - t343 * qJDD(1) + t344 * t316 + t324;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t314 * t135 - t310 * t134 - pkin(6) * (t314 * t140 + t310 * t155), t135, t325, t138, t139, t154, t171; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t310 * t135 + t314 * t134 + pkin(6) * (-t310 * t140 + t314 * t155), t134, -mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) - t316 * Ifges(3,5) + t318, t136, t145, t152, t170; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t322, t322, mrSges(3,2) * g(3) + t316 * Ifges(3,4) + Ifges(3,5) * qJDD(1) - t324, -t319, -t323, t317, -t326;];
m_new  = t1;
