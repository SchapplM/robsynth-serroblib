% Calculate vector of cutting torques with Newton-Euler for
% S6RPRPPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4,theta5]';
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
% Datum: 2019-05-05 16:26
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRPPR1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR1_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR1_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR1_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPPR1_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR1_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR1_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR1_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 16:23:55
% EndTime: 2019-05-05 16:24:18
% DurationCPUTime: 18.07s
% Computational Cost: add. (324224->341), mult. (724747->433), div. (0->0), fcn. (487345->12), ass. (0->136)
t339 = -2 * qJD(4);
t308 = sin(qJ(1));
t311 = cos(qJ(1));
t290 = t308 * g(1) - t311 * g(2);
t281 = qJDD(1) * pkin(1) + t290;
t291 = -t311 * g(1) - t308 * g(2);
t313 = qJD(1) ^ 2;
t283 = -t313 * pkin(1) + t291;
t303 = sin(pkin(9));
t305 = cos(pkin(9));
t256 = t303 * t281 + t305 * t283;
t245 = -t313 * pkin(2) + qJDD(1) * pkin(7) + t256;
t300 = -g(3) + qJDD(2);
t307 = sin(qJ(3));
t310 = cos(qJ(3));
t232 = -t307 * t245 + t310 * t300;
t333 = qJD(1) * qJD(3);
t332 = t310 * t333;
t284 = t307 * qJDD(1) + t332;
t218 = (-t284 + t332) * qJ(4) + (t307 * t310 * t313 + qJDD(3)) * pkin(3) + t232;
t233 = t310 * t245 + t307 * t300;
t285 = t310 * qJDD(1) - t307 * t333;
t336 = qJD(1) * t307;
t287 = qJD(3) * pkin(3) - qJ(4) * t336;
t299 = t310 ^ 2;
t221 = -t299 * t313 * pkin(3) + t285 * qJ(4) - qJD(3) * t287 + t233;
t302 = sin(pkin(10));
t337 = cos(pkin(10));
t270 = (t302 * t310 + t307 * t337) * qJD(1);
t200 = t218 * t337 - t302 * t221 + t270 * t339;
t335 = qJD(1) * t310;
t269 = t302 * t336 - t335 * t337;
t201 = t302 * t218 + t221 * t337 + t269 * t339;
t248 = t269 * mrSges(5,1) + t270 * mrSges(5,2);
t257 = t302 * t284 - t285 * t337;
t265 = qJD(3) * mrSges(5,1) - t270 * mrSges(5,3);
t247 = t269 * pkin(4) - t270 * qJ(5);
t312 = qJD(3) ^ 2;
t197 = -t312 * pkin(4) + qJDD(3) * qJ(5) - t269 * t247 + t201;
t255 = t305 * t281 - t303 * t283;
t324 = -qJDD(1) * pkin(2) - t255;
t222 = -t285 * pkin(3) + qJDD(4) + t287 * t336 + (-qJ(4) * t299 - pkin(7)) * t313 + t324;
t258 = t284 * t337 + t302 * t285;
t204 = (qJD(3) * t269 - t258) * qJ(5) + (qJD(3) * t270 + t257) * pkin(4) + t222;
t301 = sin(pkin(11));
t304 = cos(pkin(11));
t263 = t301 * qJD(3) + t304 * t270;
t192 = -0.2e1 * qJD(5) * t263 - t301 * t197 + t304 * t204;
t240 = t301 * qJDD(3) + t304 * t258;
t262 = t304 * qJD(3) - t301 * t270;
t190 = (t262 * t269 - t240) * pkin(8) + (t262 * t263 + t257) * pkin(5) + t192;
t193 = 0.2e1 * qJD(5) * t262 + t304 * t197 + t301 * t204;
t236 = t269 * pkin(5) - t263 * pkin(8);
t239 = t304 * qJDD(3) - t301 * t258;
t261 = t262 ^ 2;
t191 = -t261 * pkin(5) + t239 * pkin(8) - t269 * t236 + t193;
t306 = sin(qJ(6));
t309 = cos(qJ(6));
t188 = t309 * t190 - t306 * t191;
t227 = t309 * t262 - t306 * t263;
t209 = t227 * qJD(6) + t306 * t239 + t309 * t240;
t228 = t306 * t262 + t309 * t263;
t214 = -t227 * mrSges(7,1) + t228 * mrSges(7,2);
t268 = qJD(6) + t269;
t219 = -t268 * mrSges(7,2) + t227 * mrSges(7,3);
t254 = qJDD(6) + t257;
t183 = m(7) * t188 + t254 * mrSges(7,1) - t209 * mrSges(7,3) - t228 * t214 + t268 * t219;
t189 = t306 * t190 + t309 * t191;
t208 = -t228 * qJD(6) + t309 * t239 - t306 * t240;
t220 = t268 * mrSges(7,1) - t228 * mrSges(7,3);
t184 = m(7) * t189 - t254 * mrSges(7,2) + t208 * mrSges(7,3) + t227 * t214 - t268 * t220;
t175 = t309 * t183 + t306 * t184;
t230 = -t262 * mrSges(6,1) + t263 * mrSges(6,2);
t234 = -t269 * mrSges(6,2) + t262 * mrSges(6,3);
t173 = m(6) * t192 + t257 * mrSges(6,1) - t240 * mrSges(6,3) - t263 * t230 + t269 * t234 + t175;
t235 = t269 * mrSges(6,1) - t263 * mrSges(6,3);
t327 = -t306 * t183 + t309 * t184;
t174 = m(6) * t193 - t257 * mrSges(6,2) + t239 * mrSges(6,3) + t262 * t230 - t269 * t235 + t327;
t328 = -t301 * t173 + t304 * t174;
t166 = m(5) * t201 - qJDD(3) * mrSges(5,2) - t257 * mrSges(5,3) - qJD(3) * t265 - t269 * t248 + t328;
t264 = -qJD(3) * mrSges(5,2) - t269 * mrSges(5,3);
t196 = -qJDD(3) * pkin(4) - t312 * qJ(5) + t270 * t247 + qJDD(5) - t200;
t194 = -t239 * pkin(5) - t261 * pkin(8) + t263 * t236 + t196;
t322 = m(7) * t194 - t208 * mrSges(7,1) + t209 * mrSges(7,2) - t227 * t219 + t228 * t220;
t317 = -m(6) * t196 + t239 * mrSges(6,1) - t240 * mrSges(6,2) + t262 * t234 - t263 * t235 - t322;
t179 = m(5) * t200 + qJDD(3) * mrSges(5,1) - t258 * mrSges(5,3) + qJD(3) * t264 - t270 * t248 + t317;
t156 = t302 * t166 + t179 * t337;
t275 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t307 + Ifges(4,2) * t310) * qJD(1);
t276 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t307 + Ifges(4,4) * t310) * qJD(1);
t210 = Ifges(7,5) * t228 + Ifges(7,6) * t227 + Ifges(7,3) * t268;
t212 = Ifges(7,1) * t228 + Ifges(7,4) * t227 + Ifges(7,5) * t268;
t176 = -mrSges(7,1) * t194 + mrSges(7,3) * t189 + Ifges(7,4) * t209 + Ifges(7,2) * t208 + Ifges(7,6) * t254 - t228 * t210 + t268 * t212;
t211 = Ifges(7,4) * t228 + Ifges(7,2) * t227 + Ifges(7,6) * t268;
t177 = mrSges(7,2) * t194 - mrSges(7,3) * t188 + Ifges(7,1) * t209 + Ifges(7,4) * t208 + Ifges(7,5) * t254 + t227 * t210 - t268 * t211;
t223 = Ifges(6,5) * t263 + Ifges(6,6) * t262 + Ifges(6,3) * t269;
t225 = Ifges(6,1) * t263 + Ifges(6,4) * t262 + Ifges(6,5) * t269;
t158 = -mrSges(6,1) * t196 + mrSges(6,3) * t193 + Ifges(6,4) * t240 + Ifges(6,2) * t239 + Ifges(6,6) * t257 - pkin(5) * t322 + pkin(8) * t327 + t309 * t176 + t306 * t177 - t263 * t223 + t269 * t225;
t224 = Ifges(6,4) * t263 + Ifges(6,2) * t262 + Ifges(6,6) * t269;
t160 = mrSges(6,2) * t196 - mrSges(6,3) * t192 + Ifges(6,1) * t240 + Ifges(6,4) * t239 + Ifges(6,5) * t257 - pkin(8) * t175 - t306 * t176 + t309 * t177 + t262 * t223 - t269 * t224;
t242 = Ifges(5,4) * t270 - Ifges(5,2) * t269 + Ifges(5,6) * qJD(3);
t243 = Ifges(5,1) * t270 - Ifges(5,4) * t269 + Ifges(5,5) * qJD(3);
t318 = -mrSges(5,1) * t200 + mrSges(5,2) * t201 - Ifges(5,5) * t258 + Ifges(5,6) * t257 - Ifges(5,3) * qJDD(3) - pkin(4) * t317 - qJ(5) * t328 - t304 * t158 - t301 * t160 - t270 * t242 - t269 * t243;
t338 = mrSges(4,1) * t232 - mrSges(4,2) * t233 + Ifges(4,5) * t284 + Ifges(4,6) * t285 + Ifges(4,3) * qJDD(3) + pkin(3) * t156 + (t307 * t275 - t310 * t276) * qJD(1) - t318;
t282 = (-mrSges(4,1) * t310 + mrSges(4,2) * t307) * qJD(1);
t289 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t335;
t154 = m(4) * t232 + qJDD(3) * mrSges(4,1) - t284 * mrSges(4,3) + qJD(3) * t289 - t282 * t336 + t156;
t288 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t336;
t329 = t166 * t337 - t302 * t179;
t155 = m(4) * t233 - qJDD(3) * mrSges(4,2) + t285 * mrSges(4,3) - qJD(3) * t288 + t282 * t335 + t329;
t330 = -t307 * t154 + t310 * t155;
t147 = m(3) * t256 - t313 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t330;
t244 = -t313 * pkin(7) + t324;
t168 = t304 * t173 + t301 * t174;
t320 = m(5) * t222 + t257 * mrSges(5,1) + t258 * mrSges(5,2) + t269 * t264 + t270 * t265 + t168;
t316 = -m(4) * t244 + t285 * mrSges(4,1) - t284 * mrSges(4,2) - t288 * t336 + t289 * t335 - t320;
t162 = m(3) * t255 + qJDD(1) * mrSges(3,1) - t313 * mrSges(3,2) + t316;
t143 = t303 * t147 + t305 * t162;
t149 = t310 * t154 + t307 * t155;
t331 = t305 * t147 - t303 * t162;
t241 = Ifges(5,5) * t270 - Ifges(5,6) * t269 + Ifges(5,3) * qJD(3);
t144 = mrSges(5,2) * t222 - mrSges(5,3) * t200 + Ifges(5,1) * t258 - Ifges(5,4) * t257 + Ifges(5,5) * qJDD(3) - qJ(5) * t168 - qJD(3) * t242 - t301 * t158 + t304 * t160 - t269 * t241;
t321 = -mrSges(7,1) * t188 + mrSges(7,2) * t189 - Ifges(7,5) * t209 - Ifges(7,6) * t208 - Ifges(7,3) * t254 - t228 * t211 + t227 * t212;
t315 = -mrSges(6,1) * t192 + mrSges(6,2) * t193 - Ifges(6,5) * t240 - Ifges(6,6) * t239 - pkin(5) * t175 - t263 * t224 + t262 * t225 + t321;
t150 = Ifges(5,6) * qJDD(3) + t315 + (-Ifges(5,2) - Ifges(6,3)) * t257 - t270 * t241 + Ifges(5,4) * t258 + qJD(3) * t243 - mrSges(5,1) * t222 + mrSges(5,3) * t201 - pkin(4) * t168;
t274 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t307 + Ifges(4,6) * t310) * qJD(1);
t136 = -mrSges(4,1) * t244 + mrSges(4,3) * t233 + Ifges(4,4) * t284 + Ifges(4,2) * t285 + Ifges(4,6) * qJDD(3) - pkin(3) * t320 + qJ(4) * t329 + qJD(3) * t276 + t302 * t144 + t150 * t337 - t274 * t336;
t139 = mrSges(4,2) * t244 - mrSges(4,3) * t232 + Ifges(4,1) * t284 + Ifges(4,4) * t285 + Ifges(4,5) * qJDD(3) - qJ(4) * t156 - qJD(3) * t275 + t144 * t337 - t302 * t150 + t274 * t335;
t323 = mrSges(3,1) * t255 - mrSges(3,2) * t256 + Ifges(3,3) * qJDD(1) + pkin(2) * t316 + pkin(7) * t330 + t310 * t136 + t307 * t139;
t319 = mrSges(2,1) * t290 - mrSges(2,2) * t291 + Ifges(2,3) * qJDD(1) + pkin(1) * t143 + t323;
t141 = m(2) * t291 - t313 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t331;
t140 = m(2) * t290 + qJDD(1) * mrSges(2,1) - t313 * mrSges(2,2) + t143;
t137 = -mrSges(3,1) * t300 + mrSges(3,3) * t256 + t313 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t149 - t338;
t134 = mrSges(3,2) * t300 - mrSges(3,3) * t255 + Ifges(3,5) * qJDD(1) - t313 * Ifges(3,6) - pkin(7) * t149 - t307 * t136 + t310 * t139;
t133 = -mrSges(2,2) * g(3) - mrSges(2,3) * t290 + Ifges(2,5) * qJDD(1) - t313 * Ifges(2,6) - qJ(2) * t143 + t305 * t134 - t303 * t137;
t132 = Ifges(2,6) * qJDD(1) + t313 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t291 + t303 * t134 + t305 * t137 - pkin(1) * (m(3) * t300 + t149) + qJ(2) * t331;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t311 * t133 - t308 * t132 - pkin(6) * (t311 * t140 + t308 * t141), t133, t134, t139, t144, t160, t177; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t308 * t133 + t311 * t132 + pkin(6) * (-t308 * t140 + t311 * t141), t132, t137, t136, t150, t158, t176; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t319, t319, t323, t338, -t318, Ifges(6,3) * t257 - t315, -t321;];
m_new  = t1;
