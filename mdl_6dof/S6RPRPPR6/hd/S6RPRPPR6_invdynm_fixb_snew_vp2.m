% Calculate vector of cutting torques with Newton-Euler for
% S6RPRPPR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4,theta5]';
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
% Datum: 2019-05-05 17:10
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRPPR6_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR6_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR6_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR6_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR6_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR6_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR6_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR6_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR6_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:07:28
% EndTime: 2019-05-05 17:07:50
% DurationCPUTime: 12.32s
% Computational Cost: add. (208339->342), mult. (471795->422), div. (0->0), fcn. (313670->10), ass. (0->135)
t344 = -2 * qJD(4);
t308 = sin(qJ(1));
t311 = cos(qJ(1));
t287 = t308 * g(1) - t311 * g(2);
t313 = qJD(1) ^ 2;
t327 = -t313 * qJ(2) + qJDD(2) - t287;
t343 = -pkin(1) - pkin(7);
t259 = t343 * qJDD(1) + t327;
t307 = sin(qJ(3));
t310 = cos(qJ(3));
t246 = t307 * g(3) + t310 * t259;
t336 = qJD(1) * qJD(3);
t334 = t307 * t336;
t282 = t310 * qJDD(1) - t334;
t217 = (-t282 - t334) * qJ(4) + (-t307 * t310 * t313 + qJDD(3)) * pkin(3) + t246;
t247 = -t310 * g(3) + t307 * t259;
t281 = -t307 * qJDD(1) - t310 * t336;
t338 = qJD(1) * t310;
t285 = qJD(3) * pkin(3) - qJ(4) * t338;
t299 = t307 ^ 2;
t218 = -t299 * t313 * pkin(3) + t281 * qJ(4) - qJD(3) * t285 + t247;
t303 = sin(pkin(9));
t305 = cos(pkin(9));
t339 = qJD(1) * t307;
t269 = -t303 * t339 + t305 * t338;
t201 = t305 * t217 - t303 * t218 + t269 * t344;
t288 = -t311 * g(1) - t308 * g(2);
t328 = qJDD(1) * qJ(2) + 0.2e1 * qJD(2) * qJD(1) + t288;
t342 = mrSges(2,1) - mrSges(3,2);
t341 = Ifges(2,5) - Ifges(3,4);
t340 = -Ifges(2,6) + Ifges(3,5);
t268 = (t303 * t310 + t305 * t307) * qJD(1);
t202 = t303 * t217 + t305 * t218 + t268 * t344;
t237 = t268 * mrSges(5,1) + t269 * mrSges(5,2);
t243 = -t305 * t281 + t303 * t282;
t258 = qJD(3) * mrSges(5,1) - t269 * mrSges(5,3);
t236 = t268 * pkin(4) - t269 * qJ(5);
t312 = qJD(3) ^ 2;
t190 = -t312 * pkin(4) + qJDD(3) * qJ(5) - t268 * t236 + t202;
t222 = -t281 * pkin(3) + qJDD(4) + t285 * t338 + (-qJ(4) * t299 + t343) * t313 + t328;
t244 = t303 * t281 + t305 * t282;
t194 = (qJD(3) * t268 - t244) * qJ(5) + (qJD(3) * t269 + t243) * pkin(4) + t222;
t302 = sin(pkin(10));
t304 = cos(pkin(10));
t252 = t302 * qJD(3) + t304 * t269;
t185 = -0.2e1 * qJD(5) * t252 - t302 * t190 + t304 * t194;
t231 = t302 * qJDD(3) + t304 * t244;
t251 = t304 * qJD(3) - t302 * t269;
t183 = (t251 * t268 - t231) * pkin(8) + (t251 * t252 + t243) * pkin(5) + t185;
t186 = 0.2e1 * qJD(5) * t251 + t304 * t190 + t302 * t194;
t228 = t268 * pkin(5) - t252 * pkin(8);
t230 = t304 * qJDD(3) - t302 * t244;
t250 = t251 ^ 2;
t184 = -t250 * pkin(5) + t230 * pkin(8) - t268 * t228 + t186;
t306 = sin(qJ(6));
t309 = cos(qJ(6));
t181 = t309 * t183 - t306 * t184;
t220 = t309 * t251 - t306 * t252;
t199 = t220 * qJD(6) + t306 * t230 + t309 * t231;
t221 = t306 * t251 + t309 * t252;
t207 = -t220 * mrSges(7,1) + t221 * mrSges(7,2);
t266 = qJD(6) + t268;
t208 = -t266 * mrSges(7,2) + t220 * mrSges(7,3);
t242 = qJDD(6) + t243;
t176 = m(7) * t181 + t242 * mrSges(7,1) - t199 * mrSges(7,3) - t221 * t207 + t266 * t208;
t182 = t306 * t183 + t309 * t184;
t198 = -t221 * qJD(6) + t309 * t230 - t306 * t231;
t209 = t266 * mrSges(7,1) - t221 * mrSges(7,3);
t177 = m(7) * t182 - t242 * mrSges(7,2) + t198 * mrSges(7,3) + t220 * t207 - t266 * t209;
t168 = t309 * t176 + t306 * t177;
t224 = -t251 * mrSges(6,1) + t252 * mrSges(6,2);
t226 = -t268 * mrSges(6,2) + t251 * mrSges(6,3);
t166 = m(6) * t185 + t243 * mrSges(6,1) - t231 * mrSges(6,3) - t252 * t224 + t268 * t226 + t168;
t227 = t268 * mrSges(6,1) - t252 * mrSges(6,3);
t331 = -t306 * t176 + t309 * t177;
t167 = m(6) * t186 - t243 * mrSges(6,2) + t230 * mrSges(6,3) + t251 * t224 - t268 * t227 + t331;
t332 = -t302 * t166 + t304 * t167;
t159 = m(5) * t202 - qJDD(3) * mrSges(5,2) - t243 * mrSges(5,3) - qJD(3) * t258 - t268 * t237 + t332;
t257 = -qJD(3) * mrSges(5,2) - t268 * mrSges(5,3);
t189 = -qJDD(3) * pkin(4) - t312 * qJ(5) + t269 * t236 + qJDD(5) - t201;
t187 = -t230 * pkin(5) - t250 * pkin(8) + t252 * t228 + t189;
t325 = m(7) * t187 - t198 * mrSges(7,1) + t199 * mrSges(7,2) - t220 * t208 + t221 * t209;
t318 = -m(6) * t189 + t230 * mrSges(6,1) - t231 * mrSges(6,2) + t251 * t226 - t252 * t227 - t325;
t172 = m(5) * t201 + qJDD(3) * mrSges(5,1) - t244 * mrSges(5,3) + qJD(3) * t257 - t269 * t237 + t318;
t149 = t303 * t159 + t305 * t172;
t161 = t304 * t166 + t302 * t167;
t280 = (mrSges(4,1) * t307 + mrSges(4,2) * t310) * qJD(1);
t284 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t339;
t146 = m(4) * t246 + qJDD(3) * mrSges(4,1) - t282 * mrSges(4,3) + qJD(3) * t284 - t280 * t338 + t149;
t286 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t338;
t333 = t305 * t159 - t303 * t172;
t147 = m(4) * t247 - qJDD(3) * mrSges(4,2) + t281 * mrSges(4,3) - qJD(3) * t286 - t280 * t339 + t333;
t143 = -t307 * t146 + t310 * t147;
t142 = t310 * t146 + t307 * t147;
t267 = -qJDD(1) * pkin(1) + t327;
t326 = -m(3) * t267 + t313 * mrSges(3,3) - t142;
t324 = m(5) * t222 + t243 * mrSges(5,1) + t244 * mrSges(5,2) + t268 * t257 + t269 * t258 + t161;
t204 = Ifges(7,4) * t221 + Ifges(7,2) * t220 + Ifges(7,6) * t266;
t205 = Ifges(7,1) * t221 + Ifges(7,4) * t220 + Ifges(7,5) * t266;
t323 = -mrSges(7,1) * t181 + mrSges(7,2) * t182 - Ifges(7,5) * t199 - Ifges(7,6) * t198 - Ifges(7,3) * t242 - t221 * t204 + t220 * t205;
t203 = Ifges(7,5) * t221 + Ifges(7,6) * t220 + Ifges(7,3) * t266;
t169 = -mrSges(7,1) * t187 + mrSges(7,3) * t182 + Ifges(7,4) * t199 + Ifges(7,2) * t198 + Ifges(7,6) * t242 - t221 * t203 + t266 * t205;
t170 = mrSges(7,2) * t187 - mrSges(7,3) * t181 + Ifges(7,1) * t199 + Ifges(7,4) * t198 + Ifges(7,5) * t242 + t220 * t203 - t266 * t204;
t210 = Ifges(6,5) * t252 + Ifges(6,6) * t251 + Ifges(6,3) * t268;
t212 = Ifges(6,1) * t252 + Ifges(6,4) * t251 + Ifges(6,5) * t268;
t151 = -mrSges(6,1) * t189 + mrSges(6,3) * t186 + Ifges(6,4) * t231 + Ifges(6,2) * t230 + Ifges(6,6) * t243 - pkin(5) * t325 + pkin(8) * t331 + t309 * t169 + t306 * t170 - t252 * t210 + t268 * t212;
t211 = Ifges(6,4) * t252 + Ifges(6,2) * t251 + Ifges(6,6) * t268;
t153 = mrSges(6,2) * t189 - mrSges(6,3) * t185 + Ifges(6,1) * t231 + Ifges(6,4) * t230 + Ifges(6,5) * t243 - pkin(8) * t168 - t306 * t169 + t309 * t170 + t251 * t210 - t268 * t211;
t232 = Ifges(5,5) * t269 - Ifges(5,6) * t268 + Ifges(5,3) * qJD(3);
t233 = Ifges(5,4) * t269 - Ifges(5,2) * t268 + Ifges(5,6) * qJD(3);
t138 = mrSges(5,2) * t222 - mrSges(5,3) * t201 + Ifges(5,1) * t244 - Ifges(5,4) * t243 + Ifges(5,5) * qJDD(3) - qJ(5) * t161 - qJD(3) * t233 - t302 * t151 + t304 * t153 - t268 * t232;
t234 = Ifges(5,1) * t269 - Ifges(5,4) * t268 + Ifges(5,5) * qJD(3);
t315 = -mrSges(6,1) * t185 + mrSges(6,2) * t186 - Ifges(6,5) * t231 - Ifges(6,6) * t230 - pkin(5) * t168 - t252 * t211 + t251 * t212 + t323;
t144 = t315 + (-Ifges(5,2) - Ifges(6,3)) * t243 - t269 * t232 + Ifges(5,4) * t244 + qJD(3) * t234 - mrSges(5,1) * t222 + mrSges(5,3) * t202 + Ifges(5,6) * qJDD(3) - pkin(4) * t161;
t256 = t343 * t313 + t328;
t270 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t310 - Ifges(4,6) * t307) * qJD(1);
t272 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t310 - Ifges(4,4) * t307) * qJD(1);
t135 = -mrSges(4,1) * t256 + mrSges(4,3) * t247 + Ifges(4,4) * t282 + Ifges(4,2) * t281 + Ifges(4,6) * qJDD(3) - pkin(3) * t324 + qJ(4) * t333 + qJD(3) * t272 + t303 * t138 + t305 * t144 - t270 * t338;
t271 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t310 - Ifges(4,2) * t307) * qJD(1);
t137 = mrSges(4,2) * t256 - mrSges(4,3) * t246 + Ifges(4,1) * t282 + Ifges(4,4) * t281 + Ifges(4,5) * qJDD(3) - qJ(4) * t149 - qJD(3) * t271 + t305 * t138 - t303 * t144 - t270 * t339;
t262 = t313 * pkin(1) - t328;
t322 = mrSges(3,2) * t267 - mrSges(3,3) * t262 + Ifges(3,1) * qJDD(1) - pkin(7) * t142 - t307 * t135 + t310 * t137;
t156 = -m(4) * t256 + t281 * mrSges(4,1) - t282 * mrSges(4,2) - t284 * t339 - t286 * t338 - t324;
t321 = -mrSges(3,1) * t262 - pkin(2) * t156 - pkin(7) * t143 - t310 * t135 - t307 * t137;
t320 = -mrSges(5,1) * t201 + mrSges(5,2) * t202 - Ifges(5,5) * t244 + Ifges(5,6) * t243 - Ifges(5,3) * qJDD(3) - pkin(4) * t318 - qJ(5) * t332 - t304 * t151 - t302 * t153 - t269 * t233 - t268 * t234;
t317 = -m(3) * t262 + t313 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t156;
t319 = -mrSges(2,2) * t288 + pkin(1) * (-qJDD(1) * mrSges(3,2) + t326) + qJ(2) * t317 + mrSges(2,1) * t287 + Ifges(2,3) * qJDD(1) + t322;
t316 = -mrSges(4,1) * t246 + mrSges(4,2) * t247 - Ifges(4,5) * t282 - Ifges(4,6) * t281 - Ifges(4,3) * qJDD(3) - pkin(3) * t149 - t271 * t338 - t272 * t339 + t320;
t314 = -mrSges(3,1) * t267 - pkin(2) * t142 + t316;
t154 = m(2) * t288 - t313 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t317;
t141 = -m(3) * g(3) + t143;
t139 = m(2) * t287 - t313 * mrSges(2,2) + t342 * qJDD(1) + t326;
t134 = -t314 - mrSges(2,3) * t287 - qJ(2) * t141 + (-mrSges(2,2) + mrSges(3,3)) * g(3) + t341 * qJDD(1) + t340 * t313;
t133 = mrSges(2,3) * t288 - pkin(1) * t141 + t342 * g(3) - t340 * qJDD(1) + t341 * t313 + t321;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t311 * t134 - t308 * t133 - pkin(6) * (t311 * t139 + t308 * t154), t134, t322, t137, t138, t153, t170; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t308 * t134 + t311 * t133 + pkin(6) * (-t308 * t139 + t311 * t154), t133, -mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) - t313 * Ifges(3,5) + t314, t135, t144, t151, t169; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t319, t319, mrSges(3,2) * g(3) + t313 * Ifges(3,4) + Ifges(3,5) * qJDD(1) - t321, -t316, -t320, Ifges(6,3) * t243 - t315, -t323;];
m_new  = t1;
