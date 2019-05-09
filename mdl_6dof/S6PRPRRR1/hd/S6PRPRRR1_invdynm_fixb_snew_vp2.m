% Calculate vector of cutting torques with Newton-Euler for
% S6PRPRRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-05-05 00:09
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRPRRR1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR1_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR1_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR1_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_invdynm_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR1_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR1_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR1_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 00:05:13
% EndTime: 2019-05-05 00:05:36
% DurationCPUTime: 16.18s
% Computational Cost: add. (324463->298), mult. (610726->387), div. (0->0), fcn. (435782->14), ass. (0->133)
t284 = sin(qJ(5));
t285 = sin(qJ(4));
t288 = cos(qJ(5));
t289 = cos(qJ(4));
t252 = (t284 * t285 - t288 * t289) * qJD(2);
t278 = sin(pkin(11));
t281 = cos(pkin(11));
t262 = g(1) * t278 - g(2) * t281;
t263 = -g(1) * t281 - g(2) * t278;
t276 = -g(3) + qJDD(1);
t286 = sin(qJ(2));
t282 = cos(pkin(6));
t290 = cos(qJ(2));
t312 = t282 * t290;
t279 = sin(pkin(6));
t314 = t279 * t290;
t232 = t262 * t312 - t263 * t286 + t276 * t314;
t227 = qJDD(2) * pkin(2) + t232;
t313 = t282 * t286;
t315 = t279 * t286;
t233 = t262 * t313 + t290 * t263 + t276 * t315;
t291 = qJD(2) ^ 2;
t228 = -pkin(2) * t291 + t233;
t277 = sin(pkin(12));
t280 = cos(pkin(12));
t206 = t277 * t227 + t280 * t228;
t203 = -pkin(3) * t291 + qJDD(2) * pkin(8) + t206;
t245 = -t262 * t279 + t282 * t276;
t244 = qJDD(3) + t245;
t199 = -t285 * t203 + t289 * t244;
t309 = qJD(2) * qJD(4);
t307 = t289 * t309;
t259 = qJDD(2) * t285 + t307;
t196 = (-t259 + t307) * pkin(9) + (t285 * t289 * t291 + qJDD(4)) * pkin(4) + t199;
t200 = t289 * t203 + t285 * t244;
t260 = qJDD(2) * t289 - t285 * t309;
t311 = qJD(2) * t285;
t268 = qJD(4) * pkin(4) - pkin(9) * t311;
t275 = t289 ^ 2;
t197 = -pkin(4) * t275 * t291 + pkin(9) * t260 - qJD(4) * t268 + t200;
t192 = t284 * t196 + t288 * t197;
t253 = (t284 * t289 + t285 * t288) * qJD(2);
t220 = -qJD(5) * t253 - t259 * t284 + t260 * t288;
t235 = mrSges(6,1) * t252 + mrSges(6,2) * t253;
t273 = qJD(4) + qJD(5);
t243 = mrSges(6,1) * t273 - mrSges(6,3) * t253;
t272 = qJDD(4) + qJDD(5);
t236 = pkin(5) * t252 - pkin(10) * t253;
t271 = t273 ^ 2;
t189 = -pkin(5) * t271 + pkin(10) * t272 - t236 * t252 + t192;
t205 = t280 * t227 - t277 * t228;
t300 = -qJDD(2) * pkin(3) - t205;
t198 = -t260 * pkin(4) + t268 * t311 + (-pkin(9) * t275 - pkin(8)) * t291 + t300;
t221 = -qJD(5) * t252 + t259 * t288 + t260 * t284;
t193 = (t252 * t273 - t221) * pkin(10) + (t253 * t273 - t220) * pkin(5) + t198;
t283 = sin(qJ(6));
t287 = cos(qJ(6));
t186 = -t189 * t283 + t193 * t287;
t237 = -t253 * t283 + t273 * t287;
t209 = qJD(6) * t237 + t221 * t287 + t272 * t283;
t238 = t253 * t287 + t273 * t283;
t214 = -mrSges(7,1) * t237 + mrSges(7,2) * t238;
t219 = qJDD(6) - t220;
t246 = qJD(6) + t252;
t222 = -mrSges(7,2) * t246 + mrSges(7,3) * t237;
t182 = m(7) * t186 + mrSges(7,1) * t219 - mrSges(7,3) * t209 - t214 * t238 + t222 * t246;
t187 = t189 * t287 + t193 * t283;
t208 = -qJD(6) * t238 - t221 * t283 + t272 * t287;
t223 = mrSges(7,1) * t246 - mrSges(7,3) * t238;
t183 = m(7) * t187 - mrSges(7,2) * t219 + mrSges(7,3) * t208 + t214 * t237 - t223 * t246;
t303 = -t182 * t283 + t287 * t183;
t169 = m(6) * t192 - mrSges(6,2) * t272 + mrSges(6,3) * t220 - t235 * t252 - t243 * t273 + t303;
t191 = t196 * t288 - t197 * t284;
t242 = -mrSges(6,2) * t273 - mrSges(6,3) * t252;
t188 = -pkin(5) * t272 - pkin(10) * t271 + t236 * t253 - t191;
t297 = -m(7) * t188 + t208 * mrSges(7,1) - mrSges(7,2) * t209 + t237 * t222 - t223 * t238;
t178 = m(6) * t191 + mrSges(6,1) * t272 - mrSges(6,3) * t221 - t235 * t253 + t242 * t273 + t297;
t163 = t284 * t169 + t288 * t178;
t250 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t285 + Ifges(5,2) * t289) * qJD(2);
t251 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t285 + Ifges(5,4) * t289) * qJD(2);
t210 = Ifges(7,5) * t238 + Ifges(7,6) * t237 + Ifges(7,3) * t246;
t212 = Ifges(7,1) * t238 + Ifges(7,4) * t237 + Ifges(7,5) * t246;
t175 = -mrSges(7,1) * t188 + mrSges(7,3) * t187 + Ifges(7,4) * t209 + Ifges(7,2) * t208 + Ifges(7,6) * t219 - t210 * t238 + t212 * t246;
t211 = Ifges(7,4) * t238 + Ifges(7,2) * t237 + Ifges(7,6) * t246;
t176 = mrSges(7,2) * t188 - mrSges(7,3) * t186 + Ifges(7,1) * t209 + Ifges(7,4) * t208 + Ifges(7,5) * t219 + t210 * t237 - t211 * t246;
t230 = Ifges(6,4) * t253 - Ifges(6,2) * t252 + Ifges(6,6) * t273;
t231 = Ifges(6,1) * t253 - Ifges(6,4) * t252 + Ifges(6,5) * t273;
t295 = -mrSges(6,1) * t191 + mrSges(6,2) * t192 - Ifges(6,5) * t221 - Ifges(6,6) * t220 - Ifges(6,3) * t272 - pkin(5) * t297 - pkin(10) * t303 - t287 * t175 - t283 * t176 - t253 * t230 - t252 * t231;
t317 = mrSges(5,1) * t199 - mrSges(5,2) * t200 + Ifges(5,5) * t259 + Ifges(5,6) * t260 + Ifges(5,3) * qJDD(4) + pkin(4) * t163 + (t250 * t285 - t251 * t289) * qJD(2) - t295;
t258 = (-mrSges(5,1) * t289 + mrSges(5,2) * t285) * qJD(2);
t310 = qJD(2) * t289;
t265 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t310;
t161 = m(5) * t199 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t259 + qJD(4) * t265 - t258 * t311 + t163;
t264 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t311;
t304 = t288 * t169 - t178 * t284;
t162 = m(5) * t200 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t260 - qJD(4) * t264 + t258 * t310 + t304;
t305 = -t161 * t285 + t289 * t162;
t152 = m(4) * t206 - mrSges(4,1) * t291 - qJDD(2) * mrSges(4,2) + t305;
t202 = -t291 * pkin(8) + t300;
t171 = t287 * t182 + t283 * t183;
t296 = m(6) * t198 - t220 * mrSges(6,1) + mrSges(6,2) * t221 + t252 * t242 + t243 * t253 + t171;
t293 = -m(5) * t202 + t260 * mrSges(5,1) - mrSges(5,2) * t259 - t264 * t311 + t265 * t310 - t296;
t165 = m(4) * t205 + qJDD(2) * mrSges(4,1) - mrSges(4,2) * t291 + t293;
t149 = t277 * t152 + t280 * t165;
t147 = m(3) * t232 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t291 + t149;
t306 = t280 * t152 - t165 * t277;
t148 = m(3) * t233 - mrSges(3,1) * t291 - qJDD(2) * mrSges(3,2) + t306;
t138 = -t147 * t286 + t290 * t148;
t316 = pkin(7) * t138;
t155 = t289 * t161 + t285 * t162;
t308 = m(4) * t244 + t155;
t153 = m(3) * t245 + t308;
t135 = t147 * t312 + t148 * t313 - t153 * t279;
t229 = Ifges(6,5) * t253 - Ifges(6,6) * t252 + Ifges(6,3) * t273;
t156 = mrSges(6,2) * t198 - mrSges(6,3) * t191 + Ifges(6,1) * t221 + Ifges(6,4) * t220 + Ifges(6,5) * t272 - pkin(10) * t171 - t175 * t283 + t176 * t287 - t229 * t252 - t230 * t273;
t294 = mrSges(7,1) * t186 - mrSges(7,2) * t187 + Ifges(7,5) * t209 + Ifges(7,6) * t208 + Ifges(7,3) * t219 + t211 * t238 - t212 * t237;
t157 = -mrSges(6,1) * t198 + mrSges(6,3) * t192 + Ifges(6,4) * t221 + Ifges(6,2) * t220 + Ifges(6,6) * t272 - pkin(5) * t171 - t229 * t253 + t231 * t273 - t294;
t249 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t285 + Ifges(5,6) * t289) * qJD(2);
t141 = -mrSges(5,1) * t202 + mrSges(5,3) * t200 + Ifges(5,4) * t259 + Ifges(5,2) * t260 + Ifges(5,6) * qJDD(4) - pkin(4) * t296 + pkin(9) * t304 + qJD(4) * t251 + t284 * t156 + t288 * t157 - t249 * t311;
t143 = mrSges(5,2) * t202 - mrSges(5,3) * t199 + Ifges(5,1) * t259 + Ifges(5,4) * t260 + Ifges(5,5) * qJDD(4) - pkin(9) * t163 - qJD(4) * t250 + t156 * t288 - t157 * t284 + t249 * t310;
t131 = mrSges(4,2) * t244 - mrSges(4,3) * t205 + Ifges(4,5) * qJDD(2) - Ifges(4,6) * t291 - pkin(8) * t155 - t141 * t285 + t143 * t289;
t139 = -mrSges(4,1) * t244 + mrSges(4,3) * t206 + t291 * Ifges(4,5) + Ifges(4,6) * qJDD(2) - pkin(3) * t155 - t317;
t126 = -mrSges(3,1) * t245 + mrSges(3,3) * t233 + t291 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t308 + qJ(3) * t306 + t277 * t131 + t280 * t139;
t128 = mrSges(3,2) * t245 - mrSges(3,3) * t232 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t291 - qJ(3) * t149 + t131 * t280 - t139 * t277;
t298 = mrSges(4,1) * t205 - mrSges(4,2) * t206 + Ifges(4,3) * qJDD(2) + pkin(3) * t293 + pkin(8) * t305 + t289 * t141 + t285 * t143;
t130 = mrSges(3,1) * t232 - mrSges(3,2) * t233 + Ifges(3,3) * qJDD(2) + pkin(2) * t149 + t298;
t299 = mrSges(2,1) * t262 - mrSges(2,2) * t263 + pkin(1) * t135 + t126 * t314 + t128 * t315 + t282 * t130 + t279 * t316;
t136 = m(2) * t263 + t138;
t134 = t282 * t153 + (t147 * t290 + t148 * t286) * t279;
t132 = m(2) * t262 + t135;
t124 = mrSges(2,2) * t276 - mrSges(2,3) * t262 - t286 * t126 + t290 * t128 + (-t134 * t279 - t135 * t282) * pkin(7);
t123 = -mrSges(2,1) * t276 + mrSges(2,3) * t263 - pkin(1) * t134 - t279 * t130 + (t126 * t290 + t128 * t286 + t316) * t282;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t281 * t124 - t278 * t123 - qJ(1) * (t132 * t281 + t136 * t278), t124, t128, t131, t143, t156, t176; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t278 * t124 + t281 * t123 + qJ(1) * (-t132 * t278 + t136 * t281), t123, t126, t139, t141, t157, t175; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t299, t299, t130, t298, t317, -t295, t294;];
m_new  = t1;
