% Calculate vector of cutting torques with Newton-Euler for
% S6RPRPPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
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
% Datum: 2019-05-05 16:36
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRPPR2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR2_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR2_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR2_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR2_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR2_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR2_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR2_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 16:32:48
% EndTime: 2019-05-05 16:33:01
% DurationCPUTime: 7.55s
% Computational Cost: add. (117877->345), mult. (256940->414), div. (0->0), fcn. (161023->10), ass. (0->133)
t306 = sin(qJ(1));
t309 = cos(qJ(1));
t288 = t306 * g(1) - g(2) * t309;
t279 = qJDD(1) * pkin(1) + t288;
t289 = -g(1) * t309 - g(2) * t306;
t311 = qJD(1) ^ 2;
t281 = -pkin(1) * t311 + t289;
t302 = sin(pkin(9));
t303 = cos(pkin(9));
t250 = t302 * t279 + t303 * t281;
t231 = -pkin(2) * t311 + qJDD(1) * pkin(7) + t250;
t300 = -g(3) + qJDD(2);
t305 = sin(qJ(3));
t308 = cos(qJ(3));
t219 = -t305 * t231 + t308 * t300;
t334 = qJD(1) * qJD(3);
t333 = t308 * t334;
t282 = qJDD(1) * t305 + t333;
t202 = (-t282 + t333) * qJ(4) + (t305 * t308 * t311 + qJDD(3)) * pkin(3) + t219;
t220 = t308 * t231 + t305 * t300;
t283 = qJDD(1) * t308 - t305 * t334;
t338 = qJD(1) * t305;
t285 = qJD(3) * pkin(3) - qJ(4) * t338;
t299 = t308 ^ 2;
t203 = -pkin(3) * t299 * t311 + qJ(4) * t283 - qJD(3) * t285 + t220;
t301 = sin(pkin(10));
t342 = cos(pkin(10));
t268 = (t301 * t308 + t342 * t305) * qJD(1);
t347 = -2 * qJD(4);
t196 = t342 * t202 - t301 * t203 + t268 * t347;
t337 = qJD(1) * t308;
t267 = t301 * t338 - t342 * t337;
t263 = t267 * t347;
t341 = t301 * t202 + t342 * t203;
t197 = t263 + t341;
t227 = Ifges(5,4) * t268 - Ifges(5,2) * t267 + Ifges(5,6) * qJD(3);
t237 = -mrSges(6,2) * t267 - mrSges(6,3) * t268;
t251 = t301 * t282 - t342 * t283;
t252 = t342 * t282 + t301 * t283;
t258 = mrSges(6,1) * t267 - qJD(3) * mrSges(6,3);
t235 = pkin(4) * t267 - qJ(5) * t268;
t310 = qJD(3) ^ 2;
t192 = -qJDD(3) * pkin(4) - t310 * qJ(5) + t268 * t235 + qJDD(5) - t196;
t336 = qJD(3) * t267;
t186 = (t267 * t268 - qJDD(3)) * pkin(8) + (t252 + t336) * pkin(5) + t192;
t260 = pkin(5) * t268 - qJD(3) * pkin(8);
t266 = t267 ^ 2;
t249 = t303 * t279 - t302 * t281;
t326 = -qJDD(1) * pkin(2) - t249;
t204 = -t283 * pkin(3) + qJDD(4) + t285 * t338 + (-qJ(4) * t299 - pkin(7)) * t311 + t326;
t344 = -2 * qJD(5);
t314 = (-t252 + t336) * qJ(5) + t204 + (qJD(3) * pkin(4) + t344) * t268;
t189 = -t266 * pkin(5) - t268 * t260 + (pkin(4) + pkin(8)) * t251 + t314;
t304 = sin(qJ(6));
t307 = cos(qJ(6));
t183 = t186 * t307 - t189 * t304;
t253 = -qJD(3) * t304 + t267 * t307;
t213 = qJD(6) * t253 + qJDD(3) * t307 + t251 * t304;
t254 = qJD(3) * t307 + t267 * t304;
t216 = -mrSges(7,1) * t253 + mrSges(7,2) * t254;
t265 = qJD(6) + t268;
t221 = -mrSges(7,2) * t265 + mrSges(7,3) * t253;
t248 = qJDD(6) + t252;
t180 = m(7) * t183 + mrSges(7,1) * t248 - mrSges(7,3) * t213 - t216 * t254 + t221 * t265;
t184 = t186 * t304 + t189 * t307;
t212 = -qJD(6) * t254 - qJDD(3) * t304 + t251 * t307;
t222 = mrSges(7,1) * t265 - mrSges(7,3) * t254;
t181 = m(7) * t184 - mrSges(7,2) * t248 + mrSges(7,3) * t212 + t216 * t253 - t222 * t265;
t168 = t307 * t180 + t304 * t181;
t325 = t310 * pkin(4) - qJDD(3) * qJ(5) - t341;
t188 = -t251 * pkin(5) - t266 * pkin(8) - t267 * t235 + t263 + ((2 * qJD(5)) + t260) * qJD(3) - t325;
t205 = Ifges(7,5) * t254 + Ifges(7,6) * t253 + Ifges(7,3) * t265;
t207 = Ifges(7,1) * t254 + Ifges(7,4) * t253 + Ifges(7,5) * t265;
t171 = -mrSges(7,1) * t188 + mrSges(7,3) * t184 + Ifges(7,4) * t213 + Ifges(7,2) * t212 + Ifges(7,6) * t248 - t205 * t254 + t207 * t265;
t206 = Ifges(7,4) * t254 + Ifges(7,2) * t253 + Ifges(7,6) * t265;
t172 = mrSges(7,2) * t188 - mrSges(7,3) * t183 + Ifges(7,1) * t213 + Ifges(7,4) * t212 + Ifges(7,5) * t248 + t205 * t253 - t206 * t265;
t190 = qJD(3) * t344 + ((2 * qJD(4)) + t235) * t267 + t325;
t224 = Ifges(6,5) * qJD(3) - Ifges(6,6) * t268 + Ifges(6,3) * t267;
t320 = -mrSges(6,2) * t192 + mrSges(6,3) * t190 - Ifges(6,1) * qJDD(3) + Ifges(6,4) * t252 - Ifges(6,5) * t251 + pkin(8) * t168 + t304 * t171 - t307 * t172 + t268 * t224;
t185 = -m(7) * t188 + t212 * mrSges(7,1) - t213 * mrSges(7,2) + t253 * t221 - t254 * t222;
t259 = mrSges(6,1) * t268 + qJD(3) * mrSges(6,2);
t321 = -m(6) * t190 + qJDD(3) * mrSges(6,3) + qJD(3) * t259 - t185;
t323 = -m(6) * t192 - t252 * mrSges(6,1) - t268 * t237 - t168;
t226 = Ifges(6,4) * qJD(3) - Ifges(6,2) * t268 + Ifges(6,6) * t267;
t339 = Ifges(5,1) * t268 - Ifges(5,4) * t267 + Ifges(5,5) * qJD(3) - t226;
t348 = -mrSges(5,2) * t197 + pkin(4) * (-qJDD(3) * mrSges(6,2) - qJD(3) * t258 + t323) + qJ(5) * (-t251 * mrSges(6,1) - t267 * t237 + t321) + mrSges(5,1) * t196 + t268 * t227 - Ifges(5,6) * t251 + Ifges(5,5) * t252 + Ifges(5,3) * qJDD(3) - t320 + t339 * t267;
t236 = mrSges(5,1) * t267 + mrSges(5,2) * t268;
t256 = -qJD(3) * mrSges(5,2) - mrSges(5,3) * t267;
t163 = m(5) * t196 - t252 * mrSges(5,3) - t268 * t236 + (mrSges(5,1) - mrSges(6,2)) * qJDD(3) + (t256 - t258) * qJD(3) + t323;
t257 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t268;
t175 = m(5) * t197 - qJDD(3) * mrSges(5,2) - qJD(3) * t257 + (-t236 - t237) * t267 + (-mrSges(5,3) - mrSges(6,1)) * t251 + t321;
t159 = t342 * t163 + t301 * t175;
t273 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t305 + Ifges(4,2) * t308) * qJD(1);
t274 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t305 + Ifges(4,4) * t308) * qJD(1);
t346 = mrSges(4,1) * t219 - mrSges(4,2) * t220 + Ifges(4,5) * t282 + Ifges(4,6) * t283 + Ifges(4,3) * qJDD(3) + pkin(3) * t159 + (t273 * t305 - t274 * t308) * qJD(1) + t348;
t343 = Ifges(5,4) + Ifges(6,6);
t280 = (-mrSges(4,1) * t308 + mrSges(4,2) * t305) * qJD(1);
t287 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t337;
t157 = m(4) * t219 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t282 + qJD(3) * t287 - t280 * t338 + t159;
t286 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t338;
t329 = -t163 * t301 + t342 * t175;
t158 = m(4) * t220 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t283 - qJD(3) * t286 + t280 * t337 + t329;
t330 = -t157 * t305 + t308 * t158;
t150 = m(3) * t250 - mrSges(3,1) * t311 - qJDD(1) * mrSges(3,2) + t330;
t230 = -t311 * pkin(7) + t326;
t169 = -t304 * t180 + t307 * t181;
t194 = t251 * pkin(4) + t314;
t167 = m(6) * t194 - t251 * mrSges(6,2) - t252 * mrSges(6,3) - t267 * t258 - t268 * t259 + t169;
t317 = m(5) * t204 + t251 * mrSges(5,1) + mrSges(5,2) * t252 + t267 * t256 + t257 * t268 + t167;
t313 = -m(4) * t230 + t283 * mrSges(4,1) - mrSges(4,2) * t282 - t286 * t338 + t287 * t337 - t317;
t161 = m(3) * t249 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t311 + t313;
t146 = t302 * t150 + t303 * t161;
t152 = t308 * t157 + t305 * t158;
t228 = Ifges(6,1) * qJD(3) - Ifges(6,4) * t268 + Ifges(6,5) * t267;
t340 = -Ifges(5,5) * t268 + Ifges(5,6) * t267 - Ifges(5,3) * qJD(3) - t228;
t331 = t303 * t150 - t161 * t302;
t318 = -mrSges(6,1) * t190 + mrSges(6,2) * t194 - pkin(5) * t185 - pkin(8) * t169 - t307 * t171 - t304 * t172;
t147 = -mrSges(5,1) * t204 + mrSges(5,3) * t197 - pkin(4) * t167 + t340 * t268 + t343 * t252 + (-Ifges(5,2) - Ifges(6,3)) * t251 + (Ifges(5,6) - Ifges(6,5)) * qJDD(3) + t339 * qJD(3) + t318;
t322 = mrSges(7,1) * t183 - mrSges(7,2) * t184 + Ifges(7,5) * t213 + Ifges(7,6) * t212 + Ifges(7,3) * t248 + t254 * t206 - t253 * t207;
t316 = mrSges(6,1) * t192 - mrSges(6,3) * t194 + pkin(5) * t168 + t322;
t153 = -qJ(5) * t167 + t316 + t340 * t267 + (Ifges(5,1) + Ifges(6,2)) * t252 - t343 * t251 + (Ifges(5,5) - Ifges(6,4)) * qJDD(3) + (-t227 + t224) * qJD(3) - mrSges(5,3) * t196 + mrSges(5,2) * t204;
t272 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t305 + Ifges(4,6) * t308) * qJD(1);
t139 = -mrSges(4,1) * t230 + mrSges(4,3) * t220 + Ifges(4,4) * t282 + Ifges(4,2) * t283 + Ifges(4,6) * qJDD(3) - pkin(3) * t317 + qJ(4) * t329 + qJD(3) * t274 + t342 * t147 + t301 * t153 - t272 * t338;
t142 = mrSges(4,2) * t230 - mrSges(4,3) * t219 + Ifges(4,1) * t282 + Ifges(4,4) * t283 + Ifges(4,5) * qJDD(3) - qJ(4) * t159 - qJD(3) * t273 - t301 * t147 + t342 * t153 + t272 * t337;
t324 = mrSges(3,1) * t249 - mrSges(3,2) * t250 + Ifges(3,3) * qJDD(1) + pkin(2) * t313 + pkin(7) * t330 + t308 * t139 + t305 * t142;
t319 = mrSges(2,1) * t288 - mrSges(2,2) * t289 + Ifges(2,3) * qJDD(1) + pkin(1) * t146 + t324;
t144 = m(2) * t289 - mrSges(2,1) * t311 - qJDD(1) * mrSges(2,2) + t331;
t143 = m(2) * t288 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t311 + t146;
t140 = -mrSges(3,1) * t300 + mrSges(3,3) * t250 + t311 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t152 - t346;
t137 = mrSges(3,2) * t300 - mrSges(3,3) * t249 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t311 - pkin(7) * t152 - t139 * t305 + t142 * t308;
t136 = -mrSges(2,2) * g(3) - mrSges(2,3) * t288 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t311 - qJ(2) * t146 + t137 * t303 - t140 * t302;
t135 = Ifges(2,6) * qJDD(1) + t311 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t289 + t302 * t137 + t303 * t140 - pkin(1) * (m(3) * t300 + t152) + qJ(2) * t331;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t309 * t136 - t306 * t135 - pkin(6) * (t143 * t309 + t144 * t306), t136, t137, t142, t153, -t267 * t226 - t320, t172; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t306 * t136 + t309 * t135 + pkin(6) * (-t143 * t306 + t144 * t309), t135, t140, t139, t147, Ifges(6,4) * qJDD(3) - Ifges(6,2) * t252 + Ifges(6,6) * t251 - qJD(3) * t224 + t267 * t228 - t316, t171; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t319, t319, t324, t346, t348, Ifges(6,5) * qJDD(3) - Ifges(6,6) * t252 + Ifges(6,3) * t251 + qJD(3) * t226 + t268 * t228 - t318, t322;];
m_new  = t1;
