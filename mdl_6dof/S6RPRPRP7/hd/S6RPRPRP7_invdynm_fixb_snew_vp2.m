% Calculate vector of cutting torques with Newton-Euler for
% S6RPRPRP7
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
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
% Datum: 2019-05-05 17:59
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRPRP7_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP7_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP7_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP7_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP7_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP7_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP7_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP7_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP7_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:56:50
% EndTime: 2019-05-05 17:57:02
% DurationCPUTime: 6.08s
% Computational Cost: add. (97470->337), mult. (211658->402), div. (0->0), fcn. (135740->8), ass. (0->127)
t349 = -2 * qJD(4);
t305 = sin(qJ(1));
t308 = cos(qJ(1));
t286 = t305 * g(1) - t308 * g(2);
t310 = qJD(1) ^ 2;
t325 = -t310 * qJ(2) + qJDD(2) - t286;
t347 = -pkin(1) - pkin(7);
t259 = t347 * qJDD(1) + t325;
t304 = sin(qJ(3));
t307 = cos(qJ(3));
t248 = t304 * g(3) + t307 * t259;
t337 = qJD(1) * qJD(3);
t333 = t304 * t337;
t281 = t307 * qJDD(1) - t333;
t218 = (-t281 - t333) * qJ(4) + (-t304 * t307 * t310 + qJDD(3)) * pkin(3) + t248;
t249 = -t307 * g(3) + t304 * t259;
t280 = -t304 * qJDD(1) - t307 * t337;
t339 = qJD(1) * t307;
t284 = qJD(3) * pkin(3) - qJ(4) * t339;
t298 = t304 ^ 2;
t219 = -t298 * t310 * pkin(3) + t280 * qJ(4) - qJD(3) * t284 + t249;
t301 = sin(pkin(9));
t302 = cos(pkin(9));
t340 = qJD(1) * t304;
t269 = -t301 * t340 + t302 * t339;
t189 = t302 * t218 - t301 * t219 + t269 * t349;
t287 = -t308 * g(1) - t305 * g(2);
t326 = qJDD(1) * qJ(2) + 0.2e1 * qJD(2) * qJD(1) + t287;
t268 = (t301 * t307 + t302 * t304) * qJD(1);
t246 = t301 * t280 + t302 * t281;
t303 = sin(qJ(5));
t306 = cos(qJ(5));
t251 = t306 * qJD(3) - t303 * t269;
t213 = t251 * qJD(5) + t303 * qJDD(3) + t306 * t246;
t252 = t303 * qJD(3) + t306 * t269;
t223 = -t251 * mrSges(7,1) + t252 * mrSges(7,2);
t190 = t301 * t218 + t302 * t219 + t268 * t349;
t236 = t268 * pkin(4) - t269 * pkin(8);
t309 = qJD(3) ^ 2;
t184 = -t309 * pkin(4) + qJDD(3) * pkin(8) - t268 * t236 + t190;
t221 = -t280 * pkin(3) + qJDD(4) + t284 * t339 + (-qJ(4) * t298 + t347) * t310 + t326;
t245 = t302 * t280 - t301 * t281;
t187 = (qJD(3) * t268 - t246) * pkin(8) + (qJD(3) * t269 - t245) * pkin(4) + t221;
t180 = -t303 * t184 + t306 * t187;
t244 = qJDD(5) - t245;
t266 = qJD(5) + t268;
t174 = -0.2e1 * qJD(6) * t252 + (t251 * t266 - t213) * qJ(6) + (t251 * t252 + t244) * pkin(5) + t180;
t226 = -t266 * mrSges(7,2) + t251 * mrSges(7,3);
t335 = m(7) * t174 + t244 * mrSges(7,1) + t266 * t226;
t171 = -t213 * mrSges(7,3) - t252 * t223 + t335;
t181 = t306 * t184 + t303 * t187;
t198 = Ifges(6,4) * t252 + Ifges(6,2) * t251 + Ifges(6,6) * t266;
t199 = Ifges(7,1) * t252 + Ifges(7,4) * t251 + Ifges(7,5) * t266;
t200 = Ifges(6,1) * t252 + Ifges(6,4) * t251 + Ifges(6,5) * t266;
t212 = -t252 * qJD(5) + t306 * qJDD(3) - t303 * t246;
t228 = t266 * pkin(5) - t252 * qJ(6);
t250 = t251 ^ 2;
t177 = -t250 * pkin(5) + t212 * qJ(6) + 0.2e1 * qJD(6) * t251 - t266 * t228 + t181;
t197 = Ifges(7,4) * t252 + Ifges(7,2) * t251 + Ifges(7,6) * t266;
t323 = -mrSges(7,1) * t174 + mrSges(7,2) * t177 - Ifges(7,5) * t213 - Ifges(7,6) * t212 - Ifges(7,3) * t244 - t252 * t197;
t348 = mrSges(6,1) * t180 - mrSges(6,2) * t181 + Ifges(6,5) * t213 + Ifges(6,6) * t212 + Ifges(6,3) * t244 + pkin(5) * t171 + t252 * t198 - (t200 + t199) * t251 - t323;
t346 = mrSges(2,1) - mrSges(3,2);
t345 = -mrSges(6,2) - mrSges(7,2);
t344 = Ifges(2,5) - Ifges(3,4);
t343 = -Ifges(2,6) + Ifges(3,5);
t235 = t268 * mrSges(5,1) + t269 * mrSges(5,2);
t258 = qJD(3) * mrSges(5,1) - t269 * mrSges(5,3);
t224 = -t251 * mrSges(6,1) + t252 * mrSges(6,2);
t227 = -t266 * mrSges(6,2) + t251 * mrSges(6,3);
t163 = m(6) * t180 + t244 * mrSges(6,1) + t266 * t227 + (-t223 - t224) * t252 + (-mrSges(6,3) - mrSges(7,3)) * t213 + t335;
t334 = m(7) * t177 + t212 * mrSges(7,3) + t251 * t223;
t229 = t266 * mrSges(7,1) - t252 * mrSges(7,3);
t341 = -t266 * mrSges(6,1) + t252 * mrSges(6,3) - t229;
t166 = m(6) * t181 + t212 * mrSges(6,3) + t251 * t224 + t345 * t244 + t341 * t266 + t334;
t331 = -t303 * t163 + t306 * t166;
t156 = m(5) * t190 - qJDD(3) * mrSges(5,2) + t245 * mrSges(5,3) - qJD(3) * t258 - t268 * t235 + t331;
t257 = -qJD(3) * mrSges(5,2) - t268 * mrSges(5,3);
t183 = -qJDD(3) * pkin(4) - t309 * pkin(8) + t269 * t236 - t189;
t179 = -t212 * pkin(5) - t250 * qJ(6) + t252 * t228 + qJDD(6) + t183;
t329 = -m(7) * t179 + t212 * mrSges(7,1) + t251 * t226;
t316 = -m(6) * t183 + t212 * mrSges(6,1) + t345 * t213 + t251 * t227 + t341 * t252 + t329;
t168 = m(5) * t189 + qJDD(3) * mrSges(5,1) - t246 * mrSges(5,3) + qJD(3) * t257 - t269 * t235 + t316;
t148 = t301 * t156 + t302 * t168;
t160 = t306 * t163 + t303 * t166;
t279 = (mrSges(4,1) * t304 + mrSges(4,2) * t307) * qJD(1);
t283 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t340;
t145 = m(4) * t248 + qJDD(3) * mrSges(4,1) - t281 * mrSges(4,3) + qJD(3) * t283 - t279 * t339 + t148;
t285 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t339;
t332 = t302 * t156 - t301 * t168;
t146 = m(4) * t249 - qJDD(3) * mrSges(4,2) + t280 * mrSges(4,3) - qJD(3) * t285 - t279 * t340 + t332;
t142 = -t304 * t145 + t307 * t146;
t141 = t307 * t145 + t304 * t146;
t324 = -mrSges(7,1) * t179 + mrSges(7,3) * t177 + Ifges(7,4) * t213 + Ifges(7,2) * t212 + Ifges(7,6) * t244 + t266 * t199;
t195 = Ifges(7,5) * t252 + Ifges(7,6) * t251 + Ifges(7,3) * t266;
t322 = mrSges(7,2) * t179 - mrSges(7,3) * t174 + Ifges(7,1) * t213 + Ifges(7,4) * t212 + Ifges(7,5) * t244 + t251 * t195;
t267 = -qJDD(1) * pkin(1) + t325;
t321 = -m(3) * t267 + t310 * mrSges(3,3) - t141;
t320 = m(5) * t221 - t245 * mrSges(5,1) + t246 * mrSges(5,2) + t268 * t257 + t269 * t258 + t160;
t196 = Ifges(6,5) * t252 + Ifges(6,6) * t251 + Ifges(6,3) * t266;
t150 = Ifges(6,4) * t213 + Ifges(6,2) * t212 + Ifges(6,6) * t244 + t266 * t200 - mrSges(6,1) * t183 + mrSges(6,3) * t181 - pkin(5) * (t213 * mrSges(7,2) - t329) + qJ(6) * (-t244 * mrSges(7,2) - t266 * t229 + t334) + (-pkin(5) * t229 - t195 - t196) * t252 + t324;
t158 = mrSges(6,2) * t183 - mrSges(6,3) * t180 + Ifges(6,1) * t213 + Ifges(6,4) * t212 + Ifges(6,5) * t244 - qJ(6) * t171 + t251 * t196 + (-t197 - t198) * t266 + t322;
t231 = Ifges(5,5) * t269 - Ifges(5,6) * t268 + Ifges(5,3) * qJD(3);
t232 = Ifges(5,4) * t269 - Ifges(5,2) * t268 + Ifges(5,6) * qJD(3);
t137 = mrSges(5,2) * t221 - mrSges(5,3) * t189 + Ifges(5,1) * t246 + Ifges(5,4) * t245 + Ifges(5,5) * qJDD(3) - pkin(8) * t160 - qJD(3) * t232 - t303 * t150 + t306 * t158 - t268 * t231;
t233 = Ifges(5,1) * t269 - Ifges(5,4) * t268 + Ifges(5,5) * qJD(3);
t143 = -mrSges(5,1) * t221 + mrSges(5,3) * t190 + Ifges(5,4) * t246 + Ifges(5,2) * t245 + Ifges(5,6) * qJDD(3) - pkin(4) * t160 + qJD(3) * t233 - t269 * t231 - t348;
t256 = t347 * t310 + t326;
t270 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t307 - Ifges(4,6) * t304) * qJD(1);
t272 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t307 - Ifges(4,4) * t304) * qJD(1);
t134 = -mrSges(4,1) * t256 + mrSges(4,3) * t249 + Ifges(4,4) * t281 + Ifges(4,2) * t280 + Ifges(4,6) * qJDD(3) - pkin(3) * t320 + qJ(4) * t332 + qJD(3) * t272 + t301 * t137 + t302 * t143 - t270 * t339;
t271 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t307 - Ifges(4,2) * t304) * qJD(1);
t136 = mrSges(4,2) * t256 - mrSges(4,3) * t248 + Ifges(4,1) * t281 + Ifges(4,4) * t280 + Ifges(4,5) * qJDD(3) - qJ(4) * t148 - qJD(3) * t271 + t302 * t137 - t301 * t143 - t270 * t340;
t262 = t310 * pkin(1) - t326;
t319 = mrSges(3,2) * t267 - mrSges(3,3) * t262 + Ifges(3,1) * qJDD(1) - pkin(7) * t141 - t304 * t134 + t307 * t136;
t153 = -m(4) * t256 + t280 * mrSges(4,1) - t281 * mrSges(4,2) - t283 * t340 - t285 * t339 - t320;
t318 = -mrSges(3,1) * t262 - pkin(2) * t153 - pkin(7) * t142 - t307 * t134 - t304 * t136;
t317 = -mrSges(5,1) * t189 + mrSges(5,2) * t190 - Ifges(5,5) * t246 - Ifges(5,6) * t245 - Ifges(5,3) * qJDD(3) - pkin(4) * t316 - pkin(8) * t331 - t306 * t150 - t303 * t158 - t269 * t232 - t268 * t233;
t314 = -m(3) * t262 + t310 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t153;
t315 = -mrSges(2,2) * t287 + pkin(1) * (-qJDD(1) * mrSges(3,2) + t321) + qJ(2) * t314 + mrSges(2,1) * t286 + Ifges(2,3) * qJDD(1) + t319;
t313 = -mrSges(4,1) * t248 + mrSges(4,2) * t249 - Ifges(4,5) * t281 - Ifges(4,6) * t280 - Ifges(4,3) * qJDD(3) - pkin(3) * t148 - t271 * t339 - t272 * t340 + t317;
t311 = -mrSges(3,1) * t267 - pkin(2) * t141 + t313;
t151 = m(2) * t287 - t310 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t314;
t140 = -m(3) * g(3) + t142;
t138 = m(2) * t286 - t310 * mrSges(2,2) + t346 * qJDD(1) + t321;
t133 = -t311 + t343 * t310 + t344 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3) - mrSges(2,3) * t286 - qJ(2) * t140;
t132 = mrSges(2,3) * t287 - pkin(1) * t140 + t346 * g(3) - t343 * qJDD(1) + t344 * t310 + t318;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t308 * t133 - t305 * t132 - pkin(6) * (t308 * t138 + t305 * t151), t133, t319, t136, t137, t158, -t266 * t197 + t322; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t305 * t133 + t308 * t132 + pkin(6) * (-t305 * t138 + t308 * t151), t132, -mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) - t310 * Ifges(3,5) + t311, t134, t143, t150, -t252 * t195 + t324; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t315, t315, mrSges(3,2) * g(3) + t310 * Ifges(3,4) + Ifges(3,5) * qJDD(1) - t318, -t313, -t317, t348, -t251 * t199 - t323;];
m_new  = t1;
