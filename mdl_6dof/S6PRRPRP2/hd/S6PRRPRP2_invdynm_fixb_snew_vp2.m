% Calculate vector of cutting torques with Newton-Euler for
% S6PRRPRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-05-05 03:50
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRRPRP2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP2_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP2_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRP2_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP2_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP2_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP2_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRP2_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 03:43:35
% EndTime: 2019-05-05 03:43:59
% DurationCPUTime: 13.58s
% Computational Cost: add. (237930->336), mult. (507908->422), div. (0->0), fcn. (356983->12), ass. (0->134)
t354 = -2 * qJD(4);
t304 = sin(pkin(10));
t307 = cos(pkin(10));
t292 = g(1) * t304 - g(2) * t307;
t293 = -g(1) * t307 - g(2) * t304;
t302 = -g(3) + qJDD(1);
t313 = cos(qJ(2));
t308 = cos(pkin(6));
t311 = sin(qJ(2));
t344 = t308 * t311;
t305 = sin(pkin(6));
t345 = t305 * t311;
t252 = t292 * t344 + t313 * t293 + t302 * t345;
t315 = qJD(2) ^ 2;
t247 = -pkin(2) * t315 + qJDD(2) * pkin(8) + t252;
t272 = -t292 * t305 + t302 * t308;
t310 = sin(qJ(3));
t312 = cos(qJ(3));
t216 = -t310 * t247 + t312 * t272;
t337 = qJD(2) * qJD(3);
t335 = t312 * t337;
t289 = qJDD(2) * t310 + t335;
t211 = (-t289 + t335) * qJ(4) + (t310 * t312 * t315 + qJDD(3)) * pkin(3) + t216;
t217 = t312 * t247 + t310 * t272;
t290 = qJDD(2) * t312 - t310 * t337;
t340 = qJD(2) * t310;
t294 = qJD(3) * pkin(3) - qJ(4) * t340;
t301 = t312 ^ 2;
t212 = -pkin(3) * t301 * t315 + qJ(4) * t290 - qJD(3) * t294 + t217;
t303 = sin(pkin(11));
t306 = cos(pkin(11));
t277 = (t303 * t310 - t306 * t312) * qJD(2);
t205 = t303 * t211 + t306 * t212 + t277 * t354;
t278 = (t303 * t312 + t306 * t310) * qJD(2);
t254 = mrSges(5,1) * t277 + mrSges(5,2) * t278;
t264 = -t289 * t303 + t290 * t306;
t271 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t278;
t255 = pkin(4) * t277 - pkin(9) * t278;
t314 = qJD(3) ^ 2;
t202 = -pkin(4) * t314 + qJDD(3) * pkin(9) - t255 * t277 + t205;
t251 = -t311 * t293 + (t292 * t308 + t302 * t305) * t313;
t322 = -qJDD(2) * pkin(2) - t251;
t215 = -t290 * pkin(3) + qJDD(4) + t294 * t340 + (-qJ(4) * t301 - pkin(8)) * t315 + t322;
t265 = t289 * t306 + t290 * t303;
t207 = (qJD(3) * t277 - t265) * pkin(9) + (qJD(3) * t278 - t264) * pkin(4) + t215;
t309 = sin(qJ(5));
t350 = cos(qJ(5));
t199 = t350 * t202 + t309 * t207;
t267 = t309 * qJD(3) + t350 * t278;
t232 = qJD(5) * t267 - t350 * qJDD(3) + t265 * t309;
t276 = qJD(5) + t277;
t244 = mrSges(6,1) * t276 - mrSges(6,3) * t267;
t263 = qJDD(5) - t264;
t266 = -t350 * qJD(3) + t278 * t309;
t237 = pkin(5) * t266 - qJ(6) * t267;
t275 = t276 ^ 2;
t194 = -pkin(5) * t275 + qJ(6) * t263 + 0.2e1 * qJD(6) * t276 - t237 * t266 + t199;
t245 = -mrSges(7,1) * t276 + mrSges(7,2) * t267;
t336 = m(7) * t194 + t263 * mrSges(7,3) + t276 * t245;
t238 = mrSges(7,1) * t266 - mrSges(7,3) * t267;
t341 = -mrSges(6,1) * t266 - mrSges(6,2) * t267 - t238;
t348 = -mrSges(6,3) - mrSges(7,2);
t184 = m(6) * t199 - t263 * mrSges(6,2) + t348 * t232 - t276 * t244 + t341 * t266 + t336;
t198 = -t309 * t202 + t350 * t207;
t233 = -t266 * qJD(5) + t309 * qJDD(3) + t350 * t265;
t243 = -mrSges(6,2) * t276 - mrSges(6,3) * t266;
t196 = -t263 * pkin(5) - t275 * qJ(6) + t267 * t237 + qJDD(6) - t198;
t242 = -mrSges(7,2) * t266 + mrSges(7,3) * t276;
t330 = -m(7) * t196 + t263 * mrSges(7,1) + t276 * t242;
t186 = m(6) * t198 + t263 * mrSges(6,1) + t348 * t233 + t276 * t243 + t341 * t267 + t330;
t332 = t350 * t184 - t186 * t309;
t172 = m(5) * t205 - qJDD(3) * mrSges(5,2) + mrSges(5,3) * t264 - qJD(3) * t271 - t254 * t277 + t332;
t204 = t306 * t211 - t303 * t212 + t278 * t354;
t270 = -qJD(3) * mrSges(5,2) - mrSges(5,3) * t277;
t201 = -qJDD(3) * pkin(4) - t314 * pkin(9) + t278 * t255 - t204;
t197 = -0.2e1 * qJD(6) * t267 + (t266 * t276 - t233) * qJ(6) + (t267 * t276 + t232) * pkin(5) + t201;
t191 = m(7) * t197 + mrSges(7,1) * t232 - t233 * mrSges(7,3) + t242 * t266 - t267 * t245;
t319 = -m(6) * t201 - t232 * mrSges(6,1) - mrSges(6,2) * t233 - t266 * t243 - t244 * t267 - t191;
t181 = m(5) * t204 + qJDD(3) * mrSges(5,1) - mrSges(5,3) * t265 + qJD(3) * t270 - t254 * t278 + t319;
t167 = t303 * t172 + t306 * t181;
t288 = (-mrSges(4,1) * t312 + mrSges(4,2) * t310) * qJD(2);
t339 = qJD(2) * t312;
t296 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t339;
t165 = m(4) * t216 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t289 + qJD(3) * t296 - t288 * t340 + t167;
t295 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t340;
t333 = t306 * t172 - t181 * t303;
t166 = m(4) * t217 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t290 - qJD(3) * t295 + t288 * t339 + t333;
t159 = t312 * t165 + t310 * t166;
t281 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t310 + Ifges(4,2) * t312) * qJD(2);
t282 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t310 + Ifges(4,4) * t312) * qJD(2);
t222 = Ifges(7,1) * t267 + Ifges(7,4) * t276 + Ifges(7,5) * t266;
t223 = Ifges(6,1) * t267 - Ifges(6,4) * t266 + Ifges(6,5) * t276;
t329 = -mrSges(7,1) * t197 + mrSges(7,2) * t194;
t220 = Ifges(7,4) * t267 + Ifges(7,2) * t276 + Ifges(7,6) * t266;
t343 = -Ifges(6,5) * t267 + Ifges(6,6) * t266 - Ifges(6,3) * t276 - t220;
t174 = -mrSges(6,1) * t201 + mrSges(6,3) * t199 - pkin(5) * t191 + (t222 + t223) * t276 + t343 * t267 + (Ifges(6,6) - Ifges(7,6)) * t263 + (Ifges(6,4) - Ifges(7,5)) * t233 + (-Ifges(6,2) - Ifges(7,3)) * t232 + t329;
t221 = Ifges(6,4) * t267 - Ifges(6,2) * t266 + Ifges(6,6) * t276;
t218 = Ifges(7,5) * t267 + Ifges(7,6) * t276 + Ifges(7,3) * t266;
t324 = mrSges(7,2) * t196 - mrSges(7,3) * t197 + Ifges(7,1) * t233 + Ifges(7,4) * t263 + Ifges(7,5) * t232 + t276 * t218;
t176 = mrSges(6,2) * t201 - mrSges(6,3) * t198 + Ifges(6,1) * t233 - Ifges(6,4) * t232 + Ifges(6,5) * t263 - qJ(6) * t191 - t276 * t221 + t343 * t266 + t324;
t249 = Ifges(5,4) * t278 - Ifges(5,2) * t277 + Ifges(5,6) * qJD(3);
t250 = Ifges(5,1) * t278 - Ifges(5,4) * t277 + Ifges(5,5) * qJD(3);
t320 = -mrSges(5,1) * t204 + mrSges(5,2) * t205 - Ifges(5,5) * t265 - Ifges(5,6) * t264 - Ifges(5,3) * qJDD(3) - pkin(4) * t319 - pkin(9) * t332 - t350 * t174 - t309 * t176 - t278 * t249 - t277 * t250;
t351 = mrSges(4,1) * t216 - mrSges(4,2) * t217 + Ifges(4,5) * t289 + Ifges(4,6) * t290 + Ifges(4,3) * qJDD(3) + pkin(3) * t167 + (t281 * t310 - t282 * t312) * qJD(2) - t320;
t145 = -mrSges(3,1) * t272 + mrSges(3,3) * t252 + t315 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t159 - t351;
t334 = -t165 * t310 + t312 * t166;
t157 = m(3) * t252 - mrSges(3,1) * t315 - qJDD(2) * mrSges(3,2) + t334;
t246 = -t315 * pkin(8) + t322;
t178 = t309 * t184 + t350 * t186;
t321 = m(5) * t215 - t264 * mrSges(5,1) + mrSges(5,2) * t265 + t277 * t270 + t271 * t278 + t178;
t318 = -m(4) * t246 + t290 * mrSges(4,1) - mrSges(4,2) * t289 - t295 * t340 + t296 * t339 - t321;
t169 = m(3) * t251 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t315 + t318;
t154 = t313 * t157 - t169 * t311;
t353 = pkin(7) * t154 + t145 * t313;
t325 = mrSges(7,1) * t196 - mrSges(7,3) * t194 - Ifges(7,4) * t233 - Ifges(7,2) * t263 - Ifges(7,6) * t232 - t266 * t222;
t352 = -(-t221 + t218) * t267 + mrSges(6,1) * t198 - mrSges(6,2) * t199 + Ifges(6,5) * t233 - Ifges(6,6) * t232 + Ifges(6,3) * t263 + pkin(5) * (-t233 * mrSges(7,2) - t267 * t238 + t330) + qJ(6) * (-t232 * mrSges(7,2) - t266 * t238 + t336) + t266 * t223 - t325;
t346 = t169 * t313;
t158 = m(3) * t272 + t159;
t150 = t157 * t344 - t158 * t305 + t308 * t346;
t248 = Ifges(5,5) * t278 - Ifges(5,6) * t277 + Ifges(5,3) * qJD(3);
t160 = mrSges(5,2) * t215 - mrSges(5,3) * t204 + Ifges(5,1) * t265 + Ifges(5,4) * t264 + Ifges(5,5) * qJDD(3) - pkin(9) * t178 - qJD(3) * t249 - t309 * t174 + t350 * t176 - t277 * t248;
t161 = -mrSges(5,1) * t215 + mrSges(5,3) * t205 + Ifges(5,4) * t265 + Ifges(5,2) * t264 + Ifges(5,6) * qJDD(3) - pkin(4) * t178 + qJD(3) * t250 - t278 * t248 - t352;
t280 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t310 + Ifges(4,6) * t312) * qJD(2);
t146 = -mrSges(4,1) * t246 + mrSges(4,3) * t217 + Ifges(4,4) * t289 + Ifges(4,2) * t290 + Ifges(4,6) * qJDD(3) - pkin(3) * t321 + qJ(4) * t333 + qJD(3) * t282 + t303 * t160 + t306 * t161 - t280 * t340;
t151 = mrSges(4,2) * t246 - mrSges(4,3) * t216 + Ifges(4,1) * t289 + Ifges(4,4) * t290 + Ifges(4,5) * qJDD(3) - qJ(4) * t167 - qJD(3) * t281 + t160 * t306 - t161 * t303 + t280 * t339;
t141 = mrSges(3,1) * t251 - mrSges(3,2) * t252 + Ifges(3,3) * qJDD(2) + pkin(2) * t318 + pkin(8) * t334 + t312 * t146 + t310 * t151;
t143 = mrSges(3,2) * t272 - mrSges(3,3) * t251 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t315 - pkin(8) * t159 - t146 * t310 + t151 * t312;
t323 = mrSges(2,1) * t292 - mrSges(2,2) * t293 + pkin(1) * t150 + t308 * t141 + t143 * t345 + t353 * t305;
t152 = m(2) * t293 + t154;
t149 = t308 * t158 + (t157 * t311 + t346) * t305;
t147 = m(2) * t292 + t150;
t139 = mrSges(2,2) * t302 - mrSges(2,3) * t292 + t313 * t143 - t311 * t145 + (-t149 * t305 - t150 * t308) * pkin(7);
t138 = -mrSges(2,1) * t302 + mrSges(2,3) * t293 - pkin(1) * t149 - t305 * t141 + (t143 * t311 + t353) * t308;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t307 * t139 - t304 * t138 - qJ(1) * (t147 * t307 + t152 * t304), t139, t143, t151, t160, t176, -t220 * t266 + t324; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t304 * t139 + t307 * t138 + qJ(1) * (-t147 * t304 + t152 * t307), t138, t145, t146, t161, t174, -t267 * t218 - t325; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t323, t323, t141, t351, -t320, t352, Ifges(7,5) * t233 + Ifges(7,6) * t263 + Ifges(7,3) * t232 + t267 * t220 - t276 * t222 - t329;];
m_new  = t1;
