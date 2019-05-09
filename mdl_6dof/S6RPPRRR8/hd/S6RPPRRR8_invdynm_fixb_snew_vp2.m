% Calculate vector of cutting torques with Newton-Euler for
% S6RPPRRR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
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
% Datum: 2019-05-05 16:16
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPPRRR8_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR8_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR8_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR8_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR8_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR8_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR8_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR8_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR8_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 16:13:53
% EndTime: 2019-05-05 16:14:13
% DurationCPUTime: 11.67s
% Computational Cost: add. (204349->321), mult. (463465->391), div. (0->0), fcn. (330434->10), ass. (0->140)
t311 = qJD(1) ^ 2;
t300 = sin(pkin(10));
t291 = t300 ^ 2;
t301 = cos(pkin(10));
t345 = t301 ^ 2 + t291;
t336 = t345 * mrSges(4,3);
t353 = t311 * t336;
t305 = sin(qJ(1));
t309 = cos(qJ(1));
t276 = t305 * g(1) - t309 * g(2);
t327 = -t311 * qJ(2) + qJDD(2) - t276;
t347 = -pkin(1) - qJ(3);
t352 = -(2 * qJD(1) * qJD(3)) + t347 * qJDD(1) + t327;
t277 = -t309 * g(1) - t305 * g(2);
t351 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t277;
t350 = pkin(3) * t311;
t349 = mrSges(2,1) - mrSges(3,2);
t348 = -Ifges(2,6) + Ifges(3,5);
t346 = Ifges(4,6) * t300;
t246 = t300 * g(3) + t352 * t301;
t227 = (-pkin(7) * qJDD(1) - t300 * t350) * t301 + t246;
t247 = -g(3) * t301 + t352 * t300;
t340 = qJDD(1) * t300;
t228 = -pkin(7) * t340 - t291 * t350 + t247;
t304 = sin(qJ(4));
t308 = cos(qJ(4));
t209 = t304 * t227 + t308 * t228;
t343 = qJD(1) * t301;
t344 = qJD(1) * t300;
t269 = -t304 * t343 - t308 * t344;
t330 = -t300 * t304 + t301 * t308;
t270 = t330 * qJD(1);
t241 = -mrSges(5,1) * t269 + mrSges(5,2) * t270;
t264 = t270 * qJD(4);
t339 = qJDD(1) * t301;
t249 = -t304 * t339 - t308 * t340 - t264;
t259 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t270;
t248 = -pkin(4) * t269 - pkin(8) * t270;
t310 = qJD(4) ^ 2;
t198 = -pkin(4) * t310 + qJDD(4) * pkin(8) + t248 * t269 + t209;
t326 = qJDD(3) + t351;
t233 = pkin(3) * t340 + (-t345 * pkin(7) + t347) * t311 + t326;
t342 = t269 * qJD(4);
t250 = t330 * qJDD(1) + t342;
t201 = (-t250 - t342) * pkin(8) + (-t249 + t264) * pkin(4) + t233;
t303 = sin(qJ(5));
t307 = cos(qJ(5));
t188 = -t303 * t198 + t307 * t201;
t253 = qJD(4) * t307 - t270 * t303;
t218 = qJD(5) * t253 + qJDD(4) * t303 + t250 * t307;
t245 = qJDD(5) - t249;
t254 = qJD(4) * t303 + t270 * t307;
t267 = qJD(5) - t269;
t185 = (t253 * t267 - t218) * pkin(9) + (t253 * t254 + t245) * pkin(5) + t188;
t189 = t307 * t198 + t303 * t201;
t217 = -qJD(5) * t254 + qJDD(4) * t307 - t250 * t303;
t231 = pkin(5) * t267 - pkin(9) * t254;
t252 = t253 ^ 2;
t186 = -pkin(5) * t252 + pkin(9) * t217 - t231 * t267 + t189;
t302 = sin(qJ(6));
t306 = cos(qJ(6));
t183 = t185 * t306 - t186 * t302;
t219 = t253 * t306 - t254 * t302;
t194 = qJD(6) * t219 + t217 * t302 + t218 * t306;
t220 = t253 * t302 + t254 * t306;
t206 = -mrSges(7,1) * t219 + mrSges(7,2) * t220;
t262 = qJD(6) + t267;
t210 = -mrSges(7,2) * t262 + mrSges(7,3) * t219;
t239 = qJDD(6) + t245;
t178 = m(7) * t183 + mrSges(7,1) * t239 - t194 * mrSges(7,3) - t206 * t220 + t210 * t262;
t184 = t185 * t302 + t186 * t306;
t193 = -qJD(6) * t220 + t217 * t306 - t218 * t302;
t211 = mrSges(7,1) * t262 - mrSges(7,3) * t220;
t179 = m(7) * t184 - mrSges(7,2) * t239 + t193 * mrSges(7,3) + t206 * t219 - t211 * t262;
t170 = t306 * t178 + t302 * t179;
t223 = -mrSges(6,1) * t253 + mrSges(6,2) * t254;
t229 = -mrSges(6,2) * t267 + mrSges(6,3) * t253;
t168 = m(6) * t188 + mrSges(6,1) * t245 - mrSges(6,3) * t218 - t223 * t254 + t229 * t267 + t170;
t230 = mrSges(6,1) * t267 - mrSges(6,3) * t254;
t333 = -t178 * t302 + t306 * t179;
t169 = m(6) * t189 - mrSges(6,2) * t245 + mrSges(6,3) * t217 + t223 * t253 - t230 * t267 + t333;
t334 = -t168 * t303 + t307 * t169;
t161 = m(5) * t209 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t249 - qJD(4) * t259 + t241 * t269 + t334;
t208 = t227 * t308 - t304 * t228;
t258 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t269;
t197 = -qJDD(4) * pkin(4) - pkin(8) * t310 + t270 * t248 - t208;
t187 = -pkin(5) * t217 - pkin(9) * t252 + t231 * t254 + t197;
t324 = m(7) * t187 - t193 * mrSges(7,1) + t194 * mrSges(7,2) - t219 * t210 + t211 * t220;
t315 = -m(6) * t197 + t217 * mrSges(6,1) - mrSges(6,2) * t218 + t253 * t229 - t230 * t254 - t324;
t174 = m(5) * t208 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t250 + qJD(4) * t258 - t241 * t270 + t315;
t153 = t304 * t161 + t308 * t174;
t163 = t307 * t168 + t303 * t169;
t337 = Ifges(3,4) + t346;
t329 = -qJDD(1) * mrSges(4,3) - t311 * (mrSges(4,1) * t300 + mrSges(4,2) * t301);
t148 = m(4) * t246 + t329 * t301 + t153;
t335 = t308 * t161 - t304 * t174;
t149 = m(4) * t247 + t329 * t300 + t335;
t145 = -t148 * t300 + t301 * t149;
t332 = Ifges(4,1) * t301 - Ifges(4,4) * t300;
t331 = Ifges(4,4) * t301 - Ifges(4,2) * t300;
t144 = t301 * t148 + t300 * t149;
t268 = -qJDD(1) * pkin(1) + t327;
t325 = -m(3) * t268 + t311 * mrSges(3,3) - t144;
t323 = -m(5) * t233 + t249 * mrSges(5,1) - t250 * mrSges(5,2) + t269 * t258 - t270 * t259 - t163;
t203 = Ifges(7,4) * t220 + Ifges(7,2) * t219 + Ifges(7,6) * t262;
t204 = Ifges(7,1) * t220 + Ifges(7,4) * t219 + Ifges(7,5) * t262;
t322 = -mrSges(7,1) * t183 + mrSges(7,2) * t184 - Ifges(7,5) * t194 - Ifges(7,6) * t193 - Ifges(7,3) * t239 - t220 * t203 + t219 * t204;
t202 = Ifges(7,5) * t220 + Ifges(7,6) * t219 + Ifges(7,3) * t262;
t171 = -mrSges(7,1) * t187 + mrSges(7,3) * t184 + Ifges(7,4) * t194 + Ifges(7,2) * t193 + Ifges(7,6) * t239 - t202 * t220 + t204 * t262;
t172 = mrSges(7,2) * t187 - mrSges(7,3) * t183 + Ifges(7,1) * t194 + Ifges(7,4) * t193 + Ifges(7,5) * t239 + t202 * t219 - t203 * t262;
t212 = Ifges(6,5) * t254 + Ifges(6,6) * t253 + Ifges(6,3) * t267;
t214 = Ifges(6,1) * t254 + Ifges(6,4) * t253 + Ifges(6,5) * t267;
t151 = -mrSges(6,1) * t197 + mrSges(6,3) * t189 + Ifges(6,4) * t218 + Ifges(6,2) * t217 + Ifges(6,6) * t245 - pkin(5) * t324 + pkin(9) * t333 + t306 * t171 + t302 * t172 - t254 * t212 + t267 * t214;
t213 = Ifges(6,4) * t254 + Ifges(6,2) * t253 + Ifges(6,6) * t267;
t155 = mrSges(6,2) * t197 - mrSges(6,3) * t188 + Ifges(6,1) * t218 + Ifges(6,4) * t217 + Ifges(6,5) * t245 - pkin(9) * t170 - t171 * t302 + t172 * t306 + t212 * t253 - t213 * t267;
t234 = Ifges(5,5) * t270 + Ifges(5,6) * t269 + Ifges(5,3) * qJD(4);
t235 = Ifges(5,4) * t270 + Ifges(5,2) * t269 + Ifges(5,6) * qJD(4);
t140 = mrSges(5,2) * t233 - mrSges(5,3) * t208 + Ifges(5,1) * t250 + Ifges(5,4) * t249 + Ifges(5,5) * qJDD(4) - pkin(8) * t163 - qJD(4) * t235 - t151 * t303 + t155 * t307 + t234 * t269;
t236 = Ifges(5,1) * t270 + Ifges(5,4) * t269 + Ifges(5,5) * qJD(4);
t312 = mrSges(6,1) * t188 - mrSges(6,2) * t189 + Ifges(6,5) * t218 + Ifges(6,6) * t217 + Ifges(6,3) * t245 + pkin(5) * t170 + t254 * t213 - t253 * t214 - t322;
t146 = -mrSges(5,1) * t233 + mrSges(5,3) * t209 + Ifges(5,4) * t250 + Ifges(5,2) * t249 + Ifges(5,6) * qJDD(4) - pkin(4) * t163 + qJD(4) * t236 - t270 * t234 - t312;
t257 = t347 * t311 + t326;
t272 = (Ifges(4,5) * t301 - t346) * qJD(1);
t137 = -mrSges(4,1) * t257 + mrSges(4,3) * t247 + pkin(3) * t323 + pkin(7) * t335 + t331 * qJDD(1) + t304 * t140 + t308 * t146 - t272 * t343;
t139 = mrSges(4,2) * t257 - mrSges(4,3) * t246 - pkin(7) * t153 + t332 * qJDD(1) + t308 * t140 - t304 * t146 - t272 * t344;
t261 = t311 * pkin(1) - t351;
t321 = mrSges(3,2) * t268 - mrSges(3,3) * t261 + Ifges(3,1) * qJDD(1) - qJ(3) * t144 - t137 * t300 + t301 * t139;
t319 = -m(4) * t257 - mrSges(4,1) * t340 - mrSges(4,2) * t339 + t323;
t320 = -mrSges(3,1) * t261 - pkin(2) * (t319 + t353) - qJ(3) * t145 - t301 * t137 - t300 * t139;
t318 = -mrSges(5,1) * t208 + mrSges(5,2) * t209 - Ifges(5,5) * t250 - Ifges(5,6) * t249 - Ifges(5,3) * qJDD(4) - pkin(4) * t315 - pkin(8) * t334 - t307 * t151 - t303 * t155 - t270 * t235 + t269 * t236;
t316 = -m(3) * t261 + t311 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t319;
t317 = -mrSges(2,2) * t277 + pkin(1) * (-qJDD(1) * mrSges(3,2) + t325) + qJ(2) * (t316 - t353) + mrSges(2,1) * t276 + Ifges(2,3) * qJDD(1) + t321;
t314 = -mrSges(4,1) * t246 + mrSges(4,2) * t247 - Ifges(4,5) * t339 - pkin(3) * t153 + t318 + (-t300 * t332 - t301 * t331) * t311;
t313 = -mrSges(3,1) * t268 - pkin(2) * t144 + t314;
t156 = t316 + (-mrSges(2,1) - t336) * t311 + m(2) * t277 - qJDD(1) * mrSges(2,2);
t143 = -m(3) * g(3) + t145;
t141 = m(2) * t276 - t311 * mrSges(2,2) + t349 * qJDD(1) + t325;
t136 = -t313 + t348 * t311 + (Ifges(2,5) - t337) * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3) - qJ(2) * t143 - mrSges(2,3) * t276;
t135 = mrSges(2,3) * t277 - pkin(1) * t143 + (-Ifges(3,4) + Ifges(2,5)) * t311 - t348 * qJDD(1) + t349 * g(3) + t320;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t309 * t136 - t305 * t135 - pkin(6) * (t141 * t309 + t156 * t305), t136, t321, t139, t140, t155, t172; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t305 * t136 + t309 * t135 + pkin(6) * (-t141 * t305 + t156 * t309), t135, -mrSges(3,3) * g(3) - t311 * Ifges(3,5) + t337 * qJDD(1) + t313, t137, t146, t151, t171; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t317, t317, mrSges(3,2) * g(3) + t311 * Ifges(3,4) + Ifges(3,5) * qJDD(1) - t320, -Ifges(4,6) * t340 - t314, -t318, t312, -t322;];
m_new  = t1;
