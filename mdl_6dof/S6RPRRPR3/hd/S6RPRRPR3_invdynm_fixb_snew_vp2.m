% Calculate vector of cutting torques with Newton-Euler for
% S6RPRRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
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
% Datum: 2019-05-05 22:19
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRRPR3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR3_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR3_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR3_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR3_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR3_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR3_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR3_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 22:15:16
% EndTime: 2019-05-05 22:15:31
% DurationCPUTime: 7.24s
% Computational Cost: add. (132340->339), mult. (252195->408), div. (0->0), fcn. (157544->10), ass. (0->131)
t305 = sin(qJ(4));
t306 = sin(qJ(3));
t333 = qJD(1) * t306;
t339 = cos(qJ(4));
t276 = -t339 * qJD(3) + t305 * t333;
t309 = cos(qJ(3));
t331 = qJD(1) * qJD(3);
t329 = t309 * t331;
t282 = t306 * qJDD(1) + t329;
t242 = -t276 * qJD(4) + t305 * qJDD(3) + t339 * t282;
t307 = sin(qJ(1));
t310 = cos(qJ(1));
t288 = t307 * g(1) - t310 * g(2);
t278 = qJDD(1) * pkin(1) + t288;
t289 = -t310 * g(1) - t307 * g(2);
t312 = qJD(1) ^ 2;
t280 = -t312 * pkin(1) + t289;
t302 = sin(pkin(10));
t303 = cos(pkin(10));
t247 = t302 * t278 + t303 * t280;
t226 = -t312 * pkin(2) + qJDD(1) * pkin(7) + t247;
t301 = -g(3) + qJDD(2);
t216 = -t306 * t226 + t309 * t301;
t281 = (-pkin(3) * t309 - pkin(8) * t306) * qJD(1);
t311 = qJD(3) ^ 2;
t322 = qJDD(3) * pkin(3) + t311 * pkin(8) - t281 * t333 + t216;
t332 = t309 * qJD(1);
t292 = qJD(4) - t332;
t337 = t276 * t292;
t343 = (-t242 + t337) * qJ(5) - t322;
t246 = t303 * t278 - t302 * t280;
t225 = -qJDD(1) * pkin(2) - t312 * pkin(7) - t246;
t330 = t306 * t331;
t283 = t309 * qJDD(1) - t330;
t207 = (-t282 - t329) * pkin(8) + (-t283 + t330) * pkin(3) + t225;
t217 = t309 * t226 + t306 * t301;
t213 = -t311 * pkin(3) + qJDD(3) * pkin(8) + t281 * t332 + t217;
t193 = t339 * t207 - t305 * t213;
t194 = t305 * t207 + t339 * t213;
t277 = t305 * qJD(3) + t339 * t333;
t227 = Ifges(6,5) * t277 + Ifges(6,6) * t292 + Ifges(6,3) * t276;
t230 = Ifges(5,4) * t277 - Ifges(5,2) * t276 + Ifges(5,6) * t292;
t232 = Ifges(5,1) * t277 - Ifges(5,4) * t276 + Ifges(5,5) * t292;
t241 = t277 * qJD(4) - t339 * qJDD(3) + t305 * t282;
t251 = t276 * mrSges(6,1) - t277 * mrSges(6,3);
t275 = qJDD(4) - t283;
t250 = t276 * pkin(4) - t277 * qJ(5);
t291 = t292 ^ 2;
t191 = -t275 * pkin(4) - t291 * qJ(5) + t277 * t250 + qJDD(5) - t193;
t183 = (-t242 - t337) * pkin(9) + (t276 * t277 - t275) * pkin(5) + t191;
t340 = 2 * qJD(5);
t189 = -t291 * pkin(4) + t275 * qJ(5) - t276 * t250 + t292 * t340 + t194;
t257 = -t292 * pkin(5) - t277 * pkin(9);
t274 = t276 ^ 2;
t184 = -t274 * pkin(5) + t241 * pkin(9) + t292 * t257 + t189;
t304 = sin(qJ(6));
t308 = cos(qJ(6));
t181 = t308 * t183 - t304 * t184;
t244 = t308 * t276 - t304 * t277;
t202 = t244 * qJD(6) + t304 * t241 + t308 * t242;
t245 = t304 * t276 + t308 * t277;
t214 = -t244 * mrSges(7,1) + t245 * mrSges(7,2);
t290 = qJD(6) - t292;
t219 = -t290 * mrSges(7,2) + t244 * mrSges(7,3);
t271 = qJDD(6) - t275;
t176 = m(7) * t181 + t271 * mrSges(7,1) - t202 * mrSges(7,3) - t245 * t214 + t290 * t219;
t182 = t304 * t183 + t308 * t184;
t201 = -t245 * qJD(6) + t308 * t241 - t304 * t242;
t220 = t290 * mrSges(7,1) - t245 * mrSges(7,3);
t177 = m(7) * t182 - t271 * mrSges(7,2) + t201 * mrSges(7,3) + t244 * t214 - t290 * t220;
t166 = t308 * t176 + t304 * t177;
t231 = Ifges(6,1) * t277 + Ifges(6,4) * t292 + Ifges(6,5) * t276;
t205 = Ifges(7,4) * t245 + Ifges(7,2) * t244 + Ifges(7,6) * t290;
t206 = Ifges(7,1) * t245 + Ifges(7,4) * t244 + Ifges(7,5) * t290;
t323 = mrSges(7,1) * t181 - mrSges(7,2) * t182 + Ifges(7,5) * t202 + Ifges(7,6) * t201 + Ifges(7,3) * t271 + t245 * t205 - t244 * t206;
t315 = mrSges(6,1) * t191 - mrSges(6,3) * t189 - Ifges(6,4) * t242 - Ifges(6,2) * t275 - Ifges(6,6) * t241 + pkin(5) * t166 - t276 * t231 + t323;
t256 = -t276 * mrSges(6,2) + t292 * mrSges(6,3);
t320 = -m(6) * t191 + t275 * mrSges(6,1) + t292 * t256 - t166;
t167 = -t304 * t176 + t308 * t177;
t255 = -t292 * mrSges(6,1) + t277 * mrSges(6,2);
t324 = m(6) * t189 + t275 * mrSges(6,3) + t292 * t255 + t167;
t342 = (t230 - t227) * t277 + mrSges(5,1) * t193 - mrSges(5,2) * t194 + Ifges(5,5) * t242 - Ifges(5,6) * t241 + Ifges(5,3) * t275 + pkin(4) * (-t242 * mrSges(6,2) - t277 * t251 + t320) + qJ(5) * (-t241 * mrSges(6,2) - t276 * t251 + t324) + t276 * t232 - t315;
t186 = -t274 * pkin(9) + (-pkin(4) - pkin(5)) * t241 + (-pkin(4) * t292 + t257 + t340) * t277 - t343;
t178 = -m(7) * t186 + t201 * mrSges(7,1) - t202 * mrSges(7,2) + t244 * t219 - t245 * t220;
t192 = -0.2e1 * qJD(5) * t277 + (t277 * t292 + t241) * pkin(4) + t343;
t174 = m(6) * t192 + t241 * mrSges(6,1) - t242 * mrSges(6,3) - t277 * t255 + t276 * t256 + t178;
t204 = Ifges(7,5) * t245 + Ifges(7,6) * t244 + Ifges(7,3) * t290;
t169 = -mrSges(7,1) * t186 + mrSges(7,3) * t182 + Ifges(7,4) * t202 + Ifges(7,2) * t201 + Ifges(7,6) * t271 - t245 * t204 + t290 * t206;
t170 = mrSges(7,2) * t186 - mrSges(7,3) * t181 + Ifges(7,1) * t202 + Ifges(7,4) * t201 + Ifges(7,5) * t271 + t244 * t204 - t290 * t205;
t317 = -mrSges(6,1) * t192 + mrSges(6,2) * t189 - pkin(5) * t178 - pkin(9) * t167 - t308 * t169 - t304 * t170;
t229 = Ifges(6,4) * t277 + Ifges(6,2) * t292 + Ifges(6,6) * t276;
t336 = -Ifges(5,5) * t277 + Ifges(5,6) * t276 - Ifges(5,3) * t292 - t229;
t146 = mrSges(5,1) * t322 + mrSges(5,3) * t194 - pkin(4) * t174 + (t232 + t231) * t292 + t336 * t277 + (Ifges(5,6) - Ifges(6,6)) * t275 + (Ifges(5,4) - Ifges(6,5)) * t242 + (-Ifges(5,2) - Ifges(6,3)) * t241 + t317;
t319 = mrSges(6,2) * t191 - mrSges(6,3) * t192 + Ifges(6,1) * t242 + Ifges(6,4) * t275 + Ifges(6,5) * t241 - pkin(9) * t166 - t304 * t169 + t308 * t170 + t292 * t227;
t147 = -mrSges(5,2) * t322 - mrSges(5,3) * t193 + Ifges(5,1) * t242 - Ifges(5,4) * t241 + Ifges(5,5) * t275 - qJ(5) * t174 - t292 * t230 + t276 * t336 + t319;
t254 = t292 * mrSges(5,1) - t277 * mrSges(5,3);
t334 = -t276 * mrSges(5,1) - t277 * mrSges(5,2) - t251;
t338 = -mrSges(5,3) - mrSges(6,2);
t162 = m(5) * t194 - t275 * mrSges(5,2) + t338 * t241 - t292 * t254 + t334 * t276 + t324;
t253 = -t292 * mrSges(5,2) - t276 * mrSges(5,3);
t163 = m(5) * t193 + t275 * mrSges(5,1) + t338 * t242 + t292 * t253 + t334 * t277 + t320;
t160 = t339 * t162 - t305 * t163;
t173 = m(5) * t322 - t241 * mrSges(5,1) - t242 * mrSges(5,2) - t276 * t253 - t277 * t254 - t174;
t264 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t306 + Ifges(4,2) * t309) * qJD(1);
t265 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t306 + Ifges(4,4) * t309) * qJD(1);
t341 = mrSges(4,1) * t216 - mrSges(4,2) * t217 + Ifges(4,5) * t282 + Ifges(4,6) * t283 + Ifges(4,3) * qJDD(3) + pkin(3) * t173 + pkin(8) * t160 + (t264 * t306 - t265 * t309) * qJD(1) + t339 * t146 + t305 * t147;
t279 = (-mrSges(4,1) * t309 + mrSges(4,2) * t306) * qJD(1);
t285 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t333;
t158 = m(4) * t217 - qJDD(3) * mrSges(4,2) + t283 * mrSges(4,3) - qJD(3) * t285 + t279 * t332 + t160;
t286 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t332;
t172 = m(4) * t216 + qJDD(3) * mrSges(4,1) - t282 * mrSges(4,3) + qJD(3) * t286 - t279 * t333 + t173;
t327 = t309 * t158 - t306 * t172;
t150 = m(3) * t247 - t312 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t327;
t159 = t305 * t162 + t339 * t163;
t316 = -m(4) * t225 + t283 * mrSges(4,1) - t282 * mrSges(4,2) - t285 * t333 + t286 * t332 - t159;
t154 = m(3) * t246 + qJDD(1) * mrSges(3,1) - t312 * mrSges(3,2) + t316;
t145 = t302 * t150 + t303 * t154;
t152 = t306 * t158 + t309 * t172;
t328 = t303 * t150 - t302 * t154;
t263 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t306 + Ifges(4,6) * t309) * qJD(1);
t139 = mrSges(4,2) * t225 - mrSges(4,3) * t216 + Ifges(4,1) * t282 + Ifges(4,4) * t283 + Ifges(4,5) * qJDD(3) - pkin(8) * t159 - qJD(3) * t264 - t305 * t146 + t339 * t147 + t263 * t332;
t141 = -mrSges(4,1) * t225 + mrSges(4,3) * t217 + Ifges(4,4) * t282 + Ifges(4,2) * t283 + Ifges(4,6) * qJDD(3) - pkin(3) * t159 + qJD(3) * t265 - t263 * t333 - t342;
t321 = mrSges(3,1) * t246 - mrSges(3,2) * t247 + Ifges(3,3) * qJDD(1) + pkin(2) * t316 + pkin(7) * t327 + t306 * t139 + t309 * t141;
t318 = mrSges(2,1) * t288 - mrSges(2,2) * t289 + Ifges(2,3) * qJDD(1) + pkin(1) * t145 + t321;
t143 = m(2) * t289 - t312 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t328;
t142 = m(2) * t288 + qJDD(1) * mrSges(2,1) - t312 * mrSges(2,2) + t145;
t137 = -mrSges(3,1) * t301 + mrSges(3,3) * t247 + t312 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t152 - t341;
t136 = mrSges(3,2) * t301 - mrSges(3,3) * t246 + Ifges(3,5) * qJDD(1) - t312 * Ifges(3,6) - pkin(7) * t152 + t309 * t139 - t306 * t141;
t135 = -mrSges(2,2) * g(3) - mrSges(2,3) * t288 + Ifges(2,5) * qJDD(1) - t312 * Ifges(2,6) - qJ(2) * t145 + t303 * t136 - t302 * t137;
t134 = Ifges(2,6) * qJDD(1) + t312 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t289 + t302 * t136 + t303 * t137 - pkin(1) * (m(3) * t301 + t152) + qJ(2) * t328;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t310 * t135 - t307 * t134 - pkin(6) * (t310 * t142 + t307 * t143), t135, t136, t139, t147, -t276 * t229 + t319, t170; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t307 * t135 + t310 * t134 + pkin(6) * (-t307 * t142 + t310 * t143), t134, t137, t141, t146, -t277 * t227 - t315, t169; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t318, t318, t321, t341, t342, Ifges(6,5) * t242 + Ifges(6,6) * t275 + Ifges(6,3) * t241 + t277 * t229 - t292 * t231 - t317, t323;];
m_new  = t1;
