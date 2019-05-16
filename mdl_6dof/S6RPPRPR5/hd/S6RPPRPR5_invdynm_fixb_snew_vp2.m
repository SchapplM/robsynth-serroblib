% Calculate vector of cutting torques with Newton-Euler for
% S6RPPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta5]';
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
% Datum: 2019-05-05 14:22
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPPRPR5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR5_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR5_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR5_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR5_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR5_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR5_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR5_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:21:03
% EndTime: 2019-05-05 14:21:12
% DurationCPUTime: 5.00s
% Computational Cost: add. (70090->299), mult. (140737->353), div. (0->0), fcn. (79385->8), ass. (0->114)
t275 = sin(qJ(1));
t278 = cos(qJ(1));
t247 = t275 * g(1) - t278 * g(2);
t281 = qJD(1) ^ 2;
t222 = -qJDD(1) * pkin(1) - t281 * qJ(2) + qJDD(2) - t247;
t210 = -qJDD(1) * qJ(3) - (2 * qJD(3) * qJD(1)) + t222;
t248 = -t278 * g(1) - t275 * g(2);
t314 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t248;
t313 = mrSges(3,2) - mrSges(4,3);
t312 = Ifges(3,4) - Ifges(4,5);
t311 = Ifges(2,6) - Ifges(3,5);
t310 = mrSges(4,3) * t281;
t211 = qJDD(3) + (-pkin(1) - qJ(3)) * t281 + t314;
t206 = -qJDD(1) * pkin(7) + t211;
t274 = sin(qJ(4));
t277 = cos(qJ(4));
t199 = -g(3) * t277 + t274 * t206;
t241 = (mrSges(5,1) * t274 + mrSges(5,2) * t277) * qJD(1);
t307 = qJD(1) * qJD(4);
t303 = t277 * t307;
t242 = qJDD(1) * t274 + t303;
t308 = qJD(1) * t277;
t246 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t308;
t205 = -t281 * pkin(7) - t210;
t304 = t274 * t307;
t243 = qJDD(1) * t277 - t304;
t186 = (-t243 + t304) * qJ(5) + (t242 + t303) * pkin(4) + t205;
t240 = (pkin(4) * t274 - qJ(5) * t277) * qJD(1);
t280 = qJD(4) ^ 2;
t309 = qJD(1) * t274;
t189 = -pkin(4) * t280 + qJDD(4) * qJ(5) - t240 * t309 + t199;
t271 = sin(pkin(9));
t272 = cos(pkin(9));
t234 = qJD(4) * t271 + t272 * t308;
t170 = -0.2e1 * qJD(5) * t234 + t272 * t186 - t271 * t189;
t216 = qJDD(4) * t271 + t243 * t272;
t233 = qJD(4) * t272 - t271 * t308;
t168 = (t233 * t309 - t216) * pkin(8) + (t233 * t234 + t242) * pkin(5) + t170;
t171 = 0.2e1 * qJD(5) * t233 + t271 * t186 + t272 * t189;
t215 = qJDD(4) * t272 - t243 * t271;
t217 = pkin(5) * t309 - pkin(8) * t234;
t232 = t233 ^ 2;
t169 = -pkin(5) * t232 + pkin(8) * t215 - t217 * t309 + t171;
t273 = sin(qJ(6));
t276 = cos(qJ(6));
t166 = t168 * t276 - t169 * t273;
t200 = t233 * t276 - t234 * t273;
t178 = qJD(6) * t200 + t215 * t273 + t216 * t276;
t201 = t233 * t273 + t234 * t276;
t183 = -mrSges(7,1) * t200 + mrSges(7,2) * t201;
t249 = qJD(6) + t309;
t190 = -mrSges(7,2) * t249 + mrSges(7,3) * t200;
t239 = qJDD(6) + t242;
t161 = m(7) * t166 + mrSges(7,1) * t239 - t178 * mrSges(7,3) - t183 * t201 + t190 * t249;
t167 = t168 * t273 + t169 * t276;
t177 = -qJD(6) * t201 + t215 * t276 - t216 * t273;
t191 = mrSges(7,1) * t249 - mrSges(7,3) * t201;
t162 = m(7) * t167 - mrSges(7,2) * t239 + t177 * mrSges(7,3) + t183 * t200 - t191 * t249;
t153 = t276 * t161 + t273 * t162;
t202 = -mrSges(6,1) * t233 + mrSges(6,2) * t234;
t213 = -mrSges(6,2) * t309 + mrSges(6,3) * t233;
t151 = m(6) * t170 + mrSges(6,1) * t242 - mrSges(6,3) * t216 - t202 * t234 + t213 * t309 + t153;
t214 = mrSges(6,1) * t309 - mrSges(6,3) * t234;
t300 = -t161 * t273 + t276 * t162;
t152 = m(6) * t171 - mrSges(6,2) * t242 + mrSges(6,3) * t215 + t202 * t233 - t214 * t309 + t300;
t301 = -t151 * t271 + t272 * t152;
t144 = m(5) * t199 - qJDD(4) * mrSges(5,2) - mrSges(5,3) * t242 - qJD(4) * t246 - t241 * t309 + t301;
t198 = g(3) * t274 + t206 * t277;
t245 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t309;
t188 = -qJDD(4) * pkin(4) - qJ(5) * t280 + t240 * t308 + qJDD(5) - t198;
t172 = -pkin(5) * t215 - pkin(8) * t232 + t217 * t234 + t188;
t293 = m(7) * t172 - t177 * mrSges(7,1) + t178 * mrSges(7,2) - t200 * t190 + t191 * t201;
t285 = -m(6) * t188 + t215 * mrSges(6,1) - mrSges(6,2) * t216 + t233 * t213 - t214 * t234 - t293;
t157 = m(5) * t198 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t243 + qJD(4) * t245 - t241 * t308 + t285;
t133 = t274 * t144 + t277 * t157;
t146 = t272 * t151 + t271 * t152;
t302 = t277 * t144 - t274 * t157;
t299 = m(4) * t211 + qJDD(1) * mrSges(4,2) + t133;
t296 = -m(5) * t205 - t242 * mrSges(5,1) - t243 * mrSges(5,2) - t245 * t309 - t246 * t308 - t146;
t179 = Ifges(7,5) * t201 + Ifges(7,6) * t200 + Ifges(7,3) * t249;
t181 = Ifges(7,1) * t201 + Ifges(7,4) * t200 + Ifges(7,5) * t249;
t154 = -mrSges(7,1) * t172 + mrSges(7,3) * t167 + Ifges(7,4) * t178 + Ifges(7,2) * t177 + Ifges(7,6) * t239 - t179 * t201 + t181 * t249;
t180 = Ifges(7,4) * t201 + Ifges(7,2) * t200 + Ifges(7,6) * t249;
t155 = mrSges(7,2) * t172 - mrSges(7,3) * t166 + Ifges(7,1) * t178 + Ifges(7,4) * t177 + Ifges(7,5) * t239 + t179 * t200 - t180 * t249;
t194 = Ifges(6,5) * t234 + Ifges(6,6) * t233 + Ifges(6,3) * t309;
t196 = Ifges(6,1) * t234 + Ifges(6,4) * t233 + Ifges(6,5) * t309;
t129 = -mrSges(6,1) * t188 + mrSges(6,3) * t171 + Ifges(6,4) * t216 + Ifges(6,2) * t215 + Ifges(6,6) * t242 - pkin(5) * t293 + pkin(8) * t300 + t276 * t154 + t273 * t155 - t234 * t194 + t196 * t309;
t195 = Ifges(6,4) * t234 + Ifges(6,2) * t233 + Ifges(6,6) * t309;
t136 = mrSges(6,2) * t188 - mrSges(6,3) * t170 + Ifges(6,1) * t216 + Ifges(6,4) * t215 + Ifges(6,5) * t242 - pkin(8) * t153 - t154 * t273 + t155 * t276 + t194 * t233 - t195 * t309;
t226 = (Ifges(5,3) * qJD(4)) + (Ifges(5,5) * t277 - Ifges(5,6) * t274) * qJD(1);
t227 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t277 - Ifges(5,2) * t274) * qJD(1);
t122 = mrSges(5,2) * t205 - mrSges(5,3) * t198 + Ifges(5,1) * t243 - Ifges(5,4) * t242 + Ifges(5,5) * qJDD(4) - qJ(5) * t146 - qJD(4) * t227 - t129 * t271 + t136 * t272 - t226 * t309;
t228 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t277 - Ifges(5,4) * t274) * qJD(1);
t292 = -mrSges(7,1) * t166 + mrSges(7,2) * t167 - Ifges(7,5) * t178 - Ifges(7,6) * t177 - Ifges(7,3) * t239 - t201 * t180 + t200 * t181;
t282 = -mrSges(6,1) * t170 + mrSges(6,2) * t171 - Ifges(6,5) * t216 - Ifges(6,6) * t215 - pkin(5) * t153 - t234 * t195 + t233 * t196 + t292;
t124 = t282 + (-Ifges(5,2) - Ifges(6,3)) * t242 + Ifges(5,6) * qJDD(4) + Ifges(5,4) * t243 + qJD(4) * t228 - mrSges(5,1) * t205 + mrSges(5,3) * t199 - pkin(4) * t146 - t226 * t308;
t295 = mrSges(4,1) * t210 + mrSges(4,2) * g(3) + t281 * Ifges(4,4) + Ifges(4,5) * qJDD(1) + pkin(3) * t296 + pkin(7) * t302 + t274 * t122 + t277 * t124;
t220 = pkin(1) * t281 - t314;
t294 = -m(3) * t220 + t281 * mrSges(3,2) + qJDD(1) * mrSges(3,3) + t299;
t291 = mrSges(4,2) * t211 - mrSges(4,3) * t210 + Ifges(4,1) * qJDD(1) - pkin(7) * t133 + t277 * t122 - t124 * t274;
t290 = mrSges(5,1) * t198 - mrSges(5,2) * t199 + Ifges(5,5) * t243 - Ifges(5,6) * t242 + Ifges(5,3) * qJDD(4) + pkin(4) * t285 + qJ(5) * t301 + t272 * t129 + t271 * t136 + t227 * t308 + t228 * t309;
t139 = m(4) * t210 - t281 * mrSges(4,2) - qJDD(1) * mrSges(4,3) + t296;
t289 = mrSges(3,1) * t222 + pkin(2) * t139 + t295;
t288 = -m(3) * t222 + t281 * mrSges(3,3) - t139;
t287 = mrSges(4,1) * t211 - Ifges(4,4) * qJDD(1) + pkin(3) * t133 + t290;
t286 = mrSges(3,2) * t222 - mrSges(3,3) * t220 + Ifges(3,1) * qJDD(1) - qJ(3) * t139 + t291;
t284 = -mrSges(2,2) * t248 + qJ(2) * (t294 - t310) + pkin(1) * (-qJDD(1) * mrSges(3,2) + t288) + mrSges(2,1) * t247 + Ifges(2,3) * qJDD(1) + t286;
t283 = mrSges(3,1) * t220 + pkin(2) * (-t299 + t310) + qJ(3) * (-m(4) * g(3) + t302) - t287;
t137 = t288 - t281 * mrSges(2,2) + m(2) * t247 + (mrSges(2,1) - mrSges(3,2)) * qJDD(1);
t130 = (-m(3) - m(4)) * g(3) + t302;
t125 = m(2) * t248 - qJDD(1) * mrSges(2,2) + (-mrSges(2,1) - mrSges(4,3)) * t281 + t294;
t119 = -t283 + (Ifges(2,5) - t312) * t281 + (mrSges(2,1) - t313) * g(3) + t311 * qJDD(1) + mrSges(2,3) * t248 - pkin(1) * t130;
t118 = t289 - mrSges(2,3) * t247 - qJ(2) * t130 + (-mrSges(2,2) + mrSges(3,3)) * g(3) + (-Ifges(3,4) + Ifges(2,5)) * qJDD(1) - t311 * t281;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t278 * t118 - t275 * t119 - pkin(6) * (t125 * t275 + t137 * t278), t118, t286, t291, t122, t136, t155; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t275 * t118 + t278 * t119 + pkin(6) * (t125 * t278 - t137 * t275), t119, -mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) - t281 * Ifges(3,5) - t289, -mrSges(4,3) * g(3) - t281 * Ifges(4,5) - t287, t124, t129, t154; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t284, t284, Ifges(3,5) * qJDD(1) + t313 * g(3) + t312 * t281 + t283, t295, t290, Ifges(6,3) * t242 - t282, -t292;];
m_new  = t1;
