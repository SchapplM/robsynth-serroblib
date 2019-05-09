% Calculate vector of cutting torques with Newton-Euler for
% S6RPRRPP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
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
% Datum: 2019-05-05 21:24
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRRPP2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP2_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP2_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP2_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP2_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP2_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP2_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP2_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:20:42
% EndTime: 2019-05-05 21:20:53
% DurationCPUTime: 5.07s
% Computational Cost: add. (65676->338), mult. (123909->390), div. (0->0), fcn. (72492->8), ass. (0->122)
t305 = sin(qJ(1));
t307 = cos(qJ(1));
t287 = t305 * g(1) - g(2) * t307;
t277 = qJDD(1) * pkin(1) + t287;
t288 = -g(1) * t307 - g(2) * t305;
t309 = qJD(1) ^ 2;
t279 = -pkin(1) * t309 + t288;
t301 = sin(pkin(9));
t302 = cos(pkin(9));
t241 = t277 * t302 - t301 * t279;
t207 = -qJDD(1) * pkin(2) - pkin(7) * t309 - t241;
t304 = sin(qJ(3));
t306 = cos(qJ(3));
t329 = qJD(1) * qJD(3);
t327 = t306 * t329;
t281 = qJDD(1) * t304 + t327;
t294 = t304 * t329;
t282 = qJDD(1) * t306 - t294;
t188 = (-t281 - t327) * pkin(8) + (-t282 + t294) * pkin(3) + t207;
t242 = t301 * t277 + t302 * t279;
t208 = -pkin(2) * t309 + qJDD(1) * pkin(7) + t242;
t300 = -g(3) + qJDD(2);
t198 = t306 * t208 + t304 * t300;
t280 = (-pkin(3) * t306 - pkin(8) * t304) * qJD(1);
t308 = qJD(3) ^ 2;
t330 = qJD(1) * t306;
t192 = -pkin(3) * t308 + qJDD(3) * pkin(8) + t280 * t330 + t198;
t303 = sin(qJ(4));
t339 = cos(qJ(4));
t186 = t303 * t188 + t339 * t192;
t331 = qJD(1) * t304;
t275 = -t339 * qJD(3) + t303 * t331;
t276 = t303 * qJD(3) + t339 * t331;
t245 = pkin(4) * t275 - qJ(5) * t276;
t274 = -qJDD(4) + t282;
t290 = -qJD(4) + t330;
t289 = t290 ^ 2;
t340 = -2 * qJD(5);
t181 = -pkin(4) * t289 - t274 * qJ(5) - t275 * t245 + t290 * t340 + t186;
t254 = mrSges(6,1) * t290 + mrSges(6,2) * t276;
t345 = m(6) * t181 - t274 * mrSges(6,3) - t290 * t254;
t185 = t339 * t188 - t303 * t192;
t183 = t274 * pkin(4) - t289 * qJ(5) + t276 * t245 + qJDD(5) - t185;
t255 = -mrSges(6,2) * t275 - mrSges(6,3) * t290;
t344 = -m(6) * t183 - t274 * mrSges(6,1) - t290 * t255;
t239 = -t275 * qJD(4) + t303 * qJDD(3) + t339 * t281;
t197 = -t304 * t208 + t306 * t300;
t318 = qJDD(3) * pkin(3) + pkin(8) * t308 - t280 * t331 + t197;
t336 = t275 * t290;
t343 = -(t239 + t336) * qJ(5) - t318;
t249 = -mrSges(7,2) * t290 + mrSges(7,3) * t275;
t171 = -0.2e1 * qJD(6) * t276 + (-t239 + t336) * qJ(6) + (t275 * t276 + t274) * pkin(5) + t183;
t247 = -mrSges(7,1) * t275 + mrSges(7,2) * t276;
t323 = -m(7) * t171 + t239 * mrSges(7,3) + t276 * t247;
t167 = mrSges(7,1) * t274 + t249 * t290 - t323;
t252 = mrSges(7,1) * t290 - mrSges(7,3) * t276;
t238 = qJD(4) * t276 - t339 * qJDD(3) + t281 * t303;
t251 = pkin(5) * t290 - qJ(6) * t276;
t273 = t275 ^ 2;
t175 = -pkin(5) * t273 + qJ(6) * t238 + 0.2e1 * qJD(6) * t275 - t251 * t290 + t181;
t328 = m(7) * t175 + t238 * mrSges(7,3) + t275 * t247;
t168 = -mrSges(7,2) * t274 - t252 * t290 + t328;
t210 = Ifges(6,5) * t276 - Ifges(6,6) * t290 + Ifges(6,3) * t275;
t214 = Ifges(5,4) * t276 - Ifges(5,2) * t275 - Ifges(5,6) * t290;
t217 = Ifges(5,1) * t276 - Ifges(5,4) * t275 - Ifges(5,5) * t290;
t246 = mrSges(6,1) * t275 - mrSges(6,3) * t276;
t216 = Ifges(6,1) * t276 - Ifges(6,4) * t290 + Ifges(6,5) * t275;
t212 = Ifges(7,4) * t276 + Ifges(7,2) * t275 + Ifges(7,6) * t290;
t215 = Ifges(7,1) * t276 + Ifges(7,4) * t275 + Ifges(7,5) * t290;
t321 = mrSges(7,1) * t171 - mrSges(7,2) * t175 + Ifges(7,5) * t239 + Ifges(7,6) * t238 + Ifges(7,3) * t274 + t276 * t212 - t275 * t215;
t314 = mrSges(6,1) * t183 - mrSges(6,3) * t181 - Ifges(6,4) * t239 + Ifges(6,2) * t274 - Ifges(6,6) * t238 + pkin(5) * t167 - t275 * t216 + t321;
t342 = (t214 - t210) * t276 + mrSges(5,1) * t185 - mrSges(5,2) * t186 + Ifges(5,5) * t239 - Ifges(5,6) * t238 - Ifges(5,3) * t274 + pkin(4) * (-mrSges(6,2) * t239 - t246 * t276 - t167 + t344) + qJ(5) * (-mrSges(6,2) * t238 - t246 * t275 + t168 + t345) + t275 * t217 - t314;
t178 = -qJ(6) * t273 + qJDD(6) + (-pkin(4) - pkin(5)) * t238 + (pkin(4) * t290 + (2 * qJD(5)) + t251) * t276 - t343;
t169 = -m(7) * t178 + t238 * mrSges(7,1) - t239 * mrSges(7,2) + t275 * t249 - t276 * t252;
t184 = t276 * t340 + (-t276 * t290 + t238) * pkin(4) + t343;
t166 = m(6) * t184 + t238 * mrSges(6,1) - t239 * mrSges(6,3) - t276 * t254 + t275 * t255 + t169;
t209 = Ifges(7,5) * t276 + Ifges(7,6) * t275 + Ifges(7,3) * t290;
t320 = mrSges(7,1) * t178 - mrSges(7,3) * t175 - Ifges(7,4) * t239 - Ifges(7,2) * t238 - Ifges(7,6) * t274 - t290 * t215;
t313 = mrSges(6,1) * t184 - mrSges(6,2) * t181 + pkin(5) * t169 + qJ(6) * t168 - t320;
t213 = Ifges(6,4) * t276 - Ifges(6,2) * t290 + Ifges(6,6) * t275;
t334 = -Ifges(5,5) * t276 + Ifges(5,6) * t275 + Ifges(5,3) * t290 - t213;
t143 = (-t217 - t216) * t290 + (t209 + t334) * t276 + (-Ifges(5,6) + Ifges(6,6)) * t274 + (Ifges(5,4) - Ifges(6,5)) * t239 + (-Ifges(5,2) - Ifges(6,3)) * t238 + mrSges(5,1) * t318 + mrSges(5,3) * t186 - pkin(4) * t166 - t313;
t319 = mrSges(7,2) * t178 - mrSges(7,3) * t171 + Ifges(7,1) * t239 + Ifges(7,4) * t238 + Ifges(7,5) * t274 + t275 * t209;
t312 = mrSges(6,2) * t183 - mrSges(6,3) * t184 + Ifges(6,1) * t239 - Ifges(6,4) * t274 + Ifges(6,5) * t238 - qJ(6) * t167 - t290 * t210 + t319;
t149 = (t214 - t212) * t290 + t334 * t275 - Ifges(5,5) * t274 - Ifges(5,4) * t238 + Ifges(5,1) * t239 - mrSges(5,2) * t318 - mrSges(5,3) * t185 - qJ(5) * t166 + t312;
t250 = mrSges(5,2) * t290 - mrSges(5,3) * t275;
t332 = -mrSges(5,1) * t275 - mrSges(5,2) * t276 - t246;
t337 = -mrSges(5,3) - mrSges(6,2);
t161 = m(5) * t185 + (-t249 - t250) * t290 + t332 * t276 + (-mrSges(5,1) - mrSges(7,1)) * t274 + t337 * t239 + t323 + t344;
t253 = -mrSges(5,1) * t290 - mrSges(5,3) * t276;
t162 = m(5) * t186 + (-t252 + t253) * t290 + t332 * t275 + (mrSges(5,2) - mrSges(7,2)) * t274 + t337 * t238 + t328 + t345;
t157 = -t161 * t303 + t339 * t162;
t163 = m(5) * t318 - t238 * mrSges(5,1) - t239 * mrSges(5,2) - t275 * t250 - t276 * t253 - t166;
t261 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t304 + Ifges(4,2) * t306) * qJD(1);
t262 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t304 + Ifges(4,4) * t306) * qJD(1);
t341 = mrSges(4,1) * t197 - mrSges(4,2) * t198 + Ifges(4,5) * t281 + Ifges(4,6) * t282 + Ifges(4,3) * qJDD(3) + pkin(3) * t163 + pkin(8) * t157 + (t261 * t304 - t262 * t306) * qJD(1) + t339 * t143 + t303 * t149;
t335 = t290 * t212;
t278 = (-mrSges(4,1) * t306 + mrSges(4,2) * t304) * qJD(1);
t284 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t331;
t155 = m(4) * t198 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t282 - qJD(3) * t284 + t278 * t330 + t157;
t285 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t330;
t159 = m(4) * t197 + qJDD(3) * mrSges(4,1) - t281 * mrSges(4,3) + qJD(3) * t285 - t278 * t331 + t163;
t325 = t306 * t155 - t159 * t304;
t146 = m(3) * t242 - mrSges(3,1) * t309 - qJDD(1) * mrSges(3,2) + t325;
t156 = t339 * t161 + t303 * t162;
t315 = -m(4) * t207 + t282 * mrSges(4,1) - t281 * mrSges(4,2) - t284 * t331 + t285 * t330 - t156;
t151 = m(3) * t241 + qJDD(1) * mrSges(3,1) - t309 * mrSges(3,2) + t315;
t142 = t301 * t146 + t302 * t151;
t148 = t304 * t155 + t306 * t159;
t326 = t302 * t146 - t151 * t301;
t260 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t304 + Ifges(4,6) * t306) * qJD(1);
t136 = mrSges(4,2) * t207 - mrSges(4,3) * t197 + Ifges(4,1) * t281 + Ifges(4,4) * t282 + Ifges(4,5) * qJDD(3) - pkin(8) * t156 - qJD(3) * t261 - t303 * t143 + t339 * t149 + t260 * t330;
t138 = -mrSges(4,1) * t207 + mrSges(4,3) * t198 + Ifges(4,4) * t281 + Ifges(4,2) * t282 + Ifges(4,6) * qJDD(3) - pkin(3) * t156 + qJD(3) * t262 - t260 * t331 - t342;
t317 = mrSges(3,1) * t241 - mrSges(3,2) * t242 + Ifges(3,3) * qJDD(1) + pkin(2) * t315 + pkin(7) * t325 + t304 * t136 + t306 * t138;
t316 = mrSges(2,1) * t287 - mrSges(2,2) * t288 + Ifges(2,3) * qJDD(1) + pkin(1) * t142 + t317;
t140 = m(2) * t288 - mrSges(2,1) * t309 - qJDD(1) * mrSges(2,2) + t326;
t139 = m(2) * t287 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t309 + t142;
t134 = -mrSges(3,1) * t300 + mrSges(3,3) * t242 + t309 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t148 - t341;
t133 = mrSges(3,2) * t300 - mrSges(3,3) * t241 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t309 - pkin(7) * t148 + t136 * t306 - t138 * t304;
t132 = -mrSges(2,2) * g(3) - mrSges(2,3) * t287 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t309 - qJ(2) * t142 + t133 * t302 - t134 * t301;
t131 = Ifges(2,6) * qJDD(1) + t309 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t288 + t301 * t133 + t302 * t134 - pkin(1) * (m(3) * t300 + t148) + qJ(2) * t326;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t307 * t132 - t305 * t131 - pkin(6) * (t139 * t307 + t140 * t305), t132, t133, t136, t149, -t275 * t213 + t312 - t335, t319 - t335; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t305 * t132 + t307 * t131 + pkin(6) * (-t139 * t305 + t140 * t307), t131, t134, t138, t143, -t276 * t210 - t314, -t276 * t209 - t320; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t316, t316, t317, t341, t342, t290 * t216 - Ifges(6,6) * t274 + Ifges(6,3) * t238 + Ifges(6,5) * t239 + t313 + (-t209 + t213) * t276, t321;];
m_new  = t1;
