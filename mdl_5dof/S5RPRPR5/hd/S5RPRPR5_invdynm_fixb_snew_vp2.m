% Calculate vector of cutting torques with Newton-Euler for
% S5RPRPR5
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% m [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% m [3x6]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:26
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPRPR5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR5_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR5_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR5_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR5_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR5_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:25:04
% EndTime: 2022-01-23 09:25:17
% DurationCPUTime: 8.24s
% Computational Cost: add. (113962->288), mult. (293736->378), div. (0->0), fcn. (200705->10), ass. (0->125)
t264 = sin(qJ(1));
t267 = cos(qJ(1));
t244 = -t267 * g(1) - t264 * g(2);
t268 = qJD(1) ^ 2;
t302 = -t268 * pkin(1) + qJDD(1) * qJ(2) + (2 * qJD(1) * qJD(2)) + t244;
t259 = sin(pkin(8));
t261 = cos(pkin(8));
t210 = -t261 * g(3) - t302 * t259;
t297 = qJD(1) * t259;
t295 = t261 * qJD(1);
t211 = -t259 * g(3) + t302 * t261;
t284 = -pkin(2) * t261 - pkin(6) * t259;
t238 = t284 * qJD(1);
t198 = t238 * t295 + t211;
t243 = t264 * g(1) - t267 * g(2);
t277 = -t268 * qJ(2) + qJDD(2) - t243;
t212 = (-pkin(1) + t284) * qJDD(1) + t277;
t266 = cos(qJ(3));
t209 = t266 * t212;
t263 = sin(qJ(3));
t293 = qJD(1) * qJD(3);
t230 = (qJDD(1) * t266 - t263 * t293) * t259;
t292 = t261 * qJDD(1);
t246 = qJDD(3) - t292;
t247 = qJD(3) - t295;
t299 = t259 ^ 2 * t268;
t176 = t246 * pkin(3) - t230 * qJ(4) + t209 + (-pkin(3) * t266 * t299 - qJ(4) * t247 * t297 - t198) * t263;
t186 = t266 * t198 + t263 * t212;
t289 = t266 * t297;
t226 = t247 * pkin(3) - qJ(4) * t289;
t229 = (-qJDD(1) * t263 - t266 * t293) * t259;
t291 = t263 ^ 2 * t299;
t177 = -pkin(3) * t291 + t229 * qJ(4) - t247 * t226 + t186;
t258 = sin(pkin(9));
t260 = cos(pkin(9));
t221 = (-t258 * t263 + t260 * t266) * t297;
t162 = -0.2e1 * qJD(4) * t221 + t260 * t176 - t258 * t177;
t204 = t258 * t229 + t260 * t230;
t220 = (-t258 * t266 - t260 * t263) * t297;
t159 = (t220 * t247 - t204) * pkin(7) + (t220 * t221 + t246) * pkin(4) + t162;
t163 = 0.2e1 * qJD(4) * t220 + t258 * t176 + t260 * t177;
t203 = t260 * t229 - t258 * t230;
t207 = t247 * pkin(4) - t221 * pkin(7);
t219 = t220 ^ 2;
t160 = -t219 * pkin(4) + t203 * pkin(7) - t247 * t207 + t163;
t262 = sin(qJ(5));
t265 = cos(qJ(5));
t158 = t262 * t159 + t265 * t160;
t197 = t238 * t297 - t210;
t184 = -t229 * pkin(3) - qJ(4) * t291 + t226 * t289 + qJDD(4) + t197;
t165 = -t203 * pkin(4) - t219 * pkin(7) + t221 * t207 + t184;
t196 = t262 * t220 + t265 * t221;
t171 = -t196 * qJD(5) + t265 * t203 - t262 * t204;
t195 = t265 * t220 - t262 * t221;
t172 = t195 * qJD(5) + t262 * t203 + t265 * t204;
t245 = qJD(5) + t247;
t178 = Ifges(6,5) * t196 + Ifges(6,6) * t195 + Ifges(6,3) * t245;
t180 = Ifges(6,1) * t196 + Ifges(6,4) * t195 + Ifges(6,5) * t245;
t242 = qJDD(5) + t246;
t147 = -mrSges(6,1) * t165 + mrSges(6,3) * t158 + Ifges(6,4) * t172 + Ifges(6,2) * t171 + Ifges(6,6) * t242 - t196 * t178 + t245 * t180;
t157 = t265 * t159 - t262 * t160;
t179 = Ifges(6,4) * t196 + Ifges(6,2) * t195 + Ifges(6,6) * t245;
t148 = mrSges(6,2) * t165 - mrSges(6,3) * t157 + Ifges(6,1) * t172 + Ifges(6,4) * t171 + Ifges(6,5) * t242 + t195 * t178 - t245 * t179;
t190 = Ifges(5,5) * t221 + Ifges(5,6) * t220 + Ifges(5,3) * t247;
t192 = Ifges(5,1) * t221 + Ifges(5,4) * t220 + Ifges(5,5) * t247;
t188 = -t245 * mrSges(6,2) + t195 * mrSges(6,3);
t189 = t245 * mrSges(6,1) - t196 * mrSges(6,3);
t281 = m(6) * t165 - t171 * mrSges(6,1) + t172 * mrSges(6,2) - t195 * t188 + t196 * t189;
t183 = -t195 * mrSges(6,1) + t196 * mrSges(6,2);
t153 = m(6) * t157 + t242 * mrSges(6,1) - t172 * mrSges(6,3) - t196 * t183 + t245 * t188;
t154 = m(6) * t158 - t242 * mrSges(6,2) + t171 * mrSges(6,3) + t195 * t183 - t245 * t189;
t285 = -t262 * t153 + t265 * t154;
t134 = -mrSges(5,1) * t184 + mrSges(5,3) * t163 + Ifges(5,4) * t204 + Ifges(5,2) * t203 + Ifges(5,6) * t246 - pkin(4) * t281 + pkin(7) * t285 + t265 * t147 + t262 * t148 - t221 * t190 + t247 * t192;
t146 = t265 * t153 + t262 * t154;
t191 = Ifges(5,4) * t221 + Ifges(5,2) * t220 + Ifges(5,6) * t247;
t135 = mrSges(5,2) * t184 - mrSges(5,3) * t162 + Ifges(5,1) * t204 + Ifges(5,4) * t203 + Ifges(5,5) * t246 - pkin(7) * t146 - t262 * t147 + t265 * t148 + t220 * t190 - t247 * t191;
t213 = Ifges(4,3) * t247 + (Ifges(4,5) * t266 - Ifges(4,6) * t263) * t297;
t215 = Ifges(4,5) * t247 + (Ifges(4,1) * t266 - Ifges(4,4) * t263) * t297;
t205 = -t247 * mrSges(5,2) + t220 * mrSges(5,3);
t206 = t247 * mrSges(5,1) - t221 * mrSges(5,3);
t273 = m(5) * t184 - t203 * mrSges(5,1) + t204 * mrSges(5,2) - t220 * t205 + t221 * t206 + t281;
t199 = -t220 * mrSges(5,1) + t221 * mrSges(5,2);
t143 = m(5) * t162 + t246 * mrSges(5,1) - t204 * mrSges(5,3) - t221 * t199 + t247 * t205 + t146;
t144 = m(5) * t163 - t246 * mrSges(5,2) + t203 * mrSges(5,3) + t220 * t199 - t247 * t206 + t285;
t286 = -t258 * t143 + t260 * t144;
t123 = -mrSges(4,1) * t197 + mrSges(4,3) * t186 + Ifges(4,4) * t230 + Ifges(4,2) * t229 + Ifges(4,6) * t246 - pkin(3) * t273 + qJ(4) * t286 + t260 * t134 + t258 * t135 - t213 * t289 + t247 * t215;
t139 = t260 * t143 + t258 * t144;
t185 = -t263 * t198 + t209;
t214 = Ifges(4,6) * t247 + (Ifges(4,4) * t266 - Ifges(4,2) * t263) * t297;
t290 = t263 * t297;
t124 = mrSges(4,2) * t197 - mrSges(4,3) * t185 + Ifges(4,1) * t230 + Ifges(4,4) * t229 + Ifges(4,5) * t246 - qJ(4) * t139 - t258 * t134 + t260 * t135 - t213 * t290 - t247 * t214;
t225 = -t247 * mrSges(4,2) - mrSges(4,3) * t290;
t228 = (mrSges(4,1) * t263 + mrSges(4,2) * t266) * t297;
t137 = m(4) * t185 + t246 * mrSges(4,1) - t230 * mrSges(4,3) + t247 * t225 - t228 * t289 + t139;
t227 = t247 * mrSges(4,1) - mrSges(4,3) * t289;
t138 = m(4) * t186 - t246 * mrSges(4,2) + t229 * mrSges(4,3) - t247 * t227 - t228 * t290 + t286;
t133 = -t263 * t137 + t266 * t138;
t271 = -m(4) * t197 + t229 * mrSges(4,1) - t230 * mrSges(4,2) - t273;
t279 = -t225 * t263 - t227 * t266;
t283 = Ifges(3,1) * t259 + Ifges(3,4) * t261;
t301 = -((Ifges(3,4) * t259 + Ifges(3,2) * t261) * t297 - t283 * t295) * qJD(1) - mrSges(3,1) * t210 + mrSges(3,2) * t211 - pkin(2) * (t279 * t297 + t271) - pkin(6) * t133 - t266 * t123 - t263 * t124;
t300 = mrSges(3,2) * t259;
t296 = qJDD(1) * mrSges(3,3);
t234 = (-mrSges(3,1) * t261 + t300) * qJD(1);
t130 = m(3) * t211 + (qJD(1) * t234 + t296) * t261 + t133;
t152 = t271 + (-t296 + (-t234 + t279) * qJD(1)) * t259 + m(3) * t210;
t287 = t261 * t130 - t259 * t152;
t282 = Ifges(3,5) * t259 + Ifges(3,6) * t261;
t132 = t266 * t137 + t263 * t138;
t280 = t214 * t266 + t215 * t263;
t232 = -qJDD(1) * pkin(1) + t277;
t235 = t282 * qJD(1);
t120 = mrSges(3,2) * t232 - mrSges(3,3) * t210 - pkin(6) * t132 + t283 * qJDD(1) - t263 * t123 + t266 * t124 + t235 * t295;
t275 = -mrSges(6,1) * t157 + mrSges(6,2) * t158 - Ifges(6,5) * t172 - Ifges(6,6) * t171 - Ifges(6,3) * t242 - t196 * t179 + t195 * t180;
t270 = -mrSges(5,1) * t162 + mrSges(5,2) * t163 - Ifges(5,5) * t204 - Ifges(5,6) * t203 - Ifges(5,3) * t246 - pkin(4) * t146 - t221 * t191 + t220 * t192 + t275;
t269 = mrSges(4,1) * t185 - mrSges(4,2) * t186 + Ifges(4,5) * t230 + Ifges(4,6) * t229 + Ifges(4,3) * t246 + pkin(3) * t139 - t270;
t122 = -t269 + (Ifges(3,4) * qJDD(1) + (-t235 - t280) * qJD(1)) * t259 - mrSges(3,1) * t232 + mrSges(3,3) * t211 - pkin(2) * t132 + Ifges(3,2) * t292;
t274 = -m(3) * t232 + mrSges(3,1) * t292 - t132 + (t261 ^ 2 * t268 + t299) * mrSges(3,3);
t276 = -mrSges(2,2) * t244 + qJ(2) * t287 + t259 * t120 + t261 * t122 + pkin(1) * (-qJDD(1) * t300 + t274) + mrSges(2,1) * t243 + Ifges(2,3) * qJDD(1);
t128 = m(2) * t243 - t268 * mrSges(2,2) + (mrSges(2,1) - t300) * qJDD(1) + t274;
t127 = t259 * t130 + t261 * t152;
t125 = m(2) * t244 - t268 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t287;
t118 = mrSges(2,1) * g(3) + mrSges(2,3) * t244 + t268 * Ifges(2,5) - pkin(1) * t127 + (Ifges(2,6) - t282) * qJDD(1) + t301;
t117 = -mrSges(2,2) * g(3) - mrSges(2,3) * t243 + Ifges(2,5) * qJDD(1) - t268 * Ifges(2,6) - qJ(2) * t127 + t261 * t120 - t259 * t122;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t267 * t117 - t264 * t118 - pkin(5) * (t264 * t125 + t267 * t128), t117, t120, t124, t135, t148; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t264 * t117 + t267 * t118 + pkin(5) * (t267 * t125 - t264 * t128), t118, t122, t123, t134, t147; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t276, t276, t282 * qJDD(1) - t301, t280 * t297 + t269, -t270, -t275;];
m_new = t1;
