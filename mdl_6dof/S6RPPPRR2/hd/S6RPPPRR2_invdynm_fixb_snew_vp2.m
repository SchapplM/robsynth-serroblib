% Calculate vector of cutting torques with Newton-Euler for
% S6RPPPRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2,theta4]';
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
% Datum: 2019-05-05 13:37
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPPPRR2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR2_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR2_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPPRR2_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR2_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR2_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR2_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR2_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 13:35:52
% EndTime: 2019-05-05 13:36:01
% DurationCPUTime: 6.62s
% Computational Cost: add. (91248->277), mult. (190042->335), div. (0->0), fcn. (120917->10), ass. (0->125)
t269 = qJD(1) ^ 2;
t258 = sin(pkin(10));
t247 = t258 ^ 2;
t260 = cos(pkin(10));
t303 = t260 ^ 2 + t247;
t293 = t303 * mrSges(5,3);
t312 = t269 * t293;
t264 = sin(qJ(1));
t267 = cos(qJ(1));
t235 = t264 * g(1) - g(2) * t267;
t232 = qJDD(1) * pkin(1) + t235;
t236 = -g(1) * t267 - g(2) * t264;
t233 = -pkin(1) * t269 + t236;
t259 = sin(pkin(9));
t261 = cos(pkin(9));
t214 = t261 * t232 - t259 * t233;
t283 = -t269 * qJ(3) + qJDD(3) - t214;
t306 = -pkin(2) - qJ(4);
t311 = -(2 * qJD(1) * qJD(4)) + t306 * qJDD(1) + t283;
t215 = t259 * t232 + t261 * t233;
t310 = qJDD(1) * qJ(3) + (2 * qJD(3) * qJD(1)) + t215;
t263 = sin(qJ(5));
t266 = cos(qJ(5));
t287 = t258 * t266 + t260 * t263;
t223 = t287 * qJD(1);
t286 = -t258 * t263 + t260 * t266;
t224 = t286 * qJD(1);
t300 = t224 * qJD(5);
t211 = -t287 * qJDD(1) - t300;
t309 = pkin(4) * t269;
t308 = mrSges(3,1) - mrSges(4,2);
t307 = -Ifges(3,6) + Ifges(4,5);
t305 = Ifges(5,6) * t258;
t255 = -g(3) + qJDD(2);
t297 = qJDD(1) * t260;
t304 = t311 * t260;
t175 = -pkin(7) * t297 + (-t260 * t309 - t255) * t258 + t304;
t188 = t260 * t255 + t311 * t258;
t298 = qJDD(1) * t258;
t176 = -pkin(7) * t298 - t247 * t309 + t188;
t172 = t263 * t175 + t266 * t176;
t205 = mrSges(6,1) * t223 + mrSges(6,2) * t224;
t219 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t224;
t210 = pkin(5) * t223 - pkin(8) * t224;
t268 = qJD(5) ^ 2;
t168 = -pkin(5) * t268 + qJDD(5) * pkin(8) - t210 * t223 + t172;
t285 = qJDD(4) + t310;
t185 = pkin(4) * t298 + (-t303 * pkin(7) + t306) * t269 + t285;
t301 = t223 * qJD(5);
t212 = t286 * qJDD(1) - t301;
t169 = (-t212 + t301) * pkin(8) + (-t211 + t300) * pkin(5) + t185;
t262 = sin(qJ(6));
t265 = cos(qJ(6));
t165 = -t168 * t262 + t169 * t265;
t216 = qJD(5) * t265 - t224 * t262;
t183 = qJD(6) * t216 + qJDD(5) * t262 + t212 * t265;
t217 = qJD(5) * t262 + t224 * t265;
t190 = -mrSges(7,1) * t216 + mrSges(7,2) * t217;
t222 = qJD(6) + t223;
t195 = -mrSges(7,2) * t222 + mrSges(7,3) * t216;
t209 = qJDD(6) - t211;
t161 = m(7) * t165 + mrSges(7,1) * t209 - t183 * mrSges(7,3) - t190 * t217 + t195 * t222;
t166 = t168 * t265 + t169 * t262;
t182 = -qJD(6) * t217 + qJDD(5) * t265 - t212 * t262;
t196 = mrSges(7,1) * t222 - mrSges(7,3) * t217;
t162 = m(7) * t166 - mrSges(7,2) * t209 + t182 * mrSges(7,3) + t190 * t216 - t196 * t222;
t290 = -t161 * t262 + t265 * t162;
t148 = m(6) * t172 - qJDD(5) * mrSges(6,2) + mrSges(6,3) * t211 - qJD(5) * t219 - t205 * t223 + t290;
t171 = t175 * t266 - t176 * t263;
t218 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t223;
t167 = -qJDD(5) * pkin(5) - pkin(8) * t268 + t210 * t224 - t171;
t281 = -m(7) * t167 + t182 * mrSges(7,1) - t183 * mrSges(7,2) + t216 * t195 - t196 * t217;
t157 = m(6) * t171 + qJDD(5) * mrSges(6,1) - mrSges(6,3) * t212 + qJD(5) * t218 - t205 * t224 + t281;
t140 = t263 * t148 + t266 * t157;
t187 = -t258 * t255 + t304;
t284 = -qJDD(1) * mrSges(5,3) - t269 * (mrSges(5,1) * t258 + mrSges(5,2) * t260);
t137 = m(5) * t187 + t284 * t260 + t140;
t291 = t266 * t148 - t263 * t157;
t138 = m(5) * t188 + t284 * t258 + t291;
t132 = t260 * t137 + t258 * t138;
t200 = -qJDD(1) * pkin(2) + t283;
t282 = -m(4) * t200 + t269 * mrSges(4,3) - t132;
t129 = m(3) * t214 - t269 * mrSges(3,2) + t308 * qJDD(1) + t282;
t198 = t269 * pkin(2) - t310;
t194 = t306 * t269 + t285;
t150 = t265 * t161 + t262 * t162;
t280 = -m(6) * t185 + t211 * mrSges(6,1) - t212 * mrSges(6,2) - t223 * t218 - t224 * t219 - t150;
t277 = -m(5) * t194 - mrSges(5,1) * t298 - mrSges(5,2) * t297 + t280;
t273 = -m(4) * t198 + t269 * mrSges(4,2) + qJDD(1) * mrSges(4,3) - t277;
t143 = t273 - qJDD(1) * mrSges(3,2) + m(3) * t215 + (-mrSges(3,1) - t293) * t269;
t127 = t261 * t129 + t259 * t143;
t302 = t269 * (Ifges(5,5) * t260 - t305);
t294 = Ifges(4,4) + t305;
t292 = -t129 * t259 + t261 * t143;
t133 = -t137 * t258 + t260 * t138;
t131 = m(4) * t255 + t133;
t289 = Ifges(5,1) * t260 - Ifges(5,4) * t258;
t288 = Ifges(5,4) * t260 - Ifges(5,2) * t258;
t177 = Ifges(7,5) * t217 + Ifges(7,6) * t216 + Ifges(7,3) * t222;
t179 = Ifges(7,1) * t217 + Ifges(7,4) * t216 + Ifges(7,5) * t222;
t154 = -mrSges(7,1) * t167 + mrSges(7,3) * t166 + Ifges(7,4) * t183 + Ifges(7,2) * t182 + Ifges(7,6) * t209 - t177 * t217 + t179 * t222;
t178 = Ifges(7,4) * t217 + Ifges(7,2) * t216 + Ifges(7,6) * t222;
t155 = mrSges(7,2) * t167 - mrSges(7,3) * t165 + Ifges(7,1) * t183 + Ifges(7,4) * t182 + Ifges(7,5) * t209 + t177 * t216 - t178 * t222;
t201 = Ifges(6,5) * t224 - Ifges(6,6) * t223 + Ifges(6,3) * qJD(5);
t202 = Ifges(6,4) * t224 - Ifges(6,2) * t223 + Ifges(6,6) * qJD(5);
t134 = mrSges(6,2) * t185 - mrSges(6,3) * t171 + Ifges(6,1) * t212 + Ifges(6,4) * t211 + Ifges(6,5) * qJDD(5) - pkin(8) * t150 - qJD(5) * t202 - t154 * t262 + t155 * t265 - t201 * t223;
t203 = Ifges(6,1) * t224 - Ifges(6,4) * t223 + Ifges(6,5) * qJD(5);
t274 = mrSges(7,1) * t165 - mrSges(7,2) * t166 + Ifges(7,5) * t183 + Ifges(7,6) * t182 + Ifges(7,3) * t209 + t178 * t217 - t179 * t216;
t135 = -mrSges(6,1) * t185 + mrSges(6,3) * t172 + Ifges(6,4) * t212 + Ifges(6,2) * t211 + Ifges(6,6) * qJDD(5) - pkin(5) * t150 + qJD(5) * t203 - t201 * t224 - t274;
t121 = -mrSges(5,1) * t194 + mrSges(5,3) * t188 + pkin(4) * t280 + pkin(7) * t291 + t288 * qJDD(1) + t263 * t134 + t266 * t135 - t260 * t302;
t123 = mrSges(5,2) * t194 - mrSges(5,3) * t187 - pkin(7) * t140 + t289 * qJDD(1) + t266 * t134 - t263 * t135 - t258 * t302;
t279 = mrSges(4,2) * t200 - mrSges(4,3) * t198 + Ifges(4,1) * qJDD(1) - qJ(4) * t132 - t121 * t258 + t260 * t123;
t278 = -mrSges(4,1) * t198 - pkin(3) * (t277 + t312) - qJ(4) * t133 - t260 * t121 - t258 * t123;
t276 = -mrSges(6,1) * t171 + mrSges(6,2) * t172 - Ifges(6,5) * t212 - Ifges(6,6) * t211 - Ifges(6,3) * qJDD(5) - pkin(5) * t281 - pkin(8) * t290 - t265 * t154 - t262 * t155 - t224 * t202 - t223 * t203;
t275 = -mrSges(3,2) * t215 + pkin(2) * (-qJDD(1) * mrSges(4,2) + t282) + qJ(3) * (t273 - t312) + mrSges(3,1) * t214 + Ifges(3,3) * qJDD(1) + t279;
t272 = -mrSges(5,1) * t187 + mrSges(5,2) * t188 - Ifges(5,5) * t297 - pkin(4) * t140 + t276 + (-t258 * t289 - t260 * t288) * t269;
t271 = mrSges(2,1) * t235 - mrSges(2,2) * t236 + Ifges(2,3) * qJDD(1) + pkin(1) * t127 + t275;
t270 = -mrSges(4,1) * t200 - pkin(3) * t132 + t272;
t125 = m(2) * t236 - mrSges(2,1) * t269 - qJDD(1) * mrSges(2,2) + t292;
t124 = m(2) * t235 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t269 + t127;
t120 = -t270 - mrSges(3,3) * t214 - qJ(3) * t131 + (Ifges(3,5) - t294) * qJDD(1) + (mrSges(3,2) - mrSges(4,3)) * t255 + t307 * t269;
t119 = mrSges(3,3) * t215 - pkin(2) * t131 + (-Ifges(4,4) + Ifges(3,5)) * t269 - t308 * t255 - t307 * qJDD(1) + t278;
t118 = -mrSges(2,2) * g(3) - mrSges(2,3) * t235 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t269 - qJ(2) * t127 - t119 * t259 + t120 * t261;
t117 = Ifges(2,6) * qJDD(1) + t269 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t236 + t259 * t120 + t261 * t119 - pkin(1) * (m(3) * t255 + t131) + qJ(2) * t292;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t267 * t118 - t264 * t117 - pkin(6) * (t124 * t267 + t125 * t264), t118, t120, t279, t123, t134, t155; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t264 * t118 + t267 * t117 + pkin(6) * (-t124 * t264 + t125 * t267), t117, t119, mrSges(4,3) * t255 - t269 * Ifges(4,5) + t294 * qJDD(1) + t270, t121, t135, t154; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t271, t271, t275, -mrSges(4,2) * t255 + t269 * Ifges(4,4) + Ifges(4,5) * qJDD(1) - t278, -Ifges(5,6) * t298 - t272, -t276, t274;];
m_new  = t1;
