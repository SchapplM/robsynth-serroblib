% Calculate vector of cutting torques with Newton-Euler for
% S6PPRRRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-05-04 20:37
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PPRRRP2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP2_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRP2_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRRP2_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRP2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP2_invdynm_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRP2_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRP2_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRRP2_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 20:32:07
% EndTime: 2019-05-04 20:32:26
% DurationCPUTime: 14.78s
% Computational Cost: add. (288069->287), mult. (514411->364), div. (0->0), fcn. (392106->14), ass. (0->125)
t258 = sin(pkin(11));
t262 = cos(pkin(11));
t249 = -t262 * g(1) - t258 * g(2);
t257 = sin(pkin(12));
t261 = cos(pkin(12));
t248 = t258 * g(1) - t262 * g(2);
t256 = -g(3) + qJDD(1);
t260 = sin(pkin(6));
t264 = cos(pkin(6));
t279 = t248 * t264 + t256 * t260;
t201 = -t257 * t249 + t261 * t279;
t202 = t261 * t249 + t257 * t279;
t229 = -t260 * t248 + t264 * t256 + qJDD(2);
t269 = cos(qJ(3));
t263 = cos(pkin(7));
t267 = sin(qJ(3));
t294 = t263 * t267;
t259 = sin(pkin(7));
t295 = t259 * t267;
t193 = t201 * t294 + t269 * t202 + t229 * t295;
t271 = qJD(3) ^ 2;
t191 = -t271 * pkin(3) + qJDD(3) * pkin(9) + t193;
t195 = -t259 * t201 + t263 * t229;
t266 = sin(qJ(4));
t268 = cos(qJ(4));
t185 = t268 * t191 + t266 * t195;
t244 = (-pkin(4) * t268 - pkin(10) * t266) * qJD(3);
t270 = qJD(4) ^ 2;
t289 = t268 * qJD(3);
t183 = -t270 * pkin(4) + qJDD(4) * pkin(10) + t244 * t289 + t185;
t192 = -t267 * t202 + (t201 * t263 + t229 * t259) * t269;
t190 = -qJDD(3) * pkin(3) - t271 * pkin(9) - t192;
t288 = qJD(3) * qJD(4);
t285 = t268 * t288;
t245 = t266 * qJDD(3) + t285;
t286 = t266 * t288;
t246 = t268 * qJDD(3) - t286;
t187 = (-t245 - t285) * pkin(10) + (-t246 + t286) * pkin(4) + t190;
t265 = sin(qJ(5));
t303 = cos(qJ(5));
t179 = t183 * t303 + t265 * t187;
t290 = qJD(3) * t266;
t242 = t265 * qJD(4) + t290 * t303;
t217 = t242 * qJD(5) - qJDD(4) * t303 + t265 * t245;
t253 = qJD(5) - t289;
t226 = t253 * mrSges(6,1) - t242 * mrSges(6,3);
t240 = qJDD(5) - t246;
t241 = -qJD(4) * t303 + t265 * t290;
t221 = t241 * pkin(5) - t242 * qJ(6);
t252 = t253 ^ 2;
t175 = -pkin(5) * t252 + qJ(6) * t240 + 0.2e1 * qJD(6) * t253 - t221 * t241 + t179;
t227 = -t253 * mrSges(7,1) + t242 * mrSges(7,2);
t287 = m(7) * t175 + t240 * mrSges(7,3) + t253 * t227;
t222 = t241 * mrSges(7,1) - t242 * mrSges(7,3);
t291 = -t241 * mrSges(6,1) - t242 * mrSges(6,2) - t222;
t302 = -mrSges(6,3) - mrSges(7,2);
t167 = m(6) * t179 - t240 * mrSges(6,2) + t217 * t302 - t253 * t226 + t241 * t291 + t287;
t178 = -t265 * t183 + t187 * t303;
t218 = -t241 * qJD(5) + t265 * qJDD(4) + t245 * t303;
t225 = -t253 * mrSges(6,2) - t241 * mrSges(6,3);
t177 = -t240 * pkin(5) - t252 * qJ(6) + t242 * t221 + qJDD(6) - t178;
t228 = -t241 * mrSges(7,2) + t253 * mrSges(7,3);
t283 = -m(7) * t177 + t240 * mrSges(7,1) + t253 * t228;
t168 = m(6) * t178 + t240 * mrSges(6,1) + t218 * t302 + t253 * t225 + t242 * t291 + t283;
t163 = t167 * t303 - t168 * t265;
t243 = (-mrSges(5,1) * t268 + mrSges(5,2) * t266) * qJD(3);
t250 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t290;
t159 = m(5) * t185 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t246 - qJD(4) * t250 + t243 * t289 + t163;
t184 = -t266 * t191 + t268 * t195;
t182 = -qJDD(4) * pkin(4) - t270 * pkin(10) + t244 * t290 - t184;
t180 = -0.2e1 * qJD(6) * t242 + (t241 * t253 - t218) * qJ(6) + (t242 * t253 + t217) * pkin(5) + t182;
t172 = m(7) * t180 + mrSges(7,1) * t217 - t218 * mrSges(7,3) - t242 * t227 + t228 * t241;
t169 = -m(6) * t182 - t217 * mrSges(6,1) - mrSges(6,2) * t218 - t241 * t225 - t226 * t242 - t172;
t251 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t289;
t165 = m(5) * t184 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t245 + qJD(4) * t251 - t243 * t290 + t169;
t284 = t268 * t159 - t165 * t266;
t150 = m(4) * t193 - mrSges(4,1) * t271 - qJDD(3) * mrSges(4,2) + t284;
t153 = t266 * t159 + t268 * t165;
t152 = m(4) * t195 + t153;
t162 = t265 * t167 + t168 * t303;
t274 = -m(5) * t190 + t246 * mrSges(5,1) - t245 * mrSges(5,2) - t250 * t290 + t251 * t289 - t162;
t156 = m(4) * t192 + qJDD(3) * mrSges(4,1) - t271 * mrSges(4,2) + t274;
t296 = t156 * t269;
t140 = t150 * t294 - t152 * t259 + t263 * t296;
t137 = m(3) * t201 + t140;
t145 = t269 * t150 - t156 * t267;
t144 = m(3) * t202 + t145;
t307 = t137 * t261 + t144 * t257;
t207 = Ifges(7,1) * t242 + Ifges(7,4) * t253 + Ifges(7,5) * t241;
t208 = Ifges(6,1) * t242 - Ifges(6,4) * t241 + Ifges(6,5) * t253;
t282 = -mrSges(7,1) * t180 + mrSges(7,2) * t175;
t205 = Ifges(7,4) * t242 + Ifges(7,2) * t253 + Ifges(7,6) * t241;
t293 = -Ifges(6,5) * t242 + Ifges(6,6) * t241 - Ifges(6,3) * t253 - t205;
t160 = -mrSges(6,1) * t182 + mrSges(6,3) * t179 - pkin(5) * t172 + (t207 + t208) * t253 + t293 * t242 + (Ifges(6,6) - Ifges(7,6)) * t240 + (Ifges(6,4) - Ifges(7,5)) * t218 + (-Ifges(6,2) - Ifges(7,3)) * t217 + t282;
t206 = Ifges(6,4) * t242 - Ifges(6,2) * t241 + Ifges(6,6) * t253;
t203 = Ifges(7,5) * t242 + Ifges(7,6) * t253 + Ifges(7,3) * t241;
t276 = mrSges(7,2) * t177 - mrSges(7,3) * t180 + Ifges(7,1) * t218 + Ifges(7,4) * t240 + Ifges(7,5) * t217 + t253 * t203;
t161 = mrSges(6,2) * t182 - mrSges(6,3) * t178 + Ifges(6,1) * t218 - Ifges(6,4) * t217 + Ifges(6,5) * t240 - qJ(6) * t172 - t253 * t206 + t241 * t293 + t276;
t231 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t266 + Ifges(5,6) * t268) * qJD(3);
t232 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t266 + Ifges(5,2) * t268) * qJD(3);
t141 = mrSges(5,2) * t190 - mrSges(5,3) * t184 + Ifges(5,1) * t245 + Ifges(5,4) * t246 + Ifges(5,5) * qJDD(4) - pkin(10) * t162 - qJD(4) * t232 - t265 * t160 + t161 * t303 + t231 * t289;
t233 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t266 + Ifges(5,4) * t268) * qJD(3);
t277 = mrSges(7,1) * t177 - mrSges(7,3) * t175 - Ifges(7,4) * t218 - Ifges(7,2) * t240 - Ifges(7,6) * t217 - t241 * t207;
t305 = -(-t206 + t203) * t242 + mrSges(6,1) * t178 - mrSges(6,2) * t179 + Ifges(6,5) * t218 - Ifges(6,6) * t217 + Ifges(6,3) * t240 + pkin(5) * (-t218 * mrSges(7,2) - t242 * t222 + t283) + qJ(6) * (-t217 * mrSges(7,2) - t241 * t222 + t287) + t241 * t208 - t277;
t146 = -mrSges(5,1) * t190 + mrSges(5,3) * t185 + Ifges(5,4) * t245 + Ifges(5,2) * t246 + Ifges(5,6) * qJDD(4) - pkin(4) * t162 + qJD(4) * t233 - t231 * t290 - t305;
t130 = mrSges(4,1) * t192 - mrSges(4,2) * t193 + Ifges(4,3) * qJDD(3) + pkin(3) * t274 + pkin(9) * t284 + t266 * t141 + t268 * t146;
t139 = t150 * t295 + t263 * t152 + t259 * t296;
t131 = mrSges(4,2) * t195 - mrSges(4,3) * t192 + Ifges(4,5) * qJDD(3) - Ifges(4,6) * t271 - pkin(9) * t153 + t141 * t268 - t146 * t266;
t304 = mrSges(5,1) * t184 - mrSges(5,2) * t185 + Ifges(5,5) * t245 + Ifges(5,6) * t246 + Ifges(5,3) * qJDD(4) + pkin(4) * t169 + pkin(10) * t163 + (t232 * t266 - t233 * t268) * qJD(3) + t160 * t303 + t265 * t161;
t135 = -mrSges(4,1) * t195 + mrSges(4,3) * t193 + t271 * Ifges(4,5) + Ifges(4,6) * qJDD(3) - pkin(3) * t153 - t304;
t278 = pkin(8) * t145 + t131 * t267 + t135 * t269;
t123 = -mrSges(3,1) * t229 + mrSges(3,3) * t202 - pkin(2) * t139 - t259 * t130 + t263 * t278;
t125 = mrSges(3,2) * t229 - mrSges(3,3) * t201 + t269 * t131 - t267 * t135 + (-t139 * t259 - t140 * t263) * pkin(8);
t134 = -t137 * t257 + t261 * t144;
t306 = qJ(2) * t134 + t123 * t261 + t125 * t257;
t138 = m(3) * t229 + t139;
t129 = -t260 * t138 + t307 * t264;
t121 = mrSges(3,1) * t201 - mrSges(3,2) * t202 + pkin(2) * t140 + t263 * t130 + t259 * t278;
t275 = mrSges(2,1) * t248 - mrSges(2,2) * t249 + pkin(1) * t129 + t264 * t121 + t306 * t260;
t132 = m(2) * t249 + t134;
t128 = t264 * t138 + t307 * t260;
t126 = m(2) * t248 + t129;
t119 = mrSges(2,2) * t256 - mrSges(2,3) * t248 - t257 * t123 + t261 * t125 + (-t128 * t260 - t129 * t264) * qJ(2);
t118 = -mrSges(2,1) * t256 + mrSges(2,3) * t249 - pkin(1) * t128 - t260 * t121 + t306 * t264;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t262 * t119 - t258 * t118 - qJ(1) * (t126 * t262 + t132 * t258), t119, t125, t131, t141, t161, -t205 * t241 + t276; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t258 * t119 + t262 * t118 + qJ(1) * (-t126 * t258 + t132 * t262), t118, t123, t135, t146, t160, -t242 * t203 - t277; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t275, t275, t121, t130, t304, t305, Ifges(7,5) * t218 + Ifges(7,6) * t240 + Ifges(7,3) * t217 + t242 * t205 - t253 * t207 - t282;];
m_new  = t1;
