% Calculate vector of cutting torques with Newton-Euler for
% S5RPPPR1
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% m_mdh [6x1]
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
% Datum: 2020-01-03 11:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPPPR1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR1_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR1_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR1_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR1_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR1_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:20:15
% EndTime: 2020-01-03 11:20:23
% DurationCPUTime: 4.79s
% Computational Cost: add. (53134->235), mult. (123562->328), div. (0->0), fcn. (74901->10), ass. (0->117)
t230 = sin(qJ(1));
t232 = cos(qJ(1));
t208 = -t232 * g(2) - t230 * g(3);
t201 = qJDD(1) * pkin(1) + t208;
t207 = -t230 * g(2) + t232 * g(3);
t233 = qJD(1) ^ 2;
t202 = -t233 * pkin(1) + t207;
t225 = sin(pkin(7));
t228 = cos(pkin(7));
t182 = t225 * t201 + t228 * t202;
t277 = -t233 * pkin(2) + qJDD(1) * qJ(3) + (2 * qJD(1) * qJD(3)) + t182;
t181 = t228 * t201 - t225 * t202;
t244 = -t233 * qJ(3) + qJDD(3) - t181;
t224 = sin(pkin(8));
t227 = cos(pkin(8));
t253 = -pkin(3) * t227 - qJ(4) * t224;
t269 = qJD(1) * t224;
t276 = (-pkin(2) + t253) * qJDD(1) + t244 - 0.2e1 * qJD(4) * t269;
t222 = -g(1) + qJDD(2);
t164 = t227 * t222 - t277 * t224;
t268 = t227 * qJD(1);
t165 = t224 * t222 + t277 * t227;
t196 = t253 * qJD(1);
t158 = t196 * t268 + t165;
t219 = t224 ^ 2;
t223 = sin(pkin(9));
t226 = cos(pkin(9));
t271 = t224 * t226;
t249 = -pkin(4) * t227 - pkin(6) * t271;
t270 = t276 * t226;
t150 = t249 * qJDD(1) + (-t158 + (-pkin(4) * t219 * t226 + pkin(6) * t224 * t227) * t233) * t223 + t270;
t153 = t226 * t158 + t276 * t223;
t195 = t249 * qJD(1);
t273 = t219 * t233;
t264 = t223 ^ 2 * t273;
t266 = qJDD(1) * t224;
t151 = -t223 * pkin(6) * t266 - pkin(4) * t264 + t195 * t268 + t153;
t229 = sin(qJ(5));
t231 = cos(qJ(5));
t149 = t229 * t150 + t231 * t151;
t157 = t196 * t269 + qJDD(4) - t164;
t154 = -pkin(6) * t264 + (pkin(4) * qJDD(1) * t223 + qJD(1) * t195 * t226) * t224 + t157;
t246 = (-t223 * t231 - t226 * t229) * t224;
t186 = qJD(1) * t246;
t245 = (-t223 * t229 + t226 * t231) * t224;
t187 = qJD(1) * t245;
t210 = qJD(5) - t268;
t159 = Ifges(6,5) * t187 + Ifges(6,6) * t186 + Ifges(6,3) * t210;
t161 = Ifges(6,1) * t187 + Ifges(6,4) * t186 + Ifges(6,5) * t210;
t172 = -t187 * qJD(5) + qJDD(1) * t246;
t173 = t186 * qJD(5) + qJDD(1) * t245;
t265 = t227 * qJDD(1);
t209 = qJDD(5) - t265;
t137 = -mrSges(6,1) * t154 + mrSges(6,3) * t149 + Ifges(6,4) * t173 + Ifges(6,2) * t172 + Ifges(6,6) * t209 - t187 * t159 + t210 * t161;
t148 = t231 * t150 - t229 * t151;
t160 = Ifges(6,4) * t187 + Ifges(6,2) * t186 + Ifges(6,6) * t210;
t138 = mrSges(6,2) * t154 - mrSges(6,3) * t148 + Ifges(6,1) * t173 + Ifges(6,4) * t172 + Ifges(6,5) * t209 + t186 * t159 - t210 * t160;
t254 = Ifges(5,5) * t226 - Ifges(5,6) * t223;
t183 = (-Ifges(5,3) * t227 + t254 * t224) * qJD(1);
t241 = -Ifges(5,5) * t227 + (Ifges(5,1) * t226 - Ifges(5,4) * t223) * t224;
t185 = t241 * qJD(1);
t240 = -Ifges(5,6) * t227 + (Ifges(5,4) * t226 - Ifges(5,2) * t223) * t224;
t179 = -t210 * mrSges(6,2) + t186 * mrSges(6,3);
t180 = t210 * mrSges(6,1) - t187 * mrSges(6,3);
t242 = m(6) * t154 - t172 * mrSges(6,1) + t173 * mrSges(6,2) - t186 * t179 + t187 * t180;
t168 = -t186 * mrSges(6,1) + t187 * mrSges(6,2);
t144 = m(6) * t148 + t209 * mrSges(6,1) - t173 * mrSges(6,3) - t187 * t168 + t210 * t179;
t145 = m(6) * t149 - t209 * mrSges(6,2) + t172 * mrSges(6,3) + t186 * t168 - t210 * t180;
t260 = -t229 * t144 + t231 * t145;
t123 = -mrSges(5,1) * t157 + mrSges(5,3) * t153 + t229 * t138 + t231 * t137 - pkin(4) * t242 + pkin(6) * t260 + (-t183 * t271 - t227 * t185) * qJD(1) + t240 * qJDD(1);
t136 = t231 * t144 + t229 * t145;
t152 = -t223 * t158 + t270;
t184 = t240 * qJD(1);
t272 = t223 * t224;
t124 = mrSges(5,2) * t157 - mrSges(5,3) * t152 - pkin(6) * t136 - t229 * t137 + t231 * t138 + (-t183 * t272 + t184 * t227) * qJD(1) + t241 * qJDD(1);
t257 = mrSges(5,1) * t223 + mrSges(5,2) * t226;
t188 = t257 * t269;
t247 = mrSges(5,2) * t227 - mrSges(5,3) * t272;
t190 = t247 * qJD(1);
t248 = -mrSges(5,1) * t227 - mrSges(5,3) * t271;
t134 = m(5) * t152 + t248 * qJDD(1) + (-t188 * t271 - t190 * t227) * qJD(1) + t136;
t191 = t248 * qJD(1);
t135 = m(5) * t153 + t247 * qJDD(1) + (-t188 * t272 + t191 * t227) * qJD(1) + t260;
t132 = -t223 * t134 + t226 * t135;
t237 = -m(5) * t157 - t242;
t251 = -t190 * t223 - t191 * t226;
t256 = Ifges(4,1) * t224 + Ifges(4,4) * t227;
t275 = -((Ifges(4,4) * t224 + Ifges(4,2) * t227) * t269 - t256 * t268) * qJD(1) - mrSges(4,1) * t164 + mrSges(4,2) * t165 - pkin(3) * ((t251 * qJD(1) - t257 * qJDD(1)) * t224 + t237) - qJ(4) * t132 - t226 * t123 - t223 * t124;
t274 = mrSges(4,2) * t224;
t197 = (-mrSges(4,1) * t227 + t274) * qJD(1);
t129 = m(4) * t165 + (qJDD(1) * mrSges(4,3) + qJD(1) * t197) * t227 + t132;
t140 = m(4) * t164 + ((-mrSges(4,3) - t257) * qJDD(1) + (-t197 + t251) * qJD(1)) * t224 + t237;
t261 = t227 * t129 - t224 * t140;
t120 = m(3) * t182 - t233 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t261;
t131 = t226 * t134 + t223 * t135;
t176 = -qJDD(1) * pkin(2) + t244;
t238 = -m(4) * t176 + mrSges(4,1) * t265 - t131 + (t227 ^ 2 * t233 + t273) * mrSges(4,3);
t126 = m(3) * t181 - t233 * mrSges(3,2) + (mrSges(3,1) - t274) * qJDD(1) + t238;
t115 = t225 * t120 + t228 * t126;
t122 = t224 * t129 + t227 * t140;
t262 = t228 * t120 - t225 * t126;
t255 = Ifges(4,5) * t224 + Ifges(4,6) * t227;
t252 = t184 * t226 + t185 * t223;
t198 = t255 * qJD(1);
t111 = mrSges(4,2) * t176 - mrSges(4,3) * t164 - qJ(4) * t131 + t256 * qJDD(1) - t223 * t123 + t226 * t124 + t198 * t268;
t239 = -mrSges(6,1) * t148 + mrSges(6,2) * t149 - Ifges(6,5) * t173 - Ifges(6,6) * t172 - Ifges(6,3) * t209 - t187 * t160 + t186 * t161;
t234 = mrSges(5,1) * t152 - mrSges(5,2) * t153 + pkin(4) * t136 - t239;
t117 = (Ifges(4,2) + Ifges(5,3)) * t265 - t234 + ((Ifges(4,4) - t254) * qJDD(1) + (-t198 - t252) * qJD(1)) * t224 - mrSges(4,1) * t176 + mrSges(4,3) * t165 - pkin(3) * t131;
t243 = -mrSges(3,2) * t182 + qJ(3) * t261 + t224 * t111 + t227 * t117 + pkin(2) * (-mrSges(4,2) * t266 + t238) + mrSges(3,1) * t181 + Ifges(3,3) * qJDD(1);
t236 = mrSges(2,1) * t208 - mrSges(2,2) * t207 + Ifges(2,3) * qJDD(1) + pkin(1) * t115 + t243;
t113 = m(2) * t208 + qJDD(1) * mrSges(2,1) - t233 * mrSges(2,2) + t115;
t112 = m(2) * t207 - t233 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t262;
t109 = -mrSges(3,1) * t222 + mrSges(3,3) * t182 + t233 * Ifges(3,5) - pkin(2) * t122 + (Ifges(3,6) - t255) * qJDD(1) + t275;
t108 = mrSges(3,2) * t222 - mrSges(3,3) * t181 + Ifges(3,5) * qJDD(1) - t233 * Ifges(3,6) - qJ(3) * t122 + t227 * t111 - t224 * t117;
t107 = -mrSges(2,2) * g(1) - mrSges(2,3) * t208 + Ifges(2,5) * qJDD(1) - t233 * Ifges(2,6) - qJ(2) * t115 + t228 * t108 - t225 * t109;
t106 = Ifges(2,6) * qJDD(1) + t233 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t207 + t225 * t108 + t228 * t109 - pkin(1) * (m(3) * t222 + t122) + qJ(2) * t262;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t236, t107, t108, t111, t124, t138; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t230 * t107 + t232 * t106 - pkin(5) * (-t232 * t112 + t230 * t113), t106, t109, t117, t123, t137; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) - t232 * t107 + t230 * t106 + pkin(5) * (t230 * t112 + t232 * t113), t236, t243, t255 * qJDD(1) - t275, -Ifges(5,3) * t265 + (t252 * qJD(1) + t254 * qJDD(1)) * t224 + t234, -t239;];
m_new = t1;
