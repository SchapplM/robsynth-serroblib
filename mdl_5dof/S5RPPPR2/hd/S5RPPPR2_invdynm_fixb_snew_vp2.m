% Calculate vector of cutting torques with Newton-Euler for
% S5RPPPR2
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
% Datum: 2020-01-03 11:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPPPR2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR2_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR2_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR2_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR2_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:22:28
% EndTime: 2020-01-03 11:22:41
% DurationCPUTime: 6.90s
% Computational Cost: add. (64644->281), mult. (181142->386), div. (0->0), fcn. (122199->10), ass. (0->133)
t247 = sin(pkin(7));
t250 = cos(pkin(7));
t252 = sin(qJ(1));
t254 = cos(qJ(1));
t228 = -g(2) * t252 + t254 * g(3);
t255 = qJD(1) ^ 2;
t302 = -pkin(1) * t255 + qJDD(1) * qJ(2) + (2 * qJD(1) * qJD(2)) + t228;
t197 = -g(1) * t247 + t302 * t250;
t272 = -pkin(2) * t250 - qJ(3) * t247;
t220 = t272 * qJD(1);
t288 = qJD(1) * t250;
t188 = t220 * t288 + t197;
t246 = sin(pkin(8));
t249 = cos(pkin(8));
t229 = -t254 * g(2) - t252 * g(3);
t265 = -qJ(2) * t255 + qJDD(2) - t229;
t289 = qJD(1) * t247;
t301 = (-pkin(1) + t272) * qJDD(1) + t265 - 0.2e1 * qJD(3) * t289;
t168 = -t246 * t188 + t301 * t249;
t196 = -t250 * g(1) - t302 * t247;
t245 = sin(pkin(9));
t248 = cos(pkin(9));
t292 = t247 * t249;
t264 = t245 * t292 + t248 * t250;
t209 = t264 * qJD(1);
t207 = t264 * qJDD(1);
t169 = t249 * t188 + t301 * t246;
t213 = (pkin(3) * t246 - qJ(4) * t249) * t289;
t283 = t246 * t289;
t285 = qJDD(1) * t250;
t294 = t250 ^ 2 * t255;
t166 = -pkin(3) * t294 - qJ(4) * t285 - t213 * t283 + t169;
t187 = t220 * t289 + qJDD(3) - t196;
t291 = t250 * t255;
t174 = ((-qJDD(1) * t249 - t246 * t291) * qJ(4) + (qJDD(1) * t246 - t249 * t291) * pkin(3)) * t247 + t187;
t299 = 2 * qJD(4);
t162 = t248 * t166 + t245 * t174 - t209 * t299;
t282 = t249 * t289;
t210 = -t245 * t288 + t248 * t282;
t190 = pkin(4) * t209 - pkin(6) * t210;
t286 = qJDD(1) * t247;
t279 = t246 * t286;
t295 = t247 ^ 2 * t255;
t284 = t246 ^ 2 * t295;
t160 = -pkin(4) * t284 + pkin(6) * t279 - t190 * t209 + t162;
t165 = pkin(3) * t285 - qJ(4) * t294 + t213 * t282 + qJDD(4) - t168;
t208 = (-t245 * t250 + t248 * t292) * qJDD(1);
t163 = (t209 * t283 - t208) * pkin(6) + (t210 * t283 + t207) * pkin(4) + t165;
t251 = sin(qJ(5));
t253 = cos(qJ(5));
t157 = -t160 * t251 + t163 * t253;
t191 = -t210 * t251 + t253 * t283;
t192 = t210 * t253 + t251 * t283;
t176 = -mrSges(6,1) * t191 + mrSges(6,2) * t192;
t178 = qJD(5) * t191 + t208 * t253 + t251 * t279;
t206 = qJD(5) + t209;
t179 = -mrSges(6,2) * t206 + mrSges(6,3) * t191;
t205 = qJDD(5) + t207;
t154 = m(6) * t157 + mrSges(6,1) * t205 - mrSges(6,3) * t178 - t176 * t192 + t179 * t206;
t158 = t160 * t253 + t163 * t251;
t177 = -qJD(5) * t192 - t208 * t251 + t253 * t279;
t180 = mrSges(6,1) * t206 - mrSges(6,3) * t192;
t155 = m(6) * t158 - mrSges(6,2) * t205 + mrSges(6,3) * t177 + t176 * t191 - t180 * t206;
t148 = t154 * t253 + t155 * t251;
t271 = t166 * t245 - t174 * t248;
t159 = -pkin(4) * t279 - pkin(6) * t284 + (t299 + t190) * t210 + t271;
t170 = Ifges(6,5) * t192 + Ifges(6,6) * t191 + Ifges(6,3) * t206;
t172 = Ifges(6,1) * t192 + Ifges(6,4) * t191 + Ifges(6,5) * t206;
t150 = -mrSges(6,1) * t159 + mrSges(6,3) * t158 + Ifges(6,4) * t178 + Ifges(6,2) * t177 + Ifges(6,6) * t205 - t170 * t192 + t172 * t206;
t171 = Ifges(6,4) * t192 + Ifges(6,2) * t191 + Ifges(6,6) * t206;
t151 = mrSges(6,2) * t159 - mrSges(6,3) * t157 + Ifges(6,1) * t178 + Ifges(6,4) * t177 + Ifges(6,5) * t205 + t170 * t191 - t171 * t206;
t161 = -0.2e1 * qJD(4) * t210 - t271;
t181 = Ifges(5,5) * t210 - Ifges(5,6) * t209 + Ifges(5,3) * t283;
t182 = Ifges(5,4) * t210 - Ifges(5,2) * t209 + Ifges(5,6) * t283;
t293 = t246 * t247;
t136 = mrSges(5,2) * t165 - mrSges(5,3) * t161 + Ifges(5,1) * t208 - Ifges(5,4) * t207 - pkin(6) * t148 - t150 * t251 + t151 * t253 - t181 * t209 + (Ifges(5,5) * qJDD(1) - qJD(1) * t182) * t293;
t183 = Ifges(5,1) * t210 - Ifges(5,4) * t209 + Ifges(5,5) * t283;
t257 = mrSges(6,1) * t157 - mrSges(6,2) * t158 + Ifges(6,5) * t178 + Ifges(6,6) * t177 + Ifges(6,3) * t205 + t192 * t171 - t191 * t172;
t137 = -mrSges(5,1) * t165 + mrSges(5,3) * t162 + Ifges(5,4) * t208 - Ifges(5,2) * t207 - pkin(4) * t148 - t210 * t181 + (Ifges(5,6) * qJDD(1) + qJD(1) * t183) * t293 - t257;
t149 = -t154 * t251 + t253 * t155;
t189 = mrSges(5,1) * t209 + mrSges(5,2) * t210;
t195 = mrSges(5,1) * t283 - mrSges(5,3) * t210;
t146 = m(5) * t162 - mrSges(5,3) * t207 - t189 * t209 + (-mrSges(5,2) * qJDD(1) - qJD(1) * t195) * t293 + t149;
t156 = -m(6) * t159 + t177 * mrSges(6,1) - mrSges(6,2) * t178 + t191 * t179 - t180 * t192;
t194 = -mrSges(5,2) * t283 - mrSges(5,3) * t209;
t152 = m(5) * t161 - mrSges(5,3) * t208 - t189 * t210 + (mrSges(5,1) * qJDD(1) + qJD(1) * t194) * t293 + t156;
t142 = t146 * t245 + t152 * t248;
t273 = Ifges(4,5) * t249 - Ifges(4,6) * t246;
t200 = (-Ifges(4,3) * t250 + t273 * t247) * qJD(1);
t296 = Ifges(4,6) * t250;
t297 = Ifges(4,4) * t249;
t201 = (-t296 + (-Ifges(4,2) * t246 + t297) * t247) * qJD(1);
t261 = -Ifges(4,5) * t250 + (Ifges(4,1) * t249 - Ifges(4,4) * t246) * t247;
t125 = mrSges(4,2) * t187 - mrSges(4,3) * t168 - qJ(4) * t142 + t136 * t248 - t137 * t245 + (-t200 * t293 + t201 * t250) * qJD(1) + t261 * qJDD(1);
t202 = t261 * qJD(1);
t256 = -mrSges(5,1) * t161 + mrSges(5,2) * t162 - Ifges(5,5) * t208 + Ifges(5,6) * t207 - pkin(4) * t156 - pkin(6) * t149 - t253 * t150 - t251 * t151 - t210 * t182 - t209 * t183;
t126 = t256 + (-t296 + (t297 + (-Ifges(4,2) - Ifges(5,3)) * t246) * t247) * qJDD(1) + (-t200 * t292 - t202 * t250) * qJD(1) - mrSges(4,1) * t187 + mrSges(4,3) * t169 - pkin(3) * t142;
t143 = t248 * t146 - t152 * t245;
t276 = mrSges(4,1) * t246 + mrSges(4,2) * t249;
t214 = t276 * t289;
t267 = -mrSges(4,1) * t250 - mrSges(4,3) * t292;
t218 = t267 * qJD(1);
t266 = mrSges(4,2) * t250 - mrSges(4,3) * t293;
t140 = m(4) * t169 + t266 * qJDD(1) + (-t214 * t293 + t218 * t250) * qJD(1) + t143;
t147 = -m(5) * t165 - t207 * mrSges(5,1) - mrSges(5,2) * t208 - t209 * t194 - t195 * t210 - t148;
t217 = t266 * qJD(1);
t144 = m(4) * t168 + t267 * qJDD(1) + (-t214 * t292 - t217 * t250) * qJD(1) + t147;
t135 = t249 * t140 - t144 * t246;
t263 = -m(4) * t187 - t142;
t269 = -t217 * t246 - t218 * t249;
t275 = Ifges(3,1) * t247 + Ifges(3,4) * t250;
t300 = -((Ifges(3,4) * t247 + Ifges(3,2) * t250) * t289 - t275 * t288) * qJD(1) - mrSges(3,1) * t196 + mrSges(3,2) * t197 - pkin(2) * ((t269 * qJD(1) - t276 * qJDD(1)) * t247 + t263) - qJ(3) * t135 - t125 * t246 - t126 * t249;
t298 = mrSges(3,2) * t247;
t221 = (-mrSges(3,1) * t250 + t298) * qJD(1);
t132 = m(3) * t197 + (qJDD(1) * mrSges(3,3) + qJD(1) * t221) * t250 + t135;
t138 = m(3) * t196 + ((-mrSges(3,3) - t276) * qJDD(1) + (-t221 + t269) * qJD(1)) * t247 + t263;
t278 = t250 * t132 - t138 * t247;
t274 = Ifges(3,5) * t247 + Ifges(3,6) * t250;
t134 = t140 * t246 + t144 * t249;
t270 = t201 * t249 + t202 * t246;
t216 = -qJDD(1) * pkin(1) + t265;
t222 = t274 * qJD(1);
t122 = mrSges(3,2) * t216 - mrSges(3,3) * t196 - qJ(3) * t134 + t275 * qJDD(1) + t125 * t249 - t126 * t246 + t222 * t288;
t258 = mrSges(4,1) * t168 - mrSges(4,2) * t169 + pkin(3) * t147 + qJ(4) * t143 + t245 * t136 + t248 * t137;
t124 = -mrSges(3,1) * t216 + mrSges(3,3) * t197 - pkin(2) * t134 + (Ifges(3,2) + Ifges(4,3)) * t285 + ((Ifges(3,4) - t273) * qJDD(1) + (-t222 - t270) * qJD(1)) * t247 - t258;
t260 = -m(3) * t216 + mrSges(3,1) * t285 - t134 + (t294 + t295) * mrSges(3,3);
t262 = -mrSges(2,2) * t228 + qJ(2) * t278 + t247 * t122 + t250 * t124 + pkin(1) * (-mrSges(3,2) * t286 + t260) + mrSges(2,1) * t229 + Ifges(2,3) * qJDD(1);
t130 = m(2) * t229 - mrSges(2,2) * t255 + (mrSges(2,1) - t298) * qJDD(1) + t260;
t129 = t132 * t247 + t138 * t250;
t127 = m(2) * t228 - mrSges(2,1) * t255 - qJDD(1) * mrSges(2,2) + t278;
t120 = mrSges(2,1) * g(1) + mrSges(2,3) * t228 + Ifges(2,5) * t255 - pkin(1) * t129 + (Ifges(2,6) - t274) * qJDD(1) + t300;
t119 = -mrSges(2,2) * g(1) - mrSges(2,3) * t229 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t255 - qJ(2) * t129 + t122 * t250 - t124 * t247;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t262, t119, t122, t125, t136, t151; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t252 * t119 + t254 * t120 - pkin(5) * (-t127 * t254 + t130 * t252), t120, t124, t126, t137, t150; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) - t254 * t119 + t252 * t120 + pkin(5) * (t127 * t252 + t130 * t254), t262, t274 * qJDD(1) - t300, -Ifges(4,3) * t285 + (t270 * qJD(1) + t273 * qJDD(1)) * t247 + t258, Ifges(5,3) * t279 - t256, t257;];
m_new = t1;
