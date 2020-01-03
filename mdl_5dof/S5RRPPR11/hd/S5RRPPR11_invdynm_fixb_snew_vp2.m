% Calculate vector of cutting torques with Newton-Euler for
% S5RRPPR11
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
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
% Datum: 2019-12-31 19:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRPPR11_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR11_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR11_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR11_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR11_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR11_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR11_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR11_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR11_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:46:31
% EndTime: 2019-12-31 19:46:39
% DurationCPUTime: 4.71s
% Computational Cost: add. (52512->315), mult. (114673->382), div. (0->0), fcn. (65263->8), ass. (0->117)
t306 = -2 * qJD(3);
t270 = sin(qJ(1));
t273 = cos(qJ(1));
t253 = -t273 * g(1) - t270 * g(2);
t275 = qJD(1) ^ 2;
t225 = -t275 * pkin(1) + qJDD(1) * pkin(6) + t253;
t269 = sin(qJ(2));
t272 = cos(qJ(2));
t200 = -t272 * g(3) - t269 * t225;
t201 = -t269 * g(3) + t272 * t225;
t220 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t269 + Ifges(3,4) * t272) * qJD(1);
t240 = (mrSges(4,2) * t272 - mrSges(4,3) * t269) * qJD(1);
t296 = qJD(1) * qJD(2);
t293 = t272 * t296;
t242 = t269 * qJDD(1) + t293;
t294 = t269 * t296;
t243 = t272 * qJDD(1) - t294;
t298 = qJD(1) * t272;
t250 = -mrSges(4,1) * t298 - qJD(2) * mrSges(4,3);
t239 = (-pkin(2) * t272 - qJ(3) * t269) * qJD(1);
t274 = qJD(2) ^ 2;
t188 = t274 * pkin(2) - qJDD(2) * qJ(3) + qJD(2) * t306 - t239 * t298 - t201;
t297 = t269 * qJD(1);
t249 = pkin(3) * t297 - qJD(2) * qJ(4);
t265 = t272 ^ 2;
t179 = -t265 * t275 * qJ(4) + t243 * pkin(3) + qJD(2) * t249 + qJDD(4) - t188;
t266 = sin(pkin(8));
t267 = cos(pkin(8));
t230 = -t266 * qJD(2) - t267 * t298;
t205 = -mrSges(5,2) * t297 + t230 * mrSges(5,3);
t231 = t267 * qJD(2) - t266 * t298;
t206 = mrSges(5,1) * t297 - t231 * mrSges(5,3);
t207 = -t266 * qJDD(2) - t267 * t243;
t208 = t267 * qJDD(2) - t266 * t243;
t209 = pkin(4) * t297 - t231 * pkin(7);
t229 = t230 ^ 2;
t167 = -t207 * pkin(4) - t229 * pkin(7) + t231 * t209 + t179;
t268 = sin(qJ(5));
t271 = cos(qJ(5));
t198 = t268 * t230 + t271 * t231;
t174 = -t198 * qJD(5) + t271 * t207 - t268 * t208;
t197 = t271 * t230 - t268 * t231;
t175 = t197 * qJD(5) + t268 * t207 + t271 * t208;
t255 = qJD(5) + t297;
t191 = -t255 * mrSges(6,2) + t197 * mrSges(6,3);
t192 = t255 * mrSges(6,1) - t198 * mrSges(6,3);
t287 = m(6) * t167 - t174 * mrSges(6,1) + t175 * mrSges(6,2) - t197 * t191 + t198 * t192;
t157 = -m(5) * t179 + t207 * mrSges(5,1) - t208 * mrSges(5,2) + t230 * t205 - t231 * t206 - t287;
t251 = mrSges(4,1) * t297 + qJD(2) * mrSges(4,2);
t279 = -m(4) * t188 + qJDD(2) * mrSges(4,3) + qJD(2) * t251 + t240 * t298 - t157;
t252 = t270 * g(1) - t273 * g(2);
t290 = -qJDD(1) * pkin(1) - t252;
t282 = pkin(2) * t294 + t297 * t306 + (-t242 - t293) * qJ(3) + t290;
t170 = -t249 * t297 + (-pkin(3) * t265 - pkin(6)) * t275 + (-pkin(2) - qJ(4)) * t243 + t282;
t190 = -qJDD(2) * pkin(2) - t274 * qJ(3) + t239 * t297 + qJDD(3) - t200;
t183 = (-t269 * t272 * t275 - qJDD(2)) * qJ(4) + (t242 - t293) * pkin(3) + t190;
t164 = -0.2e1 * qJD(4) * t231 - t266 * t170 + t267 * t183;
t161 = (t230 * t297 - t208) * pkin(7) + (t230 * t231 + t242) * pkin(4) + t164;
t165 = 0.2e1 * qJD(4) * t230 + t267 * t170 + t266 * t183;
t162 = -t229 * pkin(4) + t207 * pkin(7) - t209 * t297 + t165;
t160 = t268 * t161 + t271 * t162;
t180 = Ifges(6,5) * t198 + Ifges(6,6) * t197 + Ifges(6,3) * t255;
t182 = Ifges(6,1) * t198 + Ifges(6,4) * t197 + Ifges(6,5) * t255;
t238 = qJDD(5) + t242;
t147 = -mrSges(6,1) * t167 + mrSges(6,3) * t160 + Ifges(6,4) * t175 + Ifges(6,2) * t174 + Ifges(6,6) * t238 - t198 * t180 + t255 * t182;
t159 = t271 * t161 - t268 * t162;
t181 = Ifges(6,4) * t198 + Ifges(6,2) * t197 + Ifges(6,6) * t255;
t148 = mrSges(6,2) * t167 - mrSges(6,3) * t159 + Ifges(6,1) * t175 + Ifges(6,4) * t174 + Ifges(6,5) * t238 + t197 * t180 - t255 * t181;
t193 = Ifges(5,5) * t231 + Ifges(5,6) * t230 + Ifges(5,3) * t297;
t195 = Ifges(5,1) * t231 + Ifges(5,4) * t230 + Ifges(5,5) * t297;
t185 = -t197 * mrSges(6,1) + t198 * mrSges(6,2);
t154 = m(6) * t159 + t238 * mrSges(6,1) - t175 * mrSges(6,3) - t198 * t185 + t255 * t191;
t155 = m(6) * t160 - t238 * mrSges(6,2) + t174 * mrSges(6,3) + t197 * t185 - t255 * t192;
t291 = -t268 * t154 + t271 * t155;
t131 = -mrSges(5,1) * t179 + mrSges(5,3) * t165 + Ifges(5,4) * t208 + Ifges(5,2) * t207 + Ifges(5,6) * t242 - pkin(4) * t287 + pkin(7) * t291 + t271 * t147 + t268 * t148 - t231 * t193 + t195 * t297;
t146 = t271 * t154 + t268 * t155;
t194 = Ifges(5,4) * t231 + Ifges(5,2) * t230 + Ifges(5,6) * t297;
t133 = mrSges(5,2) * t179 - mrSges(5,3) * t164 + Ifges(5,1) * t208 + Ifges(5,4) * t207 + Ifges(5,5) * t242 - pkin(7) * t146 - t268 * t147 + t271 * t148 + t230 * t193 - t194 * t297;
t199 = -t230 * mrSges(5,1) + t231 * mrSges(5,2);
t143 = m(5) * t164 + t242 * mrSges(5,1) - t208 * mrSges(5,3) - t231 * t199 + t205 * t297 + t146;
t144 = m(5) * t165 - t242 * mrSges(5,2) + t207 * mrSges(5,3) + t230 * t199 - t206 * t297 + t291;
t139 = t267 * t143 + t266 * t144;
t222 = Ifges(4,4) * qJD(2) + (-Ifges(4,2) * t269 - Ifges(4,6) * t272) * qJD(1);
t283 = -mrSges(4,2) * t190 + mrSges(4,3) * t188 - Ifges(4,1) * qJDD(2) + Ifges(4,4) * t242 + Ifges(4,5) * t243 + qJ(4) * t139 + t266 * t131 - t267 * t133 - t222 * t298;
t286 = -m(4) * t190 - t242 * mrSges(4,1) - t139;
t221 = Ifges(4,5) * qJD(2) + (-Ifges(4,6) * t269 - Ifges(4,3) * t272) * qJD(1);
t299 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t269 + Ifges(3,2) * t272) * qJD(1) - t221;
t305 = (-t272 * t220 + t299 * t269) * qJD(1) + mrSges(3,1) * t200 - mrSges(3,2) * t201 + Ifges(3,5) * t242 + Ifges(3,6) * t243 + Ifges(3,3) * qJDD(2) + pkin(2) * (-qJDD(2) * mrSges(4,2) - qJD(2) * t250 - t240 * t297 + t286) + qJ(3) * (t243 * mrSges(4,1) + t279) - t283;
t303 = t275 * pkin(6);
t302 = Ifges(3,4) + Ifges(4,6);
t140 = -t266 * t143 + t267 * t144;
t223 = Ifges(4,1) * qJD(2) + (-Ifges(4,4) * t269 - Ifges(4,5) * t272) * qJD(1);
t300 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t269 + Ifges(3,6) * t272) * qJD(1) + t223;
t241 = (-mrSges(3,1) * t272 + mrSges(3,2) * t269) * qJD(1);
t248 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t298;
t136 = m(3) * t200 - t242 * mrSges(3,3) + (mrSges(3,1) - mrSges(4,2)) * qJDD(2) + (t248 - t250) * qJD(2) + (-t240 - t241) * t297 + t286;
t247 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t297;
t150 = t279 - qJD(2) * t247 + m(3) * t201 - qJDD(2) * mrSges(3,2) + (mrSges(3,3) + mrSges(4,1)) * t243 + t241 * t298;
t292 = -t269 * t136 + t272 * t150;
t186 = -t243 * pkin(2) + t282 - t303;
t289 = -m(4) * t186 - t243 * mrSges(4,2) + t251 * t297 - t140;
t137 = -t242 * mrSges(4,3) + t250 * t298 - t289;
t224 = t290 - t303;
t281 = -mrSges(4,1) * t188 + mrSges(4,2) * t186 - pkin(3) * t157 - qJ(4) * t140 - t267 * t131 - t266 * t133;
t125 = -mrSges(3,1) * t224 + mrSges(3,3) * t201 - pkin(2) * t137 + (Ifges(3,2) + Ifges(4,3)) * t243 + t302 * t242 + (Ifges(3,6) - Ifges(4,5)) * qJDD(2) + (t220 - t222) * qJD(2) - t300 * t297 + t281;
t284 = -mrSges(6,1) * t159 + mrSges(6,2) * t160 - Ifges(6,5) * t175 - Ifges(6,6) * t174 - Ifges(6,3) * t238 - t198 * t181 + t197 * t182;
t280 = -mrSges(5,1) * t164 + mrSges(5,2) * t165 - Ifges(5,5) * t208 - Ifges(5,6) * t207 - Ifges(5,3) * t242 - pkin(4) * t146 - t231 * t194 + t230 * t195 + t284;
t277 = -mrSges(4,1) * t190 + mrSges(4,3) * t186 - pkin(3) * t139 + t280;
t127 = -t277 + t300 * t298 + mrSges(3,2) * t224 - mrSges(3,3) * t200 - qJ(3) * t137 - t299 * qJD(2) + (Ifges(3,5) - Ifges(4,4)) * qJDD(2) + (Ifges(3,1) + Ifges(4,2)) * t242 + t302 * t243;
t278 = -m(3) * t224 + t248 * t298 + t243 * mrSges(3,1) + (-mrSges(3,2) + mrSges(4,3)) * t242 + (-t247 * t269 - t250 * t272) * qJD(1) + t289;
t285 = mrSges(2,1) * t252 - mrSges(2,2) * t253 + Ifges(2,3) * qJDD(1) + pkin(1) * t278 + pkin(6) * t292 + t272 * t125 + t269 * t127;
t134 = m(2) * t252 + qJDD(1) * mrSges(2,1) - t275 * mrSges(2,2) + t278;
t130 = t272 * t136 + t269 * t150;
t128 = m(2) * t253 - t275 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t292;
t123 = mrSges(2,1) * g(3) + mrSges(2,3) * t253 + t275 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t130 - t305;
t122 = -mrSges(2,2) * g(3) - mrSges(2,3) * t252 + Ifges(2,5) * qJDD(1) - t275 * Ifges(2,6) - pkin(6) * t130 - t269 * t125 + t272 * t127;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t273 * t122 - t270 * t123 - pkin(5) * (t270 * t128 + t273 * t134), t122, t127, -t221 * t297 - t283, t133, t148; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t270 * t122 + t273 * t123 + pkin(5) * (t273 * t128 - t270 * t134), t123, t125, Ifges(4,4) * qJDD(2) - Ifges(4,2) * t242 - Ifges(4,6) * t243 - qJD(2) * t221 - t223 * t298 + t277, t131, t147; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t285, t285, t305, Ifges(4,5) * qJDD(2) - Ifges(4,6) * t242 - Ifges(4,3) * t243 + qJD(2) * t222 + t223 * t297 - t281, -t280, -t284;];
m_new = t1;
