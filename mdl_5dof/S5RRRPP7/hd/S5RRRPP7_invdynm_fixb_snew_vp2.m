% Calculate vector of cutting torques with Newton-Euler for
% S5RRRPP7
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 21:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRRPP7_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP7_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP7_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP7_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP7_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP7_invdynm_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP7_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP7_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP7_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:04:00
% EndTime: 2019-12-31 21:04:06
% DurationCPUTime: 2.57s
% Computational Cost: add. (26562->307), mult. (51761->354), div. (0->0), fcn. (30444->6), ass. (0->108)
t267 = sin(qJ(1));
t269 = cos(qJ(1));
t254 = t267 * g(1) - t269 * g(2);
t271 = qJD(1) ^ 2;
t230 = -qJDD(1) * pkin(1) - t271 * pkin(6) - t254;
t266 = sin(qJ(2));
t268 = cos(qJ(2));
t289 = qJD(1) * qJD(2);
t287 = t268 * t289;
t248 = qJDD(1) * t266 + t287;
t259 = t266 * t289;
t249 = qJDD(1) * t268 - t259;
t163 = (-t248 - t287) * pkin(7) + (-t249 + t259) * pkin(2) + t230;
t255 = -g(1) * t269 - g(2) * t267;
t231 = -pkin(1) * t271 + qJDD(1) * pkin(6) + t255;
t223 = -g(3) * t266 + t231 * t268;
t247 = (-pkin(2) * t268 - pkin(7) * t266) * qJD(1);
t270 = qJD(2) ^ 2;
t290 = t268 * qJD(1);
t171 = -pkin(2) * t270 + qJDD(2) * pkin(7) + t247 * t290 + t223;
t265 = sin(qJ(3));
t299 = cos(qJ(3));
t161 = t163 * t265 + t171 * t299;
t291 = qJD(1) * t266;
t244 = -qJD(2) * t299 + t265 * t291;
t245 = qJD(2) * t265 + t291 * t299;
t211 = pkin(3) * t244 - qJ(4) * t245;
t243 = -qJDD(3) + t249;
t257 = -qJD(3) + t290;
t256 = t257 ^ 2;
t300 = -2 * qJD(4);
t156 = -pkin(3) * t256 - qJ(4) * t243 - t211 * t244 + t257 * t300 + t161;
t220 = mrSges(5,1) * t257 + mrSges(5,2) * t245;
t305 = m(5) * t156 - t243 * mrSges(5,3) - t257 * t220;
t160 = t163 * t299 - t171 * t265;
t159 = t243 * pkin(3) - t256 * qJ(4) + t211 * t245 + qJDD(4) - t160;
t221 = -mrSges(5,2) * t244 - mrSges(5,3) * t257;
t304 = -m(5) * t159 - t243 * mrSges(5,1) - t257 * t221;
t208 = -qJD(3) * t244 + qJDD(2) * t265 + t248 * t299;
t222 = -t268 * g(3) - t231 * t266;
t279 = qJDD(2) * pkin(2) + t270 * pkin(7) - t247 * t291 + t222;
t296 = t244 * t257;
t303 = -(t208 + t296) * qJ(4) - t279;
t215 = -mrSges(6,2) * t257 + mrSges(6,3) * t244;
t146 = -0.2e1 * qJD(5) * t245 + (-t208 + t296) * qJ(5) + (t244 * t245 + t243) * pkin(4) + t159;
t213 = -mrSges(6,1) * t244 + mrSges(6,2) * t245;
t284 = -m(6) * t146 + mrSges(6,3) * t208 + t213 * t245;
t142 = t243 * mrSges(6,1) + t257 * t215 - t284;
t218 = mrSges(6,1) * t257 - mrSges(6,3) * t245;
t207 = qJD(3) * t245 - qJDD(2) * t299 + t248 * t265;
t217 = pkin(4) * t257 - qJ(5) * t245;
t242 = t244 ^ 2;
t150 = -pkin(4) * t242 + qJ(5) * t207 + 0.2e1 * qJD(5) * t244 - t217 * t257 + t156;
t288 = m(6) * t150 + mrSges(6,3) * t207 + t213 * t244;
t144 = -t243 * mrSges(6,2) - t257 * t218 + t288;
t179 = Ifges(5,5) * t245 - Ifges(5,6) * t257 + Ifges(5,3) * t244;
t183 = Ifges(4,4) * t245 - Ifges(4,2) * t244 - Ifges(4,6) * t257;
t186 = Ifges(4,1) * t245 - Ifges(4,4) * t244 - Ifges(4,5) * t257;
t212 = mrSges(5,1) * t244 - mrSges(5,3) * t245;
t185 = Ifges(5,1) * t245 - Ifges(5,4) * t257 + Ifges(5,5) * t244;
t181 = Ifges(6,4) * t245 + Ifges(6,2) * t244 + Ifges(6,6) * t257;
t184 = Ifges(6,1) * t245 + Ifges(6,4) * t244 + Ifges(6,5) * t257;
t282 = mrSges(6,1) * t146 - mrSges(6,2) * t150 + Ifges(6,5) * t208 + Ifges(6,6) * t207 + Ifges(6,3) * t243 + t181 * t245 - t184 * t244;
t276 = mrSges(5,1) * t159 - mrSges(5,3) * t156 - Ifges(5,4) * t208 + Ifges(5,2) * t243 - Ifges(5,6) * t207 + pkin(4) * t142 - t185 * t244 + t282;
t302 = (t183 - t179) * t245 + mrSges(4,1) * t160 - mrSges(4,2) * t161 + Ifges(4,5) * t208 - Ifges(4,6) * t207 - Ifges(4,3) * t243 + pkin(3) * (-t208 * mrSges(5,2) - t245 * t212 - t142 + t304) + qJ(4) * (-t207 * mrSges(5,2) - t244 * t212 + t144 + t305) + t244 * t186 - t276;
t153 = -t242 * qJ(5) + qJDD(5) + (-pkin(3) - pkin(4)) * t207 + (pkin(3) * t257 + (2 * qJD(4)) + t217) * t245 - t303;
t143 = -m(6) * t153 + mrSges(6,1) * t207 - mrSges(6,2) * t208 + t215 * t244 - t218 * t245;
t158 = t245 * t300 + (-t245 * t257 + t207) * pkin(3) + t303;
t139 = m(5) * t158 + mrSges(5,1) * t207 - mrSges(5,3) * t208 - t220 * t245 + t221 * t244 + t143;
t178 = Ifges(6,5) * t245 + Ifges(6,6) * t244 + Ifges(6,3) * t257;
t281 = mrSges(6,1) * t153 - mrSges(6,3) * t150 - Ifges(6,4) * t208 - Ifges(6,2) * t207 - Ifges(6,6) * t243 - t184 * t257;
t275 = mrSges(5,1) * t158 - mrSges(5,2) * t156 + pkin(4) * t143 + qJ(5) * t144 - t281;
t182 = Ifges(5,4) * t245 - Ifges(5,2) * t257 + Ifges(5,6) * t244;
t294 = -Ifges(4,5) * t245 + Ifges(4,6) * t244 + Ifges(4,3) * t257 - t182;
t123 = (-t186 - t185) * t257 + (t178 + t294) * t245 + (-Ifges(4,6) + Ifges(5,6)) * t243 + (Ifges(4,4) - Ifges(5,5)) * t208 + (-Ifges(4,2) - Ifges(5,3)) * t207 + mrSges(4,3) * t161 + mrSges(4,1) * t279 - pkin(3) * t139 - t275;
t280 = mrSges(6,2) * t153 - mrSges(6,3) * t146 + Ifges(6,1) * t208 + Ifges(6,4) * t207 + Ifges(6,5) * t243 + t178 * t244;
t274 = mrSges(5,2) * t159 - mrSges(5,3) * t158 + Ifges(5,1) * t208 - Ifges(5,4) * t243 + Ifges(5,5) * t207 - qJ(5) * t142 - t179 * t257 + t280;
t127 = (t183 - t181) * t257 + t294 * t244 - Ifges(4,5) * t243 + Ifges(4,1) * t208 - Ifges(4,4) * t207 - mrSges(4,3) * t160 - mrSges(4,2) * t279 - qJ(4) * t139 + t274;
t216 = mrSges(4,2) * t257 - mrSges(4,3) * t244;
t292 = -mrSges(4,1) * t244 - mrSges(4,2) * t245 - t212;
t297 = -mrSges(4,3) - mrSges(5,2);
t136 = m(4) * t160 + (-t215 - t216) * t257 + t292 * t245 + (-mrSges(4,1) - mrSges(6,1)) * t243 + t297 * t208 + t284 + t304;
t219 = -mrSges(4,1) * t257 - mrSges(4,3) * t245;
t137 = m(4) * t161 + (-t218 + t219) * t257 + t292 * t244 + (mrSges(4,2) - mrSges(6,2)) * t243 + t297 * t207 + t288 + t305;
t133 = -t136 * t265 + t137 * t299;
t138 = m(4) * t279 - mrSges(4,1) * t207 - mrSges(4,2) * t208 - t216 * t244 - t245 * t219 - t139;
t228 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t266 + Ifges(3,2) * t268) * qJD(1);
t229 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t266 + Ifges(3,4) * t268) * qJD(1);
t301 = mrSges(3,1) * t222 - mrSges(3,2) * t223 + Ifges(3,5) * t248 + Ifges(3,6) * t249 + Ifges(3,3) * qJDD(2) + pkin(2) * t138 + pkin(7) * t133 + (t228 * t266 - t229 * t268) * qJD(1) + t299 * t123 + t265 * t127;
t295 = t257 * t181;
t246 = (-mrSges(3,1) * t268 + mrSges(3,2) * t266) * qJD(1);
t251 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t291;
t131 = m(3) * t223 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t249 - qJD(2) * t251 + t246 * t290 + t133;
t252 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t290;
t134 = m(3) * t222 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t248 + qJD(2) * t252 - t246 * t291 + t138;
t286 = t131 * t268 - t134 * t266;
t132 = t136 * t299 + t137 * t265;
t227 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t266 + Ifges(3,6) * t268) * qJD(1);
t120 = mrSges(3,2) * t230 - mrSges(3,3) * t222 + Ifges(3,1) * t248 + Ifges(3,4) * t249 + Ifges(3,5) * qJDD(2) - pkin(7) * t132 - qJD(2) * t228 - t123 * t265 + t127 * t299 + t227 * t290;
t122 = -mrSges(3,1) * t230 + mrSges(3,3) * t223 + Ifges(3,4) * t248 + Ifges(3,2) * t249 + Ifges(3,6) * qJDD(2) - pkin(2) * t132 + qJD(2) * t229 - t227 * t291 - t302;
t277 = -m(3) * t230 + t249 * mrSges(3,1) - mrSges(3,2) * t248 - t251 * t291 + t252 * t290 - t132;
t278 = mrSges(2,1) * t254 - mrSges(2,2) * t255 + Ifges(2,3) * qJDD(1) + pkin(1) * t277 + pkin(6) * t286 + t120 * t266 + t122 * t268;
t128 = m(2) * t254 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t271 + t277;
t126 = t131 * t266 + t134 * t268;
t124 = m(2) * t255 - mrSges(2,1) * t271 - qJDD(1) * mrSges(2,2) + t286;
t118 = mrSges(2,1) * g(3) + mrSges(2,3) * t255 + t271 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t126 - t301;
t117 = -mrSges(2,2) * g(3) - mrSges(2,3) * t254 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t271 - pkin(6) * t126 + t120 * t268 - t122 * t266;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t269 * t117 - t267 * t118 - pkin(5) * (t124 * t267 + t128 * t269), t117, t120, t127, -t244 * t182 + t274 - t295, t280 - t295; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t267 * t117 + t269 * t118 + pkin(5) * (t124 * t269 - t128 * t267), t118, t122, t123, -t245 * t179 - t276, -t245 * t178 - t281; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t278, t278, t301, t302, (-t178 + t182) * t245 + t257 * t185 - Ifges(5,6) * t243 + Ifges(5,5) * t208 + Ifges(5,3) * t207 + t275, t282;];
m_new = t1;
