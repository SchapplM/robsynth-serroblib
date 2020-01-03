% Calculate vector of cutting torques with Newton-Euler for
% S5RRRPR4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRRPR4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR4_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR4_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR4_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR4_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR4_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR4_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR4_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:11:04
% EndTime: 2019-12-31 21:11:10
% DurationCPUTime: 2.86s
% Computational Cost: add. (41990->269), mult. (53297->335), div. (0->0), fcn. (27361->8), ass. (0->104)
t247 = sin(qJ(1));
t251 = cos(qJ(1));
t224 = t247 * g(1) - t251 * g(2);
t216 = qJDD(1) * pkin(1) + t224;
t225 = -t251 * g(1) - t247 * g(2);
t253 = qJD(1) ^ 2;
t217 = -t253 * pkin(1) + t225;
t246 = sin(qJ(2));
t250 = cos(qJ(2));
t180 = t246 * t216 + t250 * t217;
t236 = qJD(1) + qJD(2);
t232 = t236 ^ 2;
t234 = qJDD(1) + qJDD(2);
t177 = -t232 * pkin(2) + t234 * pkin(7) + t180;
t245 = sin(qJ(3));
t249 = cos(qJ(3));
t166 = -t249 * g(3) - t245 * t177;
t167 = -t245 * g(3) + t249 * t177;
t188 = Ifges(5,6) * qJD(3) + (Ifges(5,5) * t245 - Ifges(5,3) * t249) * t236;
t191 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t245 + Ifges(4,2) * t249) * t236;
t206 = (-mrSges(5,1) * t249 - mrSges(5,3) * t245) * t236;
t269 = qJD(3) * t249;
t268 = t236 * t269;
t208 = t245 * t234 + t268;
t273 = t236 * t245;
t209 = -qJD(3) * t273 + t249 * t234;
t205 = (-pkin(3) * t249 - qJ(4) * t245) * t236;
t252 = qJD(3) ^ 2;
t272 = t236 * t249;
t276 = 2 * qJD(4);
t155 = -t252 * pkin(3) + qJDD(3) * qJ(4) + qJD(3) * t276 + t205 * t272 + t167;
t222 = -qJD(3) * pkin(4) - pkin(8) * t273;
t243 = t249 ^ 2;
t150 = -t243 * t232 * pkin(4) - t209 * pkin(8) + qJD(3) * t222 + t155;
t157 = -qJDD(3) * pkin(3) - t252 * qJ(4) + t205 * t273 + qJDD(4) - t166;
t151 = (-t208 + t268) * pkin(8) + (-t232 * t245 * t249 - qJDD(3)) * pkin(4) + t157;
t244 = sin(qJ(5));
t248 = cos(qJ(5));
t146 = -t244 * t150 + t248 * t151;
t194 = (-t244 * t245 - t248 * t249) * t236;
t165 = t194 * qJD(5) + t248 * t208 - t244 * t209;
t195 = (-t244 * t249 + t245 * t248) * t236;
t175 = -t194 * mrSges(6,1) + t195 * mrSges(6,2);
t235 = -qJD(3) + qJD(5);
t181 = -t235 * mrSges(6,2) + t194 * mrSges(6,3);
t233 = -qJDD(3) + qJDD(5);
t141 = m(6) * t146 + t233 * mrSges(6,1) - t165 * mrSges(6,3) - t195 * t175 + t235 * t181;
t147 = t248 * t150 + t244 * t151;
t164 = -t195 * qJD(5) - t244 * t208 - t248 * t209;
t182 = t235 * mrSges(6,1) - t195 * mrSges(6,3);
t142 = m(6) * t147 - t233 * mrSges(6,2) + t164 * mrSges(6,3) + t194 * t175 - t235 * t182;
t131 = t248 * t141 + t244 * t142;
t169 = Ifges(6,4) * t195 + Ifges(6,2) * t194 + Ifges(6,6) * t235;
t170 = Ifges(6,1) * t195 + Ifges(6,4) * t194 + Ifges(6,5) * t235;
t263 = mrSges(6,1) * t146 - mrSges(6,2) * t147 + Ifges(6,5) * t165 + Ifges(6,6) * t164 + Ifges(6,3) * t233 + t195 * t169 - t194 * t170;
t256 = -mrSges(5,1) * t157 + mrSges(5,3) * t155 + Ifges(5,4) * t208 + Ifges(5,2) * qJDD(3) - Ifges(5,6) * t209 - pkin(4) * t131 - t263;
t221 = mrSges(5,2) * t272 + qJD(3) * mrSges(5,3);
t260 = -m(5) * t157 + qJDD(3) * mrSges(5,1) + qJD(3) * t221 - t131;
t132 = -t244 * t141 + t248 * t142;
t219 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t273;
t262 = m(5) * t155 + qJDD(3) * mrSges(5,3) + qJD(3) * t219 + t206 * t272 + t132;
t192 = Ifges(5,4) * qJD(3) + (Ifges(5,1) * t245 - Ifges(5,5) * t249) * t236;
t270 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t245 + Ifges(4,4) * t249) * t236 + t192;
t278 = -((t188 - t191) * t245 + t270 * t249) * t236 + mrSges(4,1) * t166 - mrSges(4,2) * t167 + Ifges(4,5) * t208 + Ifges(4,6) * t209 + Ifges(4,3) * qJDD(3) + pkin(3) * (-t208 * mrSges(5,2) - t206 * t273 + t260) + qJ(4) * (t209 * mrSges(5,2) + t262) + t256;
t275 = t232 * pkin(7);
t274 = mrSges(4,3) + mrSges(5,2);
t207 = (-mrSges(4,1) * t249 + mrSges(4,2) * t245) * t236;
t218 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t273;
t127 = m(4) * t167 - qJDD(3) * mrSges(4,2) - qJD(3) * t218 + t207 * t272 + t274 * t209 + t262;
t220 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t272;
t128 = m(4) * t166 + qJDD(3) * mrSges(4,1) + qJD(3) * t220 + (-t206 - t207) * t273 - t274 * t208 + t260;
t266 = t249 * t127 - t245 * t128;
t121 = m(3) * t180 - t232 * mrSges(3,1) - t234 * mrSges(3,2) + t266;
t179 = t250 * t216 - t246 * t217;
t265 = t234 * pkin(2) + t179;
t264 = -t208 * qJ(4) - t265;
t149 = (-pkin(8) * t243 + pkin(7)) * t232 + (pkin(3) + pkin(4)) * t209 + (qJ(4) * t269 + (-pkin(3) * qJD(3) + t222 + t276) * t245) * t236 - t264;
t143 = -m(6) * t149 + t164 * mrSges(6,1) - t165 * mrSges(6,2) + t194 * t181 - t195 * t182;
t152 = -t209 * pkin(3) - t275 + (-0.2e1 * qJD(4) * t245 + (pkin(3) * t245 - qJ(4) * t249) * qJD(3)) * t236 + t264;
t139 = m(5) * t152 - t209 * mrSges(5,1) - t208 * mrSges(5,3) - t219 * t273 - t221 * t272 + t143;
t176 = -t265 - t275;
t255 = -m(4) * t176 + t209 * mrSges(4,1) - t208 * mrSges(4,2) - t218 * t273 + t220 * t272 - t139;
t134 = m(3) * t179 + t234 * mrSges(3,1) - t232 * mrSges(3,2) + t255;
t118 = t246 * t121 + t250 * t134;
t123 = t245 * t127 + t249 * t128;
t267 = t250 * t121 - t246 * t134;
t189 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t245 + Ifges(4,6) * t249) * t236;
t190 = Ifges(5,2) * qJD(3) + (Ifges(5,4) * t245 - Ifges(5,6) * t249) * t236;
t168 = Ifges(6,5) * t195 + Ifges(6,6) * t194 + Ifges(6,3) * t235;
t137 = -mrSges(6,1) * t149 + mrSges(6,3) * t147 + Ifges(6,4) * t165 + Ifges(6,2) * t164 + Ifges(6,6) * t233 - t195 * t168 + t235 * t170;
t138 = mrSges(6,2) * t149 - mrSges(6,3) * t146 + Ifges(6,1) * t165 + Ifges(6,4) * t164 + Ifges(6,5) * t233 + t194 * t168 - t235 * t169;
t257 = -mrSges(5,1) * t152 + mrSges(5,2) * t155 - pkin(4) * t143 - pkin(8) * t132 - t248 * t137 - t244 * t138;
t112 = -mrSges(4,1) * t176 + mrSges(4,3) * t167 - pkin(3) * t139 + (-t189 - t190) * t273 + (Ifges(4,2) + Ifges(5,3)) * t209 + (Ifges(4,4) - Ifges(5,5)) * t208 + (Ifges(4,6) - Ifges(5,6)) * qJDD(3) + t270 * qJD(3) + t257;
t258 = mrSges(5,2) * t157 - mrSges(5,3) * t152 + Ifges(5,1) * t208 + Ifges(5,4) * qJDD(3) - Ifges(5,5) * t209 - pkin(8) * t131 + qJD(3) * t188 - t244 * t137 + t248 * t138 + t190 * t272;
t114 = mrSges(4,2) * t176 - mrSges(4,3) * t166 + Ifges(4,1) * t208 + Ifges(4,4) * t209 + Ifges(4,5) * qJDD(3) - qJ(4) * t139 - qJD(3) * t191 + t189 * t272 + t258;
t261 = mrSges(3,1) * t179 - mrSges(3,2) * t180 + Ifges(3,3) * t234 + pkin(2) * t255 + pkin(7) * t266 + t249 * t112 + t245 * t114;
t259 = mrSges(2,1) * t224 - mrSges(2,2) * t225 + Ifges(2,3) * qJDD(1) + pkin(1) * t118 + t261;
t116 = m(2) * t225 - t253 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t267;
t115 = m(2) * t224 + qJDD(1) * mrSges(2,1) - t253 * mrSges(2,2) + t118;
t110 = mrSges(3,1) * g(3) + mrSges(3,3) * t180 + t232 * Ifges(3,5) + Ifges(3,6) * t234 - pkin(2) * t123 - t278;
t109 = -mrSges(3,2) * g(3) - mrSges(3,3) * t179 + Ifges(3,5) * t234 - t232 * Ifges(3,6) - pkin(7) * t123 - t245 * t112 + t249 * t114;
t108 = -mrSges(2,2) * g(3) - mrSges(2,3) * t224 + Ifges(2,5) * qJDD(1) - t253 * Ifges(2,6) - pkin(6) * t118 + t250 * t109 - t246 * t110;
t107 = Ifges(2,6) * qJDD(1) + t253 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t225 + t246 * t109 + t250 * t110 - pkin(1) * (-m(3) * g(3) + t123) + pkin(6) * t267;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t251 * t108 - t247 * t107 - pkin(5) * (t251 * t115 + t247 * t116), t108, t109, t114, t258, t138; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t247 * t108 + t251 * t107 + pkin(5) * (-t247 * t115 + t251 * t116), t107, t110, t112, (-t245 * t188 - t249 * t192) * t236 + t256, t137; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t259, t259, t261, t278, Ifges(5,5) * t208 + Ifges(5,6) * qJDD(3) - Ifges(5,3) * t209 - qJD(3) * t192 + t190 * t273 - t257, t263;];
m_new = t1;
