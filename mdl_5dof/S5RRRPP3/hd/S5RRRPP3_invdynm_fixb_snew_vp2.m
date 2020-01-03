% Calculate vector of cutting torques with Newton-Euler for
% S5RRRPP3
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
% Datum: 2019-12-31 20:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRRPP3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP3_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP3_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP3_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP3_invdynm_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP3_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP3_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP3_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:53:11
% EndTime: 2019-12-31 20:53:14
% DurationCPUTime: 1.80s
% Computational Cost: add. (20032->273), mult. (25165->318), div. (0->0), fcn. (10755->6), ass. (0->102)
t235 = sin(qJ(1));
t238 = cos(qJ(1));
t213 = t235 * g(1) - g(2) * t238;
t203 = qJDD(1) * pkin(1) + t213;
t214 = -g(1) * t238 - g(2) * t235;
t240 = qJD(1) ^ 2;
t204 = -pkin(1) * t240 + t214;
t234 = sin(qJ(2));
t237 = cos(qJ(2));
t157 = t234 * t203 + t237 * t204;
t223 = qJD(1) + qJD(2);
t221 = t223 ^ 2;
t222 = qJDD(1) + qJDD(2);
t154 = -pkin(2) * t221 + pkin(7) * t222 + t157;
t233 = sin(qJ(3));
t236 = cos(qJ(3));
t149 = -t236 * g(3) - t233 * t154;
t150 = -g(3) * t233 + t236 * t154;
t166 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t233 + Ifges(4,2) * t236) * t223;
t167 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t233 + Ifges(4,4) * t236) * t223;
t168 = Ifges(6,5) * qJD(3) + (-Ifges(6,6) * t236 + Ifges(6,3) * t233) * t223;
t263 = qJD(3) * t223;
t261 = t236 * t263;
t195 = t222 * t233 + t261;
t260 = t233 * t263;
t196 = t222 * t236 - t260;
t270 = t223 * t236;
t209 = -mrSges(5,1) * t270 - qJD(3) * mrSges(5,3);
t194 = (-mrSges(6,2) * t233 - mrSges(6,3) * t236) * t223;
t191 = (-pkin(3) * t236 - qJ(4) * t233) * t223;
t239 = qJD(3) ^ 2;
t271 = t223 * t233;
t148 = -qJDD(3) * pkin(3) - qJ(4) * t239 + t191 * t271 + qJDD(4) - t149;
t275 = -2 * qJD(5);
t143 = qJD(3) * t275 + (-t221 * t233 * t236 - qJDD(3)) * qJ(5) + (t195 - t261) * pkin(4) + t148;
t210 = mrSges(6,1) * t270 + qJD(3) * mrSges(6,2);
t256 = -m(6) * t143 + qJDD(3) * mrSges(6,3) + qJD(3) * t210;
t272 = mrSges(6,1) * t195;
t135 = t194 * t271 - t256 + t272;
t247 = -pkin(3) * t239 + qJDD(3) * qJ(4) + t191 * t270 + t150;
t276 = -2 * qJD(4);
t146 = qJD(3) * t276 - t247;
t171 = Ifges(5,4) * qJD(3) + (-Ifges(5,2) * t233 - Ifges(5,6) * t236) * t223;
t207 = pkin(4) * t271 - qJD(3) * qJ(5);
t232 = t236 ^ 2;
t144 = -qJ(5) * t221 * t232 + pkin(4) * t196 + qJDD(5) + ((2 * qJD(4)) + t207) * qJD(3) + t247;
t254 = mrSges(6,2) * t144 - mrSges(6,3) * t143 + Ifges(6,1) * qJDD(3) - Ifges(6,4) * t196 + Ifges(6,5) * t195;
t244 = mrSges(5,2) * t148 - mrSges(5,3) * t146 + Ifges(5,1) * qJDD(3) - Ifges(5,4) * t195 - Ifges(5,5) * t196 - qJ(5) * t135 + t171 * t270 + t254;
t192 = (mrSges(5,2) * t236 - mrSges(5,3) * t233) * t223;
t211 = mrSges(5,1) * t271 + qJD(3) * mrSges(5,2);
t208 = mrSges(6,1) * t271 - qJD(3) * mrSges(6,3);
t257 = m(6) * t144 + qJDD(3) * mrSges(6,2) + qJD(3) * t208 + t194 * t270;
t248 = -m(5) * t146 + qJDD(3) * mrSges(5,3) + qJD(3) * t211 + t192 * t270 + t257;
t252 = -m(5) * t148 - t195 * mrSges(5,1) + t256;
t264 = -t192 - t194;
t169 = Ifges(5,5) * qJD(3) + (-Ifges(5,6) * t233 - Ifges(5,3) * t236) * t223;
t170 = Ifges(6,4) * qJD(3) + (-Ifges(6,2) * t236 + Ifges(6,6) * t233) * t223;
t265 = -t169 - t170;
t279 = ((t166 + t265) * t233 - (t167 + t168) * t236) * t223 + mrSges(4,1) * t149 - mrSges(4,2) * t150 + Ifges(4,5) * t195 + Ifges(4,6) * t196 + Ifges(4,3) * qJDD(3) + pkin(3) * (-qJDD(3) * mrSges(5,2) - qJD(3) * t209 + t264 * t271 + t252 - t272) + qJ(4) * ((mrSges(5,1) + mrSges(6,1)) * t196 + t248) + t244;
t278 = pkin(3) * t260 + t271 * t276;
t274 = -mrSges(6,1) - mrSges(4,3);
t273 = Ifges(4,4) + Ifges(5,6);
t269 = t236 * t168;
t193 = (-mrSges(4,1) * t236 + mrSges(4,2) * t233) * t223;
t205 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t271;
t129 = t193 * t270 + m(4) * t150 - qJDD(3) * mrSges(4,2) - qJD(3) * t205 + (mrSges(5,1) - t274) * t196 + t248;
t206 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t270;
t130 = m(4) * t149 + t274 * t195 + (mrSges(4,1) - mrSges(5,2)) * qJDD(3) + (t206 - t209) * qJD(3) + (-t193 + t264) * t271 + t252;
t258 = t236 * t129 - t130 * t233;
t120 = m(3) * t157 - mrSges(3,1) * t221 - mrSges(3,2) * t222 + t258;
t156 = t203 * t237 - t234 * t204;
t255 = -pkin(2) * t222 - t156;
t153 = -pkin(7) * t221 + t255;
t139 = -qJ(4) * t195 + (-pkin(4) * t232 - pkin(7)) * t221 + (-pkin(3) - qJ(5)) * t196 + (-t207 * t233 + (-qJ(4) * qJD(3) + t275) * t236) * t223 + t255 + t278;
t134 = m(6) * t139 - t195 * mrSges(6,2) - t196 * mrSges(6,3) - t208 * t271 - t210 * t270;
t145 = -pkin(3) * t196 + (-t195 - t261) * qJ(4) + t153 + t278;
t249 = -m(5) * t145 - t196 * mrSges(5,2) + t211 * t271 - t134;
t242 = -m(4) * t153 + t206 * t270 + t196 * mrSges(4,1) + (-t205 * t233 - t209 * t236) * t223 + (-mrSges(4,2) + mrSges(5,3)) * t195 + t249;
t124 = m(3) * t156 + mrSges(3,1) * t222 - mrSges(3,2) * t221 + t242;
t115 = t234 * t120 + t237 * t124;
t122 = t233 * t129 + t236 * t130;
t173 = Ifges(5,1) * qJD(3) + (-Ifges(5,4) * t233 - Ifges(5,5) * t236) * t223;
t268 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t233 + Ifges(4,6) * t236) * t223 + t173;
t266 = -t168 + t171;
t259 = t237 * t120 - t234 * t124;
t172 = Ifges(6,1) * qJD(3) + (-Ifges(6,4) * t236 + Ifges(6,5) * t233) * t223;
t253 = mrSges(6,1) * t144 - mrSges(6,3) * t139 - Ifges(6,4) * qJDD(3) + Ifges(6,2) * t196 - Ifges(6,6) * t195 - t172 * t271;
t131 = -mrSges(5,3) * t195 + t209 * t270 - t249;
t243 = mrSges(5,1) * t146 - mrSges(5,2) * t145 + pkin(4) * (-mrSges(6,1) * t196 - t257) + qJ(5) * t134 - t253;
t111 = -t268 * t271 + (Ifges(4,2) + Ifges(5,3)) * t196 + t273 * t195 + (Ifges(4,6) - Ifges(5,5)) * qJDD(3) + (t167 - t266) * qJD(3) + mrSges(4,3) * t150 - mrSges(4,1) * t153 - pkin(3) * t131 - t243;
t250 = mrSges(6,1) * t143 - mrSges(6,2) * t139 + Ifges(6,5) * qJDD(3) - Ifges(6,6) * t196 + Ifges(6,3) * t195 + qJD(3) * t170 + t172 * t270;
t245 = mrSges(5,1) * t148 - mrSges(5,3) * t145 + pkin(4) * t135 + t250;
t117 = t268 * t270 + t273 * t196 + (Ifges(4,1) + Ifges(5,2)) * t195 + (Ifges(4,5) - Ifges(5,4)) * qJDD(3) + t245 + (-t166 + t169) * qJD(3) - mrSges(4,3) * t149 + mrSges(4,2) * t153 - qJ(4) * t131;
t251 = mrSges(3,1) * t156 - mrSges(3,2) * t157 + Ifges(3,3) * t222 + pkin(2) * t242 + pkin(7) * t258 + t236 * t111 + t233 * t117;
t246 = mrSges(2,1) * t213 - mrSges(2,2) * t214 + Ifges(2,3) * qJDD(1) + pkin(1) * t115 + t251;
t113 = m(2) * t214 - mrSges(2,1) * t240 - qJDD(1) * mrSges(2,2) + t259;
t112 = m(2) * t213 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t240 + t115;
t109 = mrSges(3,1) * g(3) + mrSges(3,3) * t157 + t221 * Ifges(3,5) + Ifges(3,6) * t222 - pkin(2) * t122 - t279;
t108 = -mrSges(3,2) * g(3) - mrSges(3,3) * t156 + Ifges(3,5) * t222 - Ifges(3,6) * t221 - pkin(7) * t122 - t111 * t233 + t117 * t236;
t107 = -mrSges(2,2) * g(3) - mrSges(2,3) * t213 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t240 - pkin(6) * t115 + t108 * t237 - t109 * t234;
t106 = Ifges(2,6) * qJDD(1) + t240 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t214 + t234 * t108 + t237 * t109 - pkin(1) * (-m(3) * g(3) + t122) + pkin(6) * t259;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t238 * t107 - t235 * t106 - pkin(5) * (t112 * t238 + t113 * t235), t107, t108, t117, (t265 * t233 - t269) * t223 + t244, (-t233 * t170 - t269) * t223 + t254; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t235 * t107 + t238 * t106 + pkin(5) * (-t112 * t235 + t113 * t238), t106, t109, t111, Ifges(5,4) * qJDD(3) - Ifges(5,2) * t195 - Ifges(5,6) * t196 - qJD(3) * t169 - t173 * t270 - t245, -qJD(3) * t168 - t253; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t246, t246, t251, t279, Ifges(5,5) * qJDD(3) - Ifges(5,6) * t195 - Ifges(5,3) * t196 + t266 * qJD(3) + t173 * t271 + t243, t250;];
m_new = t1;
