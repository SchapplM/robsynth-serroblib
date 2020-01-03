% Calculate vector of cutting torques with Newton-Euler for
% S5PRRPP4
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
%   pkin=[a2,a3,a4,a5,d2,d3,theta1]';
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
% Datum: 2019-12-31 17:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PRRPP4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP4_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP4_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPP4_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRRPP4_invdynm_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP4_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPP4_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPP4_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:40:49
% EndTime: 2019-12-31 17:40:51
% DurationCPUTime: 1.57s
% Computational Cost: add. (11876->257), mult. (23417->307), div. (0->0), fcn. (10839->6), ass. (0->98)
t237 = sin(pkin(7));
t238 = cos(pkin(7));
t214 = g(1) * t237 - g(2) * t238;
t215 = -g(1) * t238 - g(2) * t237;
t240 = sin(qJ(2));
t242 = cos(qJ(2));
t159 = t240 * t214 + t242 * t215;
t244 = qJD(2) ^ 2;
t156 = -pkin(2) * t244 + qJDD(2) * pkin(6) + t159;
t239 = sin(qJ(3));
t153 = t239 * t156;
t236 = -g(3) + qJDD(1);
t241 = cos(qJ(3));
t273 = t236 * t241;
t151 = -t153 + t273;
t152 = t241 * t156 + t239 * t236;
t170 = Ifges(5,6) * qJD(3) + (Ifges(5,5) * t239 - Ifges(5,3) * t241) * qJD(2);
t174 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t239 + Ifges(4,2) * t241) * qJD(2);
t204 = (-mrSges(5,1) * t241 - mrSges(5,3) * t239) * qJD(2);
t266 = qJD(2) * qJD(3);
t207 = qJDD(2) * t239 + t241 * t266;
t208 = qJDD(2) * t241 - t239 * t266;
t203 = (-pkin(3) * t241 - qJ(4) * t239) * qJD(2);
t243 = qJD(3) ^ 2;
t268 = qJD(2) * t239;
t257 = -qJ(4) * t243 + t203 * t268 + qJDD(4) + t153;
t265 = -0.2e1 * qJD(2) * qJD(5);
t278 = pkin(4) * t244;
t279 = pkin(3) + pkin(4);
t144 = t239 * t265 - qJ(5) * t207 - t279 * qJDD(3) + (qJ(5) * t266 - t239 * t278 - t236) * t241 + t257;
t205 = (mrSges(6,1) * t241 + mrSges(6,2) * t239) * qJD(2);
t267 = qJD(2) * t241;
t220 = qJD(3) * mrSges(6,2) - mrSges(6,3) * t267;
t135 = m(6) * t144 - qJDD(3) * mrSges(6,1) - t207 * mrSges(6,3) - qJD(3) * t220 - t205 * t268;
t280 = 2 * qJD(4);
t148 = -pkin(3) * t243 + qJDD(3) * qJ(4) + qJD(3) * t280 + t203 * t267 + t152;
t150 = -qJDD(3) * pkin(3) + t257 - t273;
t216 = -qJD(3) * pkin(4) - qJ(5) * t268;
t235 = t241 ^ 2;
t143 = -qJ(5) * t208 + qJD(3) * t216 - t235 * t278 + t241 * t265 + t148;
t172 = -Ifges(6,6) * qJD(3) + (Ifges(6,4) * t239 - Ifges(6,2) * t241) * qJD(2);
t175 = -Ifges(6,5) * qJD(3) + (Ifges(6,1) * t239 - Ifges(6,4) * t241) * qJD(2);
t255 = mrSges(6,1) * t144 - mrSges(6,2) * t143 + Ifges(6,5) * t207 - Ifges(6,6) * t208 - Ifges(6,3) * qJDD(3) + t172 * t268 + t175 * t267;
t248 = -mrSges(5,1) * t150 + mrSges(5,3) * t148 + Ifges(5,4) * t207 + Ifges(5,2) * qJDD(3) - Ifges(5,6) * t208 - pkin(4) * t135 - t255;
t222 = mrSges(5,2) * t267 + qJD(3) * mrSges(5,3);
t251 = -m(5) * t150 + qJDD(3) * mrSges(5,1) + qJD(3) * t222 - t135;
t219 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t268;
t217 = -qJD(3) * mrSges(6,1) - mrSges(6,3) * t268;
t260 = m(6) * t143 + qJDD(3) * mrSges(6,2) - t208 * mrSges(6,3) + qJD(3) * t217;
t253 = m(5) * t148 + qJDD(3) * mrSges(5,3) + qJD(3) * t219 + t204 * t267 + t260;
t264 = t205 * t267;
t176 = Ifges(5,4) * qJD(3) + (Ifges(5,1) * t239 - Ifges(5,5) * t241) * qJD(2);
t270 = t176 + Ifges(4,5) * qJD(3) + (Ifges(4,1) * t239 + Ifges(4,4) * t241) * qJD(2);
t283 = -((t170 - t174) * t239 + t270 * t241) * qJD(2) + mrSges(4,1) * t151 - mrSges(4,2) * t152 + Ifges(4,5) * t207 + Ifges(4,6) * t208 + Ifges(4,3) * qJDD(3) + pkin(3) * (-mrSges(5,2) * t207 - t204 * t268 + t251) + qJ(4) * (mrSges(5,2) * t208 + t253 - t264) + t248;
t169 = -Ifges(6,3) * qJD(3) + (Ifges(6,5) * t239 - Ifges(6,6) * t241) * qJD(2);
t282 = -Ifges(6,5) * qJDD(3) - t169 * t267;
t277 = pkin(6) * t244;
t276 = mrSges(4,3) + mrSges(5,2);
t275 = Ifges(5,6) - Ifges(6,6);
t274 = qJ(4) * t241;
t206 = (-mrSges(4,1) * t241 + mrSges(4,2) * t239) * qJD(2);
t218 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t268;
t128 = m(4) * t152 - qJDD(3) * mrSges(4,2) - qJD(3) * t218 + t276 * t208 + (-t205 + t206) * t267 + t253;
t221 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t267;
t129 = m(4) * t151 + qJDD(3) * mrSges(4,1) + qJD(3) * t221 - t276 * t207 + (-t204 - t206) * t268 + t251;
t261 = t241 * t128 - t129 * t239;
t119 = m(3) * t159 - mrSges(3,1) * t244 - qJDD(2) * mrSges(3,2) + t261;
t158 = t242 * t214 - t240 * t215;
t259 = qJDD(2) * pkin(2) + t158;
t254 = -qJ(4) * t207 - t259;
t138 = qJDD(5) + (-qJ(5) * t235 + pkin(6)) * t244 + t279 * t208 + (qJD(3) * t274 + (-pkin(3) * qJD(3) + t216 + t280) * t239) * qJD(2) - t254;
t133 = -m(6) * t138 - t208 * mrSges(6,1) - t207 * mrSges(6,2) - t217 * t268 - t220 * t267;
t145 = -pkin(3) * t208 - t277 + (-0.2e1 * qJD(4) * t239 + (pkin(3) * t239 - t274) * qJD(3)) * qJD(2) + t254;
t130 = m(5) * t145 - mrSges(5,1) * t208 - t207 * mrSges(5,3) - t219 * t268 - t222 * t267 + t133;
t155 = -t259 - t277;
t246 = -m(4) * t155 + t208 * mrSges(4,1) - mrSges(4,2) * t207 - t218 * t268 + t221 * t267 - t130;
t123 = m(3) * t158 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t244 + t246;
t114 = t240 * t119 + t242 * t123;
t121 = t239 * t128 + t241 * t129;
t173 = Ifges(5,2) * qJD(3) + (Ifges(5,4) * t239 - Ifges(5,6) * t241) * qJD(2);
t272 = -t169 + t173;
t262 = t242 * t119 - t240 * t123;
t258 = mrSges(6,1) * t138 - mrSges(6,3) * t143 - Ifges(6,4) * t207 + Ifges(6,2) * t208;
t256 = mrSges(6,2) * t138 - mrSges(6,3) * t144 + Ifges(6,1) * t207 - Ifges(6,4) * t208 + qJD(3) * t172;
t171 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t239 + Ifges(4,6) * t241) * qJD(2);
t249 = mrSges(5,1) * t145 - mrSges(5,2) * t148 + pkin(4) * t133 + qJ(5) * (t260 - t264) - t258;
t110 = -t249 + (Ifges(4,2) + Ifges(5,3)) * t208 + (Ifges(4,4) - Ifges(5,5)) * t207 + (Ifges(4,6) - t275) * qJDD(3) + (t175 + t270) * qJD(3) + mrSges(4,3) * t152 - mrSges(4,1) * t155 - pkin(3) * t130 + (-t171 - t272) * t268;
t247 = mrSges(5,2) * t150 - mrSges(5,3) * t145 + Ifges(5,1) * t207 + Ifges(5,4) * qJDD(3) - Ifges(5,5) * t208 - qJ(5) * t135 + qJD(3) * t170 + t173 * t267 + t256;
t116 = t247 + (Ifges(4,5) - Ifges(6,5)) * qJDD(3) + Ifges(4,1) * t207 + Ifges(4,4) * t208 - qJD(3) * t174 - mrSges(4,3) * t151 + mrSges(4,2) * t155 - qJ(4) * t130 + (-t169 + t171) * t267;
t252 = mrSges(3,1) * t158 - mrSges(3,2) * t159 + Ifges(3,3) * qJDD(2) + pkin(2) * t246 + pkin(6) * t261 + t241 * t110 + t239 * t116;
t250 = mrSges(2,1) * t214 - mrSges(2,2) * t215 + pkin(1) * t114 + t252;
t112 = m(2) * t215 + t262;
t111 = m(2) * t214 + t114;
t108 = -mrSges(3,1) * t236 + mrSges(3,3) * t159 + t244 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t121 - t283;
t107 = mrSges(3,2) * t236 - mrSges(3,3) * t158 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t244 - pkin(6) * t121 - t110 * t239 + t116 * t241;
t106 = mrSges(2,2) * t236 - mrSges(2,3) * t214 - pkin(5) * t114 + t107 * t242 - t108 * t240;
t105 = -mrSges(2,1) * t236 + mrSges(2,3) * t215 + t240 * t107 + t242 * t108 - pkin(1) * (m(3) * t236 + t121) + pkin(5) * t262;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t238 * t106 - t237 * t105 - qJ(1) * (t111 * t238 + t112 * t237), t106, t107, t116, t247 + t282, t256 + t282; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t237 * t106 + t238 * t105 + qJ(1) * (-t111 * t237 + t112 * t238), t105, t108, t110, t248 + (-t239 * t170 - t241 * t176) * qJD(2), -Ifges(6,6) * qJDD(3) - qJD(3) * t175 - t169 * t268 - t258; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t250, t250, t252, t283, Ifges(5,5) * t207 - Ifges(5,3) * t208 + t275 * qJDD(3) + (-t175 - t176) * qJD(3) + t272 * t268 + t249, t255;];
m_new = t1;
