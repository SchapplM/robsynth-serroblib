% Calculate vector of cutting torques with Newton-Euler for
% S5RRRPP2
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
% Datum: 2019-12-31 20:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRRPP2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP2_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP2_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP2_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP2_invdynm_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP2_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP2_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP2_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:51:32
% EndTime: 2019-12-31 20:51:35
% DurationCPUTime: 1.79s
% Computational Cost: add. (20116->267), mult. (25249->318), div. (0->0), fcn. (10839->6), ass. (0->97)
t245 = sin(qJ(1));
t248 = cos(qJ(1));
t223 = t245 * g(1) - g(2) * t248;
t213 = qJDD(1) * pkin(1) + t223;
t224 = -g(1) * t248 - g(2) * t245;
t250 = qJD(1) ^ 2;
t214 = -pkin(1) * t250 + t224;
t244 = sin(qJ(2));
t247 = cos(qJ(2));
t159 = t244 * t213 + t247 * t214;
t232 = qJD(1) + qJD(2);
t230 = t232 ^ 2;
t231 = qJDD(1) + qJDD(2);
t156 = -pkin(2) * t230 + pkin(7) * t231 + t159;
t243 = sin(qJ(3));
t246 = cos(qJ(3));
t151 = -t246 * g(3) - t243 * t156;
t152 = -g(3) * t243 + t246 * t156;
t170 = Ifges(5,6) * qJD(3) + (Ifges(5,5) * t243 - Ifges(5,3) * t246) * t232;
t174 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t243 + Ifges(4,2) * t246) * t232;
t200 = (-mrSges(5,1) * t246 - mrSges(5,3) * t243) * t232;
t272 = qJD(3) * t246;
t268 = t232 * t272;
t203 = t231 * t243 + t268;
t278 = t232 * t243;
t204 = -qJD(3) * t278 + t231 * t246;
t199 = (-pkin(3) * t246 - qJ(4) * t243) * t232;
t249 = qJD(3) ^ 2;
t150 = -qJDD(3) * pkin(3) - qJ(4) * t249 + t199 * t278 + qJDD(4) - t151;
t271 = -0.2e1 * qJD(5) * t232;
t144 = t243 * t271 + (-t203 + t268) * qJ(5) + (-t230 * t243 * t246 - qJDD(3)) * pkin(4) + t150;
t201 = (mrSges(6,1) * t246 + mrSges(6,2) * t243) * t232;
t277 = t232 * t246;
t219 = qJD(3) * mrSges(6,2) - mrSges(6,3) * t277;
t135 = m(6) * t144 - qJDD(3) * mrSges(6,1) - t203 * mrSges(6,3) - qJD(3) * t219 - t201 * t278;
t282 = 2 * qJD(4);
t148 = -pkin(3) * t249 + qJDD(3) * qJ(4) + qJD(3) * t282 + t199 * t277 + t152;
t215 = -qJD(3) * pkin(4) - qJ(5) * t278;
t242 = t246 ^ 2;
t143 = -t242 * t230 * pkin(4) - qJ(5) * t204 + qJD(3) * t215 + t246 * t271 + t148;
t172 = -Ifges(6,6) * qJD(3) + (Ifges(6,4) * t243 - Ifges(6,2) * t246) * t232;
t175 = -Ifges(6,5) * qJD(3) + (Ifges(6,1) * t243 - Ifges(6,4) * t246) * t232;
t260 = mrSges(6,1) * t144 - mrSges(6,2) * t143 + Ifges(6,5) * t203 - Ifges(6,6) * t204 - Ifges(6,3) * qJDD(3) + t172 * t278 + t175 * t277;
t254 = -mrSges(5,1) * t150 + mrSges(5,3) * t148 + Ifges(5,4) * t203 + Ifges(5,2) * qJDD(3) - Ifges(5,6) * t204 - pkin(4) * t135 - t260;
t221 = mrSges(5,2) * t277 + qJD(3) * mrSges(5,3);
t257 = -m(5) * t150 + qJDD(3) * mrSges(5,1) + qJD(3) * t221 - t135;
t218 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t278;
t216 = -qJD(3) * mrSges(6,1) - mrSges(6,3) * t278;
t265 = m(6) * t143 + qJDD(3) * mrSges(6,2) - t204 * mrSges(6,3) + qJD(3) * t216;
t259 = m(5) * t148 + qJDD(3) * mrSges(5,3) + qJD(3) * t218 + t200 * t277 + t265;
t270 = t201 * t277;
t176 = Ifges(5,4) * qJD(3) + (Ifges(5,1) * t243 - Ifges(5,5) * t246) * t232;
t274 = t176 + Ifges(4,5) * qJD(3) + (Ifges(4,1) * t243 + Ifges(4,4) * t246) * t232;
t285 = -(t274 * t246 + (t170 - t174) * t243) * t232 + mrSges(4,1) * t151 - mrSges(4,2) * t152 + Ifges(4,5) * t203 + Ifges(4,6) * t204 + Ifges(4,3) * qJDD(3) + pkin(3) * (-mrSges(5,2) * t203 - t200 * t278 + t257) + qJ(4) * (mrSges(5,2) * t204 + t259 - t270) + t254;
t169 = -Ifges(6,3) * qJD(3) + (Ifges(6,5) * t243 - Ifges(6,6) * t246) * t232;
t284 = -Ifges(6,5) * qJDD(3) - t169 * t277;
t281 = t230 * pkin(7);
t280 = mrSges(4,3) + mrSges(5,2);
t279 = Ifges(5,6) - Ifges(6,6);
t202 = (-mrSges(4,1) * t246 + mrSges(4,2) * t243) * t232;
t217 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t278;
t128 = m(4) * t152 - qJDD(3) * mrSges(4,2) - qJD(3) * t217 + (-t201 + t202) * t277 + t280 * t204 + t259;
t220 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t277;
t129 = m(4) * t151 + qJDD(3) * mrSges(4,1) + qJD(3) * t220 + (-t200 - t202) * t278 - t280 * t203 + t257;
t266 = t246 * t128 - t129 * t243;
t119 = m(3) * t159 - mrSges(3,1) * t230 - mrSges(3,2) * t231 + t266;
t158 = t247 * t213 - t244 * t214;
t264 = t231 * pkin(2) + t158;
t261 = -t203 * qJ(4) - t264;
t138 = qJDD(5) + (-qJ(5) * t242 + pkin(7)) * t230 + (pkin(3) + pkin(4)) * t204 + (qJ(4) * t272 + (-pkin(3) * qJD(3) + t215 + t282) * t243) * t232 - t261;
t133 = -m(6) * t138 - t204 * mrSges(6,1) - t203 * mrSges(6,2) - t216 * t278 - t219 * t277;
t145 = -t204 * pkin(3) - t281 + (-0.2e1 * qJD(4) * t243 + (pkin(3) * t243 - qJ(4) * t246) * qJD(3)) * t232 + t261;
t130 = m(5) * t145 - mrSges(5,1) * t204 - t203 * mrSges(5,3) - t218 * t278 - t221 * t277 + t133;
t155 = -t264 - t281;
t252 = -m(4) * t155 + t204 * mrSges(4,1) - mrSges(4,2) * t203 - t217 * t278 + t220 * t277 - t130;
t123 = m(3) * t158 + mrSges(3,1) * t231 - mrSges(3,2) * t230 + t252;
t114 = t244 * t119 + t247 * t123;
t121 = t243 * t128 + t246 * t129;
t173 = Ifges(5,2) * qJD(3) + (Ifges(5,4) * t243 - Ifges(5,6) * t246) * t232;
t276 = -t169 + t173;
t267 = t247 * t119 - t244 * t123;
t263 = mrSges(6,1) * t138 - mrSges(6,3) * t143 - Ifges(6,4) * t203 + Ifges(6,2) * t204;
t262 = mrSges(6,2) * t138 - mrSges(6,3) * t144 + Ifges(6,1) * t203 - Ifges(6,4) * t204 + qJD(3) * t172;
t171 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t243 + Ifges(4,6) * t246) * t232;
t255 = mrSges(5,1) * t145 - mrSges(5,2) * t148 + pkin(4) * t133 + qJ(5) * (t265 - t270) - t263;
t110 = -t255 + (-t171 - t276) * t278 - pkin(3) * t130 + mrSges(4,3) * t152 - mrSges(4,1) * t155 + (t175 + t274) * qJD(3) + (Ifges(4,6) - t279) * qJDD(3) + (Ifges(4,4) - Ifges(5,5)) * t203 + (Ifges(4,2) + Ifges(5,3)) * t204;
t253 = mrSges(5,2) * t150 - mrSges(5,3) * t145 + Ifges(5,1) * t203 + Ifges(5,4) * qJDD(3) - Ifges(5,5) * t204 - qJ(5) * t135 + qJD(3) * t170 + t173 * t277 + t262;
t116 = t253 - qJ(4) * t130 - mrSges(4,3) * t151 + mrSges(4,2) * t155 + (-t169 + t171) * t277 - qJD(3) * t174 + Ifges(4,1) * t203 + Ifges(4,4) * t204 + (Ifges(4,5) - Ifges(6,5)) * qJDD(3);
t258 = mrSges(3,1) * t158 - mrSges(3,2) * t159 + Ifges(3,3) * t231 + pkin(2) * t252 + pkin(7) * t266 + t246 * t110 + t243 * t116;
t256 = mrSges(2,1) * t223 - mrSges(2,2) * t224 + Ifges(2,3) * qJDD(1) + pkin(1) * t114 + t258;
t112 = m(2) * t224 - mrSges(2,1) * t250 - qJDD(1) * mrSges(2,2) + t267;
t111 = m(2) * t223 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t250 + t114;
t108 = mrSges(3,1) * g(3) + mrSges(3,3) * t159 + t230 * Ifges(3,5) + Ifges(3,6) * t231 - pkin(2) * t121 - t285;
t107 = -mrSges(3,2) * g(3) - mrSges(3,3) * t158 + Ifges(3,5) * t231 - Ifges(3,6) * t230 - pkin(7) * t121 - t110 * t243 + t116 * t246;
t106 = -mrSges(2,2) * g(3) - mrSges(2,3) * t223 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t250 - pkin(6) * t114 + t107 * t247 - t108 * t244;
t105 = Ifges(2,6) * qJDD(1) + t250 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t224 + t244 * t107 + t247 * t108 - pkin(1) * (-m(3) * g(3) + t121) + pkin(6) * t267;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t248 * t106 - t245 * t105 - pkin(5) * (t111 * t248 + t112 * t245), t106, t107, t116, t253 + t284, t262 + t284; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t245 * t106 + t248 * t105 + pkin(5) * (-t111 * t245 + t112 * t248), t105, t108, t110, t254 + (-t243 * t170 - t246 * t176) * t232, -Ifges(6,6) * qJDD(3) - qJD(3) * t175 - t169 * t278 - t263; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t256, t256, t258, t285, Ifges(5,5) * t203 - Ifges(5,3) * t204 + t276 * t278 + t279 * qJDD(3) + (-t175 - t176) * qJD(3) + t255, t260;];
m_new = t1;
