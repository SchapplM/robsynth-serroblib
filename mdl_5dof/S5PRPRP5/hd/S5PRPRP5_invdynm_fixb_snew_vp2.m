% Calculate vector of cutting torques with Newton-Euler for
% S5PRPRP5
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
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PRPRP5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP5_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP5_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP5_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP5_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP5_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP5_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP5_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:37:39
% EndTime: 2019-12-05 15:37:45
% DurationCPUTime: 2.59s
% Computational Cost: add. (27025->234), mult. (59781->285), div. (0->0), fcn. (38376->8), ass. (0->101)
t224 = qJD(2) ^ 2;
t217 = sin(pkin(8));
t251 = qJD(2) * t217;
t219 = cos(pkin(8));
t262 = qJD(2) * t219;
t218 = sin(pkin(7));
t256 = cos(pkin(7));
t199 = -t256 * g(1) - t218 * g(2);
t216 = -g(3) + qJDD(1);
t221 = sin(qJ(2));
t222 = cos(qJ(2));
t186 = t222 * t199 + t221 * t216;
t178 = -t224 * pkin(2) + qJDD(2) * qJ(3) + t186;
t198 = t218 * g(1) - t256 * g(2);
t248 = qJD(2) * qJD(3);
t252 = -t219 * t198 - 0.2e1 * t217 * t248;
t259 = pkin(3) * t219;
t145 = (-pkin(6) * qJDD(2) + t224 * t259 - t178) * t217 + t252;
t150 = -t217 * t198 + (t178 + 0.2e1 * t248) * t219;
t246 = qJDD(2) * t219;
t210 = t219 ^ 2;
t255 = t210 * t224;
t146 = -pkin(3) * t255 + pkin(6) * t246 + t150;
t220 = sin(qJ(4));
t260 = cos(qJ(4));
t142 = t220 * t145 + t260 * t146;
t244 = t219 * t260;
t247 = qJDD(2) * t217;
t233 = t260 * t217 + t219 * t220;
t188 = t233 * qJD(2);
t249 = t188 * qJD(4);
t174 = -qJDD(2) * t244 + t220 * t247 + t249;
t182 = qJD(4) * mrSges(5,1) - t188 * mrSges(5,3);
t187 = -qJD(2) * t244 + t220 * t251;
t162 = t187 * pkin(4) - t188 * qJ(5);
t223 = qJD(4) ^ 2;
t135 = -t223 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) - t187 * t162 + t142;
t183 = -qJD(4) * mrSges(6,1) + t188 * mrSges(6,2);
t245 = m(6) * t135 + qJDD(4) * mrSges(6,3) + qJD(4) * t183;
t163 = t187 * mrSges(6,1) - t188 * mrSges(6,3);
t253 = -t187 * mrSges(5,1) - t188 * mrSges(5,2) - t163;
t258 = -mrSges(5,3) - mrSges(6,2);
t126 = m(5) * t142 - qJDD(4) * mrSges(5,2) - qJD(4) * t182 + t258 * t174 + t253 * t187 + t245;
t141 = t260 * t145 - t220 * t146;
t250 = t187 * qJD(4);
t175 = t233 * qJDD(2) - t250;
t181 = -qJD(4) * mrSges(5,2) - t187 * mrSges(5,3);
t137 = -qJDD(4) * pkin(4) - t223 * qJ(5) + t188 * t162 + qJDD(5) - t141;
t184 = -t187 * mrSges(6,2) + qJD(4) * mrSges(6,3);
t241 = -m(6) * t137 + qJDD(4) * mrSges(6,1) + qJD(4) * t184;
t127 = m(5) * t141 + qJDD(4) * mrSges(5,1) + qJD(4) * t181 + t258 * t175 + t253 * t188 + t241;
t120 = t220 * t126 + t260 * t127;
t149 = -t217 * t178 + t252;
t156 = Ifges(5,4) * t188 - Ifges(5,2) * t187 + Ifges(5,6) * qJD(4);
t158 = Ifges(5,1) * t188 - Ifges(5,4) * t187 + Ifges(5,5) * qJD(4);
t153 = Ifges(6,5) * t188 + Ifges(6,6) * qJD(4) + Ifges(6,3) * t187;
t157 = Ifges(6,1) * t188 + Ifges(6,4) * qJD(4) + Ifges(6,5) * t187;
t230 = mrSges(6,1) * t137 - mrSges(6,3) * t135 - Ifges(6,4) * t175 - Ifges(6,2) * qJDD(4) - Ifges(6,6) * t174 + t188 * t153 - t187 * t157;
t226 = mrSges(5,2) * t142 - t187 * t158 - qJ(5) * (-t174 * mrSges(6,2) - t187 * t163 + t245) - pkin(4) * (-t175 * mrSges(6,2) - t188 * t163 + t241) - mrSges(5,1) * t141 - t188 * t156 + Ifges(5,6) * t174 - Ifges(5,5) * t175 - Ifges(5,3) * qJDD(4) + t230;
t238 = Ifges(4,4) * t217 + Ifges(4,2) * t219;
t239 = Ifges(4,1) * t217 + Ifges(4,4) * t219;
t261 = -mrSges(4,1) * t149 + mrSges(4,2) * t150 - pkin(3) * t120 - (t238 * t251 - t239 * t262) * qJD(2) + t226;
t257 = mrSges(4,2) * t217;
t155 = Ifges(6,4) * t188 + Ifges(6,2) * qJD(4) + Ifges(6,6) * t187;
t254 = -Ifges(5,5) * t188 + Ifges(5,6) * t187 - Ifges(5,3) * qJD(4) - t155;
t234 = mrSges(4,3) * qJDD(2) + t224 * (-mrSges(4,1) * t219 + t257);
t118 = m(4) * t149 - t234 * t217 + t120;
t242 = t260 * t126 - t220 * t127;
t119 = m(4) * t150 + t234 * t219 + t242;
t114 = -t217 * t118 + t219 * t119;
t110 = m(3) * t186 - t224 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t114;
t185 = -t221 * t199 + t222 * t216;
t236 = qJDD(3) - t185;
t177 = -qJDD(2) * pkin(2) - t224 * qJ(3) + t236;
t209 = t217 ^ 2;
t151 = (-pkin(2) - t259) * qJDD(2) + (-qJ(3) + (-t209 - t210) * pkin(6)) * t224 + t236;
t139 = -0.2e1 * qJD(5) * t188 + (-t175 + t250) * qJ(5) + (t174 + t249) * pkin(4) + t151;
t128 = m(6) * t139 + t174 * mrSges(6,1) - t175 * mrSges(6,3) - t188 * t183 + t187 * t184;
t229 = m(5) * t151 + t174 * mrSges(5,1) + t175 * mrSges(5,2) + t187 * t181 + t188 * t182 + t128;
t227 = -m(4) * t177 + mrSges(4,1) * t246 - t229 + (t209 * t224 + t255) * mrSges(4,3);
t121 = -t224 * mrSges(3,2) + m(3) * t185 + (mrSges(3,1) - t257) * qJDD(2) + t227;
t243 = t222 * t110 - t221 * t121;
t240 = -mrSges(6,1) * t139 + mrSges(6,2) * t135;
t237 = Ifges(4,5) * t217 + Ifges(4,6) * t219;
t113 = t219 * t118 + t217 * t119;
t115 = -mrSges(5,1) * t151 + mrSges(5,3) * t142 - pkin(4) * t128 + t254 * t188 + (Ifges(5,4) - Ifges(6,5)) * t175 + (-Ifges(5,2) - Ifges(6,3)) * t174 + (Ifges(5,6) - Ifges(6,6)) * qJDD(4) + (t157 + t158) * qJD(4) + t240;
t231 = mrSges(6,2) * t137 - mrSges(6,3) * t139 + Ifges(6,1) * t175 + Ifges(6,4) * qJDD(4) + Ifges(6,5) * t174 + qJD(4) * t153;
t116 = mrSges(5,2) * t151 - mrSges(5,3) * t141 + Ifges(5,1) * t175 - Ifges(5,4) * t174 + Ifges(5,5) * qJDD(4) - qJ(5) * t128 - qJD(4) * t156 + t254 * t187 + t231;
t194 = t237 * qJD(2);
t104 = -mrSges(4,1) * t177 + mrSges(4,3) * t150 - pkin(3) * t229 + pkin(6) * t242 + t238 * qJDD(2) + t260 * t115 + t220 * t116 - t194 * t251;
t105 = mrSges(4,2) * t177 - mrSges(4,3) * t149 - pkin(6) * t120 + t239 * qJDD(2) - t220 * t115 + t260 * t116 + t194 * t262;
t101 = -mrSges(3,2) * t198 - mrSges(3,3) * t185 + Ifges(3,5) * qJDD(2) - t224 * Ifges(3,6) - qJ(3) * t113 - t217 * t104 + t219 * t105;
t103 = (Ifges(3,6) - t237) * qJDD(2) + t224 * Ifges(3,5) + mrSges(3,1) * t198 + mrSges(3,3) * t186 - pkin(2) * t113 + t261;
t232 = -mrSges(2,2) * t199 + pkin(5) * t243 + t221 * t101 + t222 * t103 + pkin(1) * (m(3) * t198 - t113) + mrSges(2,1) * t198;
t228 = mrSges(3,1) * t185 - mrSges(3,2) * t186 + Ifges(3,3) * qJDD(2) + pkin(2) * (-mrSges(4,2) * t247 + t227) + qJ(3) * t114 + t219 * t104 + t217 * t105;
t111 = (m(2) + m(3)) * t198 - t113;
t108 = t221 * t110 + t222 * t121;
t106 = m(2) * t199 + t243;
t99 = -mrSges(2,1) * t216 + mrSges(2,3) * t199 - pkin(1) * t108 - t228;
t98 = mrSges(2,2) * t216 - mrSges(2,3) * t198 - pkin(5) * t108 + t222 * t101 - t221 * t103;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t256 * t98 - t218 * t99 - qJ(1) * (t218 * t106 + t256 * t111), t98, t101, t105, t116, -t187 * t155 + t231; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t218 * t98 + t256 * t99 + qJ(1) * (t256 * t106 - t218 * t111), t99, t103, t104, t115, -t230; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t232, t232, t228, t237 * qJDD(2) - t261, -t226, Ifges(6,5) * t175 + Ifges(6,6) * qJDD(4) + Ifges(6,3) * t174 - qJD(4) * t157 + t188 * t155 - t240;];
m_new = t1;
