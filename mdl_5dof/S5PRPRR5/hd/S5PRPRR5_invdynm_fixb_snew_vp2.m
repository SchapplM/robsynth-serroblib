% Calculate vector of cutting torques with Newton-Euler for
% S5PRPRR5
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:55
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PRPRR5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR5_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR5_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR5_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR5_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR5_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR5_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR5_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:53:36
% EndTime: 2019-12-05 15:53:48
% DurationCPUTime: 6.16s
% Computational Cost: add. (72722->237), mult. (165352->300), div. (0->0), fcn. (117821->10), ass. (0->108)
t228 = qJD(2) ^ 2;
t220 = sin(pkin(8));
t254 = cos(pkin(8));
t204 = -g(1) * t254 - g(2) * t220;
t218 = -g(3) + qJDD(1);
t224 = sin(qJ(2));
t227 = cos(qJ(2));
t190 = t227 * t204 + t224 * t218;
t185 = -pkin(2) * t228 + qJDD(2) * qJ(3) + t190;
t219 = sin(pkin(9));
t203 = g(1) * t220 - g(2) * t254;
t221 = cos(pkin(9));
t249 = qJD(2) * qJD(3);
t252 = -t221 * t203 - 0.2e1 * t219 * t249;
t256 = pkin(3) * t221;
t163 = (-pkin(6) * qJDD(2) + t228 * t256 - t185) * t219 + t252;
t167 = -t219 * t203 + (t185 + 0.2e1 * t249) * t221;
t248 = qJDD(2) * t221;
t215 = t221 ^ 2;
t253 = t215 * t228;
t164 = -pkin(3) * t253 + pkin(6) * t248 + t167;
t223 = sin(qJ(4));
t226 = cos(qJ(4));
t145 = t226 * t163 - t164 * t223;
t238 = t219 * t226 + t221 * t223;
t237 = -t219 * t223 + t221 * t226;
t192 = t237 * qJD(2);
t250 = qJD(4) * t192;
t182 = qJDD(2) * t238 + t250;
t193 = t238 * qJD(2);
t140 = (-t182 + t250) * pkin(7) + (t192 * t193 + qJDD(4)) * pkin(4) + t145;
t146 = t223 * t163 + t226 * t164;
t181 = -qJD(4) * t193 + qJDD(2) * t237;
t188 = qJD(4) * pkin(4) - pkin(7) * t193;
t191 = t192 ^ 2;
t141 = -pkin(4) * t191 + pkin(7) * t181 - qJD(4) * t188 + t146;
t222 = sin(qJ(5));
t225 = cos(qJ(5));
t138 = t140 * t225 - t141 * t222;
t174 = t192 * t225 - t193 * t222;
t153 = qJD(5) * t174 + t181 * t222 + t182 * t225;
t175 = t192 * t222 + t193 * t225;
t159 = -mrSges(6,1) * t174 + mrSges(6,2) * t175;
t216 = qJD(4) + qJD(5);
t168 = -mrSges(6,2) * t216 + mrSges(6,3) * t174;
t213 = qJDD(4) + qJDD(5);
t135 = m(6) * t138 + mrSges(6,1) * t213 - mrSges(6,3) * t153 - t159 * t175 + t168 * t216;
t139 = t140 * t222 + t141 * t225;
t152 = -qJD(5) * t175 + t181 * t225 - t182 * t222;
t169 = mrSges(6,1) * t216 - mrSges(6,3) * t175;
t136 = m(6) * t139 - mrSges(6,2) * t213 + mrSges(6,3) * t152 + t159 * t174 - t169 * t216;
t127 = t135 * t225 + t136 * t222;
t177 = -mrSges(5,1) * t192 + mrSges(5,2) * t193;
t186 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t192;
t124 = m(5) * t145 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t182 + qJD(4) * t186 - t177 * t193 + t127;
t187 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t193;
t245 = -t135 * t222 + t136 * t225;
t125 = m(5) * t146 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t181 - qJD(4) * t187 + t177 * t192 + t245;
t120 = t124 * t226 + t125 * t223;
t166 = -t185 * t219 + t252;
t172 = Ifges(5,4) * t193 + Ifges(5,2) * t192 + Ifges(5,6) * qJD(4);
t173 = Ifges(5,1) * t193 + Ifges(5,4) * t192 + Ifges(5,5) * qJD(4);
t155 = Ifges(6,4) * t175 + Ifges(6,2) * t174 + Ifges(6,6) * t216;
t156 = Ifges(6,1) * t175 + Ifges(6,4) * t174 + Ifges(6,5) * t216;
t234 = -mrSges(6,1) * t138 + mrSges(6,2) * t139 - Ifges(6,5) * t153 - Ifges(6,6) * t152 - Ifges(6,3) * t213 - t175 * t155 + t174 * t156;
t230 = -mrSges(5,1) * t145 + mrSges(5,2) * t146 - Ifges(5,5) * t182 - Ifges(5,6) * t181 - Ifges(5,3) * qJDD(4) - pkin(4) * t127 - t193 * t172 + t192 * t173 + t234;
t243 = Ifges(4,4) * t219 + Ifges(4,2) * t221;
t244 = Ifges(4,1) * t219 + Ifges(4,4) * t221;
t257 = -mrSges(4,1) * t166 + mrSges(4,2) * t167 - pkin(3) * t120 - (t219 * t243 - t221 * t244) * t228 + t230;
t255 = mrSges(4,2) * t219;
t242 = Ifges(4,5) * t219 + Ifges(4,6) * t221;
t251 = t228 * t242;
t236 = mrSges(4,3) * qJDD(2) + t228 * (-mrSges(4,1) * t221 + t255);
t118 = m(4) * t166 - t219 * t236 + t120;
t246 = -t124 * t223 + t125 * t226;
t119 = m(4) * t167 + t221 * t236 + t246;
t114 = -t118 * t219 + t119 * t221;
t110 = m(3) * t190 - mrSges(3,1) * t228 - qJDD(2) * mrSges(3,2) + t114;
t189 = -t224 * t204 + t218 * t227;
t241 = qJDD(3) - t189;
t184 = -qJDD(2) * pkin(2) - qJ(3) * t228 + t241;
t214 = t219 ^ 2;
t170 = (-pkin(2) - t256) * qJDD(2) + (-qJ(3) + (-t214 - t215) * pkin(6)) * t228 + t241;
t143 = -pkin(4) * t181 - pkin(7) * t191 + t188 * t193 + t170;
t240 = m(6) * t143 - t152 * mrSges(6,1) + t153 * mrSges(6,2) - t174 * t168 + t175 * t169;
t233 = m(5) * t170 - t181 * mrSges(5,1) + mrSges(5,2) * t182 - t192 * t186 + t187 * t193 + t240;
t231 = -m(4) * t184 + mrSges(4,1) * t248 - t233 + (t214 * t228 + t253) * mrSges(4,3);
t130 = t231 + (mrSges(3,1) - t255) * qJDD(2) + m(3) * t189 - mrSges(3,2) * t228;
t247 = t110 * t227 - t130 * t224;
t113 = t118 * t221 + t119 * t219;
t154 = Ifges(6,5) * t175 + Ifges(6,6) * t174 + Ifges(6,3) * t216;
t128 = -mrSges(6,1) * t143 + mrSges(6,3) * t139 + Ifges(6,4) * t153 + Ifges(6,2) * t152 + Ifges(6,6) * t213 - t154 * t175 + t156 * t216;
t129 = mrSges(6,2) * t143 - mrSges(6,3) * t138 + Ifges(6,1) * t153 + Ifges(6,4) * t152 + Ifges(6,5) * t213 + t154 * t174 - t155 * t216;
t171 = Ifges(5,5) * t193 + Ifges(5,6) * t192 + Ifges(5,3) * qJD(4);
t115 = -mrSges(5,1) * t170 + mrSges(5,3) * t146 + Ifges(5,4) * t182 + Ifges(5,2) * t181 + Ifges(5,6) * qJDD(4) - pkin(4) * t240 + pkin(7) * t245 + qJD(4) * t173 + t225 * t128 + t222 * t129 - t193 * t171;
t116 = mrSges(5,2) * t170 - mrSges(5,3) * t145 + Ifges(5,1) * t182 + Ifges(5,4) * t181 + Ifges(5,5) * qJDD(4) - pkin(7) * t127 - qJD(4) * t172 - t128 * t222 + t129 * t225 + t171 * t192;
t104 = -mrSges(4,1) * t184 + mrSges(4,3) * t167 - pkin(3) * t233 + pkin(6) * t246 + qJDD(2) * t243 + t226 * t115 + t223 * t116 - t219 * t251;
t105 = mrSges(4,2) * t184 - mrSges(4,3) * t166 - pkin(6) * t120 + qJDD(2) * t244 - t115 * t223 + t116 * t226 + t221 * t251;
t101 = -mrSges(3,2) * t203 - mrSges(3,3) * t189 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t228 - qJ(3) * t113 - t104 * t219 + t105 * t221;
t103 = (Ifges(3,6) - t242) * qJDD(2) + t228 * Ifges(3,5) + mrSges(3,1) * t203 + mrSges(3,3) * t190 - pkin(2) * t113 + t257;
t235 = -mrSges(2,2) * t204 + pkin(5) * t247 + t224 * t101 + t227 * t103 + pkin(1) * (m(3) * t203 - t113) + mrSges(2,1) * t203;
t232 = mrSges(3,1) * t189 - mrSges(3,2) * t190 + Ifges(3,3) * qJDD(2) + pkin(2) * (-qJDD(2) * t255 + t231) + qJ(3) * t114 + t104 * t221 + t105 * t219;
t111 = (m(2) + m(3)) * t203 - t113;
t108 = t110 * t224 + t130 * t227;
t106 = m(2) * t204 + t247;
t99 = -mrSges(2,1) * t218 + mrSges(2,3) * t204 - pkin(1) * t108 - t232;
t98 = mrSges(2,2) * t218 - mrSges(2,3) * t203 - pkin(5) * t108 + t101 * t227 - t103 * t224;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t254 * t98 - t220 * t99 - qJ(1) * (t106 * t220 + t111 * t254), t98, t101, t105, t116, t129; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t220 * t98 + t254 * t99 + qJ(1) * (t106 * t254 - t111 * t220), t99, t103, t104, t115, t128; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t235, t235, t232, qJDD(2) * t242 - t257, -t230, -t234;];
m_new = t1;
