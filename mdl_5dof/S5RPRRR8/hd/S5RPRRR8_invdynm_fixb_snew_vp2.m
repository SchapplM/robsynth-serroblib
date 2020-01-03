% Calculate vector of cutting torques with Newton-Euler for
% S5RPRRR8
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
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
% Datum: 2019-12-31 19:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPRRR8_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR8_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR8_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR8_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR8_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR8_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR8_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR8_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR8_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:05:46
% EndTime: 2019-12-31 19:05:51
% DurationCPUTime: 2.62s
% Computational Cost: add. (42237->227), mult. (54576->280), div. (0->0), fcn. (25817->8), ass. (0->95)
t197 = qJD(1) ^ 2;
t191 = sin(qJ(1));
t195 = cos(qJ(1));
t168 = -t195 * g(1) - t191 * g(2);
t208 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t168;
t218 = -pkin(1) - pkin(2);
t146 = t218 * t197 + t208;
t167 = t191 * g(1) - t195 * g(2);
t207 = -t197 * qJ(2) + qJDD(2) - t167;
t149 = t218 * qJDD(1) + t207;
t190 = sin(qJ(3));
t194 = cos(qJ(3));
t135 = t194 * t146 + t190 * t149;
t178 = -qJD(1) + qJD(3);
t174 = t178 ^ 2;
t176 = -qJDD(1) + qJDD(3);
t132 = -(t174 * pkin(3)) + t176 * pkin(7) + t135;
t189 = sin(qJ(4));
t193 = cos(qJ(4));
t123 = t193 * g(3) - t189 * t132;
t214 = qJD(4) * t178;
t213 = t193 * t214;
t161 = t189 * t176 + t213;
t119 = (-t161 + t213) * pkin(8) + (t174 * t189 * t193 + qJDD(4)) * pkin(4) + t123;
t124 = t189 * g(3) + t193 * t132;
t162 = t193 * t176 - t189 * t214;
t216 = t178 * t189;
t165 = qJD(4) * pkin(4) - pkin(8) * t216;
t185 = t193 ^ 2;
t120 = -t185 * t174 * pkin(4) + t162 * pkin(8) - qJD(4) * t165 + t124;
t188 = sin(qJ(5));
t192 = cos(qJ(5));
t117 = t192 * t119 - t188 * t120;
t153 = (-t188 * t189 + t192 * t193) * t178;
t130 = t153 * qJD(5) + t192 * t161 + t188 * t162;
t154 = (t188 * t193 + t189 * t192) * t178;
t140 = -t153 * mrSges(6,1) + t154 * mrSges(6,2);
t177 = qJD(4) + qJD(5);
t141 = -t177 * mrSges(6,2) + t153 * mrSges(6,3);
t175 = qJDD(4) + qJDD(5);
t114 = m(6) * t117 + t175 * mrSges(6,1) - t130 * mrSges(6,3) - t154 * t140 + t177 * t141;
t118 = t188 * t119 + t192 * t120;
t129 = -t154 * qJD(5) - t188 * t161 + t192 * t162;
t142 = t177 * mrSges(6,1) - t154 * mrSges(6,3);
t115 = m(6) * t118 - t175 * mrSges(6,2) + t129 * mrSges(6,3) + t153 * t140 - t177 * t142;
t106 = t192 * t114 + t188 * t115;
t151 = (Ifges(5,6) * qJD(4)) + (Ifges(5,4) * t189 + Ifges(5,2) * t193) * t178;
t152 = (Ifges(5,5) * qJD(4)) + (Ifges(5,1) * t189 + Ifges(5,4) * t193) * t178;
t137 = Ifges(6,4) * t154 + Ifges(6,2) * t153 + Ifges(6,6) * t177;
t138 = Ifges(6,1) * t154 + Ifges(6,4) * t153 + Ifges(6,5) * t177;
t203 = -mrSges(6,1) * t117 + mrSges(6,2) * t118 - Ifges(6,5) * t130 - Ifges(6,6) * t129 - Ifges(6,3) * t175 - t154 * t137 + t153 * t138;
t219 = mrSges(5,1) * t123 - mrSges(5,2) * t124 + Ifges(5,5) * t161 + Ifges(5,6) * t162 + Ifges(5,3) * qJDD(4) + pkin(4) * t106 + (t189 * t151 - t193 * t152) * t178 - t203;
t217 = mrSges(2,1) + mrSges(3,1);
t215 = t178 * t193;
t134 = -t190 * t146 + t194 * t149;
t209 = -t176 * pkin(3) - t134;
t131 = -t174 * pkin(7) + t209;
t163 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t216;
t164 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t215;
t121 = t165 * t216 - t162 * pkin(4) + (-pkin(8) * t185 - pkin(7)) * t174 + t209;
t205 = m(6) * t121 - t129 * mrSges(6,1) + t130 * mrSges(6,2) - t153 * t141 + t154 * t142;
t110 = -m(5) * t131 + t162 * mrSges(5,1) - t161 * mrSges(5,2) - t163 * t216 + t164 * t215 - t205;
t109 = m(4) * t134 + t176 * mrSges(4,1) - (t174 * mrSges(4,2)) + t110;
t160 = (-mrSges(5,1) * t193 + mrSges(5,2) * t189) * t178;
t104 = m(5) * t123 + qJDD(4) * mrSges(5,1) - t161 * mrSges(5,3) + qJD(4) * t164 - t160 * t216 + t106;
t212 = -t188 * t114 + t192 * t115;
t105 = m(5) * t124 - qJDD(4) * mrSges(5,2) + t162 * mrSges(5,3) - qJD(4) * t163 + t160 * t215 + t212;
t102 = -t189 * t104 + t193 * t105;
t98 = m(4) * t135 - t174 * mrSges(4,1) - t176 * mrSges(4,2) + t102;
t95 = -t190 * t109 + t194 * t98;
t94 = t194 * t109 + t190 * t98;
t101 = t193 * t104 + t189 * t105;
t155 = -t197 * pkin(1) + t208;
t210 = m(3) * t155 + qJDD(1) * mrSges(3,3) + t95;
t159 = -qJDD(1) * pkin(1) + t207;
t206 = -m(3) * t159 + qJDD(1) * mrSges(3,1) + t197 * mrSges(3,3) - t94;
t136 = Ifges(6,5) * t154 + Ifges(6,6) * t153 + Ifges(6,3) * t177;
t107 = -mrSges(6,1) * t121 + mrSges(6,3) * t118 + Ifges(6,4) * t130 + Ifges(6,2) * t129 + Ifges(6,6) * t175 - t154 * t136 + t177 * t138;
t108 = mrSges(6,2) * t121 - mrSges(6,3) * t117 + Ifges(6,1) * t130 + Ifges(6,4) * t129 + Ifges(6,5) * t175 + t153 * t136 - t177 * t137;
t150 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t189 + Ifges(5,6) * t193) * t178;
t89 = -mrSges(5,1) * t131 + mrSges(5,3) * t124 + Ifges(5,4) * t161 + Ifges(5,2) * t162 + Ifges(5,6) * qJDD(4) - pkin(4) * t205 + pkin(8) * t212 + qJD(4) * t152 + t192 * t107 + t188 * t108 - t150 * t216;
t96 = mrSges(5,2) * t131 - mrSges(5,3) * t123 + Ifges(5,1) * t161 + Ifges(5,4) * t162 + Ifges(5,5) * qJDD(4) - pkin(8) * t106 - qJD(4) * t151 - t188 * t107 + t192 * t108 + t150 * t215;
t87 = mrSges(4,2) * g(3) - mrSges(4,3) * t134 + Ifges(4,5) * t176 - (t174 * Ifges(4,6)) - pkin(7) * t101 - t189 * t89 + t193 * t96;
t88 = -mrSges(4,1) * g(3) + mrSges(4,3) * t135 + t174 * Ifges(4,5) + Ifges(4,6) * t176 - pkin(3) * t101 - t219;
t204 = mrSges(3,2) * t159 + mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) + t197 * Ifges(3,6) - pkin(6) * t94 - t190 * t88 + t194 * t87;
t202 = mrSges(3,2) * t155 - pkin(2) * (-m(4) * g(3) - t101) - pkin(6) * t95 - t190 * t87 - t194 * t88;
t201 = mrSges(4,1) * t134 - mrSges(4,2) * t135 + Ifges(4,3) * t176 + pkin(3) * t110 + pkin(7) * t102 + t189 * t96 + t193 * t89;
t200 = -mrSges(3,1) * t159 + mrSges(3,3) * t155 + Ifges(3,2) * qJDD(1) - pkin(2) * t94 - t201;
t198 = -mrSges(2,2) * t168 + mrSges(2,1) * t167 + Ifges(2,3) * qJDD(1) + t200 + qJ(2) * (-t197 * mrSges(3,1) + t210) + pkin(1) * t206;
t99 = (-m(3) - m(4)) * g(3) - t101;
t91 = m(2) * t167 + qJDD(1) * mrSges(2,1) - t197 * mrSges(2,2) + t206;
t90 = m(2) * t168 - qJDD(1) * mrSges(2,2) - t217 * t197 + t210;
t85 = -mrSges(2,2) * g(3) - mrSges(2,3) * t167 + Ifges(2,5) * qJDD(1) - t197 * Ifges(2,6) - qJ(2) * t99 + t204;
t84 = mrSges(2,3) * t168 - pkin(1) * t99 + (Ifges(3,4) + Ifges(2,5)) * t197 + (Ifges(2,6) - Ifges(3,6)) * qJDD(1) + t217 * g(3) + t202;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t195 * t85 - t191 * t84 - pkin(5) * (t191 * t90 + t195 * t91), t85, t204, t87, t96, t108; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t191 * t85 + t195 * t84 + pkin(5) * (-t191 * t91 + t195 * t90), t84, t200, t88, t89, t107; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t198, t198, -mrSges(3,1) * g(3) - t197 * Ifges(3,4) + Ifges(3,6) * qJDD(1) - t202, t201, t219, -t203;];
m_new = t1;
