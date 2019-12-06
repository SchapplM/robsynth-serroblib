% Calculate vector of cutting torques with Newton-Euler for
% S5PPRRR2
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
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PPRRR2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR2_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR2_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRR2_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR2_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR2_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR2_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRR2_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:14:16
% EndTime: 2019-12-05 15:14:22
% DurationCPUTime: 3.09s
% Computational Cost: add. (39167->203), mult. (71173->265), div. (0->0), fcn. (46790->10), ass. (0->92)
t182 = sin(pkin(8));
t184 = cos(pkin(8));
t171 = -t184 * g(1) - t182 * g(2);
t180 = -g(3) + qJDD(1);
t181 = sin(pkin(9));
t183 = cos(pkin(9));
t154 = -t181 * t171 + t183 * t180;
t155 = t183 * t171 + t181 * t180;
t187 = sin(qJ(3));
t190 = cos(qJ(3));
t141 = t187 * t154 + t190 * t155;
t191 = qJD(3) ^ 2;
t135 = -t191 * pkin(3) + qJDD(3) * pkin(6) + t141;
t170 = t182 * g(1) - t184 * g(2);
t169 = qJDD(2) - t170;
t186 = sin(qJ(4));
t189 = cos(qJ(4));
t130 = -t186 * t135 + t189 * t169;
t207 = qJD(3) * qJD(4);
t206 = t189 * t207;
t166 = t186 * qJDD(3) + t206;
t127 = (-t166 + t206) * pkin(7) + (t186 * t189 * t191 + qJDD(4)) * pkin(4) + t130;
t131 = t189 * t135 + t186 * t169;
t167 = t189 * qJDD(3) - t186 * t207;
t209 = qJD(3) * t186;
t174 = qJD(4) * pkin(4) - pkin(7) * t209;
t179 = t189 ^ 2;
t128 = -t179 * t191 * pkin(4) + t167 * pkin(7) - qJD(4) * t174 + t131;
t185 = sin(qJ(5));
t188 = cos(qJ(5));
t125 = t188 * t127 - t185 * t128;
t159 = (-t185 * t186 + t188 * t189) * qJD(3);
t143 = t159 * qJD(5) + t188 * t166 + t185 * t167;
t160 = (t185 * t189 + t186 * t188) * qJD(3);
t148 = -t159 * mrSges(6,1) + t160 * mrSges(6,2);
t177 = qJD(4) + qJD(5);
t152 = -t177 * mrSges(6,2) + t159 * mrSges(6,3);
t176 = qJDD(4) + qJDD(5);
t122 = m(6) * t125 + t176 * mrSges(6,1) - t143 * mrSges(6,3) - t160 * t148 + t177 * t152;
t126 = t185 * t127 + t188 * t128;
t142 = -t160 * qJD(5) - t185 * t166 + t188 * t167;
t153 = t177 * mrSges(6,1) - t160 * mrSges(6,3);
t123 = m(6) * t126 - t176 * mrSges(6,2) + t142 * mrSges(6,3) + t159 * t148 - t177 * t153;
t112 = t188 * t122 + t185 * t123;
t157 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t186 + Ifges(5,2) * t189) * qJD(3);
t158 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t186 + Ifges(5,4) * t189) * qJD(3);
t145 = Ifges(6,4) * t160 + Ifges(6,2) * t159 + Ifges(6,6) * t177;
t146 = Ifges(6,1) * t160 + Ifges(6,4) * t159 + Ifges(6,5) * t177;
t195 = -mrSges(6,1) * t125 + mrSges(6,2) * t126 - Ifges(6,5) * t143 - Ifges(6,6) * t142 - Ifges(6,3) * t176 - t160 * t145 + t159 * t146;
t210 = mrSges(5,1) * t130 - mrSges(5,2) * t131 + Ifges(5,5) * t166 + Ifges(5,6) * t167 + Ifges(5,3) * qJDD(4) + pkin(4) * t112 + (t186 * t157 - t189 * t158) * qJD(3) - t195;
t165 = (-mrSges(5,1) * t189 + mrSges(5,2) * t186) * qJD(3);
t208 = qJD(3) * t189;
t173 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t208;
t110 = m(5) * t130 + qJDD(4) * mrSges(5,1) - t166 * mrSges(5,3) + qJD(4) * t173 - t165 * t209 + t112;
t172 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t209;
t202 = -t185 * t122 + t188 * t123;
t111 = m(5) * t131 - qJDD(4) * mrSges(5,2) + t167 * mrSges(5,3) - qJD(4) * t172 + t165 * t208 + t202;
t203 = -t186 * t110 + t189 * t111;
t101 = m(4) * t141 - t191 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t203;
t140 = t190 * t154 - t187 * t155;
t199 = -qJDD(3) * pkin(3) - t140;
t134 = -t191 * pkin(6) + t199;
t129 = t174 * t209 - t167 * pkin(4) + (-pkin(7) * t179 - pkin(6)) * t191 + t199;
t196 = m(6) * t129 - t142 * mrSges(6,1) + t143 * mrSges(6,2) - t159 * t152 + t160 * t153;
t193 = -m(5) * t134 + t167 * mrSges(5,1) - t166 * mrSges(5,2) - t172 * t209 + t173 * t208 - t196;
t116 = m(4) * t140 + qJDD(3) * mrSges(4,1) - t191 * mrSges(4,2) + t193;
t96 = t187 * t101 + t190 * t116;
t105 = t189 * t110 + t186 * t111;
t94 = m(3) * t154 + t96;
t204 = t190 * t101 - t187 * t116;
t95 = m(3) * t155 + t204;
t205 = -t181 * t94 + t183 * t95;
t200 = (-m(3) - m(4)) * t169 - t105;
t144 = Ifges(6,5) * t160 + Ifges(6,6) * t159 + Ifges(6,3) * t177;
t113 = -mrSges(6,1) * t129 + mrSges(6,3) * t126 + Ifges(6,4) * t143 + Ifges(6,2) * t142 + Ifges(6,6) * t176 - t160 * t144 + t177 * t146;
t114 = mrSges(6,2) * t129 - mrSges(6,3) * t125 + Ifges(6,1) * t143 + Ifges(6,4) * t142 + Ifges(6,5) * t176 + t159 * t144 - t177 * t145;
t156 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t186 + Ifges(5,6) * t189) * qJD(3);
t92 = -mrSges(5,1) * t134 + mrSges(5,3) * t131 + Ifges(5,4) * t166 + Ifges(5,2) * t167 + Ifges(5,6) * qJDD(4) - pkin(4) * t196 + pkin(7) * t202 + qJD(4) * t158 + t188 * t113 + t185 * t114 - t156 * t209;
t98 = mrSges(5,2) * t134 - mrSges(5,3) * t130 + Ifges(5,1) * t166 + Ifges(5,4) * t167 + Ifges(5,5) * qJDD(4) - pkin(7) * t112 - qJD(4) * t157 - t185 * t113 + t188 * t114 + t156 * t208;
t86 = mrSges(4,2) * t169 - mrSges(4,3) * t140 + Ifges(4,5) * qJDD(3) - t191 * Ifges(4,6) - pkin(6) * t105 - t186 * t92 + t189 * t98;
t90 = -mrSges(4,1) * t169 + mrSges(4,3) * t141 + t191 * Ifges(4,5) + Ifges(4,6) * qJDD(3) - pkin(3) * t105 - t210;
t82 = -mrSges(3,1) * t169 + mrSges(3,3) * t155 + t187 * t86 + t190 * t90 - pkin(2) * (m(4) * t169 + t105) + pkin(5) * t204;
t85 = mrSges(3,2) * t169 - mrSges(3,3) * t154 - pkin(5) * t96 - t187 * t90 + t190 * t86;
t198 = mrSges(2,1) * t170 - mrSges(2,2) * t171 + pkin(1) * t200 + qJ(2) * t205 + t181 * t85 + t183 * t82;
t197 = mrSges(4,1) * t140 - mrSges(4,2) * t141 + Ifges(4,3) * qJDD(3) + pkin(3) * t193 + pkin(6) * t203 + t186 * t98 + t189 * t92;
t194 = mrSges(3,1) * t154 - mrSges(3,2) * t155 + pkin(2) * t96 + t197;
t102 = m(2) * t170 + t200;
t89 = t181 * t95 + t183 * t94;
t87 = m(2) * t171 + t205;
t83 = -mrSges(2,1) * t180 + mrSges(2,3) * t171 - pkin(1) * t89 - t194;
t80 = mrSges(2,2) * t180 - mrSges(2,3) * t170 - qJ(2) * t89 - t181 * t82 + t183 * t85;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t184 * t80 - t182 * t83 - qJ(1) * (t184 * t102 + t182 * t87), t80, t85, t86, t98, t114; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t182 * t80 + t184 * t83 + qJ(1) * (-t182 * t102 + t184 * t87), t83, t82, t90, t92, t113; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t198, t198, t194, t197, t210, -t195;];
m_new = t1;
