% Calculate vector of cutting torques with Newton-Euler for
% S5RPPRR9
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-31 18:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPPRR9_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR9_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR9_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR9_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR9_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR9_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR9_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR9_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR9_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:02:22
% EndTime: 2019-12-31 18:02:25
% DurationCPUTime: 2.04s
% Computational Cost: add. (26960->224), mult. (47091->275), div. (0->0), fcn. (21493->8), ass. (0->95)
t185 = qJD(1) ^ 2;
t179 = sin(qJ(1));
t182 = cos(qJ(1));
t159 = -g(1) * t182 - g(2) * t179;
t195 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t159;
t205 = -pkin(1) - pkin(2);
t136 = t205 * t185 + t195;
t158 = g(1) * t179 - t182 * g(2);
t194 = -qJ(2) * t185 + qJDD(2) - t158;
t139 = t205 * qJDD(1) + t194;
t175 = sin(pkin(8));
t176 = cos(pkin(8));
t121 = -t175 * t136 + t139 * t176;
t118 = qJDD(1) * pkin(3) - pkin(6) * t185 - t121;
t178 = sin(qJ(4));
t181 = cos(qJ(4));
t200 = qJD(1) * qJD(4);
t198 = t181 * t200;
t153 = -qJDD(1) * t178 - t198;
t199 = t178 * t200;
t154 = -qJDD(1) * t181 + t199;
t112 = (-t153 + t198) * pkin(7) + (-t154 - t199) * pkin(4) + t118;
t122 = t176 * t136 + t175 * t139;
t119 = -pkin(3) * t185 - qJDD(1) * pkin(6) + t122;
t171 = g(3) + qJDD(3);
t116 = t181 * t119 + t178 * t171;
t152 = (pkin(4) * t181 + pkin(7) * t178) * qJD(1);
t184 = qJD(4) ^ 2;
t201 = qJD(1) * t181;
t114 = -pkin(4) * t184 + qJDD(4) * pkin(7) - t152 * t201 + t116;
t177 = sin(qJ(5));
t180 = cos(qJ(5));
t110 = t112 * t180 - t114 * t177;
t202 = qJD(1) * t178;
t149 = qJD(4) * t180 + t177 * t202;
t129 = t149 * qJD(5) + qJDD(4) * t177 + t153 * t180;
t150 = qJD(4) * t177 - t180 * t202;
t130 = -t149 * mrSges(6,1) + mrSges(6,2) * t150;
t160 = qJD(5) + t201;
t134 = -mrSges(6,2) * t160 + t149 * mrSges(6,3);
t148 = qJDD(5) - t154;
t107 = m(6) * t110 + t148 * mrSges(6,1) - t129 * mrSges(6,3) - t130 * t150 + t134 * t160;
t111 = t112 * t177 + t114 * t180;
t128 = -qJD(5) * t150 + qJDD(4) * t180 - t153 * t177;
t135 = mrSges(6,1) * t160 - mrSges(6,3) * t150;
t108 = m(6) * t111 - t148 * mrSges(6,2) + t128 * mrSges(6,3) + t149 * t130 - t135 * t160;
t102 = -t107 * t177 + t180 * t108;
t203 = t171 * t181;
t113 = -qJDD(4) * pkin(4) - pkin(7) * t184 - t203 + (-qJD(1) * t152 + t119) * t178;
t123 = Ifges(6,5) * t150 + Ifges(6,6) * t149 + Ifges(6,3) * t160;
t125 = Ifges(6,1) * t150 + Ifges(6,4) * t149 + Ifges(6,5) * t160;
t103 = -mrSges(6,1) * t113 + mrSges(6,3) * t111 + Ifges(6,4) * t129 + Ifges(6,2) * t128 + Ifges(6,6) * t148 - t123 * t150 + t125 * t160;
t124 = Ifges(6,4) * t150 + Ifges(6,2) * t149 + Ifges(6,6) * t160;
t104 = mrSges(6,2) * t113 - mrSges(6,3) * t110 + Ifges(6,1) * t129 + Ifges(6,4) * t128 + Ifges(6,5) * t148 + t149 * t123 - t124 * t160;
t109 = -m(6) * t113 + t128 * mrSges(6,1) - t129 * mrSges(6,2) + t149 * t134 - t135 * t150;
t115 = -t119 * t178 + t203;
t144 = (Ifges(5,6) * qJD(4)) + (-Ifges(5,4) * t178 - Ifges(5,2) * t181) * qJD(1);
t145 = (Ifges(5,5) * qJD(4)) + (-Ifges(5,1) * t178 - Ifges(5,4) * t181) * qJD(1);
t206 = mrSges(5,1) * t115 - mrSges(5,2) * t116 + Ifges(5,5) * t153 + Ifges(5,6) * t154 + Ifges(5,3) * qJDD(4) + pkin(4) * t109 + pkin(7) * t102 - (t144 * t178 - t145 * t181) * qJD(1) + t180 * t103 + t177 * t104;
t204 = mrSges(2,1) + mrSges(3,1);
t151 = (mrSges(5,1) * t181 - mrSges(5,2) * t178) * qJD(1);
t157 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t201;
t105 = m(5) * t115 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t153 + qJD(4) * t157 + t151 * t202 + t109;
t156 = qJD(4) * mrSges(5,1) + mrSges(5,3) * t202;
t99 = m(5) * t116 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t154 - qJD(4) * t156 - t151 * t201 + t102;
t96 = -t105 * t178 + t181 * t99;
t92 = m(4) * t122 - mrSges(4,1) * t185 + qJDD(1) * mrSges(4,2) + t96;
t101 = t107 * t180 + t108 * t177;
t100 = -m(5) * t118 + t154 * mrSges(5,1) - mrSges(5,2) * t153 + t156 * t202 - t157 * t201 - t101;
t97 = m(4) * t121 - qJDD(1) * mrSges(4,1) - mrSges(4,2) * t185 + t100;
t88 = -t175 * t97 + t176 * t92;
t87 = t175 * t92 + t176 * t97;
t95 = t105 * t181 + t178 * t99;
t140 = -pkin(1) * t185 + t195;
t196 = m(3) * t140 + qJDD(1) * mrSges(3,3) + t88;
t94 = -m(4) * t171 - t95;
t142 = -qJDD(1) * pkin(1) + t194;
t193 = -m(3) * t142 + qJDD(1) * mrSges(3,1) + t185 * mrSges(3,3) - t87;
t143 = Ifges(5,3) * qJD(4) + (-Ifges(5,5) * t178 - Ifges(5,6) * t181) * qJD(1);
t89 = mrSges(5,2) * t118 - mrSges(5,3) * t115 + Ifges(5,1) * t153 + Ifges(5,4) * t154 + Ifges(5,5) * qJDD(4) - pkin(7) * t101 - qJD(4) * t144 - t103 * t177 + t104 * t180 - t143 * t201;
t189 = mrSges(6,1) * t110 - mrSges(6,2) * t111 + Ifges(6,5) * t129 + Ifges(6,6) * t128 + Ifges(6,3) * t148 + t124 * t150 - t149 * t125;
t90 = -mrSges(5,1) * t118 + mrSges(5,3) * t116 + Ifges(5,4) * t153 + Ifges(5,2) * t154 + Ifges(5,6) * qJDD(4) - pkin(4) * t101 + qJD(4) * t145 + t143 * t202 - t189;
t81 = mrSges(4,2) * t171 - mrSges(4,3) * t121 - Ifges(4,5) * qJDD(1) - Ifges(4,6) * t185 - pkin(6) * t95 - t178 * t90 + t181 * t89;
t82 = -mrSges(4,1) * t171 + mrSges(4,3) * t122 + t185 * Ifges(4,5) - Ifges(4,6) * qJDD(1) - pkin(3) * t95 - t206;
t192 = mrSges(3,2) * t142 + mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) + t185 * Ifges(3,6) - qJ(3) * t87 - t175 * t82 + t176 * t81;
t191 = mrSges(3,2) * t140 - pkin(2) * t94 - qJ(3) * t88 - t175 * t81 - t176 * t82;
t190 = mrSges(4,1) * t121 - mrSges(4,2) * t122 - Ifges(4,3) * qJDD(1) + pkin(3) * t100 + pkin(6) * t96 + t178 * t89 + t181 * t90;
t188 = -mrSges(3,1) * t142 + mrSges(3,3) * t140 + Ifges(3,2) * qJDD(1) - pkin(2) * t87 - t190;
t186 = -mrSges(2,2) * t159 + mrSges(2,1) * t158 + Ifges(2,3) * qJDD(1) + t188 + qJ(2) * (-mrSges(3,1) * t185 + t196) + pkin(1) * t193;
t93 = -m(3) * g(3) + t94;
t84 = m(2) * t158 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t185 + t193;
t83 = m(2) * t159 - qJDD(1) * mrSges(2,2) - t204 * t185 + t196;
t79 = -mrSges(2,2) * g(3) - mrSges(2,3) * t158 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t185 - qJ(2) * t93 + t192;
t78 = mrSges(2,3) * t159 - pkin(1) * t93 + (Ifges(3,4) + Ifges(2,5)) * t185 + (Ifges(2,6) - Ifges(3,6)) * qJDD(1) + t204 * g(3) + t191;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t182 * t79 - t179 * t78 - pkin(5) * (t179 * t83 + t182 * t84), t79, t192, t81, t89, t104; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t179 * t79 + t182 * t78 + pkin(5) * (-t179 * t84 + t182 * t83), t78, t188, t82, t90, t103; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t186, t186, -mrSges(3,1) * g(3) - Ifges(3,4) * t185 + Ifges(3,6) * qJDD(1) - t191, t190, t206, t189;];
m_new = t1;
