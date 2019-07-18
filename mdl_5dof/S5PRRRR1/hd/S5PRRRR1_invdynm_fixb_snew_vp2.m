% Calculate vector of cutting torques with Newton-Euler for
% S5PRRRR1
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% Datum: 2019-07-18 13:29
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PRRRR1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(2,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR1_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR1_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_invdynm_fixb_snew_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR1_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR1_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR1_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:28:23
% EndTime: 2019-07-18 13:28:26
% DurationCPUTime: 1.46s
% Computational Cost: add. (15331->220), mult. (29257->281), div. (0->0), fcn. (20780->8), ass. (0->84)
t157 = sin(qJ(4));
t158 = sin(qJ(3));
t161 = cos(qJ(4));
t162 = cos(qJ(3));
t139 = (t157 * t158 - t161 * t162) * qJD(2);
t155 = -g(3) + qJDD(1);
t159 = sin(qJ(2));
t163 = cos(qJ(2));
t147 = -t163 * g(1) + t159 * t155;
t133 = t162 * g(2) - t158 * t147;
t134 = t158 * g(2) + t162 * t147;
t137 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t158 + Ifges(4,2) * t162) * qJD(2);
t138 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t158 + Ifges(4,4) * t162) * qJD(2);
t180 = qJD(2) * qJD(3);
t144 = t158 * qJDD(2) + t162 * t180;
t179 = t158 * t180;
t145 = t162 * qJDD(2) - t179;
t166 = qJD(2) ^ 2;
t128 = (-t162 ^ 2 * t166 - qJD(3) ^ 2) * pkin(2) + t134;
t172 = (t158 * t162 * t166 + qJDD(3)) * pkin(2) + t133;
t111 = t157 * t128 - t161 * t172;
t112 = t161 * t128 + t157 * t172;
t140 = (t157 * t162 + t158 * t161) * qJD(2);
t118 = -t140 * qJD(4) - t157 * t144 + t161 * t145;
t119 = -t139 * qJD(4) + t161 * t144 + t157 * t145;
t154 = qJD(3) + qJD(4);
t123 = Ifges(5,4) * t140 - Ifges(5,2) * t139 + Ifges(5,6) * t154;
t124 = Ifges(5,1) * t140 - Ifges(5,4) * t139 + Ifges(5,5) * t154;
t153 = qJDD(3) + qJDD(4);
t156 = sin(qJ(5));
t160 = cos(qJ(5));
t146 = t159 * g(1) + t163 * t155;
t127 = (-t145 + t179) * pkin(2) - t146;
t104 = t160 * t112 + t156 * t127;
t130 = t160 * t140 + t156 * t154;
t105 = -t130 * qJD(5) - t156 * t119 + t160 * t153;
t129 = -t156 * t140 + t160 * t154;
t106 = t129 * qJD(5) + t160 * t119 + t156 * t153;
t135 = qJD(5) + t139;
t107 = Ifges(6,5) * t130 + Ifges(6,6) * t129 + Ifges(6,3) * t135;
t109 = Ifges(6,1) * t130 + Ifges(6,4) * t129 + Ifges(6,5) * t135;
t117 = qJDD(5) - t118;
t97 = -mrSges(6,1) * t111 + mrSges(6,3) * t104 + Ifges(6,4) * t106 + Ifges(6,2) * t105 + Ifges(6,6) * t117 - t130 * t107 + t135 * t109;
t103 = -t156 * t112 + t160 * t127;
t108 = Ifges(6,4) * t130 + Ifges(6,2) * t129 + Ifges(6,6) * t135;
t98 = mrSges(6,2) * t111 - mrSges(6,3) * t103 + Ifges(6,1) * t106 + Ifges(6,4) * t105 + Ifges(6,5) * t117 + t129 * t107 - t135 * t108;
t171 = mrSges(5,1) * t111 + mrSges(5,2) * t112 - Ifges(5,5) * t119 - Ifges(5,6) * t118 - Ifges(5,3) * t153 - t140 * t123 - t139 * t124 - t156 * t98 - t160 * t97;
t120 = -t135 * mrSges(6,2) + t129 * mrSges(6,3);
t121 = t135 * mrSges(6,1) - t130 * mrSges(6,3);
t125 = t139 * mrSges(5,1) + t140 * mrSges(5,2);
t131 = -t154 * mrSges(5,2) - t139 * mrSges(5,3);
t100 = t153 * mrSges(5,1) + t105 * mrSges(6,1) - t106 * mrSges(6,2) - t119 * mrSges(5,3) + t129 * t120 - t130 * t121 - t140 * t125 + t154 * t131 + (-m(5) - m(6)) * t111;
t113 = -t129 * mrSges(6,1) + t130 * mrSges(6,2);
t101 = m(6) * t103 + t117 * mrSges(6,1) - t106 * mrSges(6,3) - t130 * t113 + t135 * t120;
t102 = m(6) * t104 - t117 * mrSges(6,2) + t105 * mrSges(6,3) + t129 * t113 - t135 * t121;
t132 = t154 * mrSges(5,1) - t140 * mrSges(5,3);
t93 = m(5) * t112 - t153 * mrSges(5,2) + t118 * mrSges(5,3) - t156 * t101 + t160 * t102 - t139 * t125 - t154 * t132;
t89 = t161 * t100 + t157 * t93;
t183 = mrSges(4,1) * t133 - mrSges(4,2) * t134 + Ifges(4,5) * t144 + Ifges(4,6) * t145 + Ifges(4,3) * qJDD(3) + pkin(2) * t89 + (t158 * t137 - t162 * t138) * qJD(2) - t171;
t182 = qJD(2) * t158;
t181 = qJD(2) * t162;
t143 = (-mrSges(4,1) * t162 + mrSges(4,2) * t158) * qJD(2);
t149 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t181;
t87 = m(4) * t133 + qJDD(3) * mrSges(4,1) - t144 * mrSges(4,3) + qJD(3) * t149 - t143 * t182 + t89;
t148 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t182;
t88 = m(4) * t134 - qJDD(3) * mrSges(4,2) + t145 * mrSges(4,3) - qJD(3) * t148 - t157 * t100 + t143 * t181 + t161 * t93;
t177 = -t158 * t88 - t162 * t87;
t136 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t158 + Ifges(4,6) * t162) * qJD(2);
t168 = m(5) * t127 - t118 * mrSges(5,1) + t119 * mrSges(5,2) + t160 * t101 + t156 * t102 + t139 * t131 + t140 * t132;
t122 = Ifges(5,5) * t140 - Ifges(5,6) * t139 + Ifges(5,3) * t154;
t90 = mrSges(5,2) * t127 + mrSges(5,3) * t111 + Ifges(5,1) * t119 + Ifges(5,4) * t118 + Ifges(5,5) * t153 - t139 * t122 - t154 * t123 - t156 * t97 + t160 * t98;
t169 = mrSges(6,1) * t103 - mrSges(6,2) * t104 + Ifges(6,5) * t106 + Ifges(6,6) * t105 + Ifges(6,3) * t117 + t130 * t108 - t129 * t109;
t94 = -mrSges(5,1) * t127 + mrSges(5,3) * t112 + Ifges(5,4) * t119 + Ifges(5,2) * t118 + Ifges(5,6) * t153 - t140 * t122 + t154 * t124 - t169;
t81 = mrSges(4,1) * t146 + mrSges(4,3) * t134 + Ifges(4,4) * t144 + Ifges(4,2) * t145 + Ifges(4,6) * qJDD(3) - pkin(2) * t168 + qJD(3) * t138 - t136 * t182 + t157 * t90 + t161 * t94;
t86 = -mrSges(4,2) * t146 - mrSges(4,3) * t133 + Ifges(4,1) * t144 + Ifges(4,4) * t145 + Ifges(4,5) * qJDD(3) - qJD(3) * t137 + t136 * t181 - t157 * t94 + t161 * t90;
t79 = mrSges(3,2) * g(2) - mrSges(3,3) * t146 + Ifges(3,5) * qJDD(2) - t166 * Ifges(3,6) - t158 * t81 + t162 * t86;
t85 = -mrSges(3,1) * g(2) + mrSges(3,3) * t147 + t166 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - t183;
t178 = t159 * t79 + mrSges(2,2) * g(1) + pkin(1) * (-m(3) * g(2) + t177) + t163 * t85;
t174 = mrSges(2,2) * t155 + mrSges(2,3) * g(2) - t159 * t85 + t163 * t79;
t173 = mrSges(3,1) * t146 - mrSges(3,2) * t147 + Ifges(3,3) * qJDD(2) + t158 * t86 + t162 * t81;
t82 = m(3) * t147 - t166 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t158 * t87 + t162 * t88;
t91 = qJDD(2) * mrSges(3,1) + t145 * mrSges(4,1) - t166 * mrSges(3,2) - t144 * mrSges(4,2) + (m(3) + m(4)) * t146 + (-t148 * t158 + t149 * t162) * qJD(2) - t168;
t170 = -mrSges(2,1) * t155 - pkin(1) * (t159 * t82 + t163 * t91) - t173;
t1 = [-mrSges(1,2) * g(3) - qJ(1) * t177 + (mrSges(1,3) - qJ(1) * (-m(2) - m(3))) * g(2) + t174, t174, t79, t86, t90, t98; qJ(1) * (-t159 * t91 + t163 * t82) + mrSges(1,1) * g(3) + (-qJ(1) * m(2) - mrSges(1,3) - mrSges(2,3)) * g(1) + t170, -mrSges(2,3) * g(1) + t170, t85, t81, t94, t97; mrSges(1,2) * g(1) + (-mrSges(1,1) - mrSges(2,1)) * g(2) + t178, -mrSges(2,1) * g(2) + t178, t173, t183, -t171, t169;];
m_new  = t1;
