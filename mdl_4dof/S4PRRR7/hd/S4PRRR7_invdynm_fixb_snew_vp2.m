% Calculate vector of cutting torques with Newton-Euler for
% S4PRRR7
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% m [3x5]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4PRRR7_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(8,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR7_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR7_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR7_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR7_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4PRRR7_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR7_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR7_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR7_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:36:07
% EndTime: 2019-12-31 16:36:10
% DurationCPUTime: 2.07s
% Computational Cost: add. (24454->191), mult. (45864->255), div. (0->0), fcn. (29731->10), ass. (0->89)
t167 = cos(qJ(2));
t158 = sin(pkin(8));
t160 = cos(pkin(8));
t151 = t158 * g(1) - t160 * g(2);
t152 = -t160 * g(1) - t158 * g(2);
t157 = -g(3) + qJDD(1);
t161 = cos(pkin(4));
t164 = sin(qJ(2));
t183 = t161 * t164;
t159 = sin(pkin(4));
t184 = t159 * t164;
t122 = t151 * t183 + t167 * t152 + t157 * t184;
t134 = -t159 * t151 + t161 * t157;
t169 = qJD(2) ^ 2;
t120 = -t169 * pkin(2) + qJDD(2) * pkin(6) + t122;
t163 = sin(qJ(3));
t166 = cos(qJ(3));
t117 = t166 * t120 + t163 * t134;
t147 = (-pkin(3) * t166 - pkin(7) * t163) * qJD(2);
t168 = qJD(3) ^ 2;
t180 = t166 * qJD(2);
t114 = -t168 * pkin(3) + qJDD(3) * pkin(7) + t147 * t180 + t117;
t121 = -t164 * t152 + (t151 * t161 + t157 * t159) * t167;
t119 = -qJDD(2) * pkin(2) - t169 * pkin(6) - t121;
t179 = qJD(2) * qJD(3);
t177 = t166 * t179;
t148 = t163 * qJDD(2) + t177;
t178 = t163 * t179;
t149 = t166 * qJDD(2) - t178;
t115 = (-t148 - t177) * pkin(7) + (-t149 + t178) * pkin(3) + t119;
t162 = sin(qJ(4));
t165 = cos(qJ(4));
t111 = -t162 * t114 + t165 * t115;
t181 = qJD(2) * t163;
t144 = t165 * qJD(3) - t162 * t181;
t129 = t144 * qJD(4) + t162 * qJDD(3) + t165 * t148;
t145 = t162 * qJD(3) + t165 * t181;
t130 = -t144 * mrSges(5,1) + t145 * mrSges(5,2);
t156 = qJD(4) - t180;
t132 = -t156 * mrSges(5,2) + t144 * mrSges(5,3);
t141 = qJDD(4) - t149;
t108 = m(5) * t111 + t141 * mrSges(5,1) - t129 * mrSges(5,3) - t145 * t130 + t156 * t132;
t112 = t165 * t114 + t162 * t115;
t128 = -t145 * qJD(4) + t165 * qJDD(3) - t162 * t148;
t133 = t156 * mrSges(5,1) - t145 * mrSges(5,3);
t109 = m(5) * t112 - t141 * mrSges(5,2) + t128 * mrSges(5,3) + t144 * t130 - t156 * t133;
t102 = -t162 * t108 + t165 * t109;
t182 = t166 * t134;
t113 = -qJDD(3) * pkin(3) - t168 * pkin(7) - t182 + (qJD(2) * t147 + t120) * t163;
t123 = Ifges(5,5) * t145 + Ifges(5,6) * t144 + Ifges(5,3) * t156;
t125 = Ifges(5,1) * t145 + Ifges(5,4) * t144 + Ifges(5,5) * t156;
t103 = -mrSges(5,1) * t113 + mrSges(5,3) * t112 + Ifges(5,4) * t129 + Ifges(5,2) * t128 + Ifges(5,6) * t141 - t145 * t123 + t156 * t125;
t124 = Ifges(5,4) * t145 + Ifges(5,2) * t144 + Ifges(5,6) * t156;
t104 = mrSges(5,2) * t113 - mrSges(5,3) * t111 + Ifges(5,1) * t129 + Ifges(5,4) * t128 + Ifges(5,5) * t141 + t144 * t123 - t156 * t124;
t110 = -m(5) * t113 + t128 * mrSges(5,1) - t129 * mrSges(5,2) + t144 * t132 - t145 * t133;
t116 = -t163 * t120 + t182;
t137 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t163 + Ifges(4,2) * t166) * qJD(2);
t138 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t163 + Ifges(4,4) * t166) * qJD(2);
t188 = mrSges(4,1) * t116 - mrSges(4,2) * t117 + Ifges(4,5) * t148 + Ifges(4,6) * t149 + Ifges(4,3) * qJDD(3) + pkin(3) * t110 + pkin(7) * t102 + t165 * t103 + t162 * t104 + (t163 * t137 - t166 * t138) * qJD(2);
t146 = (-mrSges(4,1) * t166 + mrSges(4,2) * t163) * qJD(2);
t153 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t181;
t100 = m(4) * t117 - qJDD(3) * mrSges(4,2) + t149 * mrSges(4,3) - qJD(3) * t153 + t146 * t180 + t102;
t154 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t180;
t106 = m(4) * t116 + qJDD(3) * mrSges(4,1) - t148 * mrSges(4,3) + qJD(3) * t154 - t146 * t181 + t110;
t95 = t163 * t100 + t166 * t106;
t81 = -mrSges(3,1) * t134 + mrSges(3,3) * t122 + t169 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t95 - t188;
t176 = t166 * t100 - t163 * t106;
t93 = m(3) * t122 - t169 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t176;
t101 = t165 * t108 + t162 * t109;
t172 = -m(4) * t119 + t149 * mrSges(4,1) - t148 * mrSges(4,2) - t153 * t181 + t154 * t180 - t101;
t97 = m(3) * t121 + qJDD(2) * mrSges(3,1) - t169 * mrSges(3,2) + t172;
t88 = -t164 * t97 + t167 * t93;
t189 = pkin(5) * t88 + t167 * t81;
t185 = t167 * t97;
t94 = m(3) * t134 + t95;
t85 = -t159 * t94 + t161 * t185 + t93 * t183;
t136 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t163 + Ifges(4,6) * t166) * qJD(2);
t89 = mrSges(4,2) * t119 - mrSges(4,3) * t116 + Ifges(4,1) * t148 + Ifges(4,4) * t149 + Ifges(4,5) * qJDD(3) - pkin(7) * t101 - qJD(3) * t137 - t162 * t103 + t165 * t104 + t136 * t180;
t171 = mrSges(5,1) * t111 - mrSges(5,2) * t112 + Ifges(5,5) * t129 + Ifges(5,6) * t128 + Ifges(5,3) * t141 + t145 * t124 - t144 * t125;
t90 = -mrSges(4,1) * t119 + mrSges(4,3) * t117 + Ifges(4,4) * t148 + Ifges(4,2) * t149 + Ifges(4,6) * qJDD(3) - pkin(3) * t101 + qJD(3) * t138 - t136 * t181 - t171;
t77 = mrSges(3,1) * t121 - mrSges(3,2) * t122 + Ifges(3,3) * qJDD(2) + pkin(2) * t172 + pkin(6) * t176 + t163 * t89 + t166 * t90;
t79 = mrSges(3,2) * t134 - mrSges(3,3) * t121 + Ifges(3,5) * qJDD(2) - t169 * Ifges(3,6) - pkin(6) * t95 - t163 * t90 + t166 * t89;
t173 = mrSges(2,1) * t151 - mrSges(2,2) * t152 + pkin(1) * t85 + t189 * t159 + t161 * t77 + t79 * t184;
t86 = m(2) * t152 + t88;
t84 = t161 * t94 + (t164 * t93 + t185) * t159;
t82 = m(2) * t151 + t85;
t75 = mrSges(2,2) * t157 - mrSges(2,3) * t151 - t164 * t81 + t167 * t79 + (-t159 * t84 - t161 * t85) * pkin(5);
t74 = -mrSges(2,1) * t157 + mrSges(2,3) * t152 - pkin(1) * t84 - t159 * t77 + (t164 * t79 + t189) * t161;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t160 * t75 - t158 * t74 - qJ(1) * (t158 * t86 + t160 * t82), t75, t79, t89, t104; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t158 * t75 + t160 * t74 + qJ(1) * (-t158 * t82 + t160 * t86), t74, t81, t90, t103; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t173, t173, t77, t188, t171;];
m_new = t1;
