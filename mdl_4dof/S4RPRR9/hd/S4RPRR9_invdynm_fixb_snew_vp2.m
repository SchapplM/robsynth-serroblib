% Calculate vector of cutting torques with Newton-Euler for
% S4RPRR9
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% Datum: 2019-12-31 16:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4RPRR9_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR9_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR9_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR9_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR9_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR9_invdynm_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR9_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR9_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR9_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:56:10
% EndTime: 2019-12-31 16:56:12
% DurationCPUTime: 0.93s
% Computational Cost: add. (8686->192), mult. (16199->236), div. (0->0), fcn. (8256->6), ass. (0->82)
t160 = sin(qJ(1));
t163 = cos(qJ(1));
t147 = -t163 * g(1) - t160 * g(2);
t189 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t147;
t188 = (-pkin(1) - pkin(5));
t159 = sin(qJ(3));
t187 = t159 * g(3);
t186 = mrSges(2,1) - mrSges(3,2);
t185 = -Ifges(3,4) + Ifges(2,5);
t184 = (Ifges(3,5) - Ifges(2,6));
t158 = sin(qJ(4));
t161 = cos(qJ(4));
t165 = qJD(1) ^ 2;
t121 = (t188 * t165) - t189;
t162 = cos(qJ(3));
t181 = qJD(1) * qJD(3);
t178 = t162 * t181;
t141 = -t159 * qJDD(1) - t178;
t179 = t159 * t181;
t142 = t162 * qJDD(1) - t179;
t103 = (-t142 + t179) * pkin(6) + (-t141 + t178) * pkin(3) + t121;
t146 = t160 * g(1) - t163 * g(2);
t175 = -t165 * qJ(2) + qJDD(2) - t146;
t122 = t188 * qJDD(1) + t175;
t116 = -t162 * g(3) + t159 * t122;
t140 = (pkin(3) * t159 - pkin(6) * t162) * qJD(1);
t164 = qJD(3) ^ 2;
t182 = t159 * qJD(1);
t105 = -t164 * pkin(3) + qJDD(3) * pkin(6) - t140 * t182 + t116;
t101 = t161 * t103 - t158 * t105;
t183 = qJD(1) * t162;
t137 = t161 * qJD(3) - t158 * t183;
t112 = t137 * qJD(4) + t158 * qJDD(3) + t161 * t142;
t138 = t158 * qJD(3) + t161 * t183;
t114 = -t137 * mrSges(5,1) + t138 * mrSges(5,2);
t148 = qJD(4) + t182;
t117 = -t148 * mrSges(5,2) + t137 * mrSges(5,3);
t136 = qJDD(4) - t141;
t97 = m(5) * t101 + t136 * mrSges(5,1) - t112 * mrSges(5,3) - t138 * t114 + t148 * t117;
t102 = t158 * t103 + t161 * t105;
t111 = -t138 * qJD(4) + t161 * qJDD(3) - t158 * t142;
t118 = t148 * mrSges(5,1) - t138 * mrSges(5,3);
t98 = m(5) * t102 - t136 * mrSges(5,2) + t111 * mrSges(5,3) + t137 * t114 - t148 * t118;
t87 = t158 * t98 + t161 * t97;
t177 = -t158 * t97 + t161 * t98;
t139 = (mrSges(4,1) * t159 + mrSges(4,2) * t162) * qJD(1);
t145 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t183;
t85 = m(4) * t116 - qJDD(3) * mrSges(4,2) + t141 * mrSges(4,3) - qJD(3) * t145 - t139 * t182 + t177;
t115 = t162 * t122 + t187;
t144 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t182;
t104 = -qJDD(3) * pkin(3) - t164 * pkin(6) - t187 + (qJD(1) * t140 - t122) * t162;
t173 = -m(5) * t104 + t111 * mrSges(5,1) - t112 * mrSges(5,2) + t137 * t117 - t138 * t118;
t93 = m(4) * t115 + qJDD(3) * mrSges(4,1) - t142 * mrSges(4,3) + qJD(3) * t144 - t139 * t183 + t173;
t80 = -t159 * t93 + t162 * t85;
t79 = t159 * t85 + t162 * t93;
t127 = -qJDD(1) * pkin(1) + t175;
t174 = -m(3) * t127 + (t165 * mrSges(3,3)) - t79;
t83 = -m(4) * t121 + t141 * mrSges(4,1) - t142 * mrSges(4,2) - t144 * t182 - t145 * t183 - t87;
t125 = t165 * pkin(1) + t189;
t128 = (Ifges(4,3) * qJD(3)) + (Ifges(4,5) * t162 - Ifges(4,6) * t159) * qJD(1);
t129 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t162 - Ifges(4,2) * t159) * qJD(1);
t106 = Ifges(5,5) * t138 + Ifges(5,6) * t137 + Ifges(5,3) * t148;
t108 = Ifges(5,1) * t138 + Ifges(5,4) * t137 + Ifges(5,5) * t148;
t91 = -mrSges(5,1) * t104 + mrSges(5,3) * t102 + Ifges(5,4) * t112 + Ifges(5,2) * t111 + Ifges(5,6) * t136 - t138 * t106 + t148 * t108;
t107 = Ifges(5,4) * t138 + Ifges(5,2) * t137 + Ifges(5,6) * t148;
t92 = mrSges(5,2) * t104 - mrSges(5,3) * t101 + Ifges(5,1) * t112 + Ifges(5,4) * t111 + Ifges(5,5) * t136 + t137 * t106 - t148 * t107;
t74 = mrSges(4,2) * t121 - mrSges(4,3) * t115 + Ifges(4,1) * t142 + Ifges(4,4) * t141 + Ifges(4,5) * qJDD(3) - pkin(6) * t87 - qJD(3) * t129 - t128 * t182 - t158 * t91 + t161 * t92;
t130 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t162 - Ifges(4,4) * t159) * qJD(1);
t166 = mrSges(5,1) * t101 - mrSges(5,2) * t102 + Ifges(5,5) * t112 + Ifges(5,6) * t111 + Ifges(5,3) * t136 + t138 * t107 - t137 * t108;
t75 = -mrSges(4,1) * t121 + mrSges(4,3) * t116 + Ifges(4,4) * t142 + Ifges(4,2) * t141 + Ifges(4,6) * qJDD(3) - pkin(3) * t87 + qJD(3) * t130 - t128 * t183 - t166;
t172 = mrSges(3,2) * t127 - mrSges(3,3) * t125 + Ifges(3,1) * qJDD(1) - pkin(5) * t79 - t159 * t75 + t162 * t74;
t171 = -mrSges(3,1) * t125 - pkin(2) * t83 - pkin(5) * t80 - t159 * t74 - t162 * t75;
t170 = mrSges(4,1) * t115 - mrSges(4,2) * t116 + Ifges(4,5) * t142 + Ifges(4,6) * t141 + Ifges(4,3) * qJDD(3) + pkin(3) * t173 + pkin(6) * t177 + t129 * t183 + t130 * t182 + t158 * t92 + t161 * t91;
t169 = -m(3) * t125 + t165 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t83;
t168 = -mrSges(2,2) * t147 + mrSges(2,1) * t146 + Ifges(2,3) * qJDD(1) + t172 + pkin(1) * (-qJDD(1) * mrSges(3,2) + t174) + qJ(2) * t169;
t167 = mrSges(3,1) * t127 + pkin(2) * t79 + t170;
t81 = m(2) * t147 - t165 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t169;
t78 = -m(3) * g(3) + t80;
t76 = m(2) * t146 - t165 * mrSges(2,2) + t186 * qJDD(1) + t174;
t72 = (t184 * t165) + t185 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3) + t167 - qJ(2) * t78 - mrSges(2,3) * t146;
t71 = mrSges(2,3) * t147 - pkin(1) * t78 + t186 * g(3) - t184 * qJDD(1) + t185 * t165 + t171;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t163 * t72 - t160 * t71 - pkin(4) * (t160 * t81 + t163 * t76), t72, t172, t74, t92; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t160 * t72 + t163 * t71 + pkin(4) * (-t160 * t76 + t163 * t81), t71, -mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) - (t165 * Ifges(3,5)) - t167, t75, t91; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t168, t168, mrSges(3,2) * g(3) + t165 * Ifges(3,4) + Ifges(3,5) * qJDD(1) - t171, t170, t166;];
m_new = t1;
