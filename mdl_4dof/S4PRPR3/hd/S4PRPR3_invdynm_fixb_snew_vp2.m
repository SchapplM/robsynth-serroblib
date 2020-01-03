% Calculate vector of cutting torques with Newton-Euler for
% S4PRPR3
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-31 16:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4PRPR3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR3_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR3_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR3_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR3_invdynm_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR3_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR3_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR3_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:20:50
% EndTime: 2019-12-31 16:20:51
% DurationCPUTime: 1.23s
% Computational Cost: add. (12118->161), mult. (26173->211), div. (0->0), fcn. (16590->8), ass. (0->78)
t161 = qJD(2) ^ 2;
t154 = sin(pkin(6));
t156 = cos(pkin(6));
t139 = g(1) * t154 - g(2) * t156;
t140 = -g(1) * t156 - g(2) * t154;
t158 = sin(qJ(2));
t160 = cos(qJ(2));
t126 = t158 * t139 + t160 * t140;
t123 = -pkin(2) * t161 + qJDD(2) * qJ(3) + t126;
t153 = sin(pkin(7));
t152 = -g(3) + qJDD(1);
t155 = cos(pkin(7));
t180 = qJD(2) * qJD(3);
t182 = t155 * t152 - 0.2e1 * t153 * t180;
t109 = -t123 * t153 + t182;
t110 = t153 * t152 + (t123 + 0.2e1 * t180) * t155;
t185 = pkin(3) * t155;
t105 = (-pkin(5) * qJDD(2) + t161 * t185 - t123) * t153 + t182;
t179 = qJDD(2) * t155;
t149 = t155 ^ 2;
t183 = t149 * t161;
t106 = -pkin(3) * t183 + pkin(5) * t179 + t110;
t157 = sin(qJ(4));
t159 = cos(qJ(4));
t103 = t105 * t159 - t106 * t157;
t104 = t105 * t157 + t106 * t159;
t169 = -t153 * t157 + t155 * t159;
t129 = t169 * qJD(2);
t170 = t153 * t159 + t155 * t157;
t130 = t170 * qJD(2);
t112 = Ifges(5,4) * t130 + Ifges(5,2) * t129 + Ifges(5,6) * qJD(4);
t113 = Ifges(5,1) * t130 + Ifges(5,4) * t129 + Ifges(5,5) * qJD(4);
t121 = -t130 * qJD(4) + qJDD(2) * t169;
t122 = t129 * qJD(4) + qJDD(2) * t170;
t165 = -mrSges(5,1) * t103 + mrSges(5,2) * t104 - Ifges(5,5) * t122 - Ifges(5,6) * t121 - Ifges(5,3) * qJDD(4) - t130 * t112 + t129 * t113;
t174 = Ifges(4,4) * t153 + Ifges(4,2) * t155;
t175 = Ifges(4,1) * t153 + Ifges(4,4) * t155;
t116 = -mrSges(5,1) * t129 + mrSges(5,2) * t130;
t127 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t129;
t100 = m(5) * t103 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t122 + qJD(4) * t127 - t116 * t130;
t128 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t130;
t101 = m(5) * t104 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t121 - qJD(4) * t128 + t116 * t129;
t91 = t159 * t100 + t157 * t101;
t186 = -mrSges(4,1) * t109 + mrSges(4,2) * t110 - pkin(3) * t91 - (t153 * t174 - t155 * t175) * t161 + t165;
t184 = mrSges(4,2) * t153;
t168 = mrSges(4,3) * qJDD(2) + t161 * (-mrSges(4,1) * t155 + t184);
t89 = m(4) * t109 - t153 * t168 + t91;
t176 = -t157 * t100 + t159 * t101;
t90 = m(4) * t110 + t155 * t168 + t176;
t178 = -t153 * t89 + t155 * t90;
t83 = m(3) * t126 - mrSges(3,1) * t161 - qJDD(2) * mrSges(3,2) + t178;
t125 = t160 * t139 - t158 * t140;
t172 = qJDD(3) - t125;
t120 = -qJDD(2) * pkin(2) - t161 * qJ(3) + t172;
t148 = t153 ^ 2;
t108 = (-pkin(2) - t185) * qJDD(2) + (-qJ(3) + (-t148 - t149) * pkin(5)) * t161 + t172;
t166 = m(5) * t108 - t121 * mrSges(5,1) + t122 * mrSges(5,2) - t129 * t127 + t130 * t128;
t163 = -m(4) * t120 + mrSges(4,1) * t179 - t166 + (t148 * t161 + t183) * mrSges(4,3);
t95 = m(3) * t125 - t161 * mrSges(3,2) + (mrSges(3,1) - t184) * qJDD(2) + t163;
t78 = t158 * t83 + t160 * t95;
t85 = t153 * t90 + t155 * t89;
t173 = Ifges(4,5) * t153 + Ifges(4,6) * t155;
t181 = t161 * t173;
t177 = -t158 * t95 + t160 * t83;
t111 = Ifges(5,5) * t130 + Ifges(5,6) * t129 + Ifges(5,3) * qJD(4);
t92 = -mrSges(5,1) * t108 + mrSges(5,3) * t104 + Ifges(5,4) * t122 + Ifges(5,2) * t121 + Ifges(5,6) * qJDD(4) + qJD(4) * t113 - t111 * t130;
t93 = mrSges(5,2) * t108 - mrSges(5,3) * t103 + Ifges(5,1) * t122 + Ifges(5,4) * t121 + Ifges(5,5) * qJDD(4) - qJD(4) * t112 + t111 * t129;
t74 = -mrSges(4,1) * t120 + mrSges(4,3) * t110 - pkin(3) * t166 + pkin(5) * t176 + qJDD(2) * t174 - t153 * t181 + t157 * t93 + t159 * t92;
t80 = mrSges(4,2) * t120 - mrSges(4,3) * t109 - pkin(5) * t91 + qJDD(2) * t175 + t155 * t181 - t157 * t92 + t159 * t93;
t167 = -mrSges(3,2) * t126 + qJ(3) * t178 + t153 * t80 + t155 * t74 + mrSges(3,1) * t125 + Ifges(3,3) * qJDD(2) + pkin(2) * (-qJDD(2) * t184 + t163);
t164 = mrSges(2,1) * t139 - mrSges(2,2) * t140 + pkin(1) * t78 + t167;
t76 = m(2) * t140 + t177;
t75 = m(2) * t139 + t78;
t72 = (Ifges(3,6) - t173) * qJDD(2) + t161 * Ifges(3,5) - mrSges(3,1) * t152 + mrSges(3,3) * t126 - pkin(2) * t85 + t186;
t71 = mrSges(3,2) * t152 - mrSges(3,3) * t125 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t161 - qJ(3) * t85 - t153 * t74 + t155 * t80;
t70 = mrSges(2,2) * t152 - mrSges(2,3) * t139 - pkin(4) * t78 - t158 * t72 + t160 * t71;
t69 = -mrSges(2,1) * t152 + mrSges(2,3) * t140 + t158 * t71 + t160 * t72 - pkin(1) * (m(3) * t152 + t85) + pkin(4) * t177;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t156 * t70 - t154 * t69 - qJ(1) * (t154 * t76 + t156 * t75), t70, t71, t80, t93; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t154 * t70 + t156 * t69 + qJ(1) * (-t154 * t75 + t156 * t76), t69, t72, t74, t92; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t164, t164, t167, qJDD(2) * t173 - t186, -t165;];
m_new = t1;
