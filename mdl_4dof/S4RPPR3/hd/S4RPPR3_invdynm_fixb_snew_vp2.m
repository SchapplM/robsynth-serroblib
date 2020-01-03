% Calculate vector of cutting torques with Newton-Euler for
% S4RPPR3
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
%   pkin=[a2,a3,a4,d1,d4,theta2,theta3]';
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
% Datum: 2019-12-31 16:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4RPPR3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR3_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR3_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR3_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPPR3_invdynm_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR3_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR3_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR3_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:37:52
% EndTime: 2019-12-31 16:37:53
% DurationCPUTime: 1.33s
% Computational Cost: add. (13652->172), mult. (28482->222), div. (0->0), fcn. (16590->8), ass. (0->80)
t166 = qJD(1) ^ 2;
t163 = sin(qJ(1));
t165 = cos(qJ(1));
t142 = t163 * g(1) - t165 * g(2);
t139 = qJDD(1) * pkin(1) + t142;
t143 = -t165 * g(1) - t163 * g(2);
t140 = -t166 * pkin(1) + t143;
t159 = sin(pkin(6));
t161 = cos(pkin(6));
t127 = t159 * t139 + t161 * t140;
t117 = -t166 * pkin(2) + qJDD(1) * qJ(3) + t127;
t158 = sin(pkin(7));
t157 = -g(3) + qJDD(2);
t160 = cos(pkin(7));
t185 = qJD(1) * qJD(3);
t187 = t160 * t157 - 0.2e1 * t158 * t185;
t109 = -t158 * t117 + t187;
t110 = t158 * t157 + (t117 + 0.2e1 * t185) * t160;
t190 = pkin(3) * t160;
t106 = (-pkin(5) * qJDD(1) + t166 * t190 - t117) * t158 + t187;
t184 = qJDD(1) * t160;
t153 = t160 ^ 2;
t188 = t153 * t166;
t107 = -pkin(3) * t188 + pkin(5) * t184 + t110;
t162 = sin(qJ(4));
t164 = cos(qJ(4));
t104 = t164 * t106 - t162 * t107;
t105 = t162 * t106 + t164 * t107;
t174 = -t158 * t162 + t160 * t164;
t130 = t174 * qJD(1);
t175 = t158 * t164 + t160 * t162;
t131 = t175 * qJD(1);
t115 = Ifges(5,4) * t131 + Ifges(5,2) * t130 + Ifges(5,6) * qJD(4);
t116 = Ifges(5,1) * t131 + Ifges(5,4) * t130 + Ifges(5,5) * qJD(4);
t123 = -t131 * qJD(4) + t174 * qJDD(1);
t124 = t130 * qJD(4) + t175 * qJDD(1);
t170 = -mrSges(5,1) * t104 + mrSges(5,2) * t105 - Ifges(5,5) * t124 - Ifges(5,6) * t123 - Ifges(5,3) * qJDD(4) - t131 * t115 + t130 * t116;
t179 = Ifges(4,4) * t158 + Ifges(4,2) * t160;
t180 = Ifges(4,1) * t158 + Ifges(4,4) * t160;
t119 = -t130 * mrSges(5,1) + t131 * mrSges(5,2);
t128 = -qJD(4) * mrSges(5,2) + t130 * mrSges(5,3);
t101 = m(5) * t104 + qJDD(4) * mrSges(5,1) - t124 * mrSges(5,3) + qJD(4) * t128 - t131 * t119;
t129 = qJD(4) * mrSges(5,1) - t131 * mrSges(5,3);
t102 = m(5) * t105 - qJDD(4) * mrSges(5,2) + t123 * mrSges(5,3) - qJD(4) * t129 + t130 * t119;
t92 = t164 * t101 + t162 * t102;
t191 = -mrSges(4,1) * t109 + mrSges(4,2) * t110 - pkin(3) * t92 - (t158 * t179 - t160 * t180) * t166 + t170;
t189 = mrSges(4,2) * t158;
t173 = mrSges(4,3) * qJDD(1) + t166 * (-mrSges(4,1) * t160 + t189);
t90 = m(4) * t109 - t173 * t158 + t92;
t181 = -t162 * t101 + t164 * t102;
t91 = m(4) * t110 + t173 * t160 + t181;
t183 = -t158 * t90 + t160 * t91;
t84 = m(3) * t127 - t166 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t183;
t126 = t161 * t139 - t159 * t140;
t177 = qJDD(3) - t126;
t113 = -qJDD(1) * pkin(2) - t166 * qJ(3) + t177;
t152 = t158 ^ 2;
t108 = (-pkin(2) - t190) * qJDD(1) + (-qJ(3) + (-t152 - t153) * pkin(5)) * t166 + t177;
t171 = m(5) * t108 - t123 * mrSges(5,1) + t124 * mrSges(5,2) - t130 * t128 + t131 * t129;
t168 = -m(4) * t113 + mrSges(4,1) * t184 - t171 + (t152 * t166 + t188) * mrSges(4,3);
t96 = m(3) * t126 - t166 * mrSges(3,2) + (mrSges(3,1) - t189) * qJDD(1) + t168;
t79 = t159 * t84 + t161 * t96;
t86 = t158 * t91 + t160 * t90;
t178 = Ifges(4,5) * t158 + Ifges(4,6) * t160;
t186 = t166 * t178;
t182 = -t159 * t96 + t161 * t84;
t114 = Ifges(5,5) * t131 + Ifges(5,6) * t130 + Ifges(5,3) * qJD(4);
t93 = -mrSges(5,1) * t108 + mrSges(5,3) * t105 + Ifges(5,4) * t124 + Ifges(5,2) * t123 + Ifges(5,6) * qJDD(4) + qJD(4) * t116 - t131 * t114;
t94 = mrSges(5,2) * t108 - mrSges(5,3) * t104 + Ifges(5,1) * t124 + Ifges(5,4) * t123 + Ifges(5,5) * qJDD(4) - qJD(4) * t115 + t130 * t114;
t75 = -mrSges(4,1) * t113 + mrSges(4,3) * t110 - pkin(3) * t171 + pkin(5) * t181 + t179 * qJDD(1) - t158 * t186 + t162 * t94 + t164 * t93;
t81 = mrSges(4,2) * t113 - mrSges(4,3) * t109 - pkin(5) * t92 + t180 * qJDD(1) + t160 * t186 - t162 * t93 + t164 * t94;
t172 = -mrSges(3,2) * t127 + qJ(3) * t183 + t158 * t81 + t160 * t75 + mrSges(3,1) * t126 + Ifges(3,3) * qJDD(1) + pkin(2) * (-qJDD(1) * t189 + t168);
t169 = mrSges(2,1) * t142 - mrSges(2,2) * t143 + Ifges(2,3) * qJDD(1) + pkin(1) * t79 + t172;
t77 = m(2) * t143 - t166 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t182;
t76 = m(2) * t142 + qJDD(1) * mrSges(2,1) - t166 * mrSges(2,2) + t79;
t73 = (Ifges(3,6) - t178) * qJDD(1) + t166 * Ifges(3,5) - mrSges(3,1) * t157 + mrSges(3,3) * t127 - pkin(2) * t86 + t191;
t72 = mrSges(3,2) * t157 - mrSges(3,3) * t126 + Ifges(3,5) * qJDD(1) - t166 * Ifges(3,6) - qJ(3) * t86 - t158 * t75 + t160 * t81;
t71 = -mrSges(2,2) * g(3) - mrSges(2,3) * t142 + Ifges(2,5) * qJDD(1) - t166 * Ifges(2,6) - qJ(2) * t79 - t159 * t73 + t161 * t72;
t70 = Ifges(2,6) * qJDD(1) + t166 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t143 + t159 * t72 + t161 * t73 - pkin(1) * (m(3) * t157 + t86) + qJ(2) * t182;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t165 * t71 - t163 * t70 - pkin(4) * (t163 * t77 + t165 * t76), t71, t72, t81, t94; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t163 * t71 + t165 * t70 + pkin(4) * (-t163 * t76 + t165 * t77), t70, t73, t75, t93; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t169, t169, t172, t178 * qJDD(1) - t191, -t170;];
m_new = t1;
