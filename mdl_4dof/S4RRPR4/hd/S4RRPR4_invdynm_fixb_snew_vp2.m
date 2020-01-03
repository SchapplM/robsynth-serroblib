% Calculate vector of cutting torques with Newton-Euler for
% S4RRPR4
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
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 17:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4RRPR4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR4_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR4_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR4_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR4_invdynm_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR4_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR4_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR4_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:02:28
% EndTime: 2019-12-31 17:02:30
% DurationCPUTime: 1.46s
% Computational Cost: add. (20538->173), mult. (28482->223), div. (0->0), fcn. (16590->8), ass. (0->82)
t155 = qJD(1) + qJD(2);
t151 = t155 ^ 2;
t162 = sin(qJ(1));
t165 = cos(qJ(1));
t144 = t162 * g(1) - t165 * g(2);
t139 = qJDD(1) * pkin(1) + t144;
t145 = -t165 * g(1) - t162 * g(2);
t166 = qJD(1) ^ 2;
t140 = -t166 * pkin(1) + t145;
t161 = sin(qJ(2));
t164 = cos(qJ(2));
t127 = t161 * t139 + t164 * t140;
t152 = qJDD(1) + qJDD(2);
t124 = -t151 * pkin(2) + t152 * qJ(3) + t127;
t158 = sin(pkin(7));
t159 = cos(pkin(7));
t185 = qJD(3) * t155;
t184 = -t159 * g(3) - 0.2e1 * t158 * t185;
t110 = -t158 * t124 + t184;
t111 = -t158 * g(3) + (t124 + 0.2e1 * t185) * t159;
t190 = pkin(3) * t159;
t106 = (-pkin(6) * t152 + t151 * t190 - t124) * t158 + t184;
t187 = t152 * t159;
t154 = t159 ^ 2;
t188 = t151 * t154;
t107 = -pkin(3) * t188 + pkin(6) * t187 + t111;
t160 = sin(qJ(4));
t163 = cos(qJ(4));
t104 = t163 * t106 - t160 * t107;
t105 = t160 * t106 + t163 * t107;
t173 = -t158 * t160 + t159 * t163;
t130 = t173 * t155;
t174 = t158 * t163 + t159 * t160;
t131 = t174 * t155;
t113 = Ifges(5,4) * t131 + Ifges(5,2) * t130 + Ifges(5,6) * qJD(4);
t114 = Ifges(5,1) * t131 + Ifges(5,4) * t130 + Ifges(5,5) * qJD(4);
t122 = -t131 * qJD(4) + t173 * t152;
t123 = t130 * qJD(4) + t174 * t152;
t170 = -mrSges(5,1) * t104 + mrSges(5,2) * t105 - Ifges(5,5) * t123 - Ifges(5,6) * t122 - Ifges(5,3) * qJDD(4) - t131 * t113 + t130 * t114;
t179 = Ifges(4,4) * t158 + Ifges(4,2) * t159;
t180 = Ifges(4,1) * t158 + Ifges(4,4) * t159;
t120 = -t130 * mrSges(5,1) + t131 * mrSges(5,2);
t128 = -qJD(4) * mrSges(5,2) + t130 * mrSges(5,3);
t101 = m(5) * t104 + qJDD(4) * mrSges(5,1) - t123 * mrSges(5,3) + qJD(4) * t128 - t131 * t120;
t129 = qJD(4) * mrSges(5,1) - t131 * mrSges(5,3);
t102 = m(5) * t105 - qJDD(4) * mrSges(5,2) + t122 * mrSges(5,3) - qJD(4) * t129 + t130 * t120;
t92 = t163 * t101 + t160 * t102;
t191 = -mrSges(4,1) * t110 + mrSges(4,2) * t111 - pkin(3) * t92 - (t158 * t179 - t159 * t180) * t151 + t170;
t189 = mrSges(4,2) * t158;
t176 = mrSges(4,3) * t152 + (-mrSges(4,1) * t159 + t189) * t151;
t90 = m(4) * t110 - t176 * t158 + t92;
t181 = -t160 * t101 + t163 * t102;
t91 = m(4) * t111 + t176 * t159 + t181;
t183 = -t158 * t90 + t159 * t91;
t84 = m(3) * t127 - t151 * mrSges(3,1) - t152 * mrSges(3,2) + t183;
t126 = t164 * t139 - t161 * t140;
t177 = qJDD(3) - t126;
t121 = -t152 * pkin(2) - t151 * qJ(3) + t177;
t153 = t158 ^ 2;
t109 = (-pkin(2) - t190) * t152 + (-qJ(3) + (-t153 - t154) * pkin(6)) * t151 + t177;
t171 = m(5) * t109 - t122 * mrSges(5,1) + t123 * mrSges(5,2) - t130 * t128 + t131 * t129;
t168 = -m(4) * t121 + mrSges(4,1) * t187 - t171 + (t151 * t153 + t188) * mrSges(4,3);
t96 = m(3) * t126 - t151 * mrSges(3,2) + (mrSges(3,1) - t189) * t152 + t168;
t79 = t161 * t84 + t164 * t96;
t86 = t158 * t91 + t159 * t90;
t178 = Ifges(4,5) * t158 + Ifges(4,6) * t159;
t186 = t151 * t178;
t182 = -t161 * t96 + t164 * t84;
t112 = Ifges(5,5) * t131 + Ifges(5,6) * t130 + Ifges(5,3) * qJD(4);
t93 = -mrSges(5,1) * t109 + mrSges(5,3) * t105 + Ifges(5,4) * t123 + Ifges(5,2) * t122 + Ifges(5,6) * qJDD(4) + qJD(4) * t114 - t131 * t112;
t94 = mrSges(5,2) * t109 - mrSges(5,3) * t104 + Ifges(5,1) * t123 + Ifges(5,4) * t122 + Ifges(5,5) * qJDD(4) - qJD(4) * t113 + t130 * t112;
t75 = -mrSges(4,1) * t121 + mrSges(4,3) * t111 - pkin(3) * t171 + pkin(6) * t181 + t179 * t152 - t158 * t186 + t160 * t94 + t163 * t93;
t81 = mrSges(4,2) * t121 - mrSges(4,3) * t110 - pkin(6) * t92 + t180 * t152 + t159 * t186 - t160 * t93 + t163 * t94;
t172 = -mrSges(3,2) * t127 + qJ(3) * t183 + t158 * t81 + t159 * t75 + mrSges(3,1) * t126 + Ifges(3,3) * t152 + pkin(2) * (-t152 * t189 + t168);
t169 = mrSges(2,1) * t144 - mrSges(2,2) * t145 + Ifges(2,3) * qJDD(1) + pkin(1) * t79 + t172;
t77 = m(2) * t145 - t166 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t182;
t76 = m(2) * t144 + qJDD(1) * mrSges(2,1) - t166 * mrSges(2,2) + t79;
t73 = mrSges(3,1) * g(3) + (Ifges(3,6) - t178) * t152 + t151 * Ifges(3,5) + mrSges(3,3) * t127 - pkin(2) * t86 + t191;
t72 = -mrSges(3,2) * g(3) - mrSges(3,3) * t126 + Ifges(3,5) * t152 - t151 * Ifges(3,6) - qJ(3) * t86 - t158 * t75 + t159 * t81;
t71 = -mrSges(2,2) * g(3) - mrSges(2,3) * t144 + Ifges(2,5) * qJDD(1) - t166 * Ifges(2,6) - pkin(5) * t79 - t161 * t73 + t164 * t72;
t70 = Ifges(2,6) * qJDD(1) + t166 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t145 + t161 * t72 + t164 * t73 - pkin(1) * (-m(3) * g(3) + t86) + pkin(5) * t182;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t165 * t71 - t162 * t70 - pkin(4) * (t162 * t77 + t165 * t76), t71, t72, t81, t94; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t162 * t71 + t165 * t70 + pkin(4) * (-t162 * t76 + t165 * t77), t70, t73, t75, t93; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t169, t169, t172, t178 * t152 - t191, -t170;];
m_new = t1;
