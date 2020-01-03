% Calculate vector of cutting torques with Newton-Euler for
% S4PRRR4
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
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4PRRR4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR4_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR4_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR4_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR4_invdynm_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR4_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR4_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR4_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:32:32
% EndTime: 2019-12-31 16:32:34
% DurationCPUTime: 1.37s
% Computational Cost: add. (15096->182), mult. (29853->240), div. (0->0), fcn. (18128->8), ass. (0->79)
t157 = sin(pkin(7));
t158 = cos(pkin(7));
t144 = t157 * g(1) - t158 * g(2);
t145 = -t158 * g(1) - t157 * g(2);
t161 = sin(qJ(2));
t164 = cos(qJ(2));
t127 = t161 * t144 + t164 * t145;
t165 = qJD(2) ^ 2;
t124 = -t165 * pkin(2) + qJDD(2) * pkin(5) + t127;
t156 = -g(3) + qJDD(1);
t160 = sin(qJ(3));
t163 = cos(qJ(3));
t115 = -t160 * t124 + t163 * t156;
t116 = t163 * t124 + t160 * t156;
t131 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t160 + Ifges(4,2) * t163) * qJD(2);
t132 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t160 + Ifges(4,4) * t163) * qJD(2);
t178 = qJD(2) * qJD(3);
t177 = t163 * t178;
t141 = t160 * qJDD(2) + t177;
t142 = t163 * qJDD(2) - t160 * t178;
t106 = (-t141 + t177) * pkin(6) + (t160 * t163 * t165 + qJDD(3)) * pkin(3) + t115;
t180 = qJD(2) * t160;
t148 = qJD(3) * pkin(3) - pkin(6) * t180;
t155 = t163 ^ 2;
t107 = -t155 * t165 * pkin(3) + t142 * pkin(6) - qJD(3) * t148 + t116;
t159 = sin(qJ(4));
t162 = cos(qJ(4));
t104 = t162 * t106 - t159 * t107;
t105 = t159 * t106 + t162 * t107;
t134 = (t159 * t163 + t160 * t162) * qJD(2);
t113 = -t134 * qJD(4) - t159 * t141 + t162 * t142;
t133 = (-t159 * t160 + t162 * t163) * qJD(2);
t114 = t133 * qJD(4) + t162 * t141 + t159 * t142;
t153 = qJD(3) + qJD(4);
t118 = Ifges(5,4) * t134 + Ifges(5,2) * t133 + Ifges(5,6) * t153;
t119 = Ifges(5,1) * t134 + Ifges(5,4) * t133 + Ifges(5,5) * t153;
t152 = qJDD(3) + qJDD(4);
t169 = -mrSges(5,1) * t104 + mrSges(5,2) * t105 - Ifges(5,5) * t114 - Ifges(5,6) * t113 - Ifges(5,3) * t152 - t134 * t118 + t133 * t119;
t122 = -t133 * mrSges(5,1) + t134 * mrSges(5,2);
t128 = -t153 * mrSges(5,2) + t133 * mrSges(5,3);
t101 = m(5) * t104 + t152 * mrSges(5,1) - t114 * mrSges(5,3) - t134 * t122 + t153 * t128;
t129 = t153 * mrSges(5,1) - t134 * mrSges(5,3);
t102 = m(5) * t105 - t152 * mrSges(5,2) + t113 * mrSges(5,3) + t133 * t122 - t153 * t129;
t92 = t162 * t101 + t159 * t102;
t181 = mrSges(4,1) * t115 - mrSges(4,2) * t116 + Ifges(4,5) * t141 + Ifges(4,6) * t142 + Ifges(4,3) * qJDD(3) + pkin(3) * t92 + (t160 * t131 - t163 * t132) * qJD(2) - t169;
t140 = (-mrSges(4,1) * t163 + mrSges(4,2) * t160) * qJD(2);
t179 = qJD(2) * t163;
t147 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t179;
t90 = m(4) * t115 + qJDD(3) * mrSges(4,1) - t141 * mrSges(4,3) + qJD(3) * t147 - t140 * t180 + t92;
t146 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t180;
t174 = -t159 * t101 + t162 * t102;
t91 = m(4) * t116 - qJDD(3) * mrSges(4,2) + t142 * mrSges(4,3) - qJD(3) * t146 + t140 * t179 + t174;
t176 = -t160 * t90 + t163 * t91;
t84 = m(3) * t127 - t165 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t176;
t126 = t164 * t144 - t161 * t145;
t172 = -qJDD(2) * pkin(2) - t126;
t123 = -t165 * pkin(5) + t172;
t108 = t148 * t180 - t142 * pkin(3) + (-pkin(6) * t155 - pkin(5)) * t165 + t172;
t170 = m(5) * t108 - t113 * mrSges(5,1) + t114 * mrSges(5,2) - t133 * t128 + t134 * t129;
t167 = -m(4) * t123 + t142 * mrSges(4,1) - t141 * mrSges(4,2) - t146 * t180 + t147 * t179 - t170;
t96 = m(3) * t126 + qJDD(2) * mrSges(3,1) - t165 * mrSges(3,2) + t167;
t79 = t161 * t84 + t164 * t96;
t86 = t160 * t91 + t163 * t90;
t175 = -t161 * t96 + t164 * t84;
t130 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t160 + Ifges(4,6) * t163) * qJD(2);
t117 = Ifges(5,5) * t134 + Ifges(5,6) * t133 + Ifges(5,3) * t153;
t93 = -mrSges(5,1) * t108 + mrSges(5,3) * t105 + Ifges(5,4) * t114 + Ifges(5,2) * t113 + Ifges(5,6) * t152 - t134 * t117 + t153 * t119;
t94 = mrSges(5,2) * t108 - mrSges(5,3) * t104 + Ifges(5,1) * t114 + Ifges(5,4) * t113 + Ifges(5,5) * t152 + t133 * t117 - t153 * t118;
t75 = -mrSges(4,1) * t123 + mrSges(4,3) * t116 + Ifges(4,4) * t141 + Ifges(4,2) * t142 + Ifges(4,6) * qJDD(3) - pkin(3) * t170 + pkin(6) * t174 + qJD(3) * t132 - t130 * t180 + t159 * t94 + t162 * t93;
t81 = mrSges(4,2) * t123 - mrSges(4,3) * t115 + Ifges(4,1) * t141 + Ifges(4,4) * t142 + Ifges(4,5) * qJDD(3) - pkin(6) * t92 - qJD(3) * t131 + t130 * t179 - t159 * t93 + t162 * t94;
t171 = mrSges(3,1) * t126 - mrSges(3,2) * t127 + Ifges(3,3) * qJDD(2) + pkin(2) * t167 + pkin(5) * t176 + t160 * t81 + t163 * t75;
t168 = mrSges(2,1) * t144 - mrSges(2,2) * t145 + pkin(1) * t79 + t171;
t77 = m(2) * t145 + t175;
t76 = m(2) * t144 + t79;
t73 = -mrSges(3,1) * t156 + mrSges(3,3) * t127 + t165 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t86 - t181;
t72 = mrSges(3,2) * t156 - mrSges(3,3) * t126 + Ifges(3,5) * qJDD(2) - t165 * Ifges(3,6) - pkin(5) * t86 - t160 * t75 + t163 * t81;
t71 = mrSges(2,2) * t156 - mrSges(2,3) * t144 - pkin(4) * t79 - t161 * t73 + t164 * t72;
t70 = -mrSges(2,1) * t156 + mrSges(2,3) * t145 + t161 * t72 + t164 * t73 - pkin(1) * (m(3) * t156 + t86) + pkin(4) * t175;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t158 * t71 - t157 * t70 - qJ(1) * (t157 * t77 + t158 * t76), t71, t72, t81, t94; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t157 * t71 + t158 * t70 + qJ(1) * (-t157 * t76 + t158 * t77), t70, t73, t75, t93; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t168, t168, t171, t181, -t169;];
m_new = t1;
