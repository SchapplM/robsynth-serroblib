% Calculate vector of cutting torques with Newton-Euler for
% S4RRRR2
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
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 17:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4RRRR2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR2_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR2_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR2_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR2_invdynm_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR2_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR2_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR2_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:23:10
% EndTime: 2019-12-31 17:23:13
% DurationCPUTime: 1.62s
% Computational Cost: add. (25054->195), mult. (32162->253), div. (0->0), fcn. (18128->8), ass. (0->84)
t167 = cos(qJ(3));
t165 = sin(qJ(1));
t169 = cos(qJ(1));
t150 = t165 * g(1) - t169 * g(2);
t144 = qJDD(1) * pkin(1) + t150;
t151 = -t169 * g(1) - t165 * g(2);
t170 = qJD(1) ^ 2;
t145 = -t170 * pkin(1) + t151;
t164 = sin(qJ(2));
t168 = cos(qJ(2));
t128 = t164 * t144 + t168 * t145;
t159 = qJD(1) + qJD(2);
t155 = t159 ^ 2;
t157 = qJDD(1) + qJDD(2);
t125 = -t155 * pkin(2) + t157 * pkin(6) + t128;
t163 = sin(qJ(3));
t183 = t163 * t125;
t116 = -t167 * g(3) - t183;
t117 = -t163 * g(3) + t167 * t125;
t132 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t163 + Ifges(4,2) * t167) * t159;
t133 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t163 + Ifges(4,4) * t167) * t159;
t182 = qJD(3) * t159;
t139 = t163 * t157 + t167 * t182;
t140 = t167 * t157 - t163 * t182;
t186 = pkin(3) * t155;
t107 = qJDD(3) * pkin(3) - t139 * pkin(7) - t183 + (pkin(7) * t182 + t163 * t186 - g(3)) * t167;
t185 = t159 * t163;
t148 = qJD(3) * pkin(3) - pkin(7) * t185;
t161 = t167 ^ 2;
t108 = t140 * pkin(7) - qJD(3) * t148 - t161 * t186 + t117;
t162 = sin(qJ(4));
t166 = cos(qJ(4));
t105 = t166 * t107 - t162 * t108;
t106 = t162 * t107 + t166 * t108;
t135 = (t162 * t167 + t163 * t166) * t159;
t114 = -t135 * qJD(4) - t162 * t139 + t166 * t140;
t134 = (-t162 * t163 + t166 * t167) * t159;
t115 = t134 * qJD(4) + t166 * t139 + t162 * t140;
t158 = qJD(3) + qJD(4);
t119 = Ifges(5,4) * t135 + Ifges(5,2) * t134 + Ifges(5,6) * t158;
t120 = Ifges(5,1) * t135 + Ifges(5,4) * t134 + Ifges(5,5) * t158;
t156 = qJDD(3) + qJDD(4);
t174 = -mrSges(5,1) * t105 + mrSges(5,2) * t106 - Ifges(5,5) * t115 - Ifges(5,6) * t114 - Ifges(5,3) * t156 - t135 * t119 + t134 * t120;
t123 = -t134 * mrSges(5,1) + t135 * mrSges(5,2);
t129 = -t158 * mrSges(5,2) + t134 * mrSges(5,3);
t102 = m(5) * t105 + t156 * mrSges(5,1) - t115 * mrSges(5,3) - t135 * t123 + t158 * t129;
t130 = t158 * mrSges(5,1) - t135 * mrSges(5,3);
t103 = m(5) * t106 - t156 * mrSges(5,2) + t114 * mrSges(5,3) + t134 * t123 - t158 * t130;
t93 = t166 * t102 + t162 * t103;
t187 = mrSges(4,1) * t116 - mrSges(4,2) * t117 + Ifges(4,5) * t139 + Ifges(4,6) * t140 + Ifges(4,3) * qJDD(3) + pkin(3) * t93 + (t163 * t132 - t167 * t133) * t159 - t174;
t138 = (-mrSges(4,1) * t167 + mrSges(4,2) * t163) * t159;
t184 = t159 * t167;
t147 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t184;
t91 = m(4) * t116 + qJDD(3) * mrSges(4,1) - t139 * mrSges(4,3) + qJD(3) * t147 - t138 * t185 + t93;
t146 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t185;
t179 = -t162 * t102 + t166 * t103;
t92 = m(4) * t117 - qJDD(3) * mrSges(4,2) + t140 * mrSges(4,3) - qJD(3) * t146 + t138 * t184 + t179;
t181 = -t163 * t91 + t167 * t92;
t85 = m(3) * t128 - t155 * mrSges(3,1) - t157 * mrSges(3,2) + t181;
t127 = t168 * t144 - t164 * t145;
t177 = -t157 * pkin(2) - t127;
t124 = -t155 * pkin(6) + t177;
t109 = t148 * t185 - t140 * pkin(3) + (-pkin(7) * t161 - pkin(6)) * t155 + t177;
t175 = m(5) * t109 - t114 * mrSges(5,1) + t115 * mrSges(5,2) - t134 * t129 + t135 * t130;
t172 = -m(4) * t124 + t140 * mrSges(4,1) - t139 * mrSges(4,2) - t146 * t185 + t147 * t184 - t175;
t97 = m(3) * t127 + t157 * mrSges(3,1) - t155 * mrSges(3,2) + t172;
t80 = t164 * t85 + t168 * t97;
t87 = t163 * t92 + t167 * t91;
t180 = -t164 * t97 + t168 * t85;
t131 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t163 + Ifges(4,6) * t167) * t159;
t118 = Ifges(5,5) * t135 + Ifges(5,6) * t134 + Ifges(5,3) * t158;
t94 = -mrSges(5,1) * t109 + mrSges(5,3) * t106 + Ifges(5,4) * t115 + Ifges(5,2) * t114 + Ifges(5,6) * t156 - t135 * t118 + t158 * t120;
t95 = mrSges(5,2) * t109 - mrSges(5,3) * t105 + Ifges(5,1) * t115 + Ifges(5,4) * t114 + Ifges(5,5) * t156 + t134 * t118 - t158 * t119;
t76 = -mrSges(4,1) * t124 + mrSges(4,3) * t117 + Ifges(4,4) * t139 + Ifges(4,2) * t140 + Ifges(4,6) * qJDD(3) - pkin(3) * t175 + pkin(7) * t179 + qJD(3) * t133 - t131 * t185 + t162 * t95 + t166 * t94;
t82 = mrSges(4,2) * t124 - mrSges(4,3) * t116 + Ifges(4,1) * t139 + Ifges(4,4) * t140 + Ifges(4,5) * qJDD(3) - pkin(7) * t93 - qJD(3) * t132 + t131 * t184 - t162 * t94 + t166 * t95;
t176 = mrSges(3,1) * t127 - mrSges(3,2) * t128 + Ifges(3,3) * t157 + pkin(2) * t172 + pkin(6) * t181 + t163 * t82 + t167 * t76;
t173 = mrSges(2,1) * t150 - mrSges(2,2) * t151 + Ifges(2,3) * qJDD(1) + pkin(1) * t80 + t176;
t78 = m(2) * t151 - t170 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t180;
t77 = m(2) * t150 + qJDD(1) * mrSges(2,1) - t170 * mrSges(2,2) + t80;
t74 = mrSges(3,1) * g(3) + mrSges(3,3) * t128 + t155 * Ifges(3,5) + Ifges(3,6) * t157 - pkin(2) * t87 - t187;
t73 = -mrSges(3,2) * g(3) - mrSges(3,3) * t127 + Ifges(3,5) * t157 - t155 * Ifges(3,6) - pkin(6) * t87 - t163 * t76 + t167 * t82;
t72 = -mrSges(2,2) * g(3) - mrSges(2,3) * t150 + Ifges(2,5) * qJDD(1) - t170 * Ifges(2,6) - pkin(5) * t80 - t164 * t74 + t168 * t73;
t71 = Ifges(2,6) * qJDD(1) + t170 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t151 + t164 * t73 + t168 * t74 - pkin(1) * (-m(3) * g(3) + t87) + pkin(5) * t180;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t169 * t72 - t165 * t71 - pkin(4) * (t165 * t78 + t169 * t77), t72, t73, t82, t95; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t165 * t72 + t169 * t71 + pkin(4) * (-t165 * t77 + t169 * t78), t71, t74, t76, t94; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t173, t173, t176, t187, -t174;];
m_new = t1;
