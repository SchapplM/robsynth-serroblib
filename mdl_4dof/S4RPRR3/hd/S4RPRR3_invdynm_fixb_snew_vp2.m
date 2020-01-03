% Calculate vector of cutting torques with Newton-Euler for
% S4RPRR3
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
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 16:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4RPRR3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR3_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR3_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR3_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR3_invdynm_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR3_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR3_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR3_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:49:07
% EndTime: 2019-12-31 16:49:09
% DurationCPUTime: 1.44s
% Computational Cost: add. (16630->193), mult. (32162->251), div. (0->0), fcn. (18128->8), ass. (0->81)
t166 = sin(qJ(1));
t169 = cos(qJ(1));
t150 = t166 * g(1) - g(2) * t169;
t141 = qJDD(1) * pkin(1) + t150;
t151 = -g(1) * t169 - g(2) * t166;
t170 = qJD(1) ^ 2;
t143 = -pkin(1) * t170 + t151;
t162 = sin(pkin(7));
t163 = cos(pkin(7));
t128 = t162 * t141 + t163 * t143;
t124 = -pkin(2) * t170 + qJDD(1) * pkin(5) + t128;
t161 = -g(3) + qJDD(2);
t165 = sin(qJ(3));
t168 = cos(qJ(3));
t114 = -t165 * t124 + t168 * t161;
t115 = t168 * t124 + t165 * t161;
t135 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t165 + Ifges(4,2) * t168) * qJD(1);
t136 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t165 + Ifges(4,4) * t168) * qJD(1);
t183 = qJD(1) * qJD(3);
t182 = t168 * t183;
t144 = qJDD(1) * t165 + t182;
t145 = qJDD(1) * t168 - t165 * t183;
t107 = (-t144 + t182) * pkin(6) + (t165 * t168 * t170 + qJDD(3)) * pkin(3) + t114;
t185 = qJD(1) * t165;
t149 = qJD(3) * pkin(3) - pkin(6) * t185;
t160 = t168 ^ 2;
t108 = -pkin(3) * t160 * t170 + pkin(6) * t145 - qJD(3) * t149 + t115;
t164 = sin(qJ(4));
t167 = cos(qJ(4));
t105 = t107 * t167 - t108 * t164;
t106 = t107 * t164 + t108 * t167;
t138 = (t164 * t168 + t165 * t167) * qJD(1);
t116 = -qJD(4) * t138 - t144 * t164 + t145 * t167;
t137 = (-t164 * t165 + t167 * t168) * qJD(1);
t117 = qJD(4) * t137 + t144 * t167 + t145 * t164;
t157 = qJD(3) + qJD(4);
t120 = Ifges(5,4) * t138 + Ifges(5,2) * t137 + Ifges(5,6) * t157;
t121 = Ifges(5,1) * t138 + Ifges(5,4) * t137 + Ifges(5,5) * t157;
t156 = qJDD(3) + qJDD(4);
t174 = -mrSges(5,1) * t105 + mrSges(5,2) * t106 - Ifges(5,5) * t117 - Ifges(5,6) * t116 - Ifges(5,3) * t156 - t138 * t120 + t137 * t121;
t125 = -mrSges(5,1) * t137 + mrSges(5,2) * t138;
t129 = -mrSges(5,2) * t157 + mrSges(5,3) * t137;
t102 = m(5) * t105 + mrSges(5,1) * t156 - t117 * mrSges(5,3) - t125 * t138 + t129 * t157;
t130 = mrSges(5,1) * t157 - mrSges(5,3) * t138;
t103 = m(5) * t106 - mrSges(5,2) * t156 + t116 * mrSges(5,3) + t125 * t137 - t130 * t157;
t93 = t167 * t102 + t164 * t103;
t186 = mrSges(4,1) * t114 - mrSges(4,2) * t115 + Ifges(4,5) * t144 + Ifges(4,6) * t145 + Ifges(4,3) * qJDD(3) + pkin(3) * t93 + (t135 * t165 - t136 * t168) * qJD(1) - t174;
t142 = (-mrSges(4,1) * t168 + mrSges(4,2) * t165) * qJD(1);
t184 = qJD(1) * t168;
t148 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t184;
t91 = m(4) * t114 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t144 + qJD(3) * t148 - t142 * t185 + t93;
t147 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t185;
t179 = -t102 * t164 + t167 * t103;
t92 = m(4) * t115 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t145 - qJD(3) * t147 + t142 * t184 + t179;
t180 = -t165 * t91 + t168 * t92;
t85 = m(3) * t128 - mrSges(3,1) * t170 - qJDD(1) * mrSges(3,2) + t180;
t127 = t141 * t163 - t162 * t143;
t177 = -qJDD(1) * pkin(2) - t127;
t123 = -pkin(5) * t170 + t177;
t109 = t149 * t185 - pkin(3) * t145 + (-pkin(6) * t160 - pkin(5)) * t170 + t177;
t175 = m(5) * t109 - t116 * mrSges(5,1) + t117 * mrSges(5,2) - t137 * t129 + t130 * t138;
t172 = -m(4) * t123 + t145 * mrSges(4,1) - mrSges(4,2) * t144 - t147 * t185 + t148 * t184 - t175;
t97 = m(3) * t127 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t170 + t172;
t80 = t162 * t85 + t163 * t97;
t87 = t165 * t92 + t168 * t91;
t181 = -t162 * t97 + t163 * t85;
t134 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t165 + Ifges(4,6) * t168) * qJD(1);
t119 = Ifges(5,5) * t138 + Ifges(5,6) * t137 + Ifges(5,3) * t157;
t94 = -mrSges(5,1) * t109 + mrSges(5,3) * t106 + Ifges(5,4) * t117 + Ifges(5,2) * t116 + Ifges(5,6) * t156 - t119 * t138 + t121 * t157;
t95 = mrSges(5,2) * t109 - mrSges(5,3) * t105 + Ifges(5,1) * t117 + Ifges(5,4) * t116 + Ifges(5,5) * t156 + t119 * t137 - t120 * t157;
t76 = -mrSges(4,1) * t123 + mrSges(4,3) * t115 + Ifges(4,4) * t144 + Ifges(4,2) * t145 + Ifges(4,6) * qJDD(3) - pkin(3) * t175 + pkin(6) * t179 + qJD(3) * t136 - t134 * t185 + t164 * t95 + t167 * t94;
t82 = mrSges(4,2) * t123 - mrSges(4,3) * t114 + Ifges(4,1) * t144 + Ifges(4,4) * t145 + Ifges(4,5) * qJDD(3) - pkin(6) * t93 - qJD(3) * t135 + t134 * t184 - t164 * t94 + t167 * t95;
t176 = mrSges(3,1) * t127 - mrSges(3,2) * t128 + Ifges(3,3) * qJDD(1) + pkin(2) * t172 + pkin(5) * t180 + t165 * t82 + t168 * t76;
t173 = mrSges(2,1) * t150 - mrSges(2,2) * t151 + Ifges(2,3) * qJDD(1) + pkin(1) * t80 + t176;
t78 = m(2) * t151 - mrSges(2,1) * t170 - qJDD(1) * mrSges(2,2) + t181;
t77 = m(2) * t150 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t170 + t80;
t74 = -mrSges(3,1) * t161 + mrSges(3,3) * t128 + t170 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t87 - t186;
t73 = mrSges(3,2) * t161 - mrSges(3,3) * t127 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t170 - pkin(5) * t87 - t165 * t76 + t168 * t82;
t72 = -mrSges(2,2) * g(3) - mrSges(2,3) * t150 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t170 - qJ(2) * t80 - t162 * t74 + t163 * t73;
t71 = Ifges(2,6) * qJDD(1) + t170 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t151 + t162 * t73 + t163 * t74 - pkin(1) * (m(3) * t161 + t87) + qJ(2) * t181;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t169 * t72 - t166 * t71 - pkin(4) * (t166 * t78 + t169 * t77), t72, t73, t82, t95; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t166 * t72 + t169 * t71 + pkin(4) * (-t166 * t77 + t169 * t78), t71, t74, t76, t94; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t173, t173, t176, t186, -t174;];
m_new = t1;
