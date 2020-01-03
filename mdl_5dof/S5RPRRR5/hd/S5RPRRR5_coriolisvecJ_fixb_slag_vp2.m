% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRR5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR5_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR5_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:53:47
% EndTime: 2020-01-03 11:53:51
% DurationCPUTime: 1.37s
% Computational Cost: add. (2611->227), mult. (5137->333), div. (0->0), fcn. (3002->8), ass. (0->131)
t117 = qJD(4) + qJD(5);
t121 = sin(qJ(5));
t122 = sin(qJ(4));
t124 = cos(qJ(5));
t125 = cos(qJ(4));
t98 = -t121 * t122 + t124 * t125;
t60 = t117 * t98;
t190 = t60 / 0.2e1;
t118 = qJD(1) + qJD(3);
t89 = t98 * t118;
t189 = -t89 / 0.2e1;
t48 = t60 * t118;
t188 = t48 * t98;
t99 = t121 * t125 + t122 * t124;
t61 = t117 * t99;
t49 = t61 * t118;
t187 = t49 * t99;
t186 = -mrSges(5,1) * t125 + mrSges(5,2) * t122 - mrSges(4,1);
t154 = qJD(4) * t122;
t115 = t125 * qJD(2);
t112 = cos(pkin(9)) * pkin(1) + pkin(2);
t104 = t112 * qJD(1);
t123 = sin(qJ(3));
t126 = cos(qJ(3));
t176 = pkin(1) * sin(pkin(9));
t151 = qJD(1) * t176;
t80 = t104 * t126 - t123 * t151;
t75 = t80 * qJD(3);
t162 = qJD(4) * t115 + t125 * t75;
t81 = t123 * t104 + t126 * t151;
t64 = t118 * pkin(7) + t81;
t29 = -t154 * t64 + t162;
t167 = t122 * t75;
t155 = qJD(2) * t122;
t56 = t125 * t64 + t155;
t30 = -qJD(4) * t56 - t167;
t185 = -t122 * t30 + t125 * t29;
t180 = -pkin(8) - pkin(7);
t149 = qJD(4) * t180;
t100 = t122 * t149;
t101 = t125 * t149;
t106 = t180 * t122;
t116 = t125 * pkin(8);
t107 = pkin(7) * t125 + t116;
t73 = t106 * t124 - t107 * t121;
t184 = qJD(5) * t73 + t100 * t124 + t101 * t121 - t98 * t80;
t74 = t106 * t121 + t107 * t124;
t183 = -qJD(5) * t74 - t100 * t121 + t101 * t124 + t99 * t80;
t142 = t112 * t126 - t123 * t176;
t156 = t123 * t112 + t126 * t176;
t147 = pkin(8) * t118 + t64;
t51 = t125 * t147 + t155;
t90 = t99 * t118;
t181 = t90 / 0.2e1;
t95 = pkin(7) + t156;
t179 = -pkin(8) - t95;
t178 = mrSges(6,3) * t89;
t177 = Ifges(6,4) * t90;
t175 = t125 * pkin(4);
t55 = -t122 * t64 + t115;
t174 = t55 * mrSges(5,3);
t173 = t56 * mrSges(5,3);
t172 = t90 * mrSges(6,3);
t171 = Ifges(5,4) * t122;
t169 = t121 * t51;
t85 = t142 * qJD(3);
t166 = t122 * t85;
t165 = t124 * t51;
t163 = t125 * t85;
t161 = Ifges(5,5) * qJD(4);
t160 = Ifges(5,6) * qJD(4);
t158 = t118 * t122;
t157 = t118 * t125;
t153 = qJD(5) * t121;
t152 = qJD(5) * t124;
t150 = pkin(4) * t154;
t114 = -pkin(3) - t175;
t148 = t118 * t154;
t146 = qJD(4) * t179;
t143 = t186 * t118;
t94 = -pkin(3) - t142;
t139 = t147 * t122;
t50 = t115 - t139;
t45 = qJD(4) * pkin(4) + t50;
t12 = t124 * t45 - t169;
t13 = t121 * t45 + t165;
t21 = -qJD(4) * t139 + t162;
t22 = -qJD(4) * t51 - t167;
t4 = -qJD(5) * t13 - t121 * t21 + t124 * t22;
t140 = -t12 * t60 - t4 * t99;
t68 = t179 * t122;
t69 = t125 * t95 + t116;
t34 = -t121 * t69 + t124 * t68;
t35 = t121 * t68 + t124 * t69;
t137 = -t122 * t56 - t125 * t55;
t136 = -t122 * t55 + t125 * t56;
t102 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t158;
t103 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t157;
t135 = -t122 * t102 + t125 * t103;
t134 = (Ifges(5,2) * t125 + t171) * t118;
t133 = (mrSges(5,1) * t122 + mrSges(5,2) * t125) * qJD(4);
t86 = t156 * qJD(3);
t132 = t118 * mrSges(4,2) - t135;
t76 = t81 * qJD(3);
t3 = qJD(5) * t12 + t121 * t22 + t124 * t21;
t43 = Ifges(6,2) * t89 + Ifges(6,6) * t117 + t177;
t82 = Ifges(6,4) * t89;
t44 = Ifges(6,1) * t90 + Ifges(6,5) * t117 + t82;
t58 = t114 * t118 - t80;
t129 = t4 * mrSges(6,1) - t3 * mrSges(6,2) + t12 * t178 + t43 * t181 - t58 * (mrSges(6,1) * t90 + mrSges(6,2) * t89) - t90 * (Ifges(6,1) * t89 - t177) / 0.2e1 - t117 * (Ifges(6,5) * t89 - Ifges(6,6) * t90) / 0.2e1 - Ifges(6,6) * t49 + Ifges(6,5) * t48 + (-Ifges(6,2) * t90 + t44 + t82) * t189;
t128 = m(5) * (t137 * qJD(4) + t185);
t59 = pkin(4) * t148 + t76;
t63 = -t118 * pkin(3) - t80;
t87 = t134 + t160;
t109 = Ifges(5,4) * t157;
t88 = Ifges(5,1) * t158 + t109 + t161;
t127 = t117 * (Ifges(6,5) * t60 - Ifges(6,6) * t61) / 0.2e1 + (Ifges(5,1) * t125 - t171) * t148 + t59 * (-mrSges(6,1) * t98 + mrSges(6,2) * t99) + t58 * (mrSges(6,1) * t61 + mrSges(6,2) * t60) - t75 * mrSges(4,2) + t44 * t190 - t61 * t43 / 0.2e1 + qJD(4) ^ 2 * (Ifges(5,5) * t125 - Ifges(5,6) * t122) / 0.2e1 + t63 * t133 + t186 * t76 - (t134 + t87) * t154 / 0.2e1 + (t61 * t189 - t49 * t98) * Ifges(6,2) + (t181 * t60 + t48 * t99) * Ifges(6,1) + (-t13 * t61 + t3 * t98) * mrSges(6,3) + t185 * mrSges(5,3) + (-t61 * t181 + t89 * t190 - t187 + t188) * Ifges(6,4) + (t88 + (0.3e1 * Ifges(5,4) * t125 + (Ifges(5,1) - 0.2e1 * Ifges(5,2)) * t122) * t118) * qJD(4) * t125 / 0.2e1;
t91 = t118 * t133;
t79 = t94 - t175;
t70 = t86 + t150;
t66 = mrSges(6,1) * t117 - t172;
t65 = -mrSges(6,2) * t117 + t178;
t53 = -mrSges(6,1) * t89 + mrSges(6,2) * t90;
t41 = t125 * t146 - t166;
t40 = t122 * t146 + t163;
t16 = mrSges(6,1) * t49 + mrSges(6,2) * t48;
t15 = t124 * t50 - t169;
t14 = -t121 * t50 - t165;
t6 = -qJD(5) * t35 - t121 * t40 + t124 * t41;
t5 = qJD(5) * t34 + t121 * t41 + t124 * t40;
t1 = [t143 * t86 - t132 * t85 + m(4) * (-t142 * t76 + t156 * t75 - t80 * t86 + t81 * t85) + m(6) * (t12 * t6 + t13 * t5 + t3 * t35 + t34 * t4 + t58 * t70 + t59 * t79) + (-t34 * t48 - t35 * t49 + t140) * mrSges(6,3) + t94 * t91 + t79 * t16 + t5 * t65 + t6 * t66 + t70 * t53 + t95 * t128 + m(5) * (t56 * t163 - t55 * t166 + t63 * t86 + t76 * t94) + t127 + ((-t102 * t125 - t103 * t122) * t95 + t137 * mrSges(5,3)) * qJD(4); t60 * t65 - t61 * t66 + (-t187 - t188) * mrSges(6,3) + m(5) * (t122 * t29 + t125 * t30) + m(6) * (-t12 * t61 + t13 * t60 + t3 * t99 + t4 * t98) + (m(5) * t136 + (-t122 ^ 2 - t125 ^ 2) * t118 * mrSges(5,3) + t135) * qJD(4); ((-pkin(7) * t102 - t174) * t125 + (pkin(4) * t53 - pkin(7) * t103 - t173) * t122) * qJD(4) + (-t143 - t53) * t81 + t132 * t80 + t183 * t66 + t184 * t65 + t114 * t16 - pkin(3) * t91 + pkin(7) * t128 + t127 + (-t73 * t48 - t74 * t49 + t140) * mrSges(6,3) + (t114 * t59 + t3 * t74 + t4 * t73 + (t150 - t81) * t58 + t184 * t13 + t183 * t12) * m(6) + (-pkin(3) * t76 - t136 * t80 - t63 * t81) * m(5); -m(6) * (t12 * t14 + t13 * t15) + t13 * t172 + t56 * t102 - t55 * t103 - t15 * t65 - t14 * t66 - t29 * mrSges(5,2) + t30 * mrSges(5,1) + t129 + (m(6) * (-t12 * t153 + t121 * t3 + t124 * t4 + t13 * t152) + t65 * t152 - t66 * t153 + (-t121 * t49 - t124 * t48) * mrSges(6,3)) * pkin(4) + ((-t63 * mrSges(5,2) + t161 / 0.2e1 - t88 / 0.2e1 - t109 / 0.2e1 + t174) * t125 + (-t63 * mrSges(5,1) - t160 / 0.2e1 + t87 / 0.2e1 + t173 + (t171 / 0.2e1 + (-Ifges(5,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t125) * t118 + (-m(6) * t58 - t53) * pkin(4)) * t122) * t118; t129 - t12 * t65 + (t66 + t172) * t13;];
tauc = t1(:);
