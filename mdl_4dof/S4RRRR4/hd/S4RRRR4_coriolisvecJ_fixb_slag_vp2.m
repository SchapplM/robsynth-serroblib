% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
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
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRRR4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR4_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR4_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR4_coriolisvecJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR4_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:25:41
% EndTime: 2019-12-31 17:25:46
% DurationCPUTime: 2.26s
% Computational Cost: add. (2668->284), mult. (6914->405), div. (0->0), fcn. (4475->6), ass. (0->149)
t118 = qJD(2) + qJD(3);
t211 = t118 * Ifges(4,6) / 0.2e1;
t124 = cos(qJ(2));
t116 = -pkin(2) * t124 - pkin(1);
t110 = qJD(1) * t116;
t167 = t118 * Ifges(4,5);
t210 = t110 * mrSges(4,2) + t167 / 0.2e1;
t120 = sin(qJ(3));
t121 = sin(qJ(2));
t123 = cos(qJ(3));
t103 = t120 * t124 + t123 * t121;
t97 = t103 * qJD(1);
t172 = t97 * Ifges(4,4);
t102 = t120 * t121 - t123 * t124;
t96 = t102 * qJD(1);
t209 = t211 + t172 / 0.2e1 - t96 * Ifges(4,2) / 0.2e1;
t122 = cos(qJ(4));
t82 = t118 * t102;
t71 = t82 * qJD(1);
t119 = sin(qJ(4));
t85 = t118 * t122 - t119 * t97;
t40 = t85 * qJD(4) - t122 * t71;
t86 = t118 * t119 + t122 * t97;
t41 = -t86 * qJD(4) + t119 * t71;
t10 = -mrSges(5,1) * t41 + mrSges(5,2) * t40;
t190 = -pkin(6) - pkin(5);
t147 = qJD(2) * t190;
t143 = qJD(1) * t147;
t101 = t121 * t143;
t131 = t124 * t143;
t111 = t190 * t121;
t105 = qJD(1) * t111;
t100 = qJD(2) * pkin(2) + t105;
t112 = t190 * t124;
t106 = qJD(1) * t112;
t156 = t123 * t106;
t78 = t100 * t120 - t156;
t36 = t78 * qJD(3) + t101 * t120 - t123 * t131;
t206 = -m(5) * t36 - t10;
t56 = t96 * pkin(3) - t97 * pkin(7) + t110;
t65 = pkin(7) * t118 + t78;
t16 = -t119 * t65 + t122 * t56;
t153 = qJD(2) * t121;
t150 = pkin(2) * t153;
t83 = t118 * t103;
t72 = t83 * qJD(1);
t21 = pkin(3) * t72 + pkin(7) * t71 + qJD(1) * t150;
t157 = t106 * t120;
t77 = t100 * t123 + t157;
t35 = t77 * qJD(3) + t123 * t101 + t120 * t131;
t2 = t16 * qJD(4) + t119 * t21 + t122 * t35;
t17 = t119 * t56 + t122 * t65;
t3 = -t17 * qJD(4) - t119 * t35 + t122 * t21;
t205 = -t3 * t119 + t122 * t2;
t204 = -t110 * mrSges(4,1) - t16 * mrSges(5,1) + t17 * mrSges(5,2) + t209;
t140 = mrSges(5,1) * t119 + mrSges(5,2) * t122;
t64 = -pkin(3) * t118 - t77;
t130 = t64 * t140;
t135 = Ifges(5,5) * t122 - Ifges(5,6) * t119;
t168 = Ifges(5,4) * t122;
t137 = -Ifges(5,2) * t119 + t168;
t169 = Ifges(5,4) * t119;
t139 = Ifges(5,1) * t122 - t169;
t185 = t122 / 0.2e1;
t186 = -t119 / 0.2e1;
t194 = t86 / 0.2e1;
t182 = Ifges(5,4) * t86;
t93 = qJD(4) + t96;
t32 = Ifges(5,2) * t85 + Ifges(5,6) * t93 + t182;
t84 = Ifges(5,4) * t85;
t33 = Ifges(5,1) * t86 + Ifges(5,5) * t93 + t84;
t203 = (-t119 * t17 - t122 * t16) * mrSges(5,3) + t93 * t135 / 0.2e1 + t139 * t194 + t85 * t137 / 0.2e1 + t130 + t33 * t185 + t32 * t186;
t202 = -t32 / 0.2e1;
t201 = -t119 * t16 + t122 * t17;
t200 = t3 * mrSges(5,1) - t2 * mrSges(5,2) + Ifges(5,5) * t40 + Ifges(5,6) * t41;
t199 = t40 / 0.2e1;
t198 = t41 / 0.2e1;
t197 = t72 / 0.2e1;
t196 = -t85 / 0.2e1;
t195 = -t86 / 0.2e1;
t193 = -t93 / 0.2e1;
t188 = pkin(1) * mrSges(3,1);
t187 = pkin(1) * mrSges(3,2);
t184 = m(4) * t110;
t183 = mrSges(4,3) * t96;
t92 = Ifges(4,4) * t96;
t132 = t123 * t111 + t112 * t120;
t179 = t36 * t132;
t178 = t85 * Ifges(5,6);
t177 = t86 * Ifges(5,5);
t176 = t93 * Ifges(5,3);
t174 = t97 * mrSges(4,3);
t173 = t97 * Ifges(4,1);
t171 = mrSges(4,1) * t118 + mrSges(5,1) * t85 - mrSges(5,2) * t86 - t174;
t170 = Ifges(3,4) * t121;
t164 = t119 * t96;
t162 = t122 * t96;
t161 = Ifges(3,5) * qJD(2);
t160 = Ifges(3,6) * qJD(2);
t159 = qJD(2) * mrSges(3,1);
t158 = qJD(2) * mrSges(3,2);
t155 = qJD(1) * t121;
t154 = qJD(1) * t124;
t152 = qJD(4) * t119;
t151 = qJD(4) * t122;
t146 = t161 / 0.2e1;
t145 = -t160 / 0.2e1;
t144 = t124 * t147;
t75 = pkin(3) * t97 + pkin(7) * t96;
t141 = mrSges(5,1) * t122 - mrSges(5,2) * t119;
t138 = Ifges(5,1) * t119 + t168;
t136 = Ifges(5,2) * t122 + t169;
t134 = Ifges(5,5) * t119 + Ifges(5,6) * t122;
t76 = t102 * pkin(3) - t103 * pkin(7) + t116;
t88 = t111 * t120 - t112 * t123;
t42 = -t119 * t88 + t122 * t76;
t43 = t119 * t76 + t122 * t88;
t11 = mrSges(5,1) * t72 - mrSges(5,3) * t40;
t12 = -mrSges(5,2) * t72 + mrSges(5,3) * t41;
t53 = -mrSges(5,2) * t93 + mrSges(5,3) * t85;
t54 = mrSges(5,1) * t93 - mrSges(5,3) * t86;
t126 = m(5) * (-t16 * t151 - t17 * t152 + t205) + t122 * t12 - t119 * t11 - t54 * t151 - t53 * t152;
t31 = t176 + t177 + t178;
t62 = t167 - t92 + t173;
t8 = t40 * Ifges(5,4) + t41 * Ifges(5,2) + t72 * Ifges(5,6);
t9 = t40 * Ifges(5,1) + t41 * Ifges(5,4) + t72 * Ifges(5,5);
t125 = t119 * t9 / 0.2e1 - Ifges(4,6) * t72 - Ifges(4,5) * t71 - t35 * mrSges(4,2) + t164 * t202 - t77 * t183 + t8 * t185 + t134 * t197 + t136 * t198 + t138 * t199 + t33 * t162 / 0.2e1 + (-t135 * t193 - t137 * t196 - t139 * t195 + t130 + t210) * t96 + (Ifges(5,5) * t195 + Ifges(5,6) * t196 + Ifges(5,3) * t193 + t204 + t211) * t97 + (-mrSges(4,1) - t141) * t36 + (-Ifges(4,2) * t97 + t62 - t92) * t96 / 0.2e1 - (-Ifges(4,1) * t96 - t172 + t31) * t97 / 0.2e1 + (-t16 * t162 - t17 * t164 + t205) * mrSges(5,3) + t203 * qJD(4);
t117 = Ifges(3,4) * t154;
t109 = mrSges(3,3) * t154 - t158;
t108 = -mrSges(3,3) * t155 + t159;
t107 = t121 * t147;
t95 = Ifges(3,1) * t155 + t117 + t161;
t94 = t160 + (t124 * Ifges(3,2) + t170) * qJD(1);
t89 = -mrSges(4,2) * t118 - t183;
t80 = t105 * t123 + t157;
t79 = t105 * t120 - t156;
t74 = mrSges(4,1) * t96 + mrSges(4,2) * t97;
t68 = Ifges(5,3) * t72;
t60 = pkin(2) * t155 + t75;
t45 = t88 * qJD(3) + t107 * t120 - t123 * t144;
t44 = t132 * qJD(3) + t123 * t107 + t120 * t144;
t39 = pkin(3) * t83 + pkin(7) * t82 + t150;
t30 = t119 * t75 + t122 * t77;
t29 = -t119 * t77 + t122 * t75;
t27 = t119 * t60 + t122 * t80;
t26 = -t119 * t80 + t122 * t60;
t5 = -t43 * qJD(4) - t119 * t44 + t122 * t39;
t4 = t42 * qJD(4) + t119 * t39 + t122 * t44;
t1 = [t116 * (mrSges(4,1) * t72 - mrSges(4,2) * t71) - t132 * t10 + t44 * t89 + t4 * t53 + t5 * t54 + t42 * t11 + t43 * t12 + (t176 / 0.2e1 + t178 / 0.2e1 + t177 / 0.2e1 + t31 / 0.2e1 - t204 - t209) * t83 - (-t92 / 0.2e1 + t173 / 0.2e1 + t62 / 0.2e1 + t203 + t210) * t82 - t171 * t45 + m(5) * (t16 * t5 + t17 * t4 + t2 * t43 + t3 * t42 + t45 * t64 - t179) + m(4) * (t35 * t88 + t44 * t78 - t45 * t77 - t179) + (t132 * t71 - t72 * t88 + t77 * t82 - t78 * t83) * mrSges(4,3) + (t95 / 0.2e1 - pkin(5) * t108 + t146 + (-0.2e1 * t187 + 0.3e1 / 0.2e1 * Ifges(3,4) * t124) * qJD(1)) * t124 * qJD(2) + (t68 / 0.2e1 + Ifges(4,4) * t71 - t35 * mrSges(4,3) + (Ifges(5,3) / 0.2e1 + Ifges(4,2)) * t72 + t200) * t102 + (t9 * t185 + t135 * t197 + t139 * t199 + t137 * t198 - Ifges(4,1) * t71 - Ifges(4,4) * t72 + t8 * t186 + (mrSges(4,3) + t140) * t36 + (-t119 * t2 - t122 * t3) * mrSges(5,3) + (-t201 * mrSges(5,3) + t122 * t202 + t134 * t193 + t136 * t196 + t138 * t195 + t64 * t141 + t33 * t186) * qJD(4)) * t103 + (-t94 / 0.2e1 - pkin(5) * t109 + t145 + (-0.2e1 * t188 - 0.3e1 / 0.2e1 * t170 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t124) * qJD(1) + (0.2e1 * t184 + t74 + qJD(1) * (mrSges(4,1) * t102 + mrSges(4,2) * t103)) * pkin(2)) * t153; t126 * (pkin(2) * t120 + pkin(7)) + ((-t117 / 0.2e1 - t95 / 0.2e1 + t146 + qJD(1) * t187 + (t108 - t159) * pkin(5)) * t124 + (t145 + t94 / 0.2e1 + (t188 + t170 / 0.2e1 + (Ifges(3,2) / 0.2e1 - Ifges(3,1) / 0.2e1) * t124) * qJD(1) + (t109 + t158) * pkin(5) + (-t74 - t184) * pkin(2)) * t121) * qJD(1) + (m(4) * (t120 * t35 - t123 * t36) + (-t120 * t72 + t123 * t71) * mrSges(4,3) + ((-m(4) * t77 + m(5) * t64 - t171) * t120 + (m(4) * t78 + m(5) * t201 - t119 * t54 + t122 * t53 + t89) * t123) * qJD(3)) * pkin(2) - m(5) * (t16 * t26 + t17 * t27 + t64 * t79) - t80 * t89 - t27 * t53 - t26 * t54 + t125 + t78 * t174 - m(4) * (-t77 * t79 + t78 * t80) + t171 * t79 - t206 * (-pkin(2) * t123 - pkin(3)); (t171 + t174) * t78 - m(5) * (t16 * t29 + t17 * t30 + t64 * t78) + t126 * pkin(7) - t77 * t89 - t30 * t53 - t29 * t54 + t125 + t206 * pkin(3); t68 - t64 * (mrSges(5,1) * t86 + mrSges(5,2) * t85) + (Ifges(5,1) * t85 - t182) * t195 + t32 * t194 + (Ifges(5,5) * t85 - Ifges(5,6) * t86) * t193 - t16 * t53 + t17 * t54 + (t16 * t85 + t17 * t86) * mrSges(5,3) + (-Ifges(5,2) * t86 + t33 + t84) * t196 + t200;];
tauc = t1(:);
