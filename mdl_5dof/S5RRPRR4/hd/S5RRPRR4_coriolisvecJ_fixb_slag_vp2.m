% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-05 18:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRR4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR4_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR4_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:31:50
% EndTime: 2019-12-05 18:31:53
% DurationCPUTime: 1.44s
% Computational Cost: add. (2738->243), mult. (5149->356), div. (0->0), fcn. (3054->8), ass. (0->146)
t126 = sin(qJ(5));
t127 = sin(qJ(4));
t129 = cos(qJ(5));
t130 = cos(qJ(4));
t105 = -t126 * t127 + t129 * t130;
t122 = qJD(4) + qJD(5);
t61 = t122 * t105;
t206 = t61 / 0.2e1;
t123 = qJD(1) + qJD(2);
t87 = t105 * t123;
t205 = -t87 / 0.2e1;
t50 = t61 * t123;
t204 = t105 * t50;
t106 = t126 * t130 + t127 * t129;
t62 = t122 * t106;
t51 = t62 * t123;
t203 = t106 * t51;
t202 = -mrSges(5,1) * t130 + mrSges(5,2) * t127 - mrSges(4,1);
t164 = qJD(4) * t127;
t120 = t130 * qJD(3);
t125 = cos(pkin(9));
t131 = cos(qJ(2));
t124 = sin(pkin(9));
t128 = sin(qJ(2));
t168 = t124 * t128;
t182 = pkin(1) * qJD(2);
t97 = (t125 * t131 - t168) * t182;
t82 = qJD(1) * t97;
t173 = qJD(4) * t120 + t130 * t82;
t181 = qJD(1) * pkin(1);
t158 = t131 * t181;
t109 = t123 * pkin(2) + t158;
t159 = t128 * t181;
t73 = t124 * t109 + t125 * t159;
t65 = pkin(7) * t123 + t73;
t33 = -t65 * t164 + t173;
t178 = t127 * t82;
t165 = qJD(3) * t127;
t56 = t130 * t65 + t165;
t34 = -t56 * qJD(4) - t178;
t201 = -t127 * t34 + t130 * t33;
t116 = pkin(2) * t124 + pkin(7);
t187 = -pkin(8) - t116;
t103 = t187 * t127;
t121 = t130 * pkin(8);
t104 = t116 * t130 + t121;
t59 = t103 * t126 + t104 * t129;
t112 = t124 * t159;
t96 = t125 * t158 - t112;
t150 = qJD(4) * t187;
t98 = t127 * t150;
t99 = t130 * t150;
t200 = -t59 * qJD(5) + t106 * t96 - t126 * t98 + t129 * t99;
t58 = t103 * t129 - t104 * t126;
t199 = t58 * qJD(5) - t105 * t96 + t126 * t99 + t129 * t98;
t198 = mrSges(3,1) * t128 + mrSges(3,2) * t131;
t154 = pkin(8) * t123 + t65;
t47 = t154 * t130 + t165;
t88 = t106 * t123;
t196 = t88 / 0.2e1;
t119 = pkin(1) * t131 + pkin(2);
t167 = t125 * t128;
t166 = pkin(1) * t167 + t124 * t119;
t93 = pkin(7) + t166;
t195 = -pkin(8) - t93;
t194 = mrSges(6,3) * t87;
t193 = Ifges(6,4) * t88;
t192 = pkin(2) * t125;
t191 = pkin(4) * t130;
t55 = -t127 * t65 + t120;
t190 = t55 * mrSges(5,3);
t189 = t56 * mrSges(5,3);
t188 = t88 * mrSges(6,3);
t184 = Ifges(5,4) * t127;
t180 = t126 * t47;
t177 = t127 * t97;
t176 = t129 * t47;
t174 = t130 * t97;
t172 = Ifges(5,5) * qJD(4);
t171 = Ifges(5,6) * qJD(4);
t170 = t123 * t127;
t169 = t123 * t130;
t163 = qJD(5) * t126;
t162 = qJD(5) * t129;
t157 = pkin(4) * t164;
t156 = -pkin(3) - t191;
t155 = t123 * t164;
t153 = qJD(4) * t195;
t149 = t202 * t123;
t72 = t109 * t125 - t112;
t148 = -pkin(1) * t168 + t119 * t125;
t92 = -pkin(3) - t148;
t145 = t154 * t127;
t46 = t120 - t145;
t41 = qJD(4) * pkin(4) + t46;
t12 = t129 * t41 - t180;
t13 = t126 * t41 + t176;
t25 = -qJD(4) * t145 + t173;
t26 = -t47 * qJD(4) - t178;
t4 = -t13 * qJD(5) - t126 * t25 + t129 * t26;
t146 = -t106 * t4 - t12 * t61;
t68 = t195 * t127;
t69 = t130 * t93 + t121;
t35 = -t126 * t69 + t129 * t68;
t36 = t126 * t68 + t129 * t69;
t143 = -t127 * t56 - t130 * t55;
t142 = -t127 * t55 + t130 * t56;
t107 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t170;
t108 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t169;
t141 = -t127 * t107 + t130 * t108;
t140 = pkin(1) * (t124 * t131 + t167);
t139 = (Ifges(5,2) * t130 + t184) * t123;
t138 = (mrSges(5,1) * t127 + mrSges(5,2) * t130) * qJD(4);
t95 = qJD(2) * t140;
t137 = t123 * mrSges(4,2) - t141;
t81 = qJD(1) * t95;
t3 = t12 * qJD(5) + t126 * t26 + t129 * t25;
t42 = Ifges(6,2) * t87 + Ifges(6,6) * t122 + t193;
t80 = Ifges(6,4) * t87;
t43 = Ifges(6,1) * t88 + Ifges(6,5) * t122 + t80;
t57 = t156 * t123 - t72;
t134 = t4 * mrSges(6,1) - t3 * mrSges(6,2) + t12 * t194 + t42 * t196 - t57 * (mrSges(6,1) * t88 + mrSges(6,2) * t87) - t88 * (Ifges(6,1) * t87 - t193) / 0.2e1 - t122 * (Ifges(6,5) * t87 - Ifges(6,6) * t88) / 0.2e1 - Ifges(6,6) * t51 + Ifges(6,5) * t50 + (-Ifges(6,2) * t88 + t43 + t80) * t205;
t133 = m(5) * (t143 * qJD(4) + t201);
t63 = pkin(4) * t155 + t81;
t64 = -pkin(3) * t123 - t72;
t85 = t139 + t171;
t113 = Ifges(5,4) * t169;
t86 = Ifges(5,1) * t170 + t113 + t172;
t132 = (Ifges(5,1) * t130 - t184) * t155 + qJD(4) ^ 2 * (Ifges(5,5) * t130 - Ifges(5,6) * t127) / 0.2e1 + t122 * (Ifges(6,5) * t61 - Ifges(6,6) * t62) / 0.2e1 + t63 * (-mrSges(6,1) * t105 + mrSges(6,2) * t106) - t82 * mrSges(4,2) - t62 * t42 / 0.2e1 + t57 * (mrSges(6,1) * t62 + mrSges(6,2) * t61) + t43 * t206 + t64 * t138 + t202 * t81 - (t139 + t85) * t164 / 0.2e1 + (-t105 * t51 + t62 * t205) * Ifges(6,2) + (t50 * t106 + t61 * t196) * Ifges(6,1) + (t3 * t105 - t13 * t62) * mrSges(6,3) + t201 * mrSges(5,3) + (-t62 * t196 + t87 * t206 - t203 + t204) * Ifges(6,4) + (t86 + (0.3e1 * Ifges(5,4) * t130 + (Ifges(5,1) - 0.2e1 * Ifges(5,2)) * t127) * t123) * qJD(4) * t130 / 0.2e1;
t117 = -pkin(3) - t192;
t110 = t156 - t192;
t94 = qJD(1) * t140;
t89 = t123 * t138;
t78 = t92 - t191;
t77 = t95 + t157;
t67 = mrSges(6,1) * t122 - t188;
t66 = -mrSges(6,2) * t122 + t194;
t53 = -mrSges(6,1) * t87 + mrSges(6,2) * t88;
t39 = t130 * t153 - t177;
t38 = t127 * t153 + t174;
t16 = mrSges(6,1) * t51 + mrSges(6,2) * t50;
t15 = t129 * t46 - t180;
t14 = -t126 * t46 - t176;
t6 = -t36 * qJD(5) - t126 * t38 + t129 * t39;
t5 = t35 * qJD(5) + t126 * t39 + t129 * t38;
t1 = [m(4) * (-t81 * t148 + t82 * t166 - t72 * t95 + t73 * t97) + m(6) * (t12 * t6 + t13 * t5 + t3 * t36 + t35 * t4 + t57 * t77 + t63 * t78) + (-t35 * t50 - t36 * t51 + t146) * mrSges(6,3) + t132 + m(5) * (t56 * t174 - t55 * t177 + t64 * t95 + t81 * t92) + ((-t107 * t130 - t108 * t127) * t93 + t143 * mrSges(5,3)) * qJD(4) - t137 * t97 + t149 * t95 + t92 * t89 + t77 * t53 + t78 * t16 + t5 * t66 + t6 * t67 + t93 * t133 + t198 * t182 * (-qJD(1) - t123); (-t58 * t50 - t59 * t51 + t146) * mrSges(6,3) + ((-t107 * t116 - t190) * t130 + (pkin(4) * t53 - t108 * t116 - t189) * t127) * qJD(4) + t116 * t133 + t132 + t117 * t89 + t137 * t96 + (-t149 - t53) * t94 + t110 * t16 + t200 * t67 + t199 * t66 + t198 * t181 * (-qJD(2) + t123) + (t110 * t63 + t3 * t59 + t4 * t58 + (t157 - t94) * t57 + t199 * t13 + t200 * t12) * m(6) + (t117 * t81 - t142 * t96 - t64 * t94) * m(5) + ((t124 * t82 - t125 * t81) * pkin(2) + t72 * t94 - t73 * t96) * m(4); t61 * t66 - t62 * t67 + (-t203 - t204) * mrSges(6,3) + m(5) * (t127 * t33 + t130 * t34) + m(6) * (t105 * t4 + t106 * t3 - t12 * t62 + t13 * t61) + (m(5) * t142 + (-t127 ^ 2 - t130 ^ 2) * t123 * mrSges(5,3) + t141) * qJD(4); -m(6) * (t12 * t14 + t13 * t15) + t56 * t107 - t55 * t108 - t15 * t66 - t14 * t67 - t33 * mrSges(5,2) + t34 * mrSges(5,1) + (t66 * t162 + m(6) * (-t12 * t163 + t126 * t3 + t129 * t4 + t13 * t162) - t67 * t163 + (-t126 * t51 - t129 * t50) * mrSges(6,3)) * pkin(4) + ((t190 + t172 / 0.2e1 - t64 * mrSges(5,2) - t113 / 0.2e1 - t86 / 0.2e1) * t130 + (t189 - t171 / 0.2e1 - t64 * mrSges(5,1) + t85 / 0.2e1 + (t184 / 0.2e1 + (Ifges(5,2) / 0.2e1 - Ifges(5,1) / 0.2e1) * t130) * t123 + (-m(6) * t57 - t53) * pkin(4)) * t127) * t123 + t13 * t188 + t134; (t67 + t188) * t13 - t12 * t66 + t134;];
tauc = t1(:);
