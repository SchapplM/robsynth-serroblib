% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRRR2
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
% Datum: 2019-12-05 18:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRR2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR2_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR2_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR2_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR2_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:11:17
% EndTime: 2019-12-05 18:11:31
% DurationCPUTime: 5.18s
% Computational Cost: add. (6590->342), mult. (18163->494), div. (0->0), fcn. (14133->8), ass. (0->155)
t153 = sin(pkin(9));
t157 = sin(qJ(3));
t154 = cos(pkin(9));
t160 = cos(qJ(3));
t180 = t154 * t160;
t139 = -t153 * t157 + t180;
t131 = t139 * qJD(1);
t140 = t153 * t160 + t154 * t157;
t132 = t140 * qJD(1);
t156 = sin(qJ(4));
t159 = cos(qJ(4));
t109 = t131 * t156 + t132 * t159;
t155 = sin(qJ(5));
t158 = cos(qJ(5));
t166 = t159 * t131 - t132 * t156;
t63 = t109 * t158 + t155 * t166;
t198 = t63 / 0.2e1;
t174 = qJD(4) * t159;
t175 = qJD(4) * t156;
t134 = t140 * qJD(3);
t125 = qJD(1) * t134;
t187 = pkin(6) + qJ(2);
t144 = t187 * t153;
t141 = qJD(1) * t144;
t145 = t187 * t154;
t142 = qJD(1) * t145;
t147 = qJD(2) * t180;
t171 = qJD(1) * qJD(2);
t176 = qJD(3) * t160;
t84 = -t141 * t176 + qJD(1) * t147 + (-qJD(3) * t142 - t153 * t171) * t157;
t72 = -pkin(7) * t125 + t84;
t133 = t139 * qJD(3);
t124 = qJD(1) * t133;
t115 = -t141 * t157 + t142 * t160;
t163 = t140 * qJD(2);
t85 = -qJD(1) * t163 - qJD(3) * t115;
t73 = -pkin(7) * t124 + t85;
t114 = -t160 * t141 - t142 * t157;
t90 = -pkin(7) * t132 + t114;
t89 = qJD(3) * pkin(3) + t90;
t91 = pkin(7) * t131 + t115;
t20 = t156 * t73 + t159 * t72 + t89 * t174 - t175 * t91;
t57 = -qJD(4) * t109 - t124 * t156 - t125 * t159;
t6 = pkin(8) * t57 + t20;
t88 = t159 * t91;
t41 = t156 * t89 + t88;
t21 = -qJD(4) * t41 - t156 * t72 + t159 * t73;
t56 = t166 * qJD(4) + t124 * t159 - t125 * t156;
t7 = -pkin(8) * t56 + t21;
t214 = pkin(8) * t166;
t35 = t41 + t214;
t184 = t155 * t35;
t152 = qJD(3) + qJD(4);
t222 = pkin(8) * t109;
t86 = t156 * t91;
t40 = t159 * t89 - t86;
t34 = t40 - t222;
t33 = pkin(4) * t152 + t34;
t8 = t158 * t33 - t184;
t2 = qJD(5) * t8 + t155 * t7 + t158 * t6;
t218 = -t109 * t155 + t158 * t166;
t22 = qJD(5) * t218 + t155 * t57 + t158 * t56;
t23 = -qJD(5) * t63 - t155 * t56 + t158 * t57;
t183 = t158 * t35;
t9 = t155 * t33 + t183;
t3 = -qJD(5) * t9 - t155 * t6 + t158 * t7;
t149 = qJD(5) + t152;
t188 = Ifges(6,4) * t63;
t30 = Ifges(6,2) * t218 + Ifges(6,6) * t149 + t188;
t55 = Ifges(6,4) * t218;
t31 = Ifges(6,1) * t63 + Ifges(6,5) * t149 + t55;
t235 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t22 + Ifges(6,6) * t23 + t30 * t198 - (-Ifges(6,2) * t63 + t31 + t55) * t218 / 0.2e1;
t103 = Ifges(5,4) * t166;
t170 = -pkin(2) * t154 - pkin(1);
t143 = qJD(1) * t170 + qJD(2);
t116 = -pkin(3) * t131 + t143;
t185 = Ifges(5,4) * t109;
t51 = t109 * Ifges(5,1) + t152 * Ifges(5,5) + t103;
t233 = t21 * mrSges(5,1) - t20 * mrSges(5,2) + Ifges(5,5) * t56 + Ifges(5,6) * t57 - (Ifges(5,5) * t166 - Ifges(5,6) * t109) * t152 / 0.2e1 + (t109 * t41 + t166 * t40) * mrSges(5,3) - (-Ifges(5,2) * t109 + t103 + t51) * t166 / 0.2e1 - t116 * (mrSges(5,1) * t109 + mrSges(5,2) * t166) - (Ifges(5,1) * t166 - t185) * t109 / 0.2e1 + t235;
t74 = -pkin(4) * t166 + t116;
t232 = -(Ifges(6,5) * t218 - Ifges(6,6) * t63) * t149 / 0.2e1 - t74 * (mrSges(6,1) * t63 + mrSges(6,2) * t218) - (Ifges(6,1) * t218 - t188) * t63 / 0.2e1;
t230 = t218 * t8 + t63 * t9;
t224 = t230 * mrSges(6,3) + t232;
t148 = pkin(3) * t159 + pkin(4);
t172 = qJD(5) * t158;
t173 = qJD(5) * t155;
t178 = t156 * t158;
t43 = -t156 * t90 - t88;
t36 = t43 - t214;
t44 = t159 * t90 - t86;
t37 = t44 - t222;
t221 = t155 * t37 - t158 * t36 - t148 * t173 + (-t156 * t172 + (-t155 * t159 - t178) * qJD(4)) * pkin(3);
t179 = t155 * t156;
t220 = -t155 * t36 - t158 * t37 + t148 * t172 + (-t156 * t173 + (t158 * t159 - t179) * qJD(4)) * pkin(3);
t50 = Ifges(5,2) * t166 + t152 * Ifges(5,6) + t185;
t215 = t50 / 0.2e1;
t118 = -t157 * t144 + t160 * t145;
t204 = (m(3) * qJ(2) + mrSges(3,3)) * (t153 ^ 2 + t154 ^ 2);
t200 = t218 / 0.2e1;
t196 = t166 / 0.2e1;
t194 = t109 / 0.2e1;
t192 = t133 / 0.2e1;
t191 = -t134 / 0.2e1;
t117 = -t160 * t144 - t145 * t157;
t101 = -pkin(7) * t140 + t117;
t102 = pkin(7) * t139 + t118;
t49 = t156 * t101 + t159 * t102;
t186 = Ifges(4,4) * t132;
t169 = -t57 * mrSges(5,1) + t56 * mrSges(5,2);
t168 = -t23 * mrSges(6,1) + t22 * mrSges(6,2);
t48 = t159 * t101 - t102 * t156;
t113 = t139 * t156 + t140 * t159;
t38 = -pkin(8) * t113 + t48;
t112 = t139 * t159 - t140 * t156;
t39 = pkin(8) * t112 + t49;
t24 = -t155 * t39 + t158 * t38;
t25 = t155 * t38 + t158 * t39;
t70 = t112 * t158 - t113 * t155;
t71 = t112 * t155 + t113 * t158;
t121 = -pkin(3) * t139 + t170;
t95 = -t144 * t176 + t147 + (-qJD(2) * t153 - qJD(3) * t145) * t157;
t78 = -pkin(7) * t134 + t95;
t96 = -t118 * qJD(3) - t163;
t79 = -pkin(7) * t133 + t96;
t28 = t101 * t174 - t102 * t175 + t156 * t79 + t159 * t78;
t29 = -qJD(4) * t49 - t156 * t78 + t159 * t79;
t129 = pkin(3) * t178 + t148 * t155;
t128 = -pkin(3) * t179 + t148 * t158;
t126 = Ifges(4,4) * t131;
t122 = t124 * mrSges(4,2);
t120 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t132;
t119 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t131;
t105 = t132 * Ifges(4,1) + Ifges(4,5) * qJD(3) + t126;
t104 = t131 * Ifges(4,2) + Ifges(4,6) * qJD(3) + t186;
t93 = mrSges(5,1) * t152 - mrSges(5,3) * t109;
t92 = -mrSges(5,2) * t152 + mrSges(5,3) * t166;
t83 = -pkin(4) * t112 + t121;
t80 = pkin(3) * t132 + pkin(4) * t109;
t69 = -qJD(4) * t113 - t133 * t156 - t134 * t159;
t68 = qJD(4) * t112 + t133 * t159 - t134 * t156;
t64 = -mrSges(5,1) * t166 + mrSges(5,2) * t109;
t47 = mrSges(6,1) * t149 - mrSges(6,3) * t63;
t46 = -mrSges(6,2) * t149 + mrSges(6,3) * t218;
t45 = pkin(3) * t134 - pkin(4) * t69;
t42 = pkin(3) * t125 - pkin(4) * t57;
t32 = -mrSges(6,1) * t218 + mrSges(6,2) * t63;
t27 = -qJD(5) * t71 - t155 * t68 + t158 * t69;
t26 = qJD(5) * t70 + t155 * t69 + t158 * t68;
t15 = -pkin(8) * t68 + t29;
t14 = pkin(8) * t69 + t28;
t11 = t158 * t34 - t184;
t10 = -t155 * t34 - t183;
t5 = -qJD(5) * t25 - t14 * t155 + t15 * t158;
t4 = qJD(5) * t24 + t14 * t158 + t15 * t155;
t1 = [t69 * t215 + t170 * (t125 * mrSges(4,1) + t122) + (t125 * (-mrSges(5,1) * t112 + mrSges(5,2) * t113) + t134 * t64) * pkin(3) + (-t114 * t133 - t115 * t134 - t117 * t124 - t118 * t125 + t139 * t84 - t140 * t85) * mrSges(4,3) + m(5) * (t20 * t49 + t21 * t48 + t28 * t41 + t29 * t40 + (t116 * t134 + t121 * t125) * pkin(3)) + (-t139 * t125 + t131 * t191) * Ifges(4,2) + (t139 * t124 - t125 * t140 + t131 * t192 + t132 * t191) * Ifges(4,4) + 0.2e1 * t204 * t171 + (t2 * t70 - t22 * t24 + t23 * t25 - t26 * t8 + t27 * t9 - t3 * t71) * mrSges(6,3) + t45 * t32 + t4 * t46 + t5 * t47 + t27 * t30 / 0.2e1 + t26 * t31 / 0.2e1 + m(6) * (t2 * t25 + t24 * t3 + t4 * t9 + t42 * t83 + t45 * t74 + t5 * t8) + m(4) * (t114 * t96 + t115 * t95 + t117 * t85 + t118 * t84) + t149 * (Ifges(6,5) * t26 + Ifges(6,6) * t27) / 0.2e1 + t152 * (Ifges(5,5) * t68 + Ifges(5,6) * t69) / 0.2e1 + t116 * (-mrSges(5,1) * t69 + mrSges(5,2) * t68) + t95 * t119 + t96 * t120 + t28 * t92 + t29 * t93 + t42 * (-mrSges(6,1) * t70 + mrSges(6,2) * t71) + t74 * (-mrSges(6,1) * t27 + mrSges(6,2) * t26) + t68 * t51 / 0.2e1 + t83 * t168 + t121 * t169 + (t112 * t20 - t113 * t21 - t40 * t68 + t41 * t69 - t48 * t56 + t49 * t57) * mrSges(5,3) + qJD(3) * (Ifges(4,5) * t133 - Ifges(4,6) * t134) / 0.2e1 + t143 * (mrSges(4,1) * t134 + mrSges(4,2) * t133) + (t124 * t140 + t132 * t192) * Ifges(4,1) + (t113 * t56 + t194 * t68) * Ifges(5,1) + (t57 * t112 + t196 * t69) * Ifges(5,2) + (t56 * t112 + t113 * t57 + t194 * t69 + t196 * t68) * Ifges(5,4) + (t198 * t26 + t71 * t22) * Ifges(6,1) + (t200 * t27 + t70 * t23) * Ifges(6,2) + (t198 * t27 + t200 * t26 + t70 * t22 + t71 * t23) * Ifges(6,4) + t104 * t191 + t105 * t192; t109 * t93 - t166 * t92 - t131 * t119 + t132 * t120 - t218 * t46 + t63 * t47 + t122 - (-m(5) * pkin(3) - mrSges(4,1)) * t125 - m(4) * (-t114 * t132 + t115 * t131) - m(5) * (-t109 * t40 + t166 * t41) + t168 + t169 - t204 * qJD(1) ^ 2 + (-t218 * t9 + t63 * t8 + t42) * m(6); t232 - (-Ifges(4,2) * t132 + t105 + t126) * t131 / 0.2e1 + t109 * t215 + t233 - m(5) * (t40 * t43 + t41 * t44) + (t114 * t131 + t115 * t132) * mrSges(4,3) + t220 * t46 + (t128 * t3 + t129 * t2 + t220 * t9 + t221 * t8 - t74 * t80) * m(6) + t221 * t47 - t143 * (mrSges(4,1) * t132 + mrSges(4,2) * t131) - qJD(3) * (Ifges(4,5) * t131 - Ifges(4,6) * t132) / 0.2e1 + t132 * t104 / 0.2e1 + Ifges(4,5) * t124 - Ifges(4,6) * t125 - t114 * t119 + t115 * t120 - t44 * t92 - t43 * t93 - t80 * t32 - t84 * mrSges(4,2) + t85 * mrSges(4,1) - t132 * (Ifges(4,1) * t131 - t186) / 0.2e1 + (-t128 * t22 + t129 * t23 + t230) * mrSges(6,3) + (-t132 * t64 + (-t156 * t93 + t159 * t92) * qJD(4) + (t156 * t57 - t159 * t56) * mrSges(5,3) + (-t116 * t132 + t156 * t20 + t159 * t21 + t174 * t41 - t175 * t40) * m(5)) * pkin(3); -t11 * t46 - t10 * t47 - t40 * t92 + t41 * t93 - m(6) * (t10 * t8 + t11 * t9) + t50 * t194 + (-t109 * t32 + (-t155 * t47 + t158 * t46) * qJD(5) + (t155 * t23 - t158 * t22) * mrSges(6,3) + (-t109 * t74 + t155 * t2 + t158 * t3 + t172 * t9 - t173 * t8) * m(6)) * pkin(4) + t224 + t233; -t8 * t46 + t9 * t47 + t224 + t235;];
tauc = t1(:);
