% Calculate time derivative of joint inertia matrix for
% S6PRRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRPP2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP2_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP2_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP2_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP2_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPP2_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPP2_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:48:47
% EndTime: 2019-03-08 22:48:53
% DurationCPUTime: 2.83s
% Computational Cost: add. (1670->432), mult. (4577->584), div. (0->0), fcn. (3470->8), ass. (0->183)
t224 = (m(6) + m(7));
t198 = -mrSges(6,2) + mrSges(7,3);
t222 = Ifges(6,4) + Ifges(5,5);
t223 = -Ifges(6,2) - Ifges(5,3);
t117 = sin(qJ(4));
t118 = sin(qJ(3));
t120 = cos(qJ(4));
t165 = qJD(4) * t120;
t121 = cos(qJ(3));
t168 = qJD(3) * t121;
t126 = t117 * t168 + t118 * t165;
t166 = qJD(4) * t118;
t151 = t117 * t166;
t152 = t120 * t168;
t127 = -t151 + t152;
t183 = Ifges(6,5) * t120;
t138 = Ifges(6,3) * t117 + t183;
t185 = Ifges(7,4) * t120;
t139 = Ifges(7,2) * t117 + t185;
t220 = (t138 + t139) * qJD(4);
t169 = qJD(3) * t118;
t219 = qJ(5) * t169 - qJD(5) * t121;
t186 = Ifges(7,4) * t117;
t141 = Ifges(7,1) * t120 + t186;
t184 = Ifges(6,5) * t117;
t142 = Ifges(6,1) * t120 + t184;
t188 = Ifges(5,4) * t117;
t143 = Ifges(5,1) * t120 - t188;
t218 = (t141 + t142 + t143) * qJD(4);
t177 = qJ(5) * t117;
t209 = pkin(4) + pkin(5);
t217 = -t209 * t120 - t177;
t215 = 2 * m(5);
t214 = 2 * m(7);
t213 = 2 * pkin(8);
t212 = -2 * mrSges(7,3);
t211 = m(5) / 0.2e1;
t210 = m(6) / 0.2e1;
t116 = sin(pkin(6));
t122 = cos(qJ(2));
t174 = t116 * t122;
t157 = t117 * t174;
t119 = sin(qJ(2));
t175 = t116 * t119;
t178 = cos(pkin(6));
t56 = t178 * t118 + t121 * t175;
t32 = t120 * t56 - t157;
t155 = qJD(2) * t175;
t170 = qJD(2) * t122;
t154 = t116 * t170;
t55 = t118 * t175 - t178 * t121;
t30 = -t55 * qJD(3) + t121 * t154;
t31 = t117 * t56 + t120 * t174;
t8 = -t31 * qJD(4) + t117 * t155 + t120 * t30;
t205 = t8 * qJ(5) + t32 * qJD(5);
t204 = pkin(8) * t117;
t29 = t56 * qJD(3) + t118 * t154;
t203 = t55 * t29;
t7 = -qJD(4) * t157 + t117 * t30 - t120 * t155 + t56 * t165;
t202 = t7 * t117;
t201 = t8 * t120;
t199 = -mrSges(6,1) - mrSges(7,1);
t197 = mrSges(7,2) + mrSges(6,3);
t196 = Ifges(5,6) + Ifges(7,6);
t195 = pkin(9) - qJ(6);
t37 = mrSges(5,1) * t169 - t127 * mrSges(5,3);
t167 = qJD(4) * t117;
t98 = mrSges(6,2) * t152;
t38 = t98 + (-mrSges(6,1) * qJD(3) - mrSges(6,2) * t167) * t118;
t194 = -t37 + t38;
t40 = -mrSges(5,2) * t169 - t126 * mrSges(5,3);
t41 = -t126 * mrSges(6,2) + mrSges(6,3) * t169;
t193 = t40 + t41;
t81 = (pkin(3) * t118 - pkin(9) * t121) * qJD(3);
t85 = -pkin(3) * t121 - pkin(9) * t118 - pkin(2);
t192 = t117 * t81 + t85 * t165;
t173 = t117 * t118;
t76 = mrSges(5,2) * t121 - mrSges(5,3) * t173;
t80 = -mrSges(6,2) * t173 - mrSges(6,3) * t121;
t191 = -t76 - t80;
t172 = t118 * t120;
t78 = -mrSges(5,1) * t121 - mrSges(5,3) * t172;
t79 = mrSges(6,1) * t121 + mrSges(6,2) * t172;
t190 = -t78 + t79;
t89 = -t120 * mrSges(5,1) + t117 * mrSges(5,2);
t189 = t89 - mrSges(4,1);
t187 = Ifges(5,4) * t120;
t182 = Ifges(7,5) * t121;
t181 = t118 * t29;
t180 = t121 * t30;
t171 = t120 * t121;
t106 = pkin(8) * t171;
t179 = qJD(4) * t106 + t85 * t167;
t44 = t117 * t85 + t106;
t176 = qJ(5) * t120;
t164 = qJD(5) * t117;
t163 = qJD(5) * t120;
t161 = qJD(6) * t120;
t160 = -Ifges(7,5) + t222;
t48 = t141 * t118 + t182;
t49 = -t121 * Ifges(6,4) + t142 * t118;
t50 = -t121 * Ifges(5,5) + t143 * t118;
t159 = -t48 - t49 - t50;
t158 = mrSges(7,3) * t171;
t156 = -pkin(4) - t204;
t90 = t195 * t120;
t105 = t121 * t204;
t43 = t120 * t85 - t105;
t91 = -Ifges(6,3) * t120 + t184;
t92 = -Ifges(7,2) * t120 + t186;
t93 = Ifges(5,2) * t120 + t188;
t149 = t91 / 0.2e1 + t92 / 0.2e1 - t93 / 0.2e1;
t94 = Ifges(7,1) * t117 - t185;
t95 = Ifges(6,1) * t117 - t183;
t96 = Ifges(5,1) * t117 + t187;
t148 = t94 / 0.2e1 + t95 / 0.2e1 + t96 / 0.2e1;
t147 = -t120 * t81 + t179;
t34 = -qJ(5) * t121 + t44;
t66 = -mrSges(7,1) * t167 + mrSges(7,2) * t165;
t145 = mrSges(5,1) * t117 + mrSges(5,2) * t120;
t87 = -t120 * mrSges(6,1) - t117 * mrSges(6,3);
t144 = mrSges(6,1) * t117 - mrSges(6,3) * t120;
t140 = -Ifges(5,2) * t117 + t187;
t137 = pkin(4) * t120 + t177;
t136 = pkin(4) * t117 - t176;
t134 = (m(6) * pkin(9) - t198) * t120;
t133 = pkin(8) + t136;
t132 = -t126 * Ifges(6,6) - t222 * t152 + t223 * t169;
t131 = t55 * t168 + t181;
t130 = -t209 * t117 + t176;
t45 = -Ifges(6,6) * t121 + t138 * t118;
t46 = t121 * Ifges(7,6) + t139 * t118;
t47 = -t121 * Ifges(5,6) + t140 * t118;
t129 = t196 * t121 + t45 + t46 - t47;
t128 = -pkin(8) + t130;
t13 = (-t120 * t169 - t121 * t167) * pkin(8) + t192;
t25 = -t126 * mrSges(7,1) + t127 * mrSges(7,2);
t115 = t121 * pkin(4);
t114 = Ifges(6,4) * t165;
t113 = Ifges(5,5) * t165;
t111 = Ifges(6,6) * t167;
t99 = mrSges(7,3) * t151;
t88 = mrSges(7,1) * t120 + mrSges(7,2) * t117;
t86 = t195 * t117;
t82 = -pkin(3) - t137;
t77 = mrSges(7,1) * t121 - mrSges(7,3) * t172;
t75 = -mrSges(7,2) * t121 + mrSges(7,3) * t173;
t71 = t140 * qJD(4);
t68 = (mrSges(4,1) * t118 + mrSges(4,2) * t121) * qJD(3);
t67 = t145 * qJD(4);
t65 = t144 * qJD(4);
t63 = pkin(3) - t217;
t60 = t145 * t118;
t59 = (-mrSges(7,1) * t117 + mrSges(7,2) * t120) * t118;
t58 = t144 * t118;
t54 = qJD(4) * t90 - qJD(6) * t117;
t53 = t136 * qJD(4) - t164;
t52 = -t195 * t167 - t161;
t51 = t133 * t118;
t42 = t130 * qJD(4) + t164;
t39 = mrSges(7,2) * t169 + t126 * mrSges(7,3);
t36 = t99 + (-mrSges(7,1) * t118 - t158) * qJD(3);
t35 = t115 - t43;
t33 = t128 * t118;
t27 = qJ(6) * t173 + t34;
t26 = t126 * mrSges(5,1) + t127 * mrSges(5,2);
t24 = t126 * mrSges(6,1) - t127 * mrSges(6,3);
t23 = pkin(5) * t121 + t105 + t115 + (-qJ(6) * t118 - t85) * t120;
t22 = -t96 * t166 + (Ifges(5,5) * t118 + t143 * t121) * qJD(3);
t21 = -t95 * t166 + (Ifges(6,4) * t118 + t142 * t121) * qJD(3);
t20 = -t94 * t166 + (-Ifges(7,5) * t118 + t141 * t121) * qJD(3);
t19 = -t93 * t166 + (Ifges(5,6) * t118 + t140 * t121) * qJD(3);
t18 = -t92 * t166 + (-Ifges(7,6) * t118 + t139 * t121) * qJD(3);
t17 = -t91 * t166 + (Ifges(6,6) * t118 + t138 * t121) * qJD(3);
t14 = t169 * t204 - t147;
t12 = (t137 * qJD(4) - t163) * t118 + t133 * t168;
t11 = t156 * t169 + t147;
t10 = t13 + t219;
t9 = (t217 * qJD(4) + t163) * t118 + t128 * t168;
t5 = pkin(9) * t201;
t4 = (-pkin(8) * qJD(3) + qJ(6) * qJD(4)) * t172 + (qJD(6) * t118 + (-pkin(8) * qJD(4) + qJ(6) * qJD(3)) * t121) * t117 + t192 + t219;
t3 = (-qJ(6) * t168 - t81) * t120 + (qJ(6) * t167 - t161 + (-pkin(5) + t156) * qJD(3)) * t118 + t179;
t1 = [0.2e1 * m(4) * (-t116 ^ 2 * t119 * t170 + t30 * t56 + t203) + 0.2e1 * (m(5) + t224) * (t31 * t7 + t32 * t8 + t203); (t75 - t191) * t8 + (t77 + t190) * t7 + (t24 + t26 - t25) * t55 + (t39 + t193) * t32 + (t36 + t194) * t31 + (t58 - t59 + t60) * t29 + (-t122 * t68 + (-t122 * mrSges(3,2) + (-mrSges(4,1) * t121 + mrSges(4,2) * t118 - mrSges(3,1)) * t119) * qJD(2)) * t116 + (t181 + t180 + (-t118 * t56 + t121 * t55) * qJD(3)) * mrSges(4,3) + m(5) * (t13 * t32 - t14 * t31 - t43 * t7 + t44 * t8) - m(4) * pkin(2) * t155 + m(6) * (t10 * t32 + t11 * t31 + t12 * t55 + t29 * t51 + t34 * t8 + t35 * t7) + m(7) * (t23 * t7 + t27 * t8 - t29 * t33 + t3 * t31 + t32 * t4 - t55 * t9) + (t131 * t211 + m(4) * (-t56 * t169 + t131 + t180) / 0.2e1) * t213; 0.2e1 * t3 * t77 + 0.2e1 * t14 * t78 + 0.2e1 * t11 * t79 + 0.2e1 * t10 * t80 + 0.2e1 * t12 * t58 + 0.2e1 * t9 * t59 - 0.2e1 * pkin(2) * t68 + 0.2e1 * t4 * t75 + 0.2e1 * t13 * t76 + 0.2e1 * t51 * t24 + 0.2e1 * t33 * t25 + 0.2e1 * t23 * t36 + 0.2e1 * t35 * t38 + 0.2e1 * t27 * t39 + 0.2e1 * t34 * t41 + 0.2e1 * t43 * t37 + 0.2e1 * t44 * t40 + 0.2e1 * m(6) * (t10 * t34 + t11 * t35 + t12 * t51) + (t23 * t3 + t27 * t4 + t33 * t9) * t214 + (t13 * t44 + t14 * t43) * t215 + ((0.2e1 * Ifges(4,4) * t121 + t60 * t213 + (-t159 + t182) * t120 + t129 * t117) * qJD(3) + t132) * t121 + (t26 * t213 + (t20 + t21 + t22) * t120 + (t17 + t18 - t19) * t117 + ((-0.2e1 * Ifges(4,4) + t160 * t120 + (Ifges(6,6) - t196) * t117) * t118 + ((pkin(8) ^ 2 * t215) + (2 * Ifges(4,1)) - (2 * Ifges(4,2)) - (2 * Ifges(7,3)) + t223) * t121) * qJD(3) + (t129 * t120 + (t160 * t121 + t159) * t117) * qJD(4)) * t118; -t30 * mrSges(4,2) + (t65 - t66 + t67) * t55 + (t87 - t88 + t189) * t29 + m(7) * (-t29 * t63 + t31 * t54 + t32 * t52 - t42 * t55 + t7 * t86 + t8 * t90) + m(5) * (-pkin(3) * t29 + t5) + m(6) * (t29 * t82 + t53 * t55 + t5) + 0.2e1 * (t210 + t211) * pkin(9) * (t31 * t165 - t32 * t167 + t202) + (-mrSges(5,3) + t198) * ((t117 * t32 - t120 * t31) * qJD(4) - t201 - t202); -pkin(3) * t26 + t12 * t87 + t82 * t24 + t63 * t25 + t33 * t66 + t86 * t36 + t90 * t39 + t42 * t59 + t51 * t65 + t52 * t75 + t53 * t58 + t54 * t77 + t9 * t88 + m(7) * (t23 * t54 + t27 * t52 + t3 * t86 + t33 * t42 + t4 * t90 + t63 * t9) + m(6) * (t12 * t82 + t51 * t53) + (t13 * mrSges(5,3) + t10 * mrSges(6,2) - t4 * mrSges(7,3) - t17 / 0.2e1 - t18 / 0.2e1 + t19 / 0.2e1) * t120 + (t11 * mrSges(6,2) - t14 * mrSges(5,3) - t3 * mrSges(7,3) + t20 / 0.2e1 + t21 / 0.2e1 + t22 / 0.2e1) * t117 + (t193 * t120 + t194 * t117 + m(5) * (-t14 * t117 + t13 * t120) + m(6) * (t10 * t120 + t11 * t117)) * pkin(9) + (-t114 / 0.2e1 - t111 / 0.2e1 - t113 / 0.2e1 + (Ifges(4,5) + t148 * t120 + t149 * t117 + (-m(5) * pkin(3) + t189) * pkin(8)) * qJD(3)) * t121 + (-qJD(3) * (Ifges(7,5) * t117 - Ifges(7,6) * t120) / 0.2e1 - Ifges(4,6) * qJD(3) - t117 * t71 / 0.2e1 + (mrSges(4,2) * qJD(3) + t67) * pkin(8) + t220 * t117 / 0.2e1 + ((Ifges(5,6) - Ifges(6,6)) * t120 + t222 * t117) * qJD(3) / 0.2e1 + t218 * t120 / 0.2e1) * t118 + ((t182 / 0.2e1 + t35 * mrSges(6,2) - t43 * mrSges(5,3) - t23 * mrSges(7,3) + t48 / 0.2e1 + t49 / 0.2e1 + t50 / 0.2e1 + t149 * t118) * t120 + (-t34 * mrSges(6,2) - t44 * mrSges(5,3) + t27 * mrSges(7,3) + t45 / 0.2e1 + t46 / 0.2e1 - t47 / 0.2e1 + (Ifges(7,6) / 0.2e1 + Ifges(5,6) / 0.2e1) * t121 - t148 * t118) * t117 + (t190 * t120 + t191 * t117 + m(5) * (-t117 * t44 - t120 * t43) + m(6) * (-t117 * t34 + t120 * t35)) * pkin(9)) * qJD(4); 0.2e1 * t82 * t65 + 0.2e1 * t42 * t88 + 0.2e1 * t63 * t66 + (t42 * t63 + t52 * t90 + t54 * t86) * t214 - 0.2e1 * pkin(3) * t67 + 0.2e1 * (m(6) * t82 + t87) * t53 + (t52 * t212 - t220 + t71) * t120 + (t54 * t212 + t218) * t117 + ((t86 * t212 + t94 + t95 + t96) * t120 + (0.2e1 * mrSges(7,3) * t90 + t91 + t92 - t93) * t117) * qJD(4); (-mrSges(5,1) + t199) * t7 + m(7) * (-t209 * t7 + t205) + m(6) * (-pkin(4) * t7 + t205) + (-mrSges(5,2) + t197) * t8; -t132 - t209 * t36 - pkin(4) * t38 + t10 * mrSges(6,3) - t11 * mrSges(6,1) - t13 * mrSges(5,2) + t14 * mrSges(5,1) - t3 * mrSges(7,1) + t4 * mrSges(7,2) + (Ifges(7,3) * qJD(3) + (-t160 * t117 - t196 * t120) * qJD(4)) * t118 + (t75 + t80) * qJD(5) + (t39 + t41) * qJ(5) + m(7) * (qJ(5) * t4 + qJD(5) * t27 - t209 * t3) + m(6) * (-pkin(4) * t11 + qJ(5) * t10 + qJD(5) * t34) + (-Ifges(7,5) * t120 - t196 * t117) * t168; m(7) * (qJ(5) * t52 + qJD(5) * t90 - t209 * t54) - t54 * mrSges(7,1) + t52 * mrSges(7,2) + t113 + t114 + t111 + qJD(5) * t134 + ((-mrSges(6,2) * pkin(4) + mrSges(7,3) * t209 - Ifges(7,5)) * t120 + (t198 * qJ(5) - t196) * t117 + (-m(6) * t137 + t87 + t89) * pkin(9)) * qJD(4); 0.2e1 * (t224 * qJ(5) + t197) * qJD(5); 0.2e1 * (m(7) / 0.2e1 + t210) * t7; -mrSges(6,2) * t151 + t98 + t99 + m(6) * t11 + m(7) * t3 + (t199 * t118 - t158) * qJD(3); m(7) * t54 + qJD(4) * t134; 0; 0; -m(7) * t29; m(7) * t9 + t25; m(7) * t42 + t66; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
