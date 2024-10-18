% Calculate time derivative of joint inertia matrix for
% S5RRRRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d1,d2,d3,d4,d5]';
% m [6x1]
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 22:28
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR15_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR15_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR15_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR15_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR15_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR15_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR15_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 22:24:03
% EndTime: 2024-09-27 22:24:07
% DurationCPUTime: 3.67s
% Computational Cost: add. (10797->383), mult. (31570->600), div. (0->0), fcn. (31073->12), ass. (0->181)
t156 = sin(pkin(5));
t231 = 0.2e1 * t156;
t166 = cos(qJ(2));
t195 = qJD(2) * t156;
t230 = t166 * t195;
t160 = sin(qJ(4));
t164 = cos(qJ(4));
t161 = sin(qJ(3));
t162 = sin(qJ(2));
t165 = cos(qJ(3));
t119 = (-t161 * t162 + t165 * t166) * t156;
t120 = (t161 * t166 + t162 * t165) * t156;
t84 = t119 * t164 - t120 * t160;
t228 = qJD(2) + qJD(3);
t88 = t228 * t119;
t89 = t228 * t120;
t53 = qJD(4) * t84 - t160 * t89 + t164 * t88;
t85 = t119 * t160 + t120 * t164;
t54 = -qJD(4) * t85 - t160 * t88 - t164 * t89;
t215 = Ifges(5,5) * t53 + Ifges(5,6) * t54;
t157 = cos(pkin(6));
t163 = cos(qJ(5));
t198 = t157 * t163;
t155 = sin(pkin(6));
t159 = sin(qJ(5));
t202 = t155 * t159;
t133 = pkin(4) * t198 - pkin(11) * t202;
t125 = t133 * qJD(5);
t201 = t155 * t163;
t138 = -mrSges(6,2) * t157 + mrSges(6,3) * t201;
t95 = t125 * t138;
t199 = t157 * t159;
t135 = pkin(4) * t199 + pkin(11) * t201;
t126 = t135 * qJD(5);
t137 = mrSges(6,1) * t157 - mrSges(6,3) * t202;
t96 = t126 * t137;
t229 = t95 - t96;
t132 = (-mrSges(6,1) * t163 + mrSges(6,2) * t159) * t155;
t192 = qJD(4) * t160;
t188 = pkin(3) * t192;
t113 = t155 * t132 * t188;
t153 = t155 * pkin(11);
t140 = pkin(3) * t160 + t153;
t151 = pkin(3) * t164 + pkin(4);
t108 = -t140 * t159 + t151 * t198;
t212 = pkin(3) * qJD(4);
t82 = t108 * qJD(5) + (-t160 * t199 + t163 * t164) * t212;
t72 = t82 * t138;
t109 = t140 * t163 + t151 * t199;
t83 = -t109 * qJD(5) + (-t159 * t164 - t160 * t198) * t212;
t73 = t83 * t137;
t227 = t113 + t72 + t73;
t158 = cos(pkin(5));
t220 = pkin(1) * t158;
t149 = t166 * t220;
t187 = t156 * (-pkin(8) - pkin(9));
t178 = t162 * t187;
t103 = pkin(2) * t158 + t149 + t178;
t146 = qJD(2) * t149;
t104 = qJD(2) * t178 + t146;
t148 = t162 * t220;
t105 = (t166 * t187 - t148) * qJD(2);
t200 = t156 * t166;
t136 = pkin(8) * t200 + t148;
t114 = pkin(9) * t200 + t136;
t193 = qJD(3) * t165;
t194 = qJD(3) * t161;
t61 = t103 * t193 + t165 * t104 + t161 * t105 - t114 * t194;
t75 = t161 * t103 + t165 * t114;
t62 = -qJD(3) * t75 - t104 * t161 + t165 * t105;
t86 = Ifges(4,6) * t89;
t87 = Ifges(4,5) * t88;
t226 = t62 * mrSges(4,1) - t61 * mrSges(4,2) - t86 + t87;
t174 = t155 * t158 + t157 * t84;
t64 = t159 * t174 + t163 * t85;
t218 = pkin(11) * t157;
t74 = t165 * t103 - t114 * t161;
t68 = pkin(3) * t158 - pkin(10) * t120 + t74;
t70 = pkin(10) * t119 + t75;
t44 = -t160 * t70 + t164 * t68;
t38 = pkin(4) * t158 - t218 * t85 + t44;
t139 = (-pkin(2) * t166 - pkin(1)) * t156;
t93 = -t119 * pkin(3) + t139;
t59 = -t84 * pkin(4) - t153 * t85 + t93;
t175 = t155 * t59 + t157 * t38;
t45 = t160 * t68 + t164 * t70;
t36 = pkin(11) * t174 + t45;
t17 = t159 * t175 + t163 * t36;
t225 = 2 * m(5);
t224 = 2 * m(6);
t223 = -2 * mrSges(3,3);
t190 = qJD(5) * t155;
t121 = (mrSges(6,1) * t159 + mrSges(6,2) * t163) * t190;
t222 = -0.2e1 * t121;
t221 = 0.2e1 * t139;
t219 = pkin(4) * t155;
t152 = pkin(2) * t165 + pkin(3);
t191 = qJD(4) * t164;
t197 = t160 * t161;
t91 = t152 * t191 + (-t161 * t192 + (t164 * t165 - t197) * qJD(3)) * pkin(2);
t216 = t91 * mrSges(5,2);
t214 = Ifges(6,4) * t159;
t213 = Ifges(6,4) * t163;
t186 = t162 * t195;
t127 = -pkin(8) * t186 + t146;
t211 = t127 * mrSges(3,2);
t128 = t136 * qJD(2);
t210 = t128 * mrSges(3,1);
t154 = t155 ^ 2;
t196 = t161 * t164;
t92 = -t152 * t192 + (-t161 * t191 + (-t160 * t165 - t196) * qJD(3)) * pkin(2);
t209 = t154 * t92;
t208 = t155 * t54;
t205 = t92 * t132;
t117 = Ifges(6,6) * t157 + (Ifges(6,2) * t163 + t214) * t155;
t204 = t117 * t159;
t203 = t151 * t155;
t131 = pkin(2) * t196 + t160 * t152;
t189 = 0.2e1 * mrSges(6,3);
t185 = t154 * t192;
t184 = t159 * t190;
t183 = t163 * t190;
t77 = pkin(2) * t186 + pkin(3) * t89;
t181 = mrSges(5,1) * t188;
t180 = pkin(3) * mrSges(5,2) * t191;
t179 = pkin(3) * t185;
t118 = Ifges(6,5) * t157 + (Ifges(6,1) * t159 + t213) * t155;
t122 = Ifges(6,5) * t183 - Ifges(6,6) * t184;
t123 = (-Ifges(6,2) * t159 + t213) * t190;
t124 = (Ifges(6,1) * t163 - t214) * t190;
t177 = t118 * t183 + t157 * t122 + t123 * t201 + t124 * t202;
t130 = -pkin(2) * t197 + t164 * t152;
t47 = -pkin(10) * t89 + t61;
t48 = -pkin(10) * t88 + t62;
t22 = -qJD(4) * t45 - t160 * t47 + t164 * t48;
t15 = -t218 * t53 + t22;
t31 = -pkin(4) * t54 - t153 * t53 + t77;
t176 = t15 * t157 + t155 * t31;
t63 = -t159 * t85 + t163 * t174;
t27 = qJD(5) * t63 + t163 * t53 + t199 * t54;
t28 = -t64 * qJD(5) - t159 * t53 + t54 * t198;
t10 = Ifges(6,5) * t27 + Ifges(6,6) * t28 - Ifges(6,3) * t208;
t116 = t153 + t131;
t129 = pkin(4) + t130;
t80 = -t116 * t159 + t129 * t198;
t81 = t116 * t163 + t129 * t199;
t21 = t160 * t48 + t164 * t47 + t68 * t191 - t192 * t70;
t171 = (-mrSges(4,1) * t161 - mrSges(4,2) * t165) * qJD(3) * pkin(2);
t16 = -t159 * t36 + t163 * t175;
t170 = (-t159 * t17 - t16 * t163) * qJD(5) * mrSges(6,3);
t55 = qJD(5) * t80 + t163 * t91 + t199 * t92;
t49 = t55 * t138;
t56 = -qJD(5) * t81 - t159 * t91 + t198 * t92;
t50 = t56 * t137;
t90 = t92 * mrSges(5,1);
t169 = t177 + t49 + t50 + t90 - t216;
t11 = Ifges(6,4) * t27 + Ifges(6,2) * t28 - Ifges(6,6) * t208;
t12 = Ifges(6,1) * t27 + Ifges(6,4) * t28 - Ifges(6,5) * t208;
t14 = t218 * t54 + t21;
t3 = qJD(5) * t16 + t14 * t163 + t159 * t176;
t30 = -t155 * t38 + t157 * t59;
t76 = -t155 * t84 + t157 * t158;
t34 = Ifges(6,4) * t64 + Ifges(6,2) * t63 + Ifges(6,6) * t76;
t35 = Ifges(6,1) * t64 + Ifges(6,4) * t63 + Ifges(6,5) * t76;
t4 = -t17 * qJD(5) - t14 * t159 + t176 * t163;
t6 = -t15 * t155 + t157 * t31;
t168 = -t21 * mrSges(5,2) + t12 * t202 / 0.2e1 + t6 * t132 + t4 * t137 + t3 * t138 + t22 * mrSges(5,1) + t215 + t27 * t118 / 0.2e1 + t28 * t117 / 0.2e1 + t30 * t121 - t34 * t184 / 0.2e1 - (Ifges(6,3) * t157 + (Ifges(6,5) * t159 + Ifges(6,6) * t163) * t155) * t208 / 0.2e1 + t63 * t123 / 0.2e1 + t64 * t124 / 0.2e1 + t76 * t122 / 0.2e1 + t157 * t10 / 0.2e1 + (qJD(5) * t35 + t11) * t201 / 0.2e1;
t167 = t155 * t170 + t168;
t144 = Ifges(3,5) * t230;
t134 = -pkin(8) * t156 * t162 + t149;
t107 = mrSges(4,1) * t158 - mrSges(4,3) * t120;
t106 = -mrSges(4,2) * t158 + mrSges(4,3) * t119;
t79 = mrSges(5,1) * t158 - mrSges(5,3) * t85;
t78 = -mrSges(5,2) * t158 + mrSges(5,3) * t84;
t40 = mrSges(6,1) * t76 - mrSges(6,3) * t64;
t39 = -mrSges(6,2) * t76 + mrSges(6,3) * t63;
t37 = -mrSges(6,1) * t63 + mrSges(6,2) * t64;
t19 = mrSges(6,2) * t208 + mrSges(6,3) * t28;
t18 = -mrSges(6,1) * t208 - mrSges(6,3) * t27;
t13 = -mrSges(6,1) * t28 + mrSges(6,2) * t27;
t1 = [(t144 - 0.2e1 * t210 - 0.2e1 * t211 + 0.2e1 * t215 - 0.2e1 * t86 + 0.2e1 * t87) * t158 + 0.2e1 * t88 * t120 * Ifges(4,1) + 0.2e1 * t53 * t85 * Ifges(5,1) + 0.2e1 * m(3) * (t127 * t136 - t128 * t134) + 0.2e1 * m(4) * (t61 * t75 + t62 * t74) + 0.2e1 * (t88 * t119 - t89 * t120) * Ifges(4,4) + 0.2e1 * (-t74 * t88 - t75 * t89) * mrSges(4,3) + (mrSges(4,1) * t89 + mrSges(4,2) * t88) * t221 - (Ifges(6,5) * t64 + Ifges(6,6) * t63 + Ifges(6,3) * t76) * t208 + 0.2e1 * (-t44 * t53 + t45 * t54) * mrSges(5,3) + 0.2e1 * (t53 * t84 + t54 * t85) * Ifges(5,4) + 0.2e1 * t61 * t106 + 0.2e1 * t62 * t107 + 0.2e1 * t93 * (-mrSges(5,1) * t54 + mrSges(5,2) * t53) + 0.2e1 * t77 * (-mrSges(5,1) * t84 + mrSges(5,2) * t85) + t76 * t10 + 0.2e1 * t21 * t78 + 0.2e1 * t22 * t79 + t64 * t12 + t63 * t11 + 0.2e1 * t6 * t37 + 0.2e1 * t3 * t39 + 0.2e1 * t4 * t40 + 0.2e1 * t30 * t13 + t28 * t34 + t27 * t35 + 0.2e1 * t17 * t19 + 0.2e1 * t16 * t18 + ((Ifges(3,5) * t158 + t134 * t223) * t166 + (m(4) * pkin(2) * t221 + t136 * t223 + 0.2e1 * pkin(2) * (-mrSges(4,1) * t119 + mrSges(4,2) * t120) - 0.2e1 * Ifges(3,6) * t158 + (-pkin(1) * mrSges(3,1) - Ifges(3,4) * t162) * t231) * t162) * t195 + (t16 * t4 + t17 * t3 + t30 * t6) * t224 + (t21 * t45 + t22 * t44 + t77 * t93) * t225 + ((t127 * t166 + t128 * t162) * mrSges(3,3) + (-pkin(1) * mrSges(3,2) + (Ifges(3,1) - Ifges(3,2)) * t162 + Ifges(3,4) * t166) * t230) * t231 + 0.2e1 * t84 * Ifges(5,2) * t54 - 0.2e1 * t119 * Ifges(4,2) * t89; (-t130 * t53 + t131 * t54) * mrSges(5,3) - Ifges(3,6) * t186 + (m(6) * (-t129 * t6 - t30 * t92) - t129 * t13 - t92 * t37 + t170) * t155 + t226 + t168 + (m(4) * (t161 * t61 + t165 * t62 + t193 * t75 - t194 * t74) + t106 * t193 - t107 * t194 + (-t161 * t89 - t165 * t88) * mrSges(4,3)) * pkin(2) + t91 * t78 + t92 * t79 + t80 * t18 + t81 * t19 + t56 * t40 + t55 * t39 + m(6) * (t16 * t56 + t17 * t55 + t3 * t81 + t4 * t80) - t210 - t211 + t144 + m(5) * (t130 * t22 + t131 * t21 + t44 * t92 + t45 * t91); -0.2e1 * t216 + 0.2e1 * t49 + 0.2e1 * t50 + 0.2e1 * t90 + 0.2e1 * t171 + (t129 * t209 + t55 * t81 + t56 * t80) * t224 + (t130 * t92 + t131 * t91) * t225 + (t129 * t222 - 0.2e1 * t205 + (-t204 + (-t159 * t81 - t163 * t80) * t189) * qJD(5)) * t155 + t177; m(6) * (t108 * t4 + t109 * t3 + t16 * t83 + t17 * t82 - t203 * t6) + (m(5) * (t160 * t21 + t164 * t22) + (t160 * t54 - t164 * t53) * mrSges(5,3) + ((m(5) * t45 + t78) * t164 + (-m(5) * t44 - t79 + (m(6) * t30 + t37) * t155) * t160) * qJD(4)) * pkin(3) - t13 * t203 + t108 * t18 + t109 * t19 + t82 * t39 + t83 * t40 + t167 + t226; m(6) * (t108 * t56 + t109 * t55 + t151 * t209 + t80 * t83 + t81 * t82) + t171 + ((-t160 * mrSges(5,1) - t164 * mrSges(5,2)) * qJD(4) + m(5) * (-t130 * t192 + t131 * t191 + t160 * t91 + t164 * t92) - m(6) * t129 * t185) * pkin(3) + (-t205 + (-t129 - t151) * t121 + (-t204 + ((-t108 - t80) * t163 + (-t109 - t81) * t159) * mrSges(6,3)) * qJD(5)) * t155 + t169 + t227; 0.2e1 * t113 - 0.2e1 * t180 - 0.2e1 * t181 + 0.2e1 * t72 + 0.2e1 * t73 + (t108 * t83 + t109 * t82 - t151 * t179) * t224 + (t151 * t222 + (-t204 + (-t108 * t163 - t109 * t159) * t189) * qJD(5)) * t155 + t177; m(6) * (t125 * t17 - t126 * t16 + t133 * t4 + t135 * t3 - t219 * t6) - t13 * t219 - t126 * t40 + t133 * t18 + t135 * t19 + t125 * t39 + t167; m(6) * (pkin(4) * t209 + t125 * t81 - t126 * t80 + t133 * t56 + t135 * t55) + (-t205 + (-pkin(4) - t129) * t121 + (-t204 + ((-t133 - t80) * t163 + (-t135 - t81) * t159) * mrSges(6,3)) * qJD(5)) * t155 + t169 + t229; -t180 - t181 + m(6) * (-pkin(4) * t179 - t108 * t126 + t109 * t125 + t133 * t83 + t135 * t82) + ((-pkin(4) - t151) * t121 + (-t204 + ((-t108 - t133) * t163 + (-t109 - t135) * t159) * mrSges(6,3)) * qJD(5)) * t155 + t177 + t227 + t229; 0.2e1 * t95 - 0.2e1 * t96 + (t125 * t135 - t126 * t133) * t224 + (pkin(4) * t222 + (-t204 + (-t133 * t163 - t135 * t159) * t189) * qJD(5)) * t155 + t177; mrSges(6,1) * t4 - mrSges(6,2) * t3 + t10; mrSges(6,1) * t56 - mrSges(6,2) * t55 + t122; mrSges(6,1) * t83 - mrSges(6,2) * t82 + t122; -mrSges(6,1) * t126 - mrSges(6,2) * t125 + t122; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
