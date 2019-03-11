% Calculate time derivative of joint inertia matrix for
% S6RPRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
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
% Datum: 2019-03-09 04:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPP4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP4_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP4_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP4_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP4_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP4_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP4_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:38:48
% EndTime: 2019-03-09 04:38:57
% DurationCPUTime: 4.27s
% Computational Cost: add. (4428->382), mult. (10321->542), div. (0->0), fcn. (9832->8), ass. (0->165)
t146 = sin(pkin(9));
t147 = cos(pkin(9));
t149 = sin(qJ(3));
t151 = cos(qJ(3));
t127 = t146 * t151 + t149 * t147;
t122 = t127 * qJD(3);
t225 = Ifges(7,4) + Ifges(6,5);
t226 = t122 * t225;
t224 = Ifges(6,6) - Ifges(7,6);
t125 = t146 * t149 - t151 * t147;
t121 = t125 * qJD(3);
t145 = sin(pkin(10));
t148 = sin(qJ(4));
t181 = qJD(4) * t148;
t174 = t127 * t181;
t150 = cos(qJ(4));
t185 = t150 * t121;
t158 = t174 + t185;
t191 = cos(pkin(10));
t167 = t191 * t150;
t168 = t191 * t148;
t182 = qJD(4) * t127;
t44 = t121 * t168 + t145 * t158 - t167 * t182;
t126 = t145 * t150 + t168;
t187 = t145 * t148;
t157 = t167 - t187;
t45 = -t121 * t157 - t126 * t182;
t223 = (Ifges(6,1) + Ifges(7,1)) * t45 + (Ifges(6,4) - Ifges(7,5)) * t44 + t226;
t180 = qJD(4) * t150;
t173 = t127 * t180;
t159 = -t148 * t121 + t173;
t222 = Ifges(7,2) + Ifges(5,3) + Ifges(6,3);
t119 = t126 * qJD(4);
t120 = t157 * qJD(4);
t221 = Ifges(5,5) * t180 - t224 * t119 + t225 * t120;
t219 = m(6) + m(7);
t203 = pkin(7) + qJ(2);
t132 = t203 * t146;
t133 = t203 * t147;
t218 = -t151 * t132 - t133 * t149;
t205 = pkin(4) * t145;
t137 = qJ(6) + t205;
t217 = m(7) * t137 + mrSges(7,3);
t216 = 2 * m(6);
t215 = 0.2e1 * m(7);
t214 = -2 * mrSges(4,3);
t212 = -0.2e1 * t218;
t211 = m(6) * pkin(4);
t160 = qJ(5) * t121 - qJD(5) * t127;
t64 = -t125 * qJD(2) + qJD(3) * t218;
t83 = pkin(3) * t122 + pkin(8) * t121;
t170 = -t148 * t64 + t150 * t83;
t176 = -pkin(2) * t147 - pkin(1);
t85 = pkin(3) * t125 - pkin(8) * t127 + t176;
t97 = -t149 * t132 + t133 * t151;
t88 = t150 * t97;
t7 = pkin(4) * t122 + t160 * t150 + (-t88 + (qJ(5) * t127 - t85) * t148) * qJD(4) + t170;
t178 = t148 * t83 + t150 * t64 + t85 * t180;
t9 = -qJ(5) * t173 + (-qJD(4) * t97 + t160) * t148 + t178;
t4 = t145 * t7 + t191 * t9;
t210 = t122 / 0.2e1;
t207 = -t127 / 0.2e1;
t197 = Ifges(5,4) * t148;
t135 = Ifges(5,2) * t150 + t197;
t206 = -t135 / 0.2e1;
t65 = qJD(2) * t127 + qJD(3) * t97;
t204 = t65 * t218;
t202 = -qJ(5) - pkin(8);
t23 = -mrSges(6,2) * t122 + mrSges(6,3) * t44;
t26 = mrSges(7,2) * t44 + mrSges(7,3) * t122;
t201 = t23 + t26;
t24 = mrSges(6,1) * t122 - mrSges(6,3) * t45;
t25 = -t122 * mrSges(7,1) + t45 * mrSges(7,2);
t200 = -t24 + t25;
t189 = t127 * t150;
t50 = -t148 * t97 + t150 * t85;
t28 = pkin(4) * t125 - qJ(5) * t189 + t50;
t190 = t127 * t148;
t51 = t148 * t85 + t88;
t43 = -qJ(5) * t190 + t51;
t13 = t145 * t28 + t191 * t43;
t71 = t126 * t127;
t56 = -mrSges(6,2) * t125 - mrSges(6,3) * t71;
t59 = -mrSges(7,2) * t71 + mrSges(7,3) * t125;
t199 = t56 + t59;
t72 = t157 * t127;
t57 = mrSges(6,1) * t125 - mrSges(6,3) * t72;
t58 = -mrSges(7,1) * t125 + mrSges(7,2) * t72;
t198 = -t57 + t58;
t196 = Ifges(5,4) * t150;
t195 = t120 * mrSges(7,2);
t194 = t122 * Ifges(5,5);
t193 = t122 * Ifges(5,6);
t192 = t125 * Ifges(5,6);
t179 = qJD(6) * t126;
t177 = pkin(4) * t181;
t141 = -pkin(4) * t150 - pkin(3);
t175 = t191 * pkin(4);
t19 = -t44 * mrSges(6,1) + t45 * mrSges(6,2);
t18 = -t44 * mrSges(7,1) - t45 * mrSges(7,3);
t169 = qJD(4) * t202;
t118 = qJD(5) * t150 + t148 * t169;
t153 = -qJD(5) * t148 + t150 * t169;
t68 = t118 * t145 - t153 * t191;
t69 = t118 * t191 + t145 * t153;
t134 = t202 * t150;
t98 = -t134 * t145 - t168 * t202;
t99 = -t134 * t191 + t187 * t202;
t172 = t68 * t98 + t99 * t69;
t171 = -Ifges(5,6) * t148 - (2 * Ifges(4,4));
t166 = t122 * mrSges(4,1) - t121 * mrSges(4,2);
t78 = t119 * mrSges(6,1) + t120 * mrSges(6,2);
t77 = t119 * mrSges(7,1) - t120 * mrSges(7,3);
t66 = pkin(4) * t190 - t218;
t163 = mrSges(5,1) * t148 + mrSges(5,2) * t150;
t162 = Ifges(5,1) * t150 - t197;
t161 = -Ifges(5,2) * t148 + t196;
t3 = -t145 * t9 + t191 * t7;
t12 = -t145 * t43 + t191 * t28;
t155 = -t77 - t78;
t154 = -Ifges(5,5) * t185 + t222 * t122 + t224 * t44 + t225 * t45;
t46 = pkin(4) * t159 + t65;
t140 = -t175 - pkin(5);
t136 = Ifges(5,1) * t148 + t196;
t131 = t162 * qJD(4);
t130 = t161 * qJD(4);
t129 = t163 * qJD(4);
t94 = Ifges(6,1) * t126 + Ifges(6,4) * t157;
t93 = Ifges(7,1) * t126 - Ifges(7,5) * t157;
t92 = Ifges(6,4) * t126 + Ifges(6,2) * t157;
t91 = Ifges(7,5) * t126 - Ifges(7,3) * t157;
t90 = -mrSges(6,1) * t157 + mrSges(6,2) * t126;
t89 = -mrSges(7,1) * t157 - mrSges(7,3) * t126;
t87 = mrSges(5,1) * t125 - mrSges(5,3) * t189;
t86 = -mrSges(5,2) * t125 - mrSges(5,3) * t190;
t84 = -pkin(5) * t157 - qJ(6) * t126 + t141;
t82 = Ifges(6,1) * t120 - Ifges(6,4) * t119;
t81 = Ifges(7,1) * t120 + Ifges(7,5) * t119;
t80 = Ifges(6,4) * t120 - Ifges(6,2) * t119;
t79 = Ifges(7,5) * t120 + Ifges(7,3) * t119;
t62 = Ifges(5,5) * t125 + t127 * t162;
t61 = t127 * t161 + t192;
t55 = pkin(5) * t119 - qJ(6) * t120 + t177 - t179;
t54 = -mrSges(5,2) * t122 - mrSges(5,3) * t159;
t53 = mrSges(5,1) * t122 + mrSges(5,3) * t158;
t49 = mrSges(5,1) * t159 - mrSges(5,2) * t158;
t48 = mrSges(6,1) * t71 + mrSges(6,2) * t72;
t47 = mrSges(7,1) * t71 - mrSges(7,3) * t72;
t34 = -Ifges(5,1) * t158 - Ifges(5,4) * t159 + t194;
t33 = -Ifges(5,4) * t158 - Ifges(5,2) * t159 + t193;
t32 = Ifges(6,1) * t72 - Ifges(6,4) * t71 + Ifges(6,5) * t125;
t31 = Ifges(7,1) * t72 + Ifges(7,4) * t125 + Ifges(7,5) * t71;
t30 = Ifges(6,4) * t72 - Ifges(6,2) * t71 + Ifges(6,6) * t125;
t29 = Ifges(7,5) * t72 + Ifges(7,6) * t125 + Ifges(7,3) * t71;
t22 = pkin(5) * t71 - qJ(6) * t72 + t66;
t21 = -t51 * qJD(4) + t170;
t20 = -t181 * t97 + t178;
t15 = Ifges(6,4) * t45 + Ifges(6,2) * t44 + t122 * Ifges(6,6);
t14 = Ifges(7,5) * t45 + t122 * Ifges(7,6) - Ifges(7,3) * t44;
t11 = -t125 * pkin(5) - t12;
t10 = qJ(6) * t125 + t13;
t5 = -t44 * pkin(5) - t45 * qJ(6) - t72 * qJD(6) + t46;
t2 = -t122 * pkin(5) - t3;
t1 = qJ(6) * t122 + qJD(6) * t125 + t4;
t6 = [(-0.2e1 * Ifges(4,1) * t121 - t148 * t33 + t150 * t34 + (Ifges(5,5) * t150 + t171) * t122 + (-t150 * t61 - t148 * t62 + t125 * (-Ifges(5,5) * t148 - Ifges(5,6) * t150)) * qJD(4) + 0.2e1 * (t163 + mrSges(4,3)) * t65) * t127 + t49 * t212 + (t1 * t10 + t11 * t2 + t22 * t5) * t215 + (t12 * t3 + t13 * t4 + t46 * t66) * t216 - (mrSges(4,3) * t212 - t148 * t61 + t150 * t62) * t121 + 0.2e1 * (m(3) * qJ(2) + mrSges(3,3)) * (t146 ^ 2 + t147 ^ 2) * qJD(2) + 0.2e1 * t176 * t166 + (t32 + t31) * t45 + 0.2e1 * m(5) * (t20 * t51 + t21 * t50 - t204) + 0.2e1 * m(4) * (t64 * t97 - t204) + (t30 - t29) * t44 + 0.2e1 * t20 * t86 + 0.2e1 * t21 * t87 + 0.2e1 * t66 * t19 + 0.2e1 * t51 * t54 + 0.2e1 * t4 * t56 + 0.2e1 * t3 * t57 + 0.2e1 * t2 * t58 + 0.2e1 * t1 * t59 + 0.2e1 * t46 * t48 + 0.2e1 * t50 * t53 + 0.2e1 * t5 * t47 + 0.2e1 * t13 * t23 + 0.2e1 * t12 * t24 + 0.2e1 * t11 * t25 + 0.2e1 * t10 * t26 + 0.2e1 * t22 * t18 + (-t122 * t224 + t14 - t15) * t71 + (t64 * t214 - t171 * t121 + ((2 * Ifges(4,2)) + t222) * t122 + t154) * t125 + (t223 + t226) * t72 + t97 * t122 * t214; t148 * t54 + t150 * t53 + t201 * t126 - t200 * t157 + t199 * t120 + t198 * t119 + (-t148 * t87 + t150 * t86) * qJD(4) + m(6) * (-t119 * t12 + t120 * t13 + t126 * t4 + t157 * t3) + m(7) * (t1 * t126 + t10 * t120 + t11 * t119 - t157 * t2) + m(5) * (t148 * t20 + t150 * t21 + (-t148 * t50 + t150 * t51) * qJD(4)) + t166; 0.2e1 * t219 * (-t119 * t157 + t126 * t120); (t32 / 0.2e1 + t31 / 0.2e1) * t120 + (t29 / 0.2e1 - t30 / 0.2e1) * t119 + (t194 / 0.2e1 + t65 * mrSges(5,2) + t34 / 0.2e1 - t21 * mrSges(5,3) + t130 * t207 - t121 * t206 + (t136 * t207 + pkin(4) * t48 - t51 * mrSges(5,3) - t192 / 0.2e1 - t61 / 0.2e1 + t66 * t211) * qJD(4) + (-m(5) * t21 - t53 + (-m(5) * t51 - t86) * qJD(4)) * pkin(8)) * t148 + (t193 / 0.2e1 - t65 * mrSges(5,1) + t33 / 0.2e1 + t20 * mrSges(5,3) + t127 * t131 / 0.2e1 - t121 * t136 / 0.2e1 + (t127 * t206 - t50 * mrSges(5,3) + t62 / 0.2e1) * qJD(4) + (-qJD(4) * t87 + t54 + m(5) * (-t50 * qJD(4) + t20)) * pkin(8)) * t150 + (t223 / 0.2e1 + t2 * mrSges(7,2) - t3 * mrSges(6,3) + t225 * t210) * t126 + t198 * t68 + t199 * t69 + t200 * t98 + t201 * t99 + (-m(5) * t65 - t49) * pkin(3) + m(6) * (-t12 * t68 + t13 * t69 + t141 * t46 - t3 * t98 + t4 * t99) - t218 * t129 + m(7) * (t1 * t99 + t10 * t69 + t11 * t68 + t2 * t98 + t22 * t55 + t5 * t84) + t141 * t19 + (t81 / 0.2e1 + t82 / 0.2e1) * t72 - Ifges(4,6) * t122 - Ifges(4,5) * t121 + (-t80 / 0.2e1 + t79 / 0.2e1) * t71 + t84 * t18 + t5 * t89 + t46 * t90 + t22 * t77 + t66 * t78 - t64 * mrSges(4,2) - t65 * mrSges(4,1) + t55 * t47 + (t1 * mrSges(7,2) + t4 * mrSges(6,3) - t14 / 0.2e1 + t15 / 0.2e1 + t224 * t210) * t157 + (-t10 * t119 + t11 * t120) * mrSges(7,2) + t221 * t125 / 0.2e1 + (-t119 * t13 - t12 * t120) * mrSges(6,3) + (t93 / 0.2e1 + t94 / 0.2e1) * t45 + (-t91 / 0.2e1 + t92 / 0.2e1) * t44; t219 * (t119 * t98 + t99 * t120 + t69 * t126 - t157 * t68); -0.2e1 * pkin(3) * t129 + t150 * t130 + t148 * t131 + 0.2e1 * t141 * t78 + 0.2e1 * t55 * t89 + 0.2e1 * t84 * t77 + (t81 + t82) * t126 - (t79 - t80) * t157 + (t93 + t94) * t120 + (t91 - t92) * t119 + (t150 * t136 + (0.2e1 * pkin(4) * t90 - t135) * t148) * qJD(4) + (t55 * t84 + t172) * t215 + (t141 * t177 + t172) * t216 + 0.2e1 * (mrSges(6,3) + mrSges(7,2)) * (-t119 * t99 + t120 * t98 + t126 * t68 + t157 * t69); t154 + (t145 * t4 + t191 * t3) * t211 + m(7) * (qJD(6) * t10 + t1 * t137 + t140 * t2) - Ifges(5,5) * t174 + t140 * t25 + t137 * t26 + qJD(6) * t59 + t21 * mrSges(5,1) - t20 * mrSges(5,2) - t4 * mrSges(6,2) - t2 * mrSges(7,1) + t3 * mrSges(6,1) + t1 * mrSges(7,3) + t24 * t175 + t23 * t205 - t159 * Ifges(5,6); -mrSges(5,2) * t180 - mrSges(5,1) * t181 + (-t119 * t191 + t120 * t145) * t211 + m(7) * (t119 * t140 + t120 * t137 + t179) + t155; m(7) * qJD(6) * t99 - Ifges(5,6) * t181 + t140 * t195 + (t145 * t211 - mrSges(6,2) + t217) * t69 + (m(7) * t140 - t191 * t211 - mrSges(6,1) - mrSges(7,1)) * t68 + (-mrSges(5,1) * t180 + mrSges(5,2) * t181) * pkin(8) + (-t119 * t205 - t120 * t175) * mrSges(6,3) + (qJD(6) * t157 - t119 * t137) * mrSges(7,2) + t221; 0.2e1 * t217 * qJD(6); m(6) * t46 + m(7) * t5 + t18 + t19; 0; m(6) * t177 + m(7) * t55 - t155; 0; 0; m(7) * t2 + t25; m(7) * t119; m(7) * t68 + t195; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t6(1) t6(2) t6(4) t6(7) t6(11) t6(16); t6(2) t6(3) t6(5) t6(8) t6(12) t6(17); t6(4) t6(5) t6(6) t6(9) t6(13) t6(18); t6(7) t6(8) t6(9) t6(10) t6(14) t6(19); t6(11) t6(12) t6(13) t6(14) t6(15) t6(20); t6(16) t6(17) t6(18) t6(19) t6(20) t6(21);];
Mq  = res;
