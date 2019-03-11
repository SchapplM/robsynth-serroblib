% Calculate time derivative of joint inertia matrix for
% S6PRRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 22:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPRR2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR2_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR2_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR2_inertiaDJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR2_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR2_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR2_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:56:03
% EndTime: 2019-03-08 21:56:10
% DurationCPUTime: 3.22s
% Computational Cost: add. (5000->407), mult. (11995->624), div. (0->0), fcn. (11847->12), ass. (0->182)
t130 = sin(qJ(5));
t134 = cos(qJ(5));
t149 = mrSges(6,1) * t130 + mrSges(6,2) * t134;
t111 = t149 * qJD(5);
t129 = sin(qJ(6));
t133 = cos(qJ(6));
t144 = t129 * t130 - t133 * t134;
t213 = qJD(5) + qJD(6);
t81 = t213 * t144;
t110 = t129 * t134 + t130 * t133;
t82 = t213 * t110;
t43 = t82 * mrSges(7,1) - t81 * mrSges(7,2);
t217 = -t111 - t43;
t135 = cos(qJ(3));
t126 = sin(pkin(12));
t131 = sin(qJ(3));
t182 = t126 * t131;
t185 = cos(pkin(12));
t139 = t135 * t185 - t182;
t100 = t139 * qJD(3);
t178 = t134 * t100;
t152 = t185 * t131;
t107 = t126 * t135 + t152;
t99 = t107 * qJD(3);
t216 = Ifges(6,5) * t178 + Ifges(6,3) * t99;
t65 = t144 * t107;
t172 = qJD(5) * t134;
t141 = t130 * t100 + t107 * t172;
t127 = sin(pkin(6));
t132 = sin(qJ(2));
t177 = qJD(2) * t132;
t161 = t127 * t177;
t128 = cos(pkin(6));
t181 = t127 * t132;
t102 = t128 * t131 + t135 * t181;
t136 = cos(qJ(2));
t176 = qJD(2) * t136;
t160 = t127 * t176;
t86 = -qJD(3) * t102 - t131 * t160;
t101 = t128 * t135 - t131 * t181;
t87 = qJD(3) * t101 + t135 * t160;
t42 = t126 * t86 + t185 * t87;
t180 = t127 * t136;
t63 = t126 * t101 + t102 * t185;
t50 = -t130 * t63 - t134 * t180;
t18 = qJD(5) * t50 + t130 * t161 + t134 * t42;
t142 = t130 * t180 - t134 * t63;
t19 = qJD(5) * t142 - t130 * t42 + t134 * t161;
t215 = -t130 * t19 + t18 * t134;
t173 = qJD(5) * t130;
t196 = -qJ(4) - pkin(8);
t154 = qJD(3) * t196;
t138 = -qJD(4) * t131 + t135 * t154;
t98 = qJD(4) * t135 + t131 * t154;
t59 = t126 * t138 + t185 * t98;
t175 = qJD(3) * t131;
t167 = pkin(3) * t175;
t60 = pkin(4) * t99 - pkin(9) * t100 + t167;
t123 = -pkin(3) * t135 - pkin(2);
t72 = -pkin(4) * t139 - pkin(9) * t107 + t123;
t117 = t196 * t135;
t89 = -t117 * t185 + t182 * t196;
t14 = t130 * t60 + t134 * t59 + t72 * t172 - t173 * t89;
t155 = -t130 * t59 + t134 * t60;
t80 = t134 * t89;
t40 = t130 * t72 + t80;
t15 = -t40 * qJD(5) + t155;
t214 = -t130 * t15 + t134 * t14;
t206 = m(5) * pkin(3);
t212 = t126 * t206 - mrSges(5,2);
t116 = -mrSges(6,1) * t134 + mrSges(6,2) * t130;
t162 = t185 * pkin(3);
t122 = -t162 - pkin(4);
t211 = m(6) * t122 - t185 * t206 - mrSges(5,1) + t116;
t210 = 2 * m(7);
t209 = -2 * mrSges(5,3);
t58 = t126 * t98 - t185 * t138;
t208 = 0.2e1 * t58;
t88 = -t117 * t126 - t196 * t152;
t207 = 0.2e1 * t88;
t205 = m(7) * pkin(5);
t201 = pkin(3) * t126;
t200 = t58 * t88;
t41 = t126 * t87 - t185 * t86;
t62 = -t101 * t185 + t102 * t126;
t26 = t62 * t41;
t199 = t81 * mrSges(7,3);
t198 = t82 * mrSges(7,3);
t121 = pkin(9) + t201;
t197 = pkin(10) + t121;
t195 = -Ifges(7,5) * t81 - Ifges(7,6) * t82;
t194 = Ifges(6,4) * t130;
t193 = Ifges(6,4) * t134;
t192 = Ifges(6,6) * t130;
t191 = t144 * mrSges(7,3);
t190 = t110 * mrSges(7,3);
t184 = t107 * t130;
t183 = t107 * t134;
t174 = qJD(3) * t135;
t171 = qJD(6) * t129;
t170 = qJD(6) * t133;
t169 = 0.2e1 * t131;
t22 = -t100 * t144 - t107 * t82;
t23 = -t110 * t100 + t213 * t65;
t168 = Ifges(7,5) * t22 + Ifges(7,6) * t23 + Ifges(7,3) * t99;
t166 = pkin(5) * t173;
t24 = t129 * t142 + t133 * t50;
t5 = qJD(6) * t24 + t129 * t19 + t133 * t18;
t25 = t129 * t50 - t133 * t142;
t6 = -qJD(6) * t25 - t129 * t18 + t133 * t19;
t165 = t6 * mrSges(7,1) - t5 * mrSges(7,2);
t164 = mrSges(6,3) * t172;
t163 = mrSges(6,3) * t173;
t159 = t107 * t173;
t68 = t99 * mrSges(5,1) + t100 * mrSges(5,2);
t157 = -t173 / 0.2e1;
t156 = -(2 * Ifges(5,4)) - t192;
t39 = -t130 * t89 + t134 * t72;
t153 = qJD(5) * t197;
t151 = t127 ^ 2 * t132 * t176;
t150 = t88 * t41 + t58 * t62;
t148 = Ifges(6,1) * t134 - t194;
t147 = -Ifges(6,2) * t130 + t193;
t146 = Ifges(6,5) * t130 + Ifges(6,6) * t134;
t27 = -pkin(5) * t139 - pkin(10) * t183 + t39;
t32 = -pkin(10) * t184 + t40;
t12 = -t129 * t32 + t133 * t27;
t13 = t129 * t27 + t133 * t32;
t103 = t197 * t130;
t104 = t197 * t134;
t70 = -t103 * t133 - t104 * t129;
t96 = t130 * t153;
t97 = t134 * t153;
t37 = qJD(6) * t70 - t129 * t97 - t133 * t96;
t71 = -t103 * t129 + t104 * t133;
t38 = -qJD(6) * t71 + t129 * t96 - t133 * t97;
t145 = t38 * mrSges(7,1) - t37 * mrSges(7,2) + t195;
t10 = -pkin(10) * t178 + pkin(5) * t99 + (-t80 + (pkin(10) * t107 - t72) * t130) * qJD(5) + t155;
t11 = -pkin(10) * t141 + t14;
t2 = qJD(6) * t12 + t10 * t129 + t11 * t133;
t3 = -qJD(6) * t13 + t10 * t133 - t11 * t129;
t143 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + t168;
t140 = t159 - t178;
t137 = -t86 * t131 + t87 * t135 + (-t101 * t135 - t102 * t131) * qJD(3);
t124 = Ifges(6,5) * t172;
t119 = Ifges(6,1) * t130 + t193;
t118 = Ifges(6,2) * t134 + t194;
t115 = -t134 * pkin(5) + t122;
t114 = t148 * qJD(5);
t113 = t147 * qJD(5);
t112 = (mrSges(4,1) * t131 + mrSges(4,2) * t135) * qJD(3);
t105 = (-mrSges(7,1) * t129 - mrSges(7,2) * t133) * qJD(6) * pkin(5);
t85 = Ifges(7,1) * t110 - Ifges(7,4) * t144;
t84 = Ifges(7,4) * t110 - Ifges(7,2) * t144;
t83 = mrSges(7,1) * t144 + mrSges(7,2) * t110;
t79 = -mrSges(5,1) * t139 + mrSges(5,2) * t107;
t74 = -mrSges(6,1) * t139 - mrSges(6,3) * t183;
t73 = mrSges(6,2) * t139 - mrSges(6,3) * t184;
t69 = t149 * t107;
t64 = t110 * t107;
t57 = pkin(5) * t184 + t88;
t53 = -Ifges(6,5) * t139 + t107 * t148;
t52 = -Ifges(6,6) * t139 + t107 * t147;
t49 = -mrSges(7,1) * t139 + mrSges(7,3) * t65;
t48 = mrSges(7,2) * t139 - mrSges(7,3) * t64;
t47 = -mrSges(6,2) * t99 - mrSges(6,3) * t141;
t46 = mrSges(6,1) * t99 + mrSges(6,3) * t140;
t45 = -Ifges(7,1) * t81 - Ifges(7,4) * t82;
t44 = -Ifges(7,4) * t81 - Ifges(7,2) * t82;
t36 = mrSges(6,1) * t141 - mrSges(6,2) * t140;
t35 = pkin(5) * t141 + t58;
t33 = mrSges(7,1) * t64 - mrSges(7,2) * t65;
t31 = -Ifges(7,1) * t65 - Ifges(7,4) * t64 - Ifges(7,5) * t139;
t30 = -Ifges(7,4) * t65 - Ifges(7,2) * t64 - Ifges(7,6) * t139;
t29 = -Ifges(6,1) * t140 - Ifges(6,4) * t141 + t99 * Ifges(6,5);
t28 = -Ifges(6,4) * t140 - Ifges(6,2) * t141 + t99 * Ifges(6,6);
t17 = -mrSges(7,2) * t99 + mrSges(7,3) * t23;
t16 = mrSges(7,1) * t99 - mrSges(7,3) * t22;
t9 = -mrSges(7,1) * t23 + mrSges(7,2) * t22;
t8 = Ifges(7,1) * t22 + Ifges(7,4) * t23 + t99 * Ifges(7,5);
t7 = Ifges(7,4) * t22 + Ifges(7,2) * t23 + t99 * Ifges(7,6);
t1 = [0.2e1 * m(7) * (t24 * t6 + t25 * t5 + t26) + 0.2e1 * m(6) * (-t142 * t18 + t19 * t50 + t26) + 0.2e1 * m(5) * (t63 * t42 - t151 + t26) + 0.2e1 * m(4) * (t101 * t86 + t102 * t87 - t151); t24 * t16 + t25 * t17 + t18 * t73 + t19 * t74 + t50 * t46 - t142 * t47 + t5 * t48 + t6 * t49 + (t36 + t9) * t62 + (t33 + t69) * t41 + (t100 * t62 + t107 * t41 + t139 * t42 - t63 * t99) * mrSges(5,3) + t137 * mrSges(4,3) + ((-t112 - t68) * t136 + (-t136 * mrSges(3,2) + (-mrSges(4,1) * t135 + mrSges(4,2) * t131 - mrSges(3,1) + t79) * t132) * qJD(2)) * t127 + m(5) * (t89 * t42 + t59 * t63 + (t123 * t177 - t136 * t167) * t127 + t150) + m(6) * (-t14 * t142 + t15 * t50 + t18 * t40 + t19 * t39 + t150) + m(7) * (t12 * t6 + t13 * t5 + t2 * t25 + t24 * t3 + t35 * t62 + t41 * t57) + (-pkin(2) * t161 + pkin(8) * t137) * m(4); (mrSges(5,3) * t207 - t130 * t52 + t134 * t53) * t100 + 0.2e1 * m(5) * (t123 * t167 + t59 * t89 + t200) + 0.2e1 * m(6) * (t14 * t40 + t15 * t39 + t200) + (-Ifges(7,5) * t65 - Ifges(7,6) * t64 + t89 * t209) * t99 - 0.2e1 * pkin(2) * t112 + 0.2e1 * t123 * t68 + 0.2e1 * t14 * t73 + 0.2e1 * t15 * t74 - t64 * t7 - t65 * t8 + 0.2e1 * t57 * t9 + 0.2e1 * t39 * t46 + 0.2e1 * t40 * t47 + 0.2e1 * t2 * t48 + 0.2e1 * t3 * t49 + t23 * t30 + t22 * t31 + 0.2e1 * t35 * t33 + 0.2e1 * t13 * t17 + 0.2e1 * t12 * t16 + (mrSges(5,3) * t208 + 0.2e1 * Ifges(5,1) * t100 - t130 * t28 + t134 * t29 + (Ifges(6,5) * t134 + t156) * t99 + (-t130 * t53 - t134 * t52 + t139 * t146) * qJD(5)) * t107 - (t59 * t209 + t156 * t100 + ((2 * Ifges(5,2)) + Ifges(6,3) + Ifges(7,3)) * t99 + t168 + t216) * t139 + (-Ifges(4,4) * t131 + pkin(3) * t79) * qJD(3) * t169 + (0.2e1 * Ifges(4,4) * t135 + (Ifges(4,1) - Ifges(4,2)) * t169) * t174 + t36 * t207 + t69 * t208 + (t12 * t3 + t13 * t2 + t35 * t57) * t210; -t5 * t191 - t25 * t198 - t6 * t190 + t24 * t199 + m(7) * (t24 * t38 + t25 * t37 + t5 * t71 + t6 * t70) + m(6) * ((t130 * t142 - t134 * t50) * qJD(5) + t215) * t121 + t142 * t163 - t50 * t164 - t87 * mrSges(4,2) + t86 * mrSges(4,1) + (m(7) * t166 - t217) * t62 + t212 * t42 + t215 * mrSges(6,3) + (m(7) * t115 + t211 + t83) * t41; t12 * t199 - t3 * t190 + t114 * t183 / 0.2e1 - t113 * t184 / 0.2e1 - Ifges(4,6) * t175 + t53 * t172 / 0.2e1 + m(7) * (t115 * t35 + t12 * t38 + t13 * t37 + t57 * t166 + t2 * t71 + t3 * t70) - t40 * t163 - t39 * t164 + t33 * t166 + t134 * t28 / 0.2e1 + t130 * t29 / 0.2e1 + (t178 / 0.2e1 + t107 * t157) * t119 + t115 * t9 + t122 * t36 + t110 * t8 / 0.2e1 + t88 * t111 - Ifges(5,6) * t99 + Ifges(5,5) * t100 - t82 * t30 / 0.2e1 + t35 * t83 + t23 * t84 / 0.2e1 + t22 * t85 / 0.2e1 - t81 * t31 / 0.2e1 - t64 * t44 / 0.2e1 - t65 * t45 / 0.2e1 + t70 * t16 + t71 * t17 + t57 * t43 + t37 * t48 + t38 * t49 + (-mrSges(4,1) * t174 + mrSges(4,2) * t175) * pkin(8) + (-t100 * t162 - t99 * t201) * mrSges(5,3) + (-t74 * t172 - t73 * t173 + m(6) * ((-t130 * t40 - t134 * t39) * qJD(5) + t214) + t134 * t47 - t130 * t46) * t121 + t214 * mrSges(6,3) - t141 * t118 / 0.2e1 + t212 * t59 + t211 * t58 + (Ifges(7,5) * t110 - Ifges(7,6) * t144 + t146) * t99 / 0.2e1 - t144 * t7 / 0.2e1 - (-Ifges(6,6) * t173 + t124 + t195) * t139 / 0.2e1 + t52 * t157 + Ifges(4,5) * t174 - t2 * t191 - t13 * t198; (t115 * t166 + t37 * t71 + t38 * t70) * t210 - t82 * t84 - t144 * t44 - t81 * t85 + t110 * t45 + 0.2e1 * t83 * t166 + 0.2e1 * t115 * t43 - t118 * t173 + 0.2e1 * t122 * t111 + t130 * t114 + (qJD(5) * t119 + t113) * t134 + 0.2e1 * (-t110 * t38 - t144 * t37 + t70 * t81 - t71 * t82) * mrSges(7,3); m(5) * t161 + m(7) * (t110 * t5 - t144 * t6 - t24 * t82 - t25 * t81) + m(6) * (t130 * t18 + t134 * t19 + (-t130 * t50 - t134 * t142) * qJD(5)); m(5) * t167 - t144 * t16 + t110 * t17 + t130 * t47 + t134 * t46 - t81 * t48 - t82 * t49 + (-t130 * t74 + t134 * t73) * qJD(5) + m(7) * (t110 * t2 - t12 * t82 - t13 * t81 - t144 * t3) + m(6) * (t130 * t14 + t134 * t15 + (-t130 * t39 + t134 * t40) * qJD(5)) + t68; m(7) * (t110 * t37 - t144 * t38 - t70 * t82 - t71 * t81); (-t110 * t81 + t144 * t82) * t210; t19 * mrSges(6,1) - t18 * mrSges(6,2) + (t129 * t5 + t133 * t6 + (-t129 * t24 + t133 * t25) * qJD(6)) * t205 + t165; -Ifges(6,5) * t159 + t15 * mrSges(6,1) - t14 * mrSges(6,2) - t141 * Ifges(6,6) + (m(7) * (-t12 * t171 + t129 * t2 + t13 * t170 + t133 * t3) + t48 * t170 + t129 * t17 - t49 * t171 + t133 * t16) * pkin(5) + t143 + t216; t124 + (t116 * t121 - t192) * qJD(5) + (m(7) * (t129 * t37 + t133 * t38 + (-t129 * t70 + t133 * t71) * qJD(6)) + (-t129 * t82 + t133 * t81 + (t110 * t129 - t133 * t144) * qJD(6)) * mrSges(7,3)) * pkin(5) + t145; (-t129 * t81 - t133 * t82 + (t110 * t133 + t129 * t144) * qJD(6)) * t205 + t217; 0.2e1 * t105; t165; t143; t145; -t43; t105; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
