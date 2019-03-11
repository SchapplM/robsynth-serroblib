% Calculate time derivative of joint inertia matrix for
% S6RRPPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
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
% Datum: 2019-03-09 09:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRR8_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR8_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR8_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR8_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR8_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR8_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR8_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:22:54
% EndTime: 2019-03-09 09:23:03
% DurationCPUTime: 3.62s
% Computational Cost: add. (4623->455), mult. (10648->677), div. (0->0), fcn. (9510->8), ass. (0->180)
t165 = sin(pkin(10));
t166 = cos(pkin(10));
t226 = qJD(3) * (t165 ^ 2 + t166 ^ 2);
t225 = Ifges(4,1) + Ifges(5,1);
t168 = sin(qJ(5));
t171 = cos(qJ(5));
t174 = t165 * t168 + t166 * t171;
t117 = t174 * qJD(5);
t132 = t165 * t171 - t166 * t168;
t169 = sin(qJ(2));
t172 = cos(qJ(2));
t192 = qJD(2) * t172;
t68 = -t117 * t169 + t132 * t192;
t110 = t132 * t169;
t69 = qJD(5) * t110 + t174 * t192;
t224 = Ifges(6,5) * t69 + Ifges(6,6) * t68;
t210 = -pkin(8) + qJ(3);
t141 = t210 * t165;
t143 = t210 * t166;
t95 = t168 * t141 + t171 * t143;
t193 = qJD(2) * t169;
t223 = qJ(4) * t193 - qJD(4) * t172;
t222 = qJD(5) + qJD(6);
t221 = t165 * (pkin(3) + pkin(4)) + pkin(7);
t220 = 2 * m(4);
t219 = 2 * m(5);
t218 = 2 * m(6);
t217 = 2 * m(7);
t216 = -2 * pkin(1);
t215 = 0.2e1 * pkin(7);
t181 = pkin(3) * t165 + pkin(7);
t179 = t166 * t192;
t199 = t166 * t169;
t196 = -qJ(4) * t179 - qJD(4) * t199;
t85 = t181 * t192 + t196;
t213 = m(5) * t85;
t167 = sin(qJ(6));
t170 = cos(qJ(6));
t111 = t174 * t169;
t61 = t110 * t170 - t111 * t167;
t23 = qJD(6) * t61 + t167 * t68 + t170 * t69;
t62 = t110 * t167 + t111 * t170;
t24 = -qJD(6) * t62 - t167 * t69 + t170 * t68;
t209 = Ifges(7,5) * t23 + Ifges(7,6) * t24;
t118 = t132 * qJD(5);
t81 = -t132 * t167 - t170 * t174;
t40 = qJD(6) * t81 - t117 * t170 - t118 * t167;
t82 = t132 * t170 - t167 * t174;
t41 = -qJD(6) * t82 + t117 * t167 - t118 * t170;
t208 = Ifges(7,5) * t40 + Ifges(7,6) * t41;
t140 = -pkin(2) * t172 - t169 * qJ(3) - pkin(1);
t201 = t165 * t172;
t154 = pkin(7) * t201;
t162 = t172 * pkin(3);
t75 = pkin(4) * t172 + t154 + t162 + (-pkin(8) * t169 - t140) * t166;
t202 = t165 * t169;
t198 = t166 * t172;
t107 = pkin(7) * t198 + t165 * t140;
t97 = -qJ(4) * t172 + t107;
t83 = pkin(8) * t202 + t97;
t37 = t168 * t75 + t171 * t83;
t207 = mrSges(5,3) * t166;
t206 = Ifges(4,4) * t165;
t205 = Ifges(4,4) * t166;
t204 = Ifges(5,5) * t165;
t203 = Ifges(5,5) * t166;
t116 = -t169 * qJD(3) + (pkin(2) * t169 - qJ(3) * t172) * qJD(2);
t200 = t166 * t116;
t197 = -Ifges(6,5) * t117 - Ifges(6,6) * t118;
t180 = t165 * t192;
t113 = mrSges(4,1) * t180 + mrSges(4,2) * t179;
t195 = qJ(3) * t226;
t191 = qJD(3) * t168;
t190 = qJD(3) * t171;
t189 = qJD(4) * t165;
t187 = qJD(5) * t168;
t186 = qJD(5) * t171;
t185 = qJD(6) * t167;
t184 = qJD(6) * t170;
t139 = -t166 * pkin(3) - t165 * qJ(4) - pkin(2);
t183 = pkin(7) * t193;
t182 = -pkin(7) * t165 - pkin(3);
t133 = -t167 * t168 + t170 * t171;
t90 = t222 * t133;
t134 = t167 * t171 + t168 * t170;
t91 = t222 * t134;
t178 = -t91 * mrSges(7,1) - t90 * mrSges(7,2);
t36 = -t168 * t83 + t171 * t75;
t94 = t171 * t141 - t143 * t168;
t106 = t166 * t140 - t154;
t119 = t166 * pkin(4) - t139;
t123 = -mrSges(5,1) * t193 + mrSges(5,2) * t179;
t34 = -t68 * mrSges(6,1) + t69 * mrSges(6,2);
t6 = -t24 * mrSges(7,1) + t23 * mrSges(7,2);
t16 = -t41 * mrSges(7,1) + t40 * mrSges(7,2);
t25 = pkin(5) * t172 - t111 * pkin(9) + t36;
t26 = pkin(9) * t110 + t37;
t12 = -t167 * t26 + t170 * t25;
t54 = -t200 + (-pkin(8) * t198 + (-pkin(4) + t182) * t169) * qJD(2);
t108 = t165 * t116;
t55 = t108 + (-pkin(7) * t199 + pkin(8) * t201) * qJD(2) + t223;
t15 = -qJD(5) * t37 - t168 * t55 + t171 * t54;
t7 = -pkin(5) * t193 - pkin(9) * t69 + t15;
t14 = t168 * t54 + t171 * t55 + t75 * t186 - t187 * t83;
t9 = pkin(9) * t68 + t14;
t2 = qJD(6) * t12 + t167 * t7 + t170 * t9;
t13 = t167 * t25 + t170 * t26;
t3 = -qJD(6) * t13 - t167 * t9 + t170 * t7;
t176 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + t209;
t65 = -pkin(9) * t132 + t94;
t66 = -pkin(9) * t174 + t95;
t32 = -t167 * t66 + t170 * t65;
t57 = t141 * t186 - t143 * t187 + t165 * t191 + t166 * t190;
t45 = -pkin(9) * t118 + t57;
t58 = -qJD(5) * t95 + t165 * t190 - t166 * t191;
t46 = pkin(9) * t117 + t58;
t10 = qJD(6) * t32 + t167 * t46 + t170 * t45;
t33 = t167 * t65 + t170 * t66;
t11 = -qJD(6) * t33 - t167 * t45 + t170 * t46;
t175 = t11 * mrSges(7,1) - t10 * mrSges(7,2) + t208;
t77 = t118 * mrSges(6,1) - t117 * mrSges(6,2);
t93 = -t166 * t183 + t108;
t153 = qJ(4) * t199;
t96 = -t169 * t221 + t153;
t74 = t192 * t221 + t196;
t173 = (-Ifges(5,4) - Ifges(4,5)) * t166 + (Ifges(4,6) - Ifges(5,6)) * t165;
t145 = mrSges(5,1) * t180;
t142 = -mrSges(5,1) * t166 - mrSges(5,3) * t165;
t138 = -mrSges(5,2) * t202 - mrSges(5,3) * t172;
t137 = mrSges(5,1) * t172 + mrSges(5,2) * t199;
t136 = -mrSges(4,1) * t172 - mrSges(4,3) * t199;
t135 = mrSges(4,2) * t172 - mrSges(4,3) * t202;
t126 = (-mrSges(7,1) * t167 - mrSges(7,2) * t170) * qJD(6) * pkin(5);
t124 = (-mrSges(5,2) * t201 + mrSges(5,3) * t169) * qJD(2);
t122 = (mrSges(4,1) * t169 - mrSges(4,3) * t198) * qJD(2);
t121 = (-mrSges(4,2) * t169 - mrSges(4,3) * t201) * qJD(2);
t120 = (mrSges(5,1) * t165 - t207) * t169;
t112 = -mrSges(5,3) * t179 + t145;
t109 = t169 * t181 - t153;
t105 = pkin(5) * t118 + t189;
t104 = (t169 * Ifges(4,5) + (t166 * Ifges(4,1) - t206) * t172) * qJD(2);
t103 = (t169 * Ifges(5,4) + (t166 * Ifges(5,1) + t204) * t172) * qJD(2);
t102 = (t169 * Ifges(4,6) + (-t165 * Ifges(4,2) + t205) * t172) * qJD(2);
t101 = (t169 * Ifges(5,6) + (t165 * Ifges(5,3) + t203) * t172) * qJD(2);
t100 = mrSges(6,1) * t172 - t111 * mrSges(6,3);
t99 = -mrSges(6,2) * t172 + t110 * mrSges(6,3);
t98 = -t106 + t162;
t92 = t165 * t183 + t200;
t89 = pkin(5) * t174 + t119;
t88 = Ifges(6,1) * t132 - Ifges(6,4) * t174;
t87 = Ifges(6,4) * t132 - Ifges(6,2) * t174;
t86 = mrSges(6,1) * t174 + mrSges(6,2) * t132;
t80 = t182 * t193 - t200;
t79 = -Ifges(6,1) * t117 - Ifges(6,4) * t118;
t78 = -Ifges(6,4) * t117 - Ifges(6,2) * t118;
t73 = t93 + t223;
t70 = -mrSges(6,1) * t110 + mrSges(6,2) * t111;
t60 = Ifges(6,1) * t111 + Ifges(6,4) * t110 + Ifges(6,5) * t172;
t59 = Ifges(6,4) * t111 + Ifges(6,2) * t110 + Ifges(6,6) * t172;
t56 = -pkin(5) * t110 + t96;
t53 = -mrSges(6,1) * t193 - mrSges(6,3) * t69;
t52 = mrSges(6,2) * t193 + mrSges(6,3) * t68;
t51 = mrSges(7,1) * t172 - t62 * mrSges(7,3);
t50 = -mrSges(7,2) * t172 + t61 * mrSges(7,3);
t44 = Ifges(7,1) * t82 + Ifges(7,4) * t81;
t43 = Ifges(7,4) * t82 + Ifges(7,2) * t81;
t42 = -mrSges(7,1) * t81 + mrSges(7,2) * t82;
t35 = t68 * pkin(5) + t74;
t31 = -mrSges(7,1) * t61 + mrSges(7,2) * t62;
t30 = Ifges(6,1) * t69 + Ifges(6,4) * t68 - Ifges(6,5) * t193;
t29 = Ifges(6,4) * t69 + Ifges(6,2) * t68 - Ifges(6,6) * t193;
t28 = Ifges(7,1) * t62 + Ifges(7,4) * t61 + Ifges(7,5) * t172;
t27 = Ifges(7,4) * t62 + Ifges(7,2) * t61 + Ifges(7,6) * t172;
t20 = mrSges(7,2) * t193 + mrSges(7,3) * t24;
t19 = -mrSges(7,1) * t193 - mrSges(7,3) * t23;
t18 = Ifges(7,1) * t40 + Ifges(7,4) * t41;
t17 = Ifges(7,4) * t40 + Ifges(7,2) * t41;
t5 = Ifges(7,1) * t23 + Ifges(7,4) * t24 - Ifges(7,5) * t193;
t4 = Ifges(7,4) * t23 + Ifges(7,2) * t24 - Ifges(7,6) * t193;
t1 = [(((mrSges(3,1) * t216) - Ifges(6,5) * t111 - Ifges(7,5) * t62 - Ifges(6,6) * t110 - Ifges(7,6) * t61 + (-(2 * Ifges(3,4)) - t173) * t169) * t169 + ((mrSges(3,2) * t216) + 0.2e1 * (Ifges(3,4) + t173) * t172 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) - (2 * Ifges(4,3)) - (2 * Ifges(5,2)) - (2 * Ifges(6,3)) - (2 * Ifges(7,3)) + pkin(7) ^ 2 * t220 + (mrSges(4,2) * t215 + t225 * t166) * t166 + (mrSges(4,1) * t215 + (Ifges(5,3) + Ifges(4,2)) * t165 + 0.2e1 * (-Ifges(4,4) + Ifges(5,5)) * t166) * t165) * t169) * t172) * qJD(2) + (t209 + t224) * t172 + (t113 * t215 + (t103 + t104) * t166 + (t101 - t102) * t165) * t169 + 0.2e1 * t93 * t135 + 0.2e1 * t92 * t136 + 0.2e1 * t80 * t137 + 0.2e1 * t73 * t138 + 0.2e1 * t85 * t120 + 0.2e1 * t107 * t121 + 0.2e1 * t106 * t122 + 0.2e1 * t98 * t123 + 0.2e1 * t97 * t124 + t110 * t29 + t111 * t30 + 0.2e1 * t109 * t112 + 0.2e1 * t15 * t100 + 0.2e1 * t96 * t34 + 0.2e1 * t14 * t99 + t68 * t59 + t69 * t60 - 0.2e1 * t74 * t70 + t61 * t4 + t62 * t5 + 0.2e1 * t2 * t50 + 0.2e1 * t3 * t51 + 0.2e1 * t37 * t52 + 0.2e1 * t36 * t53 + 0.2e1 * t56 * t6 - 0.2e1 * t35 * t31 + t24 * t27 + t23 * t28 + 0.2e1 * t12 * t19 + 0.2e1 * t13 * t20 + (t12 * t3 + t13 * t2 - t35 * t56) * t217 + (t14 * t37 + t15 * t36 - t74 * t96) * t218 + (t109 * t85 + t73 * t97 + t80 * t98) * t219 + (t106 * t92 + t107 * t93) * t220; (t112 + t213) * t139 + ((pkin(7) * mrSges(3,2) - Ifges(3,6) - Ifges(7,5) * t82 / 0.2e1 - Ifges(7,6) * t81 / 0.2e1 - Ifges(6,5) * t132 / 0.2e1 + Ifges(6,6) * t174 / 0.2e1 + (Ifges(4,6) / 0.2e1 - Ifges(5,6) / 0.2e1) * t166 + (Ifges(4,5) / 0.2e1 + Ifges(5,4) / 0.2e1) * t165) * t169 + (-t165 * (Ifges(4,2) * t166 + t206) / 0.2e1 + t165 * (-Ifges(5,3) * t166 + t204) / 0.2e1 + Ifges(3,5) + (-m(4) * pkin(2) - mrSges(4,1) * t166 + mrSges(4,2) * t165 - mrSges(3,1)) * pkin(7) + (t225 * t165 - t203 + t205) * t166 / 0.2e1) * t172) * qJD(2) + (t117 * t36 - t118 * t37 - t132 * t15 - t14 * t174) * mrSges(6,3) - t174 * t29 / 0.2e1 + (t208 + t197) * t172 / 0.2e1 + (-t12 * t40 + t13 * t41 + t2 * t81 - t3 * t82) * mrSges(7,3) + m(6) * (-t119 * t74 + t14 * t95 + t15 * t94 + t36 * t58 + t37 * t57) + (-t101 / 0.2e1 + t102 / 0.2e1 + t73 * mrSges(5,2) + t93 * mrSges(4,3) + (t135 + t138) * qJD(3) + (t121 + t124) * qJ(3) + m(4) * (qJ(3) * t93 + qJD(3) * t107) + m(5) * (qJ(3) * t73 + qJD(3) * t97)) * t166 + t132 * t30 / 0.2e1 + t85 * t142 - t117 * t60 / 0.2e1 - t118 * t59 / 0.2e1 + t119 * t34 + t110 * t78 / 0.2e1 + t111 * t79 / 0.2e1 - pkin(2) * t113 + t58 * t100 + t105 * t31 - t74 * t86 + t68 * t87 / 0.2e1 + t69 * t88 / 0.2e1 + t89 * t6 + t94 * t53 + t95 * t52 + t96 * t77 + t57 * t99 + t81 * t4 / 0.2e1 + t82 * t5 / 0.2e1 + t61 * t17 / 0.2e1 + t62 * t18 / 0.2e1 + t10 * t50 + t11 * t51 + t56 * t16 - t35 * t42 + t24 * t43 / 0.2e1 + t23 * t44 / 0.2e1 + t40 * t28 / 0.2e1 + t41 * t27 / 0.2e1 + t32 * t19 + t33 * t20 + m(7) * (t10 * t13 + t105 * t56 + t11 * t12 + t2 * t33 + t3 * t32 - t35 * t89) + (t103 / 0.2e1 + t104 / 0.2e1 + t80 * mrSges(5,2) - t92 * mrSges(4,3) + (-t136 + t137) * qJD(3) + (-t122 + t123) * qJ(3) + m(4) * (-qJ(3) * t92 - qJD(3) * t106) + m(5) * (qJ(3) * t80 + qJD(3) * t98) + (-m(5) * t109 + m(6) * t96 - t120 + t70) * qJD(4)) * t165; 0.2e1 * t105 * t42 - t117 * t88 - t118 * t87 + 0.2e1 * t119 * t77 - t174 * t78 + t132 * t79 + 0.2e1 * t89 * t16 + t81 * t17 + t82 * t18 + t40 * t44 + t41 * t43 + (t10 * t33 + t105 * t89 + t11 * t32) * t217 + t195 * t220 + (t119 * t189 + t57 * t95 + t58 * t94) * t218 + (-t139 * t189 + t195) * t219 + 0.2e1 * (mrSges(5,2) + mrSges(4,3)) * t226 + 0.2e1 * (-t142 + t86) * t189 + 0.2e1 * (t10 * t81 - t11 * t82 - t32 * t40 + t33 * t41) * mrSges(7,3) + 0.2e1 * (t117 * t94 - t118 * t95 - t132 * t58 - t174 * t57) * mrSges(6,3); t145 + (m(4) * pkin(7) - t207) * t192 + t213 + m(6) * t74 + m(7) * t35 - t6 - t34 + t113; -m(7) * t105 + (-m(5) - m(6)) * t189 - t77 - t16; 0; t133 * t19 + t134 * t20 + t168 * t52 + t171 * t53 + t90 * t50 - t91 * t51 + (-t100 * t168 + t171 * t99) * qJD(5) + m(7) * (-t12 * t91 + t13 * t90 + t133 * t3 + t134 * t2) + m(6) * (t14 * t168 + t15 * t171 + (-t168 * t36 + t171 * t37) * qJD(5)) + m(5) * t80 + t123; m(5) * t165 * qJD(3) + m(7) * (t10 * t134 + t11 * t133 - t32 * t91 + t33 * t90) + m(6) * (t168 * t57 + t171 * t58 + (-t168 * t94 + t171 * t95) * qJD(5)) + (-t133 * t40 + t134 * t41 + t81 * t90 + t82 * t91) * mrSges(7,3) + (t171 * t117 - t168 * t118 + (t132 * t168 - t171 * t174) * qJD(5)) * mrSges(6,3); 0; (-t133 * t91 + t134 * t90) * t217; t15 * mrSges(6,1) - t14 * mrSges(6,2) + (-Ifges(6,3) - Ifges(7,3)) * t193 + (m(7) * (-t12 * t185 + t13 * t184 + t167 * t2 + t170 * t3) + t50 * t184 + t167 * t20 - t51 * t185 + t170 * t19) * pkin(5) + t176 + t224; t58 * mrSges(6,1) - t57 * mrSges(6,2) + (m(7) * (t10 * t167 + t11 * t170 + (-t167 * t32 + t170 * t33) * qJD(6)) + (t167 * t41 - t170 * t40 + (t167 * t82 + t170 * t81) * qJD(6)) * mrSges(7,3)) * pkin(5) + t175 + t197; 0; (-mrSges(6,1) * t168 - mrSges(6,2) * t171) * qJD(5) + m(7) * (t167 * t90 - t170 * t91 + (-t133 * t167 + t134 * t170) * qJD(6)) * pkin(5) + t178; 0.2e1 * t126; -Ifges(7,3) * t193 + t176; t175; 0; t178; t126; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
