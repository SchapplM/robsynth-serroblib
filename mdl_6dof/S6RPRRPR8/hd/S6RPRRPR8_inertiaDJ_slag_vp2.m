% Calculate time derivative of joint inertia matrix for
% S6RPRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 05:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPR8_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR8_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR8_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR8_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR8_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR8_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR8_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:22:42
% EndTime: 2019-03-09 05:22:49
% DurationCPUTime: 3.66s
% Computational Cost: add. (5102->479), mult. (11637->713), div. (0->0), fcn. (10296->8), ass. (0->198)
t231 = 2 * qJD(2);
t238 = Ifges(5,3) + Ifges(6,3);
t170 = sin(qJ(4));
t171 = sin(qJ(3));
t208 = qJD(3) * t171;
t192 = t170 * t208;
t173 = cos(qJ(4));
t174 = cos(qJ(3));
t203 = qJD(4) * t174;
t193 = t173 * t203;
t177 = t192 - t193;
t167 = sin(pkin(10));
t168 = cos(pkin(10));
t137 = t167 * t173 + t168 * t170;
t181 = t167 * t170 - t168 * t173;
t234 = qJD(4) * t181;
t73 = t137 * t208 + t174 * t234;
t129 = t137 * qJD(4);
t75 = -t129 * t174 + t181 * t208;
t38 = -t73 * mrSges(6,1) + t75 * mrSges(6,2);
t169 = sin(qJ(6));
t172 = cos(qJ(6));
t117 = t137 * t174;
t119 = t181 * t174;
t66 = -t117 * t172 + t119 * t169;
t26 = qJD(6) * t66 + t169 * t73 + t172 * t75;
t68 = -t117 * t169 - t119 * t172;
t28 = -qJD(6) * t68 - t169 * t75 + t172 * t73;
t6 = -t28 * mrSges(7,1) + t26 * mrSges(7,2);
t237 = -t38 - t6;
t222 = pkin(8) * t174;
t224 = pkin(3) * t171;
t147 = qJ(2) - t222 + t224;
t175 = -pkin(1) - pkin(7);
t211 = t171 * t175;
t107 = t170 * t147 + t173 * t211;
t135 = qJD(2) + (pkin(3) * t174 + pkin(8) * t171) * qJD(3);
t176 = -qJD(4) * t107 + t173 * t135;
t189 = -t170 * t175 + pkin(4);
t198 = t173 * t208;
t202 = qJD(5) * t173;
t205 = qJD(4) * t170;
t41 = qJ(5) * t198 + (qJ(5) * t205 + qJD(3) * t189 - t202) * t174 + t176;
t206 = qJD(3) * t175;
t196 = t174 * t206;
t204 = qJD(4) * t173;
t199 = t170 * t135 + t147 * t204 + t173 * t196;
t46 = -qJ(5) * t193 + (-qJD(5) * t174 + (qJ(5) * qJD(3) - qJD(4) * t175) * t171) * t170 + t199;
t14 = -t167 * t46 + t168 * t41;
t207 = qJD(3) * t174;
t10 = pkin(5) * t207 - pkin(9) * t75 + t14;
t15 = t167 * t41 + t168 * t46;
t11 = pkin(9) * t73 + t15;
t134 = t173 * t147;
t210 = t173 * t174;
t89 = -qJ(5) * t210 + t171 * t189 + t134;
t212 = t170 * t174;
t97 = -qJ(5) * t212 + t107;
t49 = -t167 * t97 + t168 * t89;
t31 = pkin(5) * t171 + pkin(9) * t119 + t49;
t50 = t167 * t89 + t168 * t97;
t35 = -pkin(9) * t117 + t50;
t12 = -t169 * t35 + t172 * t31;
t2 = qJD(6) * t12 + t10 * t169 + t11 * t172;
t13 = t169 * t31 + t172 * t35;
t3 = -qJD(6) * t13 + t10 * t172 - t11 * t169;
t236 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t90 = -t137 * t169 - t172 * t181;
t47 = qJD(6) * t90 - t129 * t169 - t172 * t234;
t91 = t137 * t172 - t169 * t181;
t48 = -qJD(6) * t91 - t129 * t172 + t169 * t234;
t16 = -t48 * mrSges(7,1) + t47 * mrSges(7,2);
t85 = t129 * mrSges(6,1) - mrSges(6,2) * t234;
t235 = -t16 - t85;
t209 = t170 ^ 2 + t173 ^ 2;
t233 = 2 * m(6);
t232 = 2 * m(7);
t230 = t90 / 0.2e1;
t229 = t91 / 0.2e1;
t227 = -t181 / 0.2e1;
t226 = t137 / 0.2e1;
t225 = -t170 / 0.2e1;
t223 = pkin(4) * t167;
t221 = -qJ(5) - pkin(8);
t220 = Ifges(5,4) * t170;
t219 = Ifges(5,4) * t173;
t218 = Ifges(5,5) * t170;
t217 = Ifges(5,6) * t170;
t216 = Ifges(5,6) * t173;
t58 = -t170 * t196 + t176;
t215 = t170 * t58;
t214 = t171 * Ifges(5,6);
t150 = -mrSges(5,1) * t173 + mrSges(5,2) * t170;
t213 = -mrSges(4,1) + t150;
t188 = qJD(4) * t221;
t125 = t170 * t188 + t202;
t126 = -qJD(5) * t170 + t173 * t188;
t77 = t168 * t125 + t167 * t126;
t149 = t221 * t170;
t151 = t221 * t173;
t99 = t167 * t149 - t168 * t151;
t201 = Ifges(7,5) * t26 + Ifges(7,6) * t28 + Ifges(7,3) * t207;
t200 = pkin(4) * t205;
t160 = -pkin(4) * t173 - pkin(3);
t197 = t171 * t207;
t157 = t171 * t206;
t195 = t171 * t205;
t194 = t170 * t203;
t116 = t137 * t171;
t118 = t181 * t171;
t65 = -t116 * t172 + t118 * t169;
t72 = -qJD(3) * t117 + t171 * t234;
t74 = -qJD(3) * t119 - t129 * t171;
t25 = qJD(6) * t65 + t169 * t72 + t172 * t74;
t67 = -t116 * t169 - t118 * t172;
t27 = -qJD(6) * t67 - t169 * t74 + t172 * t72;
t191 = t27 * mrSges(7,1) - t25 * mrSges(7,2);
t190 = -Ifges(5,5) * t173 + (2 * Ifges(4,4));
t106 = -t170 * t211 + t134;
t57 = -t175 * t195 + t199;
t187 = -t106 * qJD(4) + t57;
t76 = -t125 * t167 + t168 * t126;
t98 = t168 * t149 + t151 * t167;
t138 = pkin(4) * t212 - t174 * t175;
t44 = Ifges(7,6) * t48;
t45 = Ifges(7,5) * t47;
t78 = -pkin(9) * t137 + t98;
t79 = -pkin(9) * t181 + t99;
t36 = -t169 * t79 + t172 * t78;
t55 = pkin(9) * t234 + t76;
t56 = -pkin(9) * t129 + t77;
t8 = qJD(6) * t36 + t169 * t55 + t172 * t56;
t37 = t169 * t78 + t172 * t79;
t9 = -qJD(6) * t37 - t169 * t56 + t172 * t55;
t186 = t9 * mrSges(7,1) - t8 * mrSges(7,2) + t44 + t45;
t185 = mrSges(5,1) * t170 + mrSges(5,2) * t173;
t159 = pkin(4) * t168 + pkin(5);
t123 = t159 * t172 - t169 * t223;
t112 = t123 * qJD(6);
t124 = t159 * t169 + t172 * t223;
t113 = t124 * qJD(6);
t184 = -t113 * mrSges(7,1) - t112 * mrSges(7,2);
t183 = Ifges(5,1) * t173 - t220;
t153 = Ifges(5,1) * t170 + t219;
t182 = -Ifges(5,2) * t170 + t219;
t152 = Ifges(5,2) * t173 + t220;
t102 = -pkin(4) * t177 + t157;
t179 = Ifges(6,5) * t75 + Ifges(5,6) * t192 + Ifges(6,6) * t73 + t207 * t238 + t201;
t178 = t194 + t198;
t164 = Ifges(5,5) * t204;
t146 = mrSges(5,1) * t171 - mrSges(5,3) * t210;
t145 = -mrSges(5,2) * t171 - mrSges(5,3) * t212;
t144 = t183 * qJD(4);
t143 = t182 * qJD(4);
t142 = t185 * qJD(4);
t132 = t185 * t174;
t122 = Ifges(6,5) * t234;
t121 = Ifges(6,6) * t129;
t115 = Ifges(5,5) * t171 + t174 * t183;
t114 = t174 * t182 + t214;
t108 = pkin(5) * t181 + t160;
t105 = -mrSges(5,2) * t207 + mrSges(5,3) * t177;
t104 = mrSges(5,1) * t207 + mrSges(5,3) * t178;
t103 = pkin(5) * t129 + t200;
t101 = mrSges(6,1) * t171 + mrSges(6,3) * t119;
t100 = -mrSges(6,2) * t171 - mrSges(6,3) * t117;
t96 = Ifges(6,1) * t137 - Ifges(6,4) * t181;
t95 = Ifges(6,4) * t137 - Ifges(6,2) * t181;
t94 = mrSges(6,1) * t181 + mrSges(6,2) * t137;
t93 = t117 * pkin(5) + t138;
t88 = -mrSges(5,1) * t177 - mrSges(5,2) * t178;
t87 = -Ifges(6,1) * t234 - Ifges(6,4) * t129;
t86 = -Ifges(6,4) * t234 - Ifges(6,2) * t129;
t82 = -t153 * t203 + (Ifges(5,5) * t174 - t171 * t183) * qJD(3);
t81 = -t152 * t203 + (Ifges(5,6) * t174 - t171 * t182) * qJD(3);
t80 = mrSges(6,1) * t117 - mrSges(6,2) * t119;
t64 = -Ifges(6,1) * t119 - Ifges(6,4) * t117 + Ifges(6,5) * t171;
t63 = -Ifges(6,4) * t119 - Ifges(6,2) * t117 + Ifges(6,6) * t171;
t62 = mrSges(6,1) * t207 - mrSges(6,3) * t75;
t61 = -mrSges(6,2) * t207 + mrSges(6,3) * t73;
t60 = mrSges(7,1) * t171 - mrSges(7,3) * t68;
t59 = -mrSges(7,2) * t171 + mrSges(7,3) * t66;
t54 = -pkin(5) * t73 + t102;
t53 = Ifges(7,1) * t91 + Ifges(7,4) * t90;
t52 = Ifges(7,4) * t91 + Ifges(7,2) * t90;
t51 = -mrSges(7,1) * t90 + mrSges(7,2) * t91;
t34 = -mrSges(7,1) * t66 + mrSges(7,2) * t68;
t33 = Ifges(6,1) * t75 + Ifges(6,4) * t73 + Ifges(6,5) * t207;
t32 = Ifges(6,4) * t75 + Ifges(6,2) * t73 + Ifges(6,6) * t207;
t30 = Ifges(7,1) * t68 + Ifges(7,4) * t66 + Ifges(7,5) * t171;
t29 = Ifges(7,4) * t68 + Ifges(7,2) * t66 + Ifges(7,6) * t171;
t20 = -mrSges(7,2) * t207 + mrSges(7,3) * t28;
t19 = mrSges(7,1) * t207 - mrSges(7,3) * t26;
t18 = Ifges(7,1) * t47 + Ifges(7,4) * t48;
t17 = Ifges(7,4) * t47 + Ifges(7,2) * t48;
t5 = Ifges(7,1) * t26 + Ifges(7,4) * t28 + Ifges(7,5) * t207;
t4 = Ifges(7,4) * t26 + Ifges(7,2) * t28 + Ifges(7,6) * t207;
t1 = [0.2e1 * m(5) * (t106 * t58 + t107 * t57) + (mrSges(3,3) + (m(3) + m(4)) * qJ(2)) * t231 + 0.2e1 * t57 * t145 + 0.2e1 * t58 * t146 + 0.2e1 * t138 * t38 - t119 * t33 - t117 * t32 + 0.2e1 * t15 * t100 + 0.2e1 * t14 * t101 + 0.2e1 * t102 * t80 + 0.2e1 * t106 * t104 + 0.2e1 * t107 * t105 + 0.2e1 * t93 * t6 + t73 * t63 + t75 * t64 + t66 * t4 + t68 * t5 + 0.2e1 * t2 * t59 + 0.2e1 * t3 * t60 + 0.2e1 * t50 * t61 + 0.2e1 * t49 * t62 + 0.2e1 * t54 * t34 + t28 * t29 + t26 * t30 + 0.2e1 * t12 * t19 + 0.2e1 * t13 * t20 + ((mrSges(4,2) * t231) - t170 * t81 + t173 * t82 - 0.2e1 * t175 * t88 + (-t173 * t114 - t170 * t115 + t171 * (-t216 - t218)) * qJD(4) + (0.2e1 * qJ(2) * mrSges(4,1) + Ifges(7,5) * t68 + Ifges(7,6) * t66 - Ifges(6,5) * t119 - Ifges(6,6) * t117 + (-t190 - t217) * t174 + (-0.2e1 * m(5) * t175 ^ 2 - (2 * Ifges(4,1)) + (2 * Ifges(4,2)) + Ifges(7,3) + t238) * t171) * qJD(3)) * t174 + (t12 * t3 + t13 * t2 + t54 * t93) * t232 + (t102 * t138 + t14 * t49 + t15 * t50) * t233 + (mrSges(4,1) * t231 + (-0.2e1 * qJ(2) * mrSges(4,2) + t170 * t114 - t173 * t115 + 0.2e1 * t175 * t132 + t171 * t190) * qJD(3) + t179) * t171; t74 * t100 + t72 * t101 - t116 * t62 - t118 * t61 + t65 * t19 + t67 * t20 + t25 * t59 + t27 * t60 + (-t88 + (t145 * t173 - t146 * t170) * qJD(3) + t237) * t174 + (-t170 * t104 + t173 * t105 + (-t145 * t170 - t146 * t173) * qJD(4) + (t132 + t34 + t80) * qJD(3)) * t171 + m(7) * (t12 * t27 + t13 * t25 - t174 * t54 + t2 * t67 + t208 * t93 + t3 * t65) + m(6) * (-t102 * t174 - t116 * t14 - t118 * t15 + t138 * t208 + t49 * t72 + t50 * t74) + m(5) * ((-t106 * t170 + t107 * t173) * t207 + (-0.2e1 * t196 - t215 + t173 * t57 + (-t106 * t173 - t107 * t170) * qJD(4)) * t171); 0.2e1 * m(7) * (t25 * t67 + t27 * t65 - t197) + 0.2e1 * m(6) * (-t116 * t72 - t118 * t74 - t197) + 0.2e1 * m(5) * (-0.1e1 + t209) * t197; (-t12 * t47 + t13 * t48 + t2 * t90 - t3 * t91) * mrSges(7,3) + t160 * t38 - t129 * t63 / 0.2e1 + t138 * t85 - t119 * t87 / 0.2e1 - t117 * t86 / 0.2e1 + t98 * t62 + t99 * t61 + t77 * t100 + t76 * t101 + t102 * t94 + t103 * t34 + t108 * t6 - pkin(3) * t88 + t93 * t16 + t73 * t95 / 0.2e1 + t75 * t96 / 0.2e1 + t66 * t17 / 0.2e1 + t68 * t18 / 0.2e1 + t8 * t59 + t9 * t60 + t48 * t29 / 0.2e1 + t28 * t52 / 0.2e1 + t26 * t53 / 0.2e1 + t54 * t51 + t47 * t30 / 0.2e1 + t37 * t20 + t36 * t19 + m(7) * (t103 * t93 + t108 * t54 + t12 * t9 + t13 * t8 + t2 * t37 + t3 * t36) + m(5) * (-pkin(3) * t157 + (-t107 * t205 - t215) * pkin(8)) + t33 * t226 + t32 * t227 + t5 * t229 + t4 * t230 + m(6) * (t102 * t160 + t138 * t200 + t14 * t98 + t15 * t99 + t49 * t76 + t50 * t77) + (-t153 * t208 / 0.2e1 + t81 / 0.2e1 + qJD(4) * t115 / 0.2e1 + t187 * mrSges(5,3) + (m(5) * t187 - qJD(4) * t146 + t105) * pkin(8)) * t173 - t234 * t64 / 0.2e1 + (-t129 * t50 - t137 * t14 - t15 * t181 + t234 * t49) * mrSges(6,3) + (t164 / 0.2e1 - t122 / 0.2e1 - t121 / 0.2e1 + t45 / 0.2e1 + t44 / 0.2e1 + (t175 * t213 - Ifges(4,5)) * qJD(3)) * t171 + (t152 * t208 / 0.2e1 + t82 / 0.2e1 - t58 * mrSges(5,3) - pkin(8) * t104 + (-pkin(8) * t145 + pkin(4) * t80 - t107 * mrSges(5,3) - t214 / 0.2e1 - t114 / 0.2e1) * qJD(4)) * t170 + (t143 * t225 + t173 * t144 / 0.2e1 - t175 * t142 + (-t173 * t152 / 0.2e1 + t153 * t225) * qJD(4) + (-t175 * mrSges(4,2) - Ifges(4,6) + Ifges(7,5) * t229 + Ifges(7,6) * t230 + Ifges(6,5) * t226 + Ifges(6,6) * t227 + t218 / 0.2e1 + t216 / 0.2e1) * qJD(3)) * t174; (-t142 + t235) * t174 + m(7) * (-t103 * t174 + t25 * t37 + t27 * t36 + t65 * t9 + t67 * t8) + m(6) * (-pkin(4) * t194 - t116 * t76 - t118 * t77 + t72 * t98 + t74 * t99) + (t25 * t90 - t27 * t91 - t47 * t65 + t48 * t67) * mrSges(7,3) + (-t116 * t234 + t118 * t129 - t137 * t72 - t181 * t74) * mrSges(6,3) + ((mrSges(5,3) * t209 - mrSges(4,2)) * t174 + m(5) * (t209 * t222 - t224) + (m(6) * t160 + m(7) * t108 + t213 + t51 + t94) * t171) * qJD(3); -0.2e1 * pkin(3) * t142 + 0.2e1 * t103 * t51 + 0.2e1 * t108 * t16 - t129 * t95 - t234 * t96 - t181 * t86 + t137 * t87 + t173 * t143 + t170 * t144 + 0.2e1 * t160 * t85 + t90 * t17 + t91 * t18 + t47 * t53 + t48 * t52 + (t173 * t153 + (0.2e1 * pkin(4) * t94 - t152) * t170) * qJD(4) + (t103 * t108 + t36 * t9 + t37 * t8) * t232 + (t160 * t200 + t76 * t98 + t77 * t99) * t233 + 0.2e1 * (-t36 * t47 + t37 * t48 + t8 * t90 - t9 * t91) * mrSges(7,3) + 0.2e1 * (-t129 * t99 - t137 * t76 - t181 * t77 + t234 * t98) * mrSges(6,3); -t178 * Ifges(5,5) + (t167 * t61 + t168 * t62 + m(6) * (t14 * t168 + t15 * t167)) * pkin(4) - Ifges(5,6) * t193 + t123 * t19 + t124 * t20 + t112 * t59 - t113 * t60 + t58 * mrSges(5,1) - t57 * mrSges(5,2) + t14 * mrSges(6,1) - t15 * mrSges(6,2) + m(7) * (t112 * t13 - t113 * t12 + t123 * t3 + t124 * t2) + t179 + t236; t72 * mrSges(6,1) - t74 * mrSges(6,2) + (-t173 * t207 + t195) * mrSges(5,2) + (-t170 * t207 - t171 * t204) * mrSges(5,1) + m(7) * (t112 * t67 - t113 * t65 + t123 * t27 + t124 * t25) + m(6) * (t167 * t74 + t168 * t72) * pkin(4) + t191; m(7) * (t112 * t37 - t113 * t36 + t123 * t9 + t124 * t8) - t77 * mrSges(6,2) + t76 * mrSges(6,1) + t164 - t122 - t121 + (pkin(8) * t150 - t217) * qJD(4) + (m(6) * (t167 * t77 + t168 * t76) + (-t129 * t167 + t168 * t234) * mrSges(6,3)) * pkin(4) + (t112 * t90 + t113 * t91 - t123 * t47 + t124 * t48) * mrSges(7,3) + t186; 0.2e1 * m(7) * (t112 * t124 - t113 * t123) + 0.2e1 * t184; m(6) * t102 + m(7) * t54 - t237; (m(6) + m(7)) * t208; m(6) * t200 + m(7) * t103 - t235; 0; 0; t201 + t236; t191; t186; t184; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
