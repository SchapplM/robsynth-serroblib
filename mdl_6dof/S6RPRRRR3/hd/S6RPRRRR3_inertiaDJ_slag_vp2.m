% Calculate time derivative of joint inertia matrix for
% S6RPRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 07:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRR3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR3_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR3_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR3_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR3_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR3_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR3_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:00:25
% EndTime: 2019-03-09 07:00:34
% DurationCPUTime: 4.22s
% Computational Cost: add. (7681->509), mult. (17500->751), div. (0->0), fcn. (15631->10), ass. (0->201)
t169 = sin(qJ(3));
t167 = sin(qJ(5));
t168 = sin(qJ(4));
t171 = cos(qJ(5));
t172 = cos(qJ(4));
t181 = t167 * t168 - t171 * t172;
t173 = cos(qJ(3));
t207 = qJD(3) * t173;
t132 = t167 * t172 + t168 * t171;
t240 = qJD(4) + qJD(5);
t99 = t240 * t132;
t64 = -t169 * t99 - t181 * t207;
t115 = t181 * t169;
t65 = t115 * t240 - t132 * t207;
t29 = -t65 * mrSges(6,1) + t64 * mrSges(6,2);
t166 = sin(qJ(6));
t170 = cos(qJ(6));
t114 = t132 * t169;
t76 = -t114 * t170 + t115 * t166;
t25 = qJD(6) * t76 + t166 * t65 + t170 * t64;
t77 = -t114 * t166 - t115 * t170;
t26 = -qJD(6) * t77 - t166 * t64 + t170 * t65;
t6 = -t26 * mrSges(7,1) + t25 * mrSges(7,2);
t179 = -t29 - t6;
t204 = qJD(4) * t172;
t176 = t168 * t207 + t169 * t204;
t205 = qJD(4) * t169;
t192 = t168 * t205;
t193 = t172 * t207;
t177 = -t192 + t193;
t91 = t176 * mrSges(5,1) + t177 * mrSges(5,2);
t243 = -t91 + t179;
t197 = -cos(pkin(11)) * pkin(1) - pkin(2);
t242 = 0.2e1 * t197;
t232 = -pkin(9) - pkin(8);
t150 = t232 * t168;
t151 = t232 * t172;
t105 = t167 * t150 - t171 * t151;
t127 = -pkin(3) * t173 - t169 * pkin(8) + t197;
t156 = sin(pkin(11)) * pkin(1) + pkin(7);
t211 = t172 * t173;
t141 = t156 * t211;
t97 = t168 * t127 + t141;
t208 = qJD(3) * t169;
t241 = -Ifges(5,5) * t193 - Ifges(5,3) * t208;
t209 = t168 ^ 2 + t172 ^ 2;
t239 = -Ifges(6,5) * t64 - Ifges(6,6) * t65 - Ifges(6,3) * t208;
t238 = 2 * m(5);
t237 = 2 * m(6);
t236 = 2 * m(7);
t235 = 0.2e1 * t156;
t92 = -t132 * t166 - t170 * t181;
t234 = t92 / 0.2e1;
t93 = t132 * t170 - t166 * t181;
t233 = t93 / 0.2e1;
t231 = -t181 / 0.2e1;
t230 = t132 / 0.2e1;
t224 = Ifges(5,4) * t168;
t147 = Ifges(5,2) * t172 + t224;
t229 = -t147 / 0.2e1;
t228 = -t168 / 0.2e1;
t227 = pkin(3) * t169;
t226 = pkin(8) * t173;
t157 = pkin(4) * t171 + pkin(5);
t200 = qJD(6) * t170;
t201 = qJD(6) * t166;
t215 = t166 * t167;
t89 = t157 * t200 + (-t167 * t201 + (t170 * t171 - t215) * qJD(5)) * pkin(4);
t225 = t89 * mrSges(7,2);
t117 = t172 * t127;
t212 = t169 * t172;
t78 = -pkin(9) * t212 + t117 + (-t156 * t168 - pkin(4)) * t173;
t213 = t168 * t169;
t87 = -pkin(9) * t213 + t97;
t47 = t167 * t78 + t171 * t87;
t223 = Ifges(5,4) * t172;
t222 = Ifges(5,5) * t168;
t221 = Ifges(5,6) * t168;
t220 = Ifges(5,6) * t172;
t144 = (-t226 + t227) * qJD(3);
t195 = t156 * t208;
t210 = t172 * t144 + t168 * t195;
t63 = -t97 * qJD(4) + t210;
t219 = t168 * t63;
t218 = t173 * Ifges(5,6);
t146 = -mrSges(5,1) * t172 + mrSges(5,2) * t168;
t217 = -mrSges(4,1) + t146;
t216 = t156 * t173;
t214 = t167 * t170;
t118 = pkin(4) * t213 + t169 * t156;
t206 = qJD(4) * t168;
t203 = qJD(5) * t167;
t202 = qJD(5) * t171;
t199 = -Ifges(7,5) * t25 - Ifges(7,6) * t26 - Ifges(7,3) * t208;
t198 = pkin(4) * t206;
t145 = t156 * t207;
t103 = pkin(4) * t176 + t145;
t158 = -pkin(4) * t172 - pkin(3);
t196 = qJD(4) * t232;
t194 = t169 * t207;
t191 = t173 * t206;
t90 = -t157 * t201 + (-t167 * t200 + (-t166 * t171 - t214) * qJD(5)) * pkin(4);
t86 = t90 * mrSges(7,1);
t188 = t86 - t225;
t187 = (2 * Ifges(4,4)) + t221;
t46 = -t167 * t87 + t171 * t78;
t62 = t127 * t204 + t168 * t144 + (-t172 * t208 - t191) * t156;
t96 = -t168 * t216 + t117;
t186 = -qJD(4) * t96 + t62;
t104 = t171 * t150 + t151 * t167;
t142 = t168 * t196;
t143 = t172 * t196;
t67 = t171 * t142 + t167 * t143 + t150 * t202 + t151 * t203;
t38 = -pkin(10) * t99 + t67;
t68 = -qJD(5) * t105 - t142 * t167 + t171 * t143;
t98 = t240 * t181;
t39 = pkin(10) * t98 + t68;
t83 = -pkin(10) * t132 + t104;
t84 = -pkin(10) * t181 + t105;
t44 = -t166 * t84 + t170 * t83;
t10 = qJD(6) * t44 + t166 * t39 + t170 * t38;
t45 = t166 * t83 + t170 * t84;
t11 = -qJD(6) * t45 - t166 * t38 + t170 * t39;
t35 = -qJD(6) * t93 + t166 * t98 - t170 * t99;
t32 = Ifges(7,6) * t35;
t34 = qJD(6) * t92 - t166 * t99 - t170 * t98;
t33 = Ifges(7,5) * t34;
t185 = t11 * mrSges(7,1) - t10 * mrSges(7,2) + t32 + t33;
t184 = mrSges(5,1) * t168 + mrSges(5,2) * t172;
t183 = Ifges(5,1) * t172 - t224;
t148 = Ifges(5,1) * t168 + t223;
t182 = -Ifges(5,2) * t168 + t223;
t30 = -pkin(5) * t173 + t115 * pkin(10) + t46;
t31 = -pkin(10) * t114 + t47;
t15 = -t166 * t31 + t170 * t30;
t16 = t166 * t30 + t170 * t31;
t48 = (pkin(4) * t169 - pkin(9) * t211) * qJD(3) + (-t141 + (pkin(9) * t169 - t127) * t168) * qJD(4) + t210;
t50 = -pkin(9) * t176 + t62;
t14 = -qJD(5) * t47 - t167 * t50 + t171 * t48;
t7 = pkin(5) * t208 - pkin(10) * t64 + t14;
t13 = t167 * t48 + t171 * t50 + t78 * t202 - t203 * t87;
t8 = pkin(10) * t65 + t13;
t2 = qJD(6) * t15 + t166 * t7 + t170 * t8;
t3 = -qJD(6) * t16 - t166 * t8 + t170 * t7;
t180 = t3 * mrSges(7,1) - t2 * mrSges(7,2) - t199;
t178 = (-mrSges(6,1) * t167 - mrSges(6,2) * t171) * qJD(5) * pkin(4);
t94 = Ifges(6,6) * t99;
t95 = Ifges(6,5) * t98;
t175 = t68 * mrSges(6,1) - t67 * mrSges(6,2) + t185 - t94 - t95;
t174 = t14 * mrSges(6,1) - t13 * mrSges(6,2) + t180 - t239;
t162 = Ifges(5,5) * t204;
t140 = -mrSges(5,1) * t173 - mrSges(5,3) * t212;
t139 = mrSges(5,2) * t173 - mrSges(5,3) * t213;
t138 = t183 * qJD(4);
t137 = t182 * qJD(4);
t136 = t184 * qJD(4);
t128 = (-mrSges(7,1) * t166 - mrSges(7,2) * t170) * qJD(6) * pkin(5);
t126 = t184 * t169;
t120 = pkin(4) * t214 + t157 * t166;
t119 = -pkin(4) * t215 + t157 * t170;
t113 = -Ifges(5,5) * t173 + t169 * t183;
t112 = t169 * t182 - t218;
t110 = pkin(5) * t181 + t158;
t109 = -mrSges(5,2) * t208 - mrSges(5,3) * t176;
t108 = mrSges(5,1) * t208 - mrSges(5,3) * t177;
t107 = -mrSges(6,1) * t173 + t115 * mrSges(6,3);
t106 = mrSges(6,2) * t173 - t114 * mrSges(6,3);
t102 = Ifges(6,1) * t132 - Ifges(6,4) * t181;
t101 = Ifges(6,4) * t132 - Ifges(6,2) * t181;
t100 = mrSges(6,1) * t181 + mrSges(6,2) * t132;
t88 = pkin(5) * t114 + t118;
t85 = pkin(5) * t99 + t198;
t82 = mrSges(6,1) * t114 - mrSges(6,2) * t115;
t81 = -t148 * t205 + (Ifges(5,5) * t169 + t173 * t183) * qJD(3);
t80 = -t147 * t205 + (Ifges(5,6) * t169 + t173 * t182) * qJD(3);
t75 = -Ifges(6,1) * t115 - Ifges(6,4) * t114 - Ifges(6,5) * t173;
t74 = -Ifges(6,4) * t115 - Ifges(6,2) * t114 - Ifges(6,6) * t173;
t70 = -mrSges(7,1) * t173 - t77 * mrSges(7,3);
t69 = mrSges(7,2) * t173 + t76 * mrSges(7,3);
t58 = -Ifges(6,1) * t98 - Ifges(6,4) * t99;
t57 = -Ifges(6,4) * t98 - Ifges(6,2) * t99;
t56 = mrSges(6,1) * t99 - mrSges(6,2) * t98;
t55 = -mrSges(6,2) * t208 + mrSges(6,3) * t65;
t54 = mrSges(6,1) * t208 - mrSges(6,3) * t64;
t53 = Ifges(7,1) * t93 + Ifges(7,4) * t92;
t52 = Ifges(7,4) * t93 + Ifges(7,2) * t92;
t51 = -mrSges(7,1) * t92 + mrSges(7,2) * t93;
t41 = -pkin(5) * t65 + t103;
t40 = -mrSges(7,1) * t76 + mrSges(7,2) * t77;
t37 = Ifges(7,1) * t77 + Ifges(7,4) * t76 - Ifges(7,5) * t173;
t36 = Ifges(7,4) * t77 + Ifges(7,2) * t76 - Ifges(7,6) * t173;
t28 = Ifges(6,1) * t64 + Ifges(6,4) * t65 + Ifges(6,5) * t208;
t27 = Ifges(6,4) * t64 + Ifges(6,2) * t65 + Ifges(6,6) * t208;
t21 = -mrSges(7,2) * t208 + mrSges(7,3) * t26;
t20 = mrSges(7,1) * t208 - mrSges(7,3) * t25;
t19 = Ifges(7,1) * t34 + Ifges(7,4) * t35;
t18 = Ifges(7,4) * t34 + Ifges(7,2) * t35;
t17 = -mrSges(7,1) * t35 + mrSges(7,2) * t34;
t5 = Ifges(7,1) * t25 + Ifges(7,4) * t26 + Ifges(7,5) * t208;
t4 = Ifges(7,4) * t25 + Ifges(7,2) * t26 + Ifges(7,6) * t208;
t1 = [(t91 * t235 - t168 * t80 + t172 * t81 + (-t168 * t113 - t172 * t112 - t173 * (-t220 - t222)) * qJD(4) + (-Ifges(6,5) * t115 - Ifges(6,6) * t114 + Ifges(7,5) * t77 + Ifges(7,6) * t76 + mrSges(4,1) * t242 + (Ifges(5,5) * t172 - t187) * t169 + (t156 ^ 2 * t238 + (2 * Ifges(4,1)) - (2 * Ifges(4,2)) - Ifges(5,3) - Ifges(6,3) - Ifges(7,3)) * t173) * qJD(3)) * t169 + ((mrSges(4,2) * t242 - t168 * t112 + t172 * t113 + t126 * t235 + t173 * t187) * qJD(3) + t199 + t239 + t241) * t173 + (t15 * t3 + t16 * t2 + t41 * t88) * t236 + (t103 * t118 + t13 * t47 + t14 * t46) * t237 + (t97 * t62 + t96 * t63) * t238 + 0.2e1 * t15 * t20 + 0.2e1 * t16 * t21 + t26 * t36 + t25 * t37 + 0.2e1 * t41 * t40 + 0.2e1 * t46 * t54 + 0.2e1 * t47 * t55 + 0.2e1 * t2 * t69 + 0.2e1 * t3 * t70 + t65 * t74 + t64 * t75 + t76 * t4 + t77 * t5 + 0.2e1 * t88 * t6 + 0.2e1 * t103 * t82 + 0.2e1 * t13 * t106 + 0.2e1 * t14 * t107 + 0.2e1 * t96 * t108 + 0.2e1 * t97 * t109 - t114 * t27 - t115 * t28 + 0.2e1 * t118 * t29 + 0.2e1 * t62 * t139 + 0.2e1 * t63 * t140; t64 * t106 + t65 * t107 - t114 * t54 - t115 * t55 + t76 * t20 + t77 * t21 + t25 * t69 + t26 * t70 + ((t139 * t172 - t140 * t168) * qJD(3) + t243) * t173 + (-t168 * t108 + t172 * t109 + (-t139 * t168 - t140 * t172) * qJD(4) + (t126 + t40 + t82) * qJD(3)) * t169 + m(7) * (t15 * t26 + t16 * t25 - t173 * t41 + t2 * t77 + t208 * t88 + t3 * t76) + m(6) * (-t103 * t173 - t14 * t114 - t13 * t115 + t118 * t208 + t46 * t65 + t47 * t64) + m(5) * ((-t168 * t96 + t172 * t97 - t216) * t207 + (t195 - t219 + t172 * t62 + (-t168 * t97 - t172 * t96) * qJD(4)) * t169); 0.2e1 * m(7) * (t77 * t25 + t76 * t26 - t194) + 0.2e1 * m(6) * (-t114 * t65 - t115 * t64 - t194) + 0.2e1 * m(5) * (-0.1e1 + t209) * t194; (-t13 * t181 - t132 * t14 + t46 * t98 - t47 * t99) * mrSges(6,3) + (-t162 / 0.2e1 - t33 / 0.2e1 - t32 / 0.2e1 + t95 / 0.2e1 + t94 / 0.2e1 + (t156 * t217 + Ifges(4,5)) * qJD(3)) * t173 + (-t15 * t34 + t16 * t35 + t2 * t92 - t3 * t93) * mrSges(7,3) + m(5) * (-pkin(3) * t145 + (-t206 * t97 - t219) * pkin(8)) + m(7) * (t10 * t16 + t11 * t15 + t110 * t41 + t2 * t45 + t3 * t44 + t85 * t88) + t28 * t230 + t27 * t231 + t5 * t233 + t4 * t234 + m(6) * (t103 * t158 + t104 * t14 + t105 * t13 + t118 * t198 + t46 * t68 + t47 * t67) + (qJD(4) * t113 / 0.2e1 + t148 * t207 / 0.2e1 + t80 / 0.2e1 + t186 * mrSges(5,3) + (m(5) * t186 - qJD(4) * t140 + t109) * pkin(8)) * t172 + t35 * t36 / 0.2e1 + t34 * t37 / 0.2e1 + t44 * t20 + t45 * t21 + t41 * t51 + t26 * t52 / 0.2e1 + t25 * t53 / 0.2e1 + (-pkin(8) * t108 - t63 * mrSges(5,3) + t207 * t229 + t81 / 0.2e1 + (-pkin(8) * t139 - t97 * mrSges(5,3) + pkin(4) * t82 - t112 / 0.2e1 + t218 / 0.2e1) * qJD(4)) * t168 + (t172 * t138 / 0.2e1 + t137 * t228 + t156 * t136 + (t148 * t228 + t172 * t229) * qJD(4) + (t156 * mrSges(4,2) + t222 / 0.2e1 + t220 / 0.2e1 + Ifges(6,5) * t230 + Ifges(6,6) * t231 + Ifges(7,5) * t233 + Ifges(7,6) * t234 - Ifges(4,6)) * qJD(3)) * t169 + t10 * t69 + t11 * t70 + t76 * t18 / 0.2e1 + t77 * t19 / 0.2e1 + t85 * t40 + t88 * t17 - pkin(3) * t91 - t98 * t75 / 0.2e1 - t99 * t74 / 0.2e1 + t65 * t101 / 0.2e1 + t64 * t102 / 0.2e1 + t103 * t100 + t104 * t54 + t105 * t55 + t67 * t106 + t68 * t107 + t110 * t6 - t114 * t57 / 0.2e1 - t115 * t58 / 0.2e1 + t118 * t56 + t158 * t29; (-t136 - t17 - t56) * t173 + m(7) * (t10 * t77 + t11 * t76 - t173 * t85 + t45 * t25 + t44 * t26) + m(6) * (-pkin(4) * t191 + t104 * t65 + t105 * t64 - t68 * t114 - t67 * t115) + (t25 * t92 - t26 * t93 - t34 * t76 + t35 * t77) * mrSges(7,3) + (-t114 * t98 + t115 * t99 - t132 * t65 - t181 * t64) * mrSges(6,3) + ((mrSges(5,3) * t209 - mrSges(4,2)) * t173 + m(5) * (t209 * t226 - t227) + (m(6) * t158 + m(7) * t110 + t100 + t217 + t51) * t169) * qJD(3); -0.2e1 * pkin(3) * t136 - t99 * t101 - t98 * t102 + 0.2e1 * t110 * t17 - t181 * t57 + t132 * t58 + t172 * t137 + t168 * t138 + 0.2e1 * t158 * t56 + t92 * t18 + t93 * t19 + t34 * t53 + t35 * t52 + 0.2e1 * t85 * t51 + (t172 * t148 + (0.2e1 * pkin(4) * t100 - t147) * t168) * qJD(4) + (t104 * t68 + t105 * t67 + t158 * t198) * t237 + (t10 * t45 + t11 * t44 + t110 * t85) * t236 + 0.2e1 * (t10 * t92 - t11 * t93 - t34 * t44 + t35 * t45) * mrSges(7,3) + 0.2e1 * (t104 * t98 - t105 * t99 - t132 * t68 - t181 * t67) * mrSges(6,3); t174 - Ifges(5,5) * t192 + m(7) * (t119 * t3 + t120 * t2 + t15 * t90 + t16 * t89) - t176 * Ifges(5,6) + (t106 * t202 - t107 * t203 + m(6) * (t13 * t167 + t14 * t171 + t202 * t47 - t203 * t46) + t171 * t54 + t167 * t55) * pkin(4) - t62 * mrSges(5,2) + t63 * mrSges(5,1) + t89 * t69 + t90 * t70 + t119 * t20 + t120 * t21 - t241; m(7) * (t119 * t26 + t120 * t25 + t76 * t90 + t77 * t89) + m(6) * (t167 * t64 + t171 * t65 + (t114 * t167 - t115 * t171) * qJD(5)) * pkin(4) + t243; m(7) * (t10 * t120 + t11 * t119 + t44 * t90 + t45 * t89) + t162 + (pkin(8) * t146 - t221) * qJD(4) + (-t119 * t34 + t120 * t35 + t89 * t92 - t90 * t93) * mrSges(7,3) + (m(6) * (t167 * t67 + t171 * t68 + (-t104 * t167 + t105 * t171) * qJD(5)) + (-t167 * t99 + t171 * t98 + (t132 * t167 - t171 * t181) * qJD(5)) * mrSges(6,3)) * pkin(4) + t175; (t119 * t90 + t120 * t89) * t236 - 0.2e1 * t225 + 0.2e1 * t86 + 0.2e1 * t178; (t69 * t200 + t166 * t21 - t70 * t201 + t170 * t20 + m(7) * (-t15 * t201 + t16 * t200 + t166 * t2 + t170 * t3)) * pkin(5) + t174; m(7) * (t166 * t25 + t170 * t26 + (-t166 * t76 + t170 * t77) * qJD(6)) * pkin(5) + t179; (m(7) * (t10 * t166 + t11 * t170 + (-t166 * t44 + t170 * t45) * qJD(6)) + (t166 * t35 - t170 * t34 + (t166 * t93 + t170 * t92) * qJD(6)) * mrSges(7,3)) * pkin(5) + t175; t178 + (m(7) * (-t119 * t201 + t120 * t200 + t166 * t89 + t170 * t90) - mrSges(7,2) * t200 - mrSges(7,1) * t201) * pkin(5) + t188; 0.2e1 * t128; t180; -t6; t185; t188; t128; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
