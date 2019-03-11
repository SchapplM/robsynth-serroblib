% Calculate time derivative of joint inertia matrix for
% S6RPRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
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
% Datum: 2019-03-09 07:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRR9_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR9_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR9_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR9_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR9_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR9_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR9_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:23:05
% EndTime: 2019-03-09 07:23:15
% DurationCPUTime: 4.51s
% Computational Cost: add. (7399->526), mult. (16537->785), div. (0->0), fcn. (14662->8), ass. (0->207)
t239 = 2 * qJD(2);
t171 = sin(qJ(5));
t172 = sin(qJ(4));
t175 = cos(qJ(5));
t176 = cos(qJ(4));
t187 = t171 * t172 - t175 * t176;
t243 = qJD(4) + qJD(5);
t245 = t243 * t187;
t173 = sin(qJ(3));
t215 = qJD(3) * t173;
t198 = t172 * t215;
t177 = cos(qJ(3));
t210 = qJD(4) * t177;
t199 = t176 * t210;
t182 = t198 - t199;
t236 = -pkin(9) - pkin(8);
t155 = t236 * t172;
t156 = t236 * t176;
t108 = t171 * t155 - t175 * t156;
t214 = qJD(3) * t177;
t244 = Ifges(5,6) * t198 + Ifges(5,3) * t214;
t231 = pkin(8) * t177;
t232 = pkin(3) * t173;
t150 = qJ(2) - t231 + t232;
t178 = -pkin(1) - pkin(7);
t218 = t173 * t178;
t158 = t176 * t218;
t115 = t172 * t150 + t158;
t216 = t172 ^ 2 + t176 ^ 2;
t139 = t171 * t176 + t172 * t175;
t102 = t243 * t139;
t65 = -t102 * t177 + t187 * t215;
t67 = t139 * t215 + t177 * t245;
t242 = Ifges(6,5) * t65 + Ifges(6,6) * t67 + Ifges(6,3) * t214;
t241 = 2 * m(6);
t240 = 2 * m(7);
t170 = sin(qJ(6));
t174 = cos(qJ(6));
t97 = -t139 * t170 - t174 * t187;
t238 = t97 / 0.2e1;
t98 = t139 * t174 - t170 * t187;
t237 = t98 / 0.2e1;
t235 = -t187 / 0.2e1;
t234 = t139 / 0.2e1;
t233 = -t172 / 0.2e1;
t162 = pkin(4) * t175 + pkin(5);
t206 = qJD(6) * t174;
t207 = qJD(6) * t170;
t221 = t170 * t171;
t91 = t162 * t206 + (-t171 * t207 + (t174 * t175 - t221) * qJD(5)) * pkin(4);
t230 = t91 * mrSges(7,2);
t219 = t172 * t177;
t103 = -pkin(9) * t219 + t115;
t135 = t176 * t150;
t193 = -t172 * t178 + pkin(4);
t217 = t176 * t177;
t94 = -pkin(9) * t217 + t173 * t193 + t135;
t51 = t175 * t103 + t171 * t94;
t229 = Ifges(5,4) * t172;
t228 = Ifges(5,4) * t176;
t227 = Ifges(5,5) * t172;
t226 = Ifges(5,6) * t172;
t225 = Ifges(5,6) * t176;
t136 = qJD(2) + (pkin(3) * t177 + pkin(8) * t173) * qJD(3);
t126 = t176 * t136;
t213 = qJD(3) * t178;
t201 = t177 * t213;
t72 = -qJD(4) * t115 - t172 * t201 + t126;
t224 = t172 * t72;
t223 = t173 * Ifges(5,6);
t152 = -mrSges(5,1) * t176 + mrSges(5,2) * t172;
t222 = -mrSges(4,1) + t152;
t220 = t171 * t174;
t212 = qJD(4) * t172;
t211 = qJD(4) * t176;
t209 = qJD(5) * t171;
t208 = qJD(5) * t175;
t120 = t139 * t177;
t122 = t187 * t177;
t78 = -t120 * t174 + t122 * t170;
t26 = qJD(6) * t78 + t170 * t67 + t174 * t65;
t80 = -t120 * t170 - t122 * t174;
t28 = -qJD(6) * t80 - t170 * t65 + t174 * t67;
t205 = Ifges(7,5) * t26 + Ifges(7,6) * t28 + Ifges(7,3) * t214;
t204 = pkin(4) * t212;
t163 = -pkin(4) * t176 - pkin(3);
t203 = qJD(4) * t236;
t202 = t173 * t214;
t160 = t173 * t213;
t200 = t172 * t210;
t197 = t173 * t212;
t64 = -qJD(3) * t122 - t102 * t173;
t66 = -qJD(3) * t120 + t173 * t245;
t119 = t139 * t173;
t121 = t187 * t173;
t77 = -t119 * t174 + t121 * t170;
t25 = qJD(6) * t77 + t170 * t66 + t174 * t64;
t79 = -t119 * t170 - t121 * t174;
t27 = -qJD(6) * t79 - t170 * t64 + t174 * t66;
t196 = t27 * mrSges(7,1) - t25 * mrSges(7,2);
t92 = -t162 * t207 + (-t171 * t206 + (-t170 * t175 - t220) * qJD(5)) * pkin(4);
t88 = t92 * mrSges(7,1);
t195 = t88 - t230;
t194 = -Ifges(5,5) * t176 + (2 * Ifges(4,4));
t50 = -t103 * t171 + t175 * t94;
t114 = -t172 * t218 + t135;
t71 = t172 * t136 + t150 * t211 + t176 * t201 - t178 * t197;
t192 = -qJD(4) * t114 + t71;
t107 = t175 * t155 + t156 * t171;
t137 = pkin(4) * t219 - t177 * t178;
t148 = t172 * t203;
t149 = t176 * t203;
t69 = t175 * t148 + t171 * t149 + t155 * t208 + t156 * t209;
t39 = -pkin(10) * t102 + t69;
t70 = -qJD(5) * t108 - t148 * t171 + t175 * t149;
t40 = pkin(10) * t245 + t70;
t84 = -pkin(10) * t139 + t107;
t85 = -pkin(10) * t187 + t108;
t43 = -t170 * t85 + t174 * t84;
t10 = qJD(6) * t43 + t170 * t40 + t174 * t39;
t44 = t170 * t84 + t174 * t85;
t11 = -qJD(6) * t44 - t170 * t39 + t174 * t40;
t35 = -qJD(6) * t98 - t102 * t174 + t170 * t245;
t32 = Ifges(7,6) * t35;
t34 = qJD(6) * t97 - t102 * t170 - t174 * t245;
t33 = Ifges(7,5) * t34;
t191 = t11 * mrSges(7,1) - t10 * mrSges(7,2) + t32 + t33;
t190 = mrSges(5,1) * t172 + mrSges(5,2) * t176;
t189 = Ifges(5,1) * t176 - t229;
t154 = Ifges(5,1) * t172 + t228;
t188 = -Ifges(5,2) * t172 + t228;
t153 = Ifges(5,2) * t176 + t229;
t38 = pkin(5) * t173 + pkin(10) * t122 + t50;
t42 = -pkin(10) * t120 + t51;
t18 = -t170 * t42 + t174 * t38;
t19 = t170 * t38 + t174 * t42;
t48 = t126 + (-t158 + (pkin(9) * t177 - t150) * t172) * qJD(4) + (pkin(9) * t173 * t176 + t177 * t193) * qJD(3);
t57 = pkin(9) * t182 + t71;
t14 = -qJD(5) * t51 - t171 * t57 + t175 * t48;
t7 = pkin(5) * t214 - pkin(10) * t65 + t14;
t13 = -t103 * t209 + t171 * t48 + t175 * t57 + t94 * t208;
t8 = pkin(10) * t67 + t13;
t2 = qJD(6) * t18 + t170 * t7 + t174 * t8;
t3 = -qJD(6) * t19 - t170 * t8 + t174 * t7;
t186 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + t205;
t185 = t66 * mrSges(6,1) - t64 * mrSges(6,2) + t196;
t111 = -pkin(4) * t182 + t160;
t184 = (-mrSges(6,1) * t171 - mrSges(6,2) * t175) * qJD(5) * pkin(4);
t183 = t176 * t215 + t200;
t100 = Ifges(6,5) * t245;
t99 = Ifges(6,6) * t102;
t181 = t70 * mrSges(6,1) - t69 * mrSges(6,2) - t100 + t191 - t99;
t180 = t14 * mrSges(6,1) - t13 * mrSges(6,2) + t186 + t242;
t167 = Ifges(5,5) * t211;
t147 = mrSges(5,1) * t173 - mrSges(5,3) * t217;
t146 = -mrSges(5,2) * t173 - mrSges(5,3) * t219;
t145 = t189 * qJD(4);
t144 = t188 * qJD(4);
t143 = t190 * qJD(4);
t132 = (-t170 * mrSges(7,1) - t174 * mrSges(7,2)) * qJD(6) * pkin(5);
t131 = t190 * t177;
t124 = pkin(4) * t220 + t162 * t170;
t123 = -pkin(4) * t221 + t162 * t174;
t118 = Ifges(5,5) * t173 + t177 * t189;
t117 = t177 * t188 + t223;
t116 = pkin(5) * t187 + t163;
t113 = -mrSges(5,2) * t214 + mrSges(5,3) * t182;
t112 = mrSges(5,1) * t214 + mrSges(5,3) * t183;
t110 = mrSges(6,1) * t173 + mrSges(6,3) * t122;
t109 = -mrSges(6,2) * t173 - mrSges(6,3) * t120;
t106 = Ifges(6,1) * t139 - Ifges(6,4) * t187;
t105 = Ifges(6,4) * t139 - Ifges(6,2) * t187;
t104 = mrSges(6,1) * t187 + mrSges(6,2) * t139;
t96 = t120 * pkin(5) + t137;
t93 = -mrSges(5,1) * t182 - mrSges(5,2) * t183;
t87 = pkin(5) * t102 + t204;
t83 = mrSges(6,1) * t120 - mrSges(6,2) * t122;
t82 = -t154 * t210 + (Ifges(5,5) * t177 - t173 * t189) * qJD(3);
t81 = -t153 * t210 + (Ifges(5,6) * t177 - t173 * t188) * qJD(3);
t76 = -Ifges(6,1) * t122 - Ifges(6,4) * t120 + Ifges(6,5) * t173;
t75 = -Ifges(6,4) * t122 - Ifges(6,2) * t120 + Ifges(6,6) * t173;
t74 = mrSges(7,1) * t173 - mrSges(7,3) * t80;
t73 = -mrSges(7,2) * t173 + mrSges(7,3) * t78;
t60 = -Ifges(6,1) * t245 - Ifges(6,4) * t102;
t59 = -Ifges(6,4) * t245 - Ifges(6,2) * t102;
t58 = mrSges(6,1) * t102 - mrSges(6,2) * t245;
t56 = -mrSges(6,2) * t214 + mrSges(6,3) * t67;
t55 = mrSges(6,1) * t214 - mrSges(6,3) * t65;
t54 = Ifges(7,1) * t98 + Ifges(7,4) * t97;
t53 = Ifges(7,4) * t98 + Ifges(7,2) * t97;
t52 = -mrSges(7,1) * t97 + mrSges(7,2) * t98;
t45 = -pkin(5) * t67 + t111;
t41 = -mrSges(7,1) * t78 + mrSges(7,2) * t80;
t37 = Ifges(7,1) * t80 + Ifges(7,4) * t78 + Ifges(7,5) * t173;
t36 = Ifges(7,4) * t80 + Ifges(7,2) * t78 + Ifges(7,6) * t173;
t31 = -mrSges(6,1) * t67 + mrSges(6,2) * t65;
t30 = Ifges(6,1) * t65 + Ifges(6,4) * t67 + Ifges(6,5) * t214;
t29 = Ifges(6,4) * t65 + Ifges(6,2) * t67 + Ifges(6,6) * t214;
t21 = -mrSges(7,2) * t214 + mrSges(7,3) * t28;
t20 = mrSges(7,1) * t214 - mrSges(7,3) * t26;
t17 = Ifges(7,1) * t34 + Ifges(7,4) * t35;
t16 = Ifges(7,4) * t34 + Ifges(7,2) * t35;
t15 = -mrSges(7,1) * t35 + mrSges(7,2) * t34;
t6 = -mrSges(7,1) * t28 + mrSges(7,2) * t26;
t5 = Ifges(7,1) * t26 + Ifges(7,4) * t28 + Ifges(7,5) * t214;
t4 = Ifges(7,4) * t26 + Ifges(7,2) * t28 + Ifges(7,6) * t214;
t1 = [(t18 * t3 + t19 * t2 + t45 * t96) * t240 + (t111 * t137 + t13 * t51 + t14 * t50) * t241 + 0.2e1 * m(5) * (t114 * t72 + t115 * t71) + (mrSges(4,1) * t239 + (-0.2e1 * qJ(2) * mrSges(4,2) + t172 * t117 - t176 * t118 + 0.2e1 * t178 * t131 + t173 * t194) * qJD(3) + t205 + t242 + t244) * t173 + (mrSges(3,3) + (m(3) + m(4)) * qJ(2)) * t239 + 0.2e1 * t18 * t20 + 0.2e1 * t19 * t21 + t28 * t36 + t26 * t37 + 0.2e1 * t45 * t41 + 0.2e1 * t50 * t55 + 0.2e1 * t51 * t56 + 0.2e1 * t2 * t73 + 0.2e1 * t3 * t74 + t67 * t75 + t65 * t76 + t78 * t4 + t80 * t5 + 0.2e1 * t96 * t6 + 0.2e1 * t13 * t109 + 0.2e1 * t14 * t110 + 0.2e1 * t111 * t83 + 0.2e1 * t114 * t112 + 0.2e1 * t115 * t113 - t120 * t29 - t122 * t30 + 0.2e1 * t137 * t31 + 0.2e1 * t71 * t146 + 0.2e1 * t72 * t147 + (mrSges(4,2) * t239 - t172 * t81 + t176 * t82 - 0.2e1 * t178 * t93 + (-t176 * t117 - t172 * t118 + t173 * (-t225 - t227)) * qJD(4) + (-Ifges(6,5) * t122 - Ifges(6,6) * t120 + Ifges(7,5) * t80 + Ifges(7,6) * t78 + 0.2e1 * qJ(2) * mrSges(4,1) + (-t194 - t226) * t177 + (-0.2e1 * m(5) * t178 ^ 2 - (2 * Ifges(4,1)) + (2 * Ifges(4,2)) + Ifges(5,3) + Ifges(6,3) + Ifges(7,3)) * t173) * qJD(3)) * t177; t64 * t109 + t66 * t110 - t119 * t55 - t121 * t56 + t77 * t20 + t79 * t21 + t25 * t73 + t27 * t74 + (-t31 - t6 - t93 + (t146 * t176 - t147 * t172) * qJD(3)) * t177 + (-t172 * t112 + t176 * t113 + (-t146 * t172 - t147 * t176) * qJD(4) + (t131 + t41 + t83) * qJD(3)) * t173 + m(7) * (-t177 * t45 + t18 * t27 + t19 * t25 + t2 * t79 + t215 * t96 + t3 * t77) + m(6) * (-t111 * t177 - t119 * t14 - t121 * t13 + t137 * t215 + t50 * t66 + t51 * t64) + m(5) * ((-t114 * t172 + t115 * t176) * t214 + (-0.2e1 * t201 - t224 + t176 * t71 + (-t114 * t176 - t115 * t172) * qJD(4)) * t173); 0.2e1 * m(7) * (t25 * t79 + t27 * t77 - t202) + 0.2e1 * m(6) * (-t119 * t66 - t121 * t64 - t202) + 0.2e1 * m(5) * (-0.1e1 + t216) * t202; -t245 * t76 / 0.2e1 + (t167 / 0.2e1 + t33 / 0.2e1 + t32 / 0.2e1 - t100 / 0.2e1 - t99 / 0.2e1 + (t178 * t222 - Ifges(4,5)) * qJD(3)) * t173 + (-t18 * t34 + t19 * t35 + t2 * t97 - t3 * t98) * mrSges(7,3) + t4 * t238 + t30 * t234 + t29 * t235 + t5 * t237 + (t176 * t145 / 0.2e1 + t144 * t233 - t178 * t143 + (-t176 * t153 / 0.2e1 + t154 * t233) * qJD(4) + (-t178 * mrSges(4,2) + t227 / 0.2e1 + t225 / 0.2e1 + Ifges(6,5) * t234 + Ifges(6,6) * t235 + Ifges(7,5) * t237 + Ifges(7,6) * t238 - Ifges(4,6)) * qJD(3)) * t177 + (t153 * t215 / 0.2e1 - pkin(8) * t112 - t72 * mrSges(5,3) + t82 / 0.2e1 + (-pkin(8) * t146 - t115 * mrSges(5,3) + pkin(4) * t83 - t117 / 0.2e1 - t223 / 0.2e1) * qJD(4)) * t172 + (-t154 * t215 / 0.2e1 + qJD(4) * t118 / 0.2e1 + t81 / 0.2e1 + t192 * mrSges(5,3) + (m(5) * t192 - qJD(4) * t147 + t113) * pkin(8)) * t176 + m(6) * (t107 * t14 + t108 * t13 + t111 * t163 + t137 * t204 + t50 * t70 + t51 * t69) + m(7) * (t10 * t19 + t11 * t18 + t116 * t45 + t2 * t44 + t3 * t43 + t87 * t96) + (-t102 * t51 - t13 * t187 - t139 * t14 + t245 * t50) * mrSges(6,3) + t35 * t36 / 0.2e1 + t34 * t37 / 0.2e1 + t43 * t20 + t44 * t21 + t45 * t52 + t28 * t53 / 0.2e1 + t26 * t54 / 0.2e1 + t10 * t73 + t11 * t74 + t78 * t16 / 0.2e1 + t80 * t17 / 0.2e1 + t87 * t41 - pkin(3) * t93 + t96 * t15 - t102 * t75 / 0.2e1 + t67 * t105 / 0.2e1 + t65 * t106 / 0.2e1 + t107 * t55 + t108 * t56 + t69 * t109 + t70 * t110 + t111 * t104 + t116 * t6 - t120 * t59 / 0.2e1 - t122 * t60 / 0.2e1 + t137 * t58 + t163 * t31 + m(5) * (-pkin(3) * t160 + (-t115 * t212 - t224) * pkin(8)); (-t143 - t15 - t58) * t177 + m(7) * (t10 * t79 + t11 * t77 - t177 * t87 + t25 * t44 + t27 * t43) + m(6) * (-pkin(4) * t200 + t107 * t66 + t108 * t64 - t119 * t70 - t121 * t69) + (t25 * t97 - t27 * t98 - t34 * t77 + t35 * t79) * mrSges(7,3) + (t102 * t121 - t119 * t245 - t139 * t66 - t187 * t64) * mrSges(6,3) + ((mrSges(5,3) * t216 - mrSges(4,2)) * t177 + m(5) * (t216 * t231 - t232) + (m(6) * t163 + m(7) * t116 + t104 + t222 + t52) * t173) * qJD(3); -0.2e1 * pkin(3) * t143 - t245 * t106 - t102 * t105 + 0.2e1 * t116 * t15 - t187 * t59 + t139 * t60 + t176 * t144 + t172 * t145 + t97 * t16 + 0.2e1 * t163 * t58 + t98 * t17 + t34 * t54 + t35 * t53 + 0.2e1 * t87 * t52 + (t176 * t154 + (0.2e1 * pkin(4) * t104 - t153) * t172) * qJD(4) + (t107 * t70 + t108 * t69 + t163 * t204) * t241 + (t10 * t44 + t11 * t43 + t116 * t87) * t240 + 0.2e1 * (t10 * t97 - t11 * t98 - t34 * t43 + t35 * t44) * mrSges(7,3) + 0.2e1 * (-t102 * t108 + t107 * t245 - t139 * t70 - t187 * t69) * mrSges(6,3); -Ifges(5,6) * t199 - t183 * Ifges(5,5) + m(7) * (t123 * t3 + t124 * t2 + t18 * t92 + t19 * t91) + (t109 * t208 - t110 * t209 + m(6) * (t13 * t171 + t14 * t175 + t208 * t51 - t209 * t50) + t175 * t55 + t171 * t56) * pkin(4) + t180 - t71 * mrSges(5,2) + t72 * mrSges(5,1) + t91 * t73 + t92 * t74 + t123 * t20 + t124 * t21 + t244; m(7) * (t123 * t27 + t124 * t25 + t77 * t92 + t79 * t91) + (-t176 * t214 + t197) * mrSges(5,2) + (-t172 * t214 - t173 * t211) * mrSges(5,1) + m(6) * (t171 * t64 + t175 * t66 + (t119 * t171 - t121 * t175) * qJD(5)) * pkin(4) + t185; m(7) * (t10 * t124 + t11 * t123 + t43 * t92 + t44 * t91) + t167 + (pkin(8) * t152 - t226) * qJD(4) + (-t123 * t34 + t124 * t35 + t91 * t97 - t92 * t98) * mrSges(7,3) + (m(6) * (t171 * t69 + t175 * t70 + (-t107 * t171 + t108 * t175) * qJD(5)) + (t175 * t245 - t171 * t102 + (t139 * t171 - t175 * t187) * qJD(5)) * mrSges(6,3)) * pkin(4) + t181; (t123 * t92 + t124 * t91) * t240 - 0.2e1 * t230 + 0.2e1 * t88 + 0.2e1 * t184; (-t74 * t207 + t174 * t20 + m(7) * (t170 * t2 + t174 * t3 - t18 * t207 + t19 * t206) + t73 * t206 + t170 * t21) * pkin(5) + t180; m(7) * (t170 * t25 + t174 * t27 + (-t170 * t77 + t174 * t79) * qJD(6)) * pkin(5) + t185; (m(7) * (t10 * t170 + t11 * t174 + (-t170 * t43 + t174 * t44) * qJD(6)) + (t170 * t35 - t174 * t34 + (t170 * t98 + t174 * t97) * qJD(6)) * mrSges(7,3)) * pkin(5) + t181; t184 + (m(7) * (-t123 * t207 + t124 * t206 + t170 * t91 + t174 * t92) - mrSges(7,2) * t206 - mrSges(7,1) * t207) * pkin(5) + t195; 0.2e1 * t132; t186; t196; t191; t195; t132; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
