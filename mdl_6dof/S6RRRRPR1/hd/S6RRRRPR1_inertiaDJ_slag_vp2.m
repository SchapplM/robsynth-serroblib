% Calculate time derivative of joint inertia matrix for
% S6RRRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 21:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPR1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR1_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR1_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR1_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR1_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR1_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR1_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:52:59
% EndTime: 2019-03-09 21:53:10
% DurationCPUTime: 5.03s
% Computational Cost: add. (14186->389), mult. (30848->593), div. (0->0), fcn. (31236->10), ass. (0->188)
t247 = qJD(2) + qJD(3);
t141 = sin(qJ(6));
t145 = cos(qJ(6));
t192 = qJD(6) * t145;
t140 = sin(pkin(11));
t205 = cos(pkin(11));
t142 = sin(qJ(4));
t146 = cos(qJ(4));
t143 = sin(qJ(3));
t144 = sin(qJ(2));
t147 = cos(qJ(3));
t148 = cos(qJ(2));
t116 = -t143 * t144 + t147 * t148;
t92 = t247 * t116;
t117 = t143 * t148 + t147 * t144;
t93 = t247 * t117;
t166 = t142 * t93 - t146 * t92;
t87 = t116 * t146 - t117 * t142;
t50 = qJD(4) * t87 - t166;
t88 = t116 * t142 + t117 * t146;
t51 = -qJD(4) * t88 - t142 * t92 - t146 * t93;
t39 = t140 * t51 + t205 * t50;
t68 = t140 * t87 + t205 * t88;
t163 = t141 * t39 + t68 * t192;
t235 = -pkin(8) - pkin(7);
t127 = t235 * t144;
t250 = t148 * t235;
t251 = -t143 * t127 + t147 * t250;
t95 = t127 * t147 + t143 * t250;
t194 = qJD(4) * t146;
t195 = qJD(4) * t142;
t249 = -mrSges(5,1) * t195 - mrSges(5,2) * t194;
t159 = -pkin(9) * t117 + t95;
t172 = pkin(9) * t116 - t251;
t60 = -t142 * t172 + t146 * t159;
t155 = -t88 * qJ(5) + t60;
t61 = t142 * t159 + t146 * t172;
t46 = qJ(5) * t87 + t61;
t25 = t140 * t155 + t205 * t46;
t67 = t140 * t88 - t205 * t87;
t135 = -pkin(2) * t148 - pkin(1);
t101 = -t116 * pkin(3) + t135;
t73 = -t87 * pkin(4) + t101;
t40 = t67 * pkin(5) - t68 * pkin(10) + t73;
t16 = -t141 * t25 + t145 * t40;
t17 = t141 * t40 + t145 * t25;
t209 = t145 * t17;
t248 = -t141 * t16 + t209;
t181 = t16 * t192;
t193 = qJD(6) * t141;
t38 = t140 * t50 - t205 * t51;
t82 = qJD(2) * t144 * pkin(2) + pkin(3) * t93;
t44 = -pkin(4) * t51 + t82;
t13 = pkin(5) * t38 - pkin(10) * t39 + t44;
t204 = qJD(6) * t17;
t189 = t147 * t235;
t190 = t143 * t235;
t150 = ((-t142 * t190 + t146 * t189) * t148 + (-t142 * t189 - t146 * t190) * t144) * qJD(2);
t158 = qJD(4) * t159;
t164 = qJD(4) * t172;
t151 = -t50 * qJ(5) - t88 * qJD(5) - t142 * t158 - t146 * t164;
t232 = t93 * pkin(9);
t70 = t247 * t95;
t153 = t70 - t232;
t233 = t92 * pkin(9);
t71 = t247 * t251;
t154 = t71 - t233;
t21 = (t153 + t158) * t146 + (t154 - t164) * t142;
t19 = t51 * qJ(5) + t87 * qJD(5) + t21;
t196 = qJD(3) * t147;
t197 = qJD(3) * t143;
t7 = t205 * t19 + (-t142 * (t127 * t196 + t197 * t250 - t232) + t146 * (-t127 * t197 + t196 * t250 - t233) + t151 + t150) * t140;
t3 = t13 * t145 - t141 * t7 - t204;
t227 = t3 * t141;
t246 = -t17 * t193 - t181 - t227;
t184 = t68 * t193;
t208 = t145 * t39;
t162 = t184 - t208;
t14 = mrSges(7,1) * t38 + mrSges(7,3) * t162;
t15 = -mrSges(7,2) * t38 - mrSges(7,3) * t163;
t211 = t141 * t68;
t42 = -mrSges(7,2) * t67 - mrSges(7,3) * t211;
t207 = t145 * t68;
t43 = mrSges(7,1) * t67 - mrSges(7,3) * t207;
t245 = -t141 * t14 + t145 * t15 - t43 * t192 - t42 * t193;
t174 = t205 * t142;
t217 = pkin(3) * qJD(4);
t106 = (t140 * t146 + t174) * t217;
t102 = t106 * mrSges(6,1);
t133 = pkin(3) * t146 + pkin(4);
t201 = t140 * t142;
t109 = -pkin(3) * t201 + t133 * t205;
t104 = -pkin(5) - t109;
t171 = mrSges(7,1) * t141 + mrSges(7,2) * t145;
t120 = t171 * qJD(6);
t91 = t104 * t120;
t123 = -mrSges(7,1) * t145 + mrSges(7,2) * t141;
t94 = t106 * t123;
t107 = (t146 * t205 - t201) * t217;
t138 = t141 ^ 2;
t222 = mrSges(7,3) * t138;
t97 = t107 * t222;
t139 = t145 ^ 2;
t221 = mrSges(7,3) * t139;
t98 = t107 * t221;
t244 = t91 + t94 + t97 + t98 - t102;
t243 = 2 * m(5);
t242 = 2 * m(6);
t241 = 2 * m(7);
t6 = t140 * t19 - t205 * (-t142 * t153 + t146 * t154 + t151);
t240 = 0.2e1 * t6;
t24 = t140 * t46 - t155 * t205;
t239 = 0.2e1 * t24;
t238 = 0.2e1 * t44;
t237 = 0.2e1 * t135;
t236 = m(6) * pkin(4);
t234 = t24 * t6;
t230 = pkin(4) * t140;
t2 = qJD(6) * t16 + t13 * t141 + t145 * t7;
t229 = t145 * t2;
t134 = pkin(2) * t147 + pkin(3);
t200 = t142 * t143;
t85 = t134 * t194 + (-t143 * t195 + (t146 * t147 - t200) * qJD(3)) * pkin(2);
t199 = t143 * t146;
t86 = -t134 * t195 + (-t143 * t194 + (-t142 * t147 - t199) * qJD(3)) * pkin(2);
t65 = t140 * t85 - t205 * t86;
t228 = t24 * t65;
t226 = t38 * mrSges(6,3);
t66 = t140 * t86 + t205 * t85;
t225 = t66 * mrSges(6,2);
t224 = t85 * mrSges(5,2);
t223 = Ifges(7,5) * t208 + Ifges(7,3) * t38;
t220 = Ifges(7,4) * t141;
t219 = Ifges(7,4) * t145;
t218 = Ifges(7,6) * t141;
t216 = t106 * t24;
t215 = t107 * mrSges(6,2);
t111 = -pkin(2) * t200 + t146 * t134;
t108 = pkin(4) + t111;
t112 = pkin(2) * t199 + t134 * t142;
t80 = t140 * t108 + t205 * t112;
t76 = pkin(10) + t80;
t206 = t145 * t76;
t203 = qJD(6) * t76;
t110 = pkin(3) * t174 + t140 * t133;
t198 = t138 + t139;
t191 = 0.2e1 * t148;
t179 = t205 * pkin(4);
t178 = t38 * mrSges(6,1) + t39 * mrSges(6,2);
t177 = -t193 / 0.2e1;
t176 = -(2 * Ifges(6,4)) - t218;
t175 = t198 * t66;
t173 = t198 * t107;
t170 = Ifges(7,1) * t145 - t220;
t169 = -Ifges(7,2) * t141 + t219;
t168 = Ifges(7,5) * t141 + Ifges(7,6) * t145;
t167 = -t141 * t43 + t145 * t42;
t79 = t108 * t205 - t140 * t112;
t121 = t169 * qJD(6);
t122 = t170 * qJD(6);
t124 = Ifges(7,2) * t145 + t220;
t125 = Ifges(7,1) * t141 + t219;
t161 = t145 * t121 + t141 * t122 - t124 * t193 + t125 * t192;
t160 = (-mrSges(4,1) * t143 - mrSges(4,2) * t147) * qJD(3) * pkin(2);
t157 = -t227 + (-t141 * t17 - t145 * t16) * qJD(6);
t56 = t65 * t123;
t62 = t66 * t222;
t63 = t66 * t221;
t64 = t65 * mrSges(6,1);
t75 = -pkin(5) - t79;
t72 = t75 * t120;
t83 = t86 * mrSges(5,1);
t156 = t161 + t56 + t62 + t63 - t64 + t72 + t83 - t224;
t10 = -Ifges(7,4) * t162 - Ifges(7,2) * t163 + t38 * Ifges(7,6);
t11 = -Ifges(7,1) * t162 - Ifges(7,4) * t163 + t38 * Ifges(7,5);
t136 = Ifges(7,5) * t192;
t22 = t166 * pkin(9) - t61 * qJD(4) + (-t142 * t95 + t146 * t251) * qJD(3) + t150;
t28 = t67 * Ifges(7,6) + t169 * t68;
t29 = t67 * Ifges(7,5) + t170 * t68;
t152 = -t21 * mrSges(5,2) - t7 * mrSges(6,2) + mrSges(7,3) * t229 + t22 * mrSges(5,1) + t24 * t120 + t28 * t177 + t29 * t192 / 0.2e1 + Ifges(6,5) * t39 + Ifges(5,6) * t51 + Ifges(5,5) * t50 - t121 * t211 / 0.2e1 + t122 * t207 / 0.2e1 + t67 * (-Ifges(7,6) * t193 + t136) / 0.2e1 + t141 * t11 / 0.2e1 + t145 * t10 / 0.2e1 + (-mrSges(6,1) + t123) * t6 + (t168 / 0.2e1 - Ifges(6,6)) * t38 - t163 * t124 / 0.2e1 + (t208 / 0.2e1 + t68 * t177) * t125;
t149 = t71 * mrSges(4,1) - t70 * mrSges(4,2) + Ifges(4,5) * t92 - Ifges(4,6) * t93 + t152;
t132 = -t179 - pkin(5);
t131 = pkin(10) + t230;
t105 = pkin(10) + t110;
t103 = t132 * t120;
t41 = t171 * t68;
t12 = mrSges(7,1) * t163 - mrSges(7,2) * t162;
t1 = [0.2e1 * t51 * Ifges(5,2) * t87 - 0.2e1 * t116 * Ifges(4,2) * t93 + 0.2e1 * t50 * t88 * Ifges(5,1) + 0.2e1 * t92 * t117 * Ifges(4,1) + t12 * t239 + t41 * t240 + (t16 * t3 + t17 * t2 + t234) * t241 + (t25 * t7 + t44 * t73 + t234) * t242 + (t101 * t82 + t21 * t61 + t22 * t60) * t243 - 0.2e1 * t25 * t226 + 0.2e1 * (t50 * t87 + t51 * t88) * Ifges(5,4) + 0.2e1 * (t21 * t87 - t22 * t88 - t50 * t60 + t51 * t61) * mrSges(5,3) + 0.2e1 * t73 * t178 + (mrSges(4,1) * t93 + mrSges(4,2) * t92) * t237 + 0.2e1 * (t116 * t92 - t117 * t93) * Ifges(4,4) + 0.2e1 * m(4) * (-t251 * t70 + t71 * t95) + 0.2e1 * (t116 * t70 - t117 * t71 + t251 * t93 - t92 * t95) * mrSges(4,3) + 0.2e1 * t101 * (-mrSges(5,1) * t51 + mrSges(5,2) * t50) + 0.2e1 * t82 * (-mrSges(5,1) * t87 + mrSges(5,2) * t88) + 0.2e1 * t2 * t42 + 0.2e1 * t3 * t43 + 0.2e1 * t16 * t14 + 0.2e1 * t17 * t15 + ((-mrSges(3,2) * pkin(1) + Ifges(3,4) * t148) * t191 + (0.2e1 * pkin(2) * (-mrSges(4,1) * t116 + mrSges(4,2) * t117) - 0.2e1 * pkin(1) * mrSges(3,1) + m(4) * pkin(2) * t237 - 0.2e1 * Ifges(3,4) * t144 + (Ifges(3,1) - Ifges(3,2)) * t191) * t144) * qJD(2) + (mrSges(6,1) * t238 - 0.2e1 * mrSges(6,3) * t7 + t176 * t39 + ((2 * Ifges(6,2)) + Ifges(7,3)) * t38 + t223) * t67 + (mrSges(6,3) * t239 - t141 * t28 + t145 * t29) * t39 + (mrSges(6,2) * t238 + mrSges(6,3) * t240 + 0.2e1 * Ifges(6,1) * t39 - t141 * t10 + t145 * t11 + (Ifges(7,5) * t145 + t176) * t38 + (-t141 * t29 - t145 * t28 - t168 * t67) * qJD(6)) * t68; t149 + (Ifges(3,5) * t148 - Ifges(3,6) * t144 + (-mrSges(3,1) * t148 + mrSges(3,2) * t144) * pkin(7)) * qJD(2) + (t76 * t15 + t66 * t42 + (-mrSges(7,3) * t16 - t43 * t76) * qJD(6)) * t145 + m(7) * (-t76 * t181 + t2 * t206 + t66 * t209 + t6 * t75 + t228) + m(6) * (t25 * t66 - t6 * t79 + t7 * t80 + t228) + m(5) * (t111 * t22 + t112 * t21 + t60 * t86 + t61 * t85) + t75 * t12 + t65 * t41 + (m(4) * (t143 * t70 + t147 * t71 + (-t143 * t95 - t147 * t251) * qJD(3)) + (-t143 * t93 - t147 * t92 + (t116 * t147 + t117 * t143) * qJD(3)) * mrSges(4,3)) * pkin(2) + (-t42 * t203 - t66 * t43 - t76 * t14 + m(7) * (-t16 * t66 - t17 * t203 - t3 * t76) + (-t3 - t204) * mrSges(7,3)) * t141 + (-t111 * t50 + t112 * t51 + t85 * t87 - t86 * t88) * mrSges(5,3) + (-t38 * t80 - t39 * t79 + t65 * t68 - t66 * t67) * mrSges(6,3); -0.2e1 * t224 - 0.2e1 * t225 + 0.2e1 * t56 + 0.2e1 * t62 + 0.2e1 * t63 - 0.2e1 * t64 + 0.2e1 * t72 + 0.2e1 * t83 + 0.2e1 * t160 + (t175 * t76 + t65 * t75) * t241 + (-t65 * t79 + t66 * t80) * t242 + (t111 * t86 + t112 * t85) * t243 + t161; t149 + m(7) * (t104 * t6 + t216) + m(6) * (-t109 * t6 + t110 * t7 + t216) + t104 * t12 + t106 * t41 + (m(5) * (t142 * t21 + t146 * t22 + (-t142 * t60 + t146 * t61) * qJD(4)) + (t142 * t51 - t146 * t50 + (t142 * t88 + t146 * t87) * qJD(4)) * mrSges(5,3)) * pkin(3) + (m(7) * (t229 + t246) + t245) * t105 + t157 * mrSges(7,3) + (t106 * t68 - t109 * t39 - t110 * t38) * mrSges(6,3) + (m(6) * t25 + m(7) * t248 - t67 * mrSges(6,3) + t167) * t107; m(7) * (t104 * t65 + t105 * t175 + t106 * t75 + t173 * t76) + t160 + (m(5) * (-t111 * t195 + t112 * t194 + t142 * t85 + t146 * t86) + t249) * pkin(3) + t156 + (-t107 - t66) * mrSges(6,2) + m(6) * (-t106 * t79 + t107 * t80 - t109 * t65 + t110 * t66) + t244; -0.2e1 * t215 - 0.2e1 * t102 + 0.2e1 * t91 + 0.2e1 * t94 + 0.2e1 * t97 + 0.2e1 * t98 + 0.2e1 * (-mrSges(5,1) * t142 - mrSges(5,2) * t146) * t217 + (t104 * t106 + t105 * t173) * t241 + (-t106 * t109 + t107 * t110) * t242 + t161; -t226 * t230 - t39 * mrSges(6,3) * t179 + t152 + (t140 * t7 - t205 * t6) * t236 + (m(7) * t6 + t12) * t132 + t246 * mrSges(7,3) + (m(7) * (t157 + t229) + t245) * t131; t103 - t225 + (t140 * t66 - t205 * t65) * t236 + m(7) * (t131 * t175 + t132 * t65) + t156; t103 + m(7) * (t132 * t106 + t131 * t173) + (-t106 * t205 + t107 * t140) * t236 - t215 + t161 + t249 * pkin(3) + t244; 0.2e1 * t103 + t161; t145 * t14 + t141 * t15 + t167 * qJD(6) + m(7) * (qJD(6) * t248 + t141 * t2 + t145 * t3) + m(6) * t44 + t178; 0; 0; 0; 0; mrSges(7,1) * t3 - mrSges(7,2) * t2 - Ifges(7,5) * t184 - Ifges(7,6) * t163 + t223; t136 - t171 * t66 + (-mrSges(7,1) * t206 + (mrSges(7,2) * t76 - Ifges(7,6)) * t141) * qJD(6); t136 - t171 * t107 + (t105 * t123 - t218) * qJD(6); t136 + (t123 * t131 - t218) * qJD(6); -t120; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
