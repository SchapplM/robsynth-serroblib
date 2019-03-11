% Calculate time derivative of joint inertia matrix for
% S6RRPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 11:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRP3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP3_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP3_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP3_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP3_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP3_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP3_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:47:23
% EndTime: 2019-03-09 11:47:37
% DurationCPUTime: 6.26s
% Computational Cost: add. (6680->432), mult. (14850->613), div. (0->0), fcn. (14270->8), ass. (0->177)
t242 = Ifges(6,4) + Ifges(7,4);
t238 = Ifges(6,5) + Ifges(7,5);
t237 = Ifges(6,6) + Ifges(7,6);
t159 = sin(pkin(10));
t160 = cos(pkin(10));
t163 = sin(qJ(2));
t166 = cos(qJ(2));
t141 = t159 * t166 + t160 * t163;
t136 = t141 * qJD(2);
t161 = sin(qJ(5));
t162 = sin(qJ(4));
t164 = cos(qJ(5));
t165 = cos(qJ(4));
t144 = t161 * t165 + t162 * t164;
t231 = qJD(4) + qJD(5);
t106 = t231 * t144;
t140 = t159 * t163 - t160 * t166;
t137 = t140 * qJD(2);
t174 = t161 * t162 - t164 * t165;
t35 = -t106 * t141 + t137 * t174;
t87 = t174 * t141;
t36 = t144 * t137 + t231 * t87;
t241 = (Ifges(6,2) + Ifges(7,2)) * t36 + t242 * t35 + t237 * t136;
t240 = t242 * t36 + (Ifges(6,1) + Ifges(7,1)) * t35 + t238 * t136;
t105 = t231 * t174;
t239 = -t238 * t105 - t237 * t106;
t236 = Ifges(6,3) + Ifges(7,3);
t195 = qJD(4) * t165;
t172 = -t162 * t137 + t141 * t195;
t235 = pkin(4) * t164;
t234 = -mrSges(6,1) - mrSges(7,1);
t233 = pkin(4) * qJD(5);
t216 = -qJ(3) - pkin(7);
t149 = t216 * t163;
t150 = t216 * t166;
t114 = t149 * t159 - t150 * t160;
t104 = t165 * t114;
t157 = -pkin(2) * t166 - pkin(1);
t93 = t140 * pkin(3) - t141 * pkin(8) + t157;
t61 = t162 * t93 + t104;
t199 = t165 * t137;
t232 = -Ifges(5,5) * t199 + Ifges(5,3) * t136;
t154 = pkin(2) * t159 + pkin(8);
t217 = pkin(9) + t154;
t138 = t217 * t162;
t139 = t217 * t165;
t92 = -t161 * t138 + t164 * t139;
t178 = mrSges(5,1) * t162 + mrSges(5,2) * t165;
t145 = t178 * qJD(4);
t230 = 2 * m(6);
t229 = 2 * m(7);
t228 = -2 * mrSges(4,3);
t113 = -t160 * t149 - t150 * t159;
t226 = 0.2e1 * t113;
t225 = 0.2e1 * t157;
t222 = -t141 / 0.2e1;
t211 = Ifges(5,4) * t162;
t151 = Ifges(5,2) * t165 + t211;
t219 = -t151 / 0.2e1;
t218 = pkin(5) * t106;
t27 = -mrSges(7,2) * t136 + mrSges(7,3) * t36;
t28 = -mrSges(6,2) * t136 + mrSges(6,3) * t36;
t215 = t27 + t28;
t202 = t141 * t165;
t60 = -t114 * t162 + t165 * t93;
t42 = pkin(4) * t140 - pkin(9) * t202 + t60;
t203 = t141 * t162;
t51 = -pkin(9) * t203 + t61;
t22 = t161 * t42 + t164 * t51;
t86 = t144 * t141;
t70 = -mrSges(7,2) * t140 - mrSges(7,3) * t86;
t71 = -mrSges(6,2) * t140 - mrSges(6,3) * t86;
t214 = t70 + t71;
t72 = mrSges(7,1) * t140 + mrSges(7,3) * t87;
t73 = mrSges(6,1) * t140 + mrSges(6,3) * t87;
t213 = t72 + t73;
t212 = -t106 * mrSges(7,1) + t105 * mrSges(7,2);
t210 = Ifges(5,4) * t165;
t209 = Ifges(5,6) * t162;
t182 = qJD(2) * t216;
t135 = qJD(3) * t166 + t163 * t182;
t170 = -t163 * qJD(3) + t166 * t182;
t82 = t135 * t159 - t160 * t170;
t208 = t113 * t82;
t207 = t136 * Ifges(5,5);
t206 = t136 * Ifges(5,6);
t205 = t140 * Ifges(5,6);
t193 = qJD(5) * t164;
t204 = (-t105 * t161 + t144 * t193) * pkin(4);
t196 = qJD(4) * t162;
t194 = qJD(5) * t161;
t191 = pkin(2) * qJD(2) * t163;
t190 = pkin(4) * t196;
t83 = t160 * t135 + t159 * t170;
t84 = pkin(3) * t136 + pkin(8) * t137 + t191;
t183 = -t162 * t83 + t165 * t84;
t17 = pkin(9) * t199 + pkin(4) * t136 + (-t104 + (pkin(9) * t141 - t93) * t162) * qJD(4) + t183;
t23 = -t114 * t196 + t162 * t84 + t165 * t83 + t93 * t195;
t20 = -pkin(9) * t172 + t23;
t6 = -qJD(5) * t22 - t161 * t20 + t164 * t17;
t2 = pkin(5) * t136 - qJ(6) * t35 + qJD(6) * t87 + t6;
t25 = mrSges(7,1) * t136 - mrSges(7,3) * t35;
t189 = m(7) * t2 + t25;
t155 = -pkin(2) * t160 - pkin(3);
t188 = t141 * t196;
t186 = t174 * t194;
t11 = -t36 * mrSges(7,1) + t35 * mrSges(7,2);
t185 = (-mrSges(6,2) - mrSges(7,2)) * t164;
t63 = t106 * mrSges(6,1) - t105 * mrSges(6,2);
t184 = -(2 * Ifges(4,4)) - t209;
t21 = -t161 * t51 + t164 * t42;
t181 = qJD(4) * t217;
t180 = 0.2e1 * t191;
t91 = -t164 * t138 - t139 * t161;
t133 = t162 * t181;
t134 = t165 * t181;
t59 = -qJD(5) * t92 + t133 * t161 - t164 * t134;
t38 = qJ(6) * t105 - qJD(6) * t144 + t59;
t179 = m(7) * t38 + t105 * mrSges(7,3);
t81 = pkin(4) * t203 + t113;
t177 = Ifges(5,1) * t165 - t211;
t176 = -Ifges(5,2) * t162 + t210;
t175 = -t63 + t212;
t148 = -pkin(4) * t165 + t155;
t173 = t236 * t136 + t237 * t36 + t238 * t35;
t56 = pkin(4) * t172 + t82;
t5 = t161 * t17 + t164 * t20 + t42 * t193 - t194 * t51;
t171 = t188 + t199;
t58 = -t164 * t133 - t161 * t134 - t138 * t193 - t139 * t194;
t169 = -t161 * t106 + (t144 * t161 - t164 * t174) * qJD(5);
t37 = -qJ(6) * t106 - qJD(6) * t174 + t58;
t168 = t59 * mrSges(6,1) + t38 * mrSges(7,1) - t58 * mrSges(6,2) - t37 * mrSges(7,2) + t239;
t3 = qJ(6) * t36 - qJD(6) * t86 + t5;
t167 = t6 * mrSges(6,1) + t2 * mrSges(7,1) - t5 * mrSges(6,2) - t3 * mrSges(7,2) + t173;
t158 = Ifges(5,5) * t195;
t156 = pkin(5) + t235;
t152 = Ifges(5,1) * t162 + t210;
t147 = t177 * qJD(4);
t146 = t176 * qJD(4);
t125 = t137 * mrSges(4,2);
t115 = pkin(5) * t174 + t148;
t112 = Ifges(6,1) * t144 - Ifges(6,4) * t174;
t111 = Ifges(7,1) * t144 - Ifges(7,4) * t174;
t110 = Ifges(6,4) * t144 - Ifges(6,2) * t174;
t109 = Ifges(7,4) * t144 - Ifges(7,2) * t174;
t108 = mrSges(6,1) * t174 + mrSges(6,2) * t144;
t107 = mrSges(7,1) * t174 + mrSges(7,2) * t144;
t96 = mrSges(5,1) * t140 - mrSges(5,3) * t202;
t95 = -mrSges(5,2) * t140 - mrSges(5,3) * t203;
t90 = t190 + t218;
t77 = -qJ(6) * t174 + t92;
t76 = -qJ(6) * t144 + t91;
t75 = Ifges(5,5) * t140 + t141 * t177;
t74 = t141 * t176 + t205;
t69 = -mrSges(5,2) * t136 - mrSges(5,3) * t172;
t68 = mrSges(5,1) * t136 + mrSges(5,3) * t171;
t67 = -Ifges(6,1) * t105 - Ifges(6,4) * t106;
t66 = -Ifges(7,1) * t105 - Ifges(7,4) * t106;
t65 = -Ifges(6,4) * t105 - Ifges(6,2) * t106;
t64 = -Ifges(7,4) * t105 - Ifges(7,2) * t106;
t57 = mrSges(5,1) * t172 - mrSges(5,2) * t171;
t54 = mrSges(6,1) * t86 - mrSges(6,2) * t87;
t53 = mrSges(7,1) * t86 - mrSges(7,2) * t87;
t52 = pkin(5) * t86 + t81;
t48 = -Ifges(6,1) * t87 - Ifges(6,4) * t86 + Ifges(6,5) * t140;
t47 = -Ifges(7,1) * t87 - Ifges(7,4) * t86 + Ifges(7,5) * t140;
t46 = -Ifges(6,4) * t87 - Ifges(6,2) * t86 + Ifges(6,6) * t140;
t45 = -Ifges(7,4) * t87 - Ifges(7,2) * t86 + Ifges(7,6) * t140;
t44 = -Ifges(5,1) * t171 - Ifges(5,4) * t172 + t207;
t43 = -Ifges(5,4) * t171 - Ifges(5,2) * t172 + t206;
t26 = mrSges(6,1) * t136 - mrSges(6,3) * t35;
t24 = -t61 * qJD(4) + t183;
t18 = -pkin(5) * t36 + t56;
t14 = -qJ(6) * t86 + t22;
t13 = pkin(5) * t140 + qJ(6) * t87 + t21;
t12 = -mrSges(6,1) * t36 + mrSges(6,2) * t35;
t1 = [0.2e1 * ((-t163 ^ 2 + t166 ^ 2) * Ifges(3,4) - pkin(1) * (mrSges(3,1) * t163 + mrSges(3,2) * t166) + (-Ifges(3,2) + Ifges(3,1)) * t163 * t166) * qJD(2) + (mrSges(4,1) * t180 + t83 * t228 - t184 * t137 + ((2 * Ifges(4,2)) + Ifges(5,3) + t236) * t136 + t173 + t232) * t140 + (mrSges(4,1) * t225 + t114 * t228 - t237 * t86 - t238 * t87) * t136 + (mrSges(4,2) * t180 + t165 * t44 - t162 * t43 - 0.2e1 * Ifges(4,1) * t137 + (Ifges(5,5) * t165 + t184) * t136 + (-t165 * t74 - t162 * t75 + t140 * (-Ifges(5,5) * t162 - Ifges(5,6) * t165)) * qJD(4) + 0.2e1 * (mrSges(4,3) + t178) * t82) * t141 - t241 * t86 - t240 * t87 - t125 * t225 + t57 * t226 + (t13 * t2 + t14 * t3 + t18 * t52) * t229 + (t21 * t6 + t22 * t5 + t56 * t81) * t230 + 0.2e1 * t13 * t25 + 0.2e1 * t21 * t26 + 0.2e1 * t14 * t27 + 0.2e1 * t22 * t28 + 0.2e1 * t52 * t11 + 0.2e1 * t18 * t53 + 0.2e1 * t56 * t54 + 0.2e1 * t60 * t68 + 0.2e1 * t61 * t69 + 0.2e1 * t3 * t70 + 0.2e1 * t5 * t71 + 0.2e1 * t2 * t72 + 0.2e1 * t6 * t73 + 0.2e1 * t81 * t12 + 0.2e1 * t23 * t95 + 0.2e1 * t24 * t96 + 0.2e1 * m(5) * (t23 * t61 + t24 * t60 + t208) + 0.2e1 * m(4) * (t114 * t83 + t157 * t191 + t208) + (t47 + t48) * t35 - (mrSges(4,3) * t226 - t162 * t74 + t165 * t75) * t137 + (t45 + t46) * t36; -t52 * t212 + (m(5) * t82 + t57) * t155 + ((-t136 * t159 + t137 * t160) * mrSges(4,3) + m(4) * (t159 * t83 - t160 * t82)) * pkin(2) + (t13 * t105 - t14 * t106 - t2 * t144 - t174 * t3) * mrSges(7,3) + (t105 * t21 - t106 * t22 - t144 * t6 - t174 * t5) * mrSges(6,3) + (t141 * t147 / 0.2e1 - t137 * t152 / 0.2e1 + t23 * mrSges(5,3) + t43 / 0.2e1 - t82 * mrSges(5,1) + t206 / 0.2e1 + (t141 * t219 - t60 * mrSges(5,3) + t75 / 0.2e1) * qJD(4) + (-qJD(4) * t96 + m(5) * (-t60 * qJD(4) + t23) + t69) * t154) * t165 - t241 * t174 / 0.2e1 + t240 * t144 / 0.2e1 + (t238 * t144 - t237 * t174) * t136 / 0.2e1 + (t158 + t239) * t140 / 0.2e1 + m(6) * (t148 * t56 + t21 * t59 + t22 * t58 + t5 * t92 + t6 * t91) + (Ifges(3,5) * t166 - Ifges(3,6) * t163 + (-mrSges(3,1) * t166 + mrSges(3,2) * t163) * pkin(7)) * qJD(2) + m(7) * (t115 * t18 + t13 * t38 + t14 * t37 + t2 * t76 + t3 * t77 + t52 * t90) - (t66 / 0.2e1 + t67 / 0.2e1) * t87 - (t64 / 0.2e1 + t65 / 0.2e1) * t86 + (t109 / 0.2e1 + t110 / 0.2e1) * t36 + t37 * t70 + t58 * t71 + t38 * t72 + t59 * t73 + t76 * t25 + t77 * t27 + t81 * t63 - t82 * mrSges(4,1) - t83 * mrSges(4,2) + (t146 * t222 - t137 * t219 - t24 * mrSges(5,3) + t44 / 0.2e1 + t82 * mrSges(5,2) + t207 / 0.2e1 + (t152 * t222 - t61 * mrSges(5,3) - t74 / 0.2e1 - t205 / 0.2e1 + (m(6) * t81 + t54) * pkin(4)) * qJD(4) + (-m(5) * t24 - t68 + (-m(5) * t61 - t95) * qJD(4)) * t154) * t162 + t90 * t53 + t91 * t26 + t92 * t28 + t18 * t107 + t56 * t108 + t115 * t11 - Ifges(4,6) * t136 - Ifges(4,5) * t137 + t113 * t145 + t148 * t12 - (t47 / 0.2e1 + t48 / 0.2e1) * t105 - (t45 / 0.2e1 + t46 / 0.2e1) * t106 + (t111 / 0.2e1 + t112 / 0.2e1) * t35; 0.2e1 * t90 * t107 - 0.2e1 * t115 * t212 + 0.2e1 * t155 * t145 + t165 * t146 + t162 * t147 + 0.2e1 * t148 * t63 + (t66 + t67) * t144 - (t64 + t65) * t174 - (t109 + t110) * t106 - (t111 + t112) * t105 + (t165 * t152 + (0.2e1 * pkin(4) * t108 - t151) * t162) * qJD(4) + (t148 * t190 + t58 * t92 + t59 * t91) * t230 + (t115 * t90 + t37 * t77 + t38 * t76) * t229 + 0.2e1 * (t105 * t76 - t106 * t77 - t144 * t38 - t174 * t37) * mrSges(7,3) + 0.2e1 * (t105 * t91 - t106 * t92 - t144 * t59 - t174 * t58) * mrSges(6,3); m(4) * t191 + t136 * mrSges(4,1) + t162 * t69 + t165 * t68 - t125 + t215 * t144 - (t25 + t26) * t174 - t213 * t106 - t214 * t105 + (-t162 * t96 + t165 * t95) * qJD(4) + m(7) * (-t105 * t14 - t106 * t13 + t144 * t3 - t174 * t2) + m(6) * (-t105 * t22 - t106 * t21 + t144 * t5 - t174 * t6) + m(5) * (t162 * t23 + t165 * t24 + (-t162 * t60 + t165 * t61) * qJD(4)); m(6) * (-t105 * t92 - t106 * t91 + t144 * t58 - t174 * t59) + m(7) * (-t105 * t77 - t106 * t76 + t144 * t37 - t174 * t38); 0.4e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * (-t105 * t144 + t106 * t174); t189 * t156 + (t164 * t26 + t215 * t161 + (-t161 * t213 + t164 * t214) * qJD(5) + m(7) * (-t13 * t194 + t14 * t193 + t161 * t3) + m(6) * (t161 * t5 + t164 * t6 + t193 * t22 - t194 * t21)) * pkin(4) - Ifges(5,5) * t188 - t172 * Ifges(5,6) + t167 - t23 * mrSges(5,2) + t24 * mrSges(5,1) + t232; t158 + t179 * t156 + (-t209 + (-mrSges(5,1) * t165 + mrSges(5,2) * t162) * t154) * qJD(4) + (m(6) * (t161 * t58 + t164 * t59 + t193 * t92 - t194 * t91) + m(7) * (t161 * t37 + t193 * t77 - t194 * t76) + t169 * mrSges(7,3) + (t164 * t105 + t169) * mrSges(6,3)) * pkin(4) + t168; -t145 + m(6) * ((-t106 * t164 + t186) * pkin(4) + t204) + m(7) * (pkin(4) * t186 - t106 * t156 + t204) + t175; 0.2e1 * (t185 + ((-t156 + t235) * m(7) + t234) * t161) * t233; pkin(5) * t189 + t167; pkin(5) * t179 + t168; -m(7) * t218 + t175; (t185 + (-m(7) * pkin(5) + t234) * t161) * t233; 0; m(7) * t18 + t11; m(7) * t90 - t212; 0; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
