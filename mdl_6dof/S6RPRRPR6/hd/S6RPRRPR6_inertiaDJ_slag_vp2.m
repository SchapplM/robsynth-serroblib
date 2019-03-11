% Calculate time derivative of joint inertia matrix for
% S6RPRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 05:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPR6_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR6_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR6_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR6_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR6_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR6_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR6_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:15:28
% EndTime: 2019-03-09 05:15:36
% DurationCPUTime: 3.85s
% Computational Cost: add. (8390->436), mult. (18954->652), div. (0->0), fcn. (19179->10), ass. (0->180)
t235 = Ifges(5,3) + Ifges(6,3);
t169 = sin(pkin(11));
t171 = cos(pkin(11));
t174 = sin(qJ(4));
t177 = cos(qJ(4));
t149 = t169 * t177 + t171 * t174;
t142 = t149 * qJD(4);
t184 = t169 * t174 - t171 * t177;
t143 = t184 * qJD(4);
t201 = qJD(4) * t177;
t234 = Ifges(5,5) * t201 - Ifges(6,5) * t143 - Ifges(6,6) * t142;
t170 = sin(pkin(10));
t172 = cos(pkin(10));
t175 = sin(qJ(3));
t178 = cos(qJ(3));
t148 = t170 * t175 - t178 * t172;
t144 = t148 * qJD(3);
t150 = t170 * t178 + t175 * t172;
t195 = t150 * t201;
t182 = -t174 * t144 + t195;
t173 = sin(qJ(6));
t176 = cos(qJ(6));
t197 = -pkin(2) * t172 - pkin(1);
t107 = pkin(3) * t148 - pkin(8) * t150 + t197;
t218 = pkin(7) + qJ(2);
t157 = t218 * t170;
t158 = t218 * t172;
t117 = -t175 * t157 + t158 * t178;
t112 = t177 * t117;
t145 = t150 * qJD(3);
t183 = qJ(5) * t144 - qJD(5) * t150;
t106 = pkin(3) * t145 + pkin(8) * t144;
t231 = -t178 * t157 - t158 * t175;
t89 = -t148 * qJD(2) + qJD(3) * t231;
t191 = t177 * t106 - t174 * t89;
t24 = pkin(4) * t145 + t183 * t177 + (-t112 + (qJ(5) * t150 - t107) * t174) * qJD(4) + t191;
t200 = t174 * t106 + t107 * t201 + t177 * t89;
t27 = -qJ(5) * t195 + (-qJD(4) * t117 + t183) * t174 + t200;
t11 = -t169 * t27 + t171 * t24;
t63 = -t142 * t150 + t144 * t184;
t4 = pkin(5) * t145 - pkin(9) * t63 + t11;
t12 = t169 * t24 + t171 * t27;
t62 = t143 * t150 + t144 * t149;
t7 = pkin(9) * t62 + t12;
t207 = t150 * t177;
t78 = t177 * t107 - t117 * t174;
t50 = pkin(4) * t148 - qJ(5) * t207 + t78;
t208 = t150 * t174;
t79 = t174 * t107 + t112;
t61 = -qJ(5) * t208 + t79;
t30 = -t169 * t61 + t171 * t50;
t98 = t184 * t150;
t20 = pkin(5) * t148 + pkin(9) * t98 + t30;
t31 = t169 * t50 + t171 * t61;
t97 = t149 * t150;
t21 = -pkin(9) * t97 + t31;
t9 = -t173 * t21 + t176 * t20;
t2 = qJD(6) * t9 + t173 * t4 + t176 * t7;
t10 = t173 * t20 + t176 * t21;
t3 = -qJD(6) * t10 - t173 * t7 + t176 * t4;
t232 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t189 = mrSges(5,1) * t174 + mrSges(5,2) * t177;
t154 = t189 * qJD(4);
t230 = 2 * m(6);
t229 = 2 * m(7);
t228 = -2 * mrSges(4,3);
t226 = -0.2e1 * t231;
t225 = m(6) * pkin(4);
t221 = -t150 / 0.2e1;
t215 = Ifges(5,4) * t174;
t161 = Ifges(5,2) * t177 + t215;
t220 = -t161 / 0.2e1;
t219 = pkin(4) * t169;
t217 = -qJ(5) - pkin(8);
t108 = -t149 * t173 - t176 * t184;
t73 = qJD(6) * t108 - t142 * t173 - t143 * t176;
t109 = t149 * t176 - t173 * t184;
t74 = -qJD(6) * t109 - t142 * t176 + t143 * t173;
t216 = Ifges(7,5) * t73 + Ifges(7,6) * t74;
t214 = Ifges(5,4) * t177;
t213 = Ifges(5,6) * t174;
t90 = qJD(2) * t150 + qJD(3) * t117;
t212 = t231 * t90;
t211 = t145 * Ifges(5,5);
t210 = t145 * Ifges(5,6);
t209 = t148 * Ifges(5,6);
t204 = t177 * t144;
t193 = qJD(4) * t217;
t140 = qJD(5) * t177 + t174 * t193;
t141 = -qJD(5) * t174 + t177 * t193;
t93 = t171 * t140 + t169 * t141;
t159 = t217 * t174;
t160 = t217 * t177;
t119 = t169 * t159 - t171 * t160;
t202 = qJD(4) * t174;
t64 = t173 * t98 - t176 * t97;
t18 = qJD(6) * t64 + t173 * t62 + t176 * t63;
t65 = -t173 * t97 - t176 * t98;
t19 = -qJD(6) * t65 - t173 * t63 + t176 * t62;
t199 = Ifges(7,5) * t18 + Ifges(7,6) * t19 + Ifges(7,3) * t145;
t198 = pkin(4) * t202;
t165 = -pkin(4) * t177 - pkin(3);
t196 = t150 * t202;
t36 = -t62 * mrSges(6,1) + t63 * mrSges(6,2);
t8 = -t19 * mrSges(7,1) + t18 * mrSges(7,2);
t41 = -t74 * mrSges(7,1) + t73 * mrSges(7,2);
t194 = -(2 * Ifges(4,4)) - t213;
t192 = t145 * mrSges(4,1) - t144 * mrSges(4,2);
t103 = t142 * mrSges(6,1) - t143 * mrSges(6,2);
t92 = -t140 * t169 + t171 * t141;
t118 = t171 * t159 + t160 * t169;
t91 = pkin(4) * t208 - t231;
t164 = pkin(4) * t171 + pkin(5);
t138 = t164 * t176 - t173 * t219;
t127 = t138 * qJD(6);
t139 = t164 * t173 + t176 * t219;
t128 = t139 * qJD(6);
t188 = -t128 * mrSges(7,1) - t127 * mrSges(7,2);
t187 = Ifges(5,1) * t177 - t215;
t186 = -Ifges(5,2) * t174 + t214;
t94 = -pkin(9) * t149 + t118;
t95 = -pkin(9) * t184 + t119;
t56 = -t173 * t95 + t176 * t94;
t57 = t173 * t94 + t176 * t95;
t80 = pkin(9) * t143 + t92;
t81 = -pkin(9) * t142 + t93;
t28 = qJD(6) * t56 + t173 * t80 + t176 * t81;
t29 = -qJD(6) * t57 - t173 * t81 + t176 * t80;
t185 = t29 * mrSges(7,1) - t28 * mrSges(7,2) + t216;
t181 = t196 + t204;
t180 = -t103 - t41;
t179 = -Ifges(5,5) * t204 + Ifges(6,5) * t63 + Ifges(6,6) * t62 + t235 * t145 + t199;
t67 = t182 * pkin(4) + t90;
t162 = Ifges(5,1) * t174 + t214;
t156 = t187 * qJD(4);
t155 = t186 * qJD(4);
t121 = pkin(5) * t184 + t165;
t120 = pkin(5) * t142 + t198;
t115 = Ifges(6,1) * t149 - Ifges(6,4) * t184;
t114 = Ifges(6,4) * t149 - Ifges(6,2) * t184;
t113 = mrSges(6,1) * t184 + mrSges(6,2) * t149;
t111 = mrSges(5,1) * t148 - mrSges(5,3) * t207;
t110 = -mrSges(5,2) * t148 - mrSges(5,3) * t208;
t105 = -Ifges(6,1) * t143 - Ifges(6,4) * t142;
t104 = -Ifges(6,4) * t143 - Ifges(6,2) * t142;
t87 = Ifges(5,5) * t148 + t150 * t187;
t86 = t150 * t186 + t209;
t85 = mrSges(6,1) * t148 + mrSges(6,3) * t98;
t84 = -mrSges(6,2) * t148 - mrSges(6,3) * t97;
t83 = -mrSges(5,2) * t145 - mrSges(5,3) * t182;
t82 = mrSges(5,1) * t145 + mrSges(5,3) * t181;
t77 = Ifges(7,1) * t109 + Ifges(7,4) * t108;
t76 = Ifges(7,4) * t109 + Ifges(7,2) * t108;
t75 = -mrSges(7,1) * t108 + mrSges(7,2) * t109;
t69 = mrSges(5,1) * t182 - mrSges(5,2) * t181;
t68 = mrSges(6,1) * t97 - mrSges(6,2) * t98;
t66 = pkin(5) * t97 + t91;
t54 = -Ifges(5,1) * t181 - Ifges(5,4) * t182 + t211;
t53 = -Ifges(5,4) * t181 - Ifges(5,2) * t182 + t210;
t52 = -Ifges(6,1) * t98 - Ifges(6,4) * t97 + Ifges(6,5) * t148;
t51 = -Ifges(6,4) * t98 - Ifges(6,2) * t97 + Ifges(6,6) * t148;
t49 = mrSges(7,1) * t148 - mrSges(7,3) * t65;
t48 = -mrSges(7,2) * t148 + mrSges(7,3) * t64;
t45 = mrSges(6,1) * t145 - mrSges(6,3) * t63;
t44 = -mrSges(6,2) * t145 + mrSges(6,3) * t62;
t43 = Ifges(7,1) * t73 + Ifges(7,4) * t74;
t42 = Ifges(7,4) * t73 + Ifges(7,2) * t74;
t40 = -qJD(4) * t79 + t191;
t39 = -t117 * t202 + t200;
t38 = -t62 * pkin(5) + t67;
t37 = -mrSges(7,1) * t64 + mrSges(7,2) * t65;
t35 = Ifges(7,1) * t65 + Ifges(7,4) * t64 + Ifges(7,5) * t148;
t34 = Ifges(7,4) * t65 + Ifges(7,2) * t64 + Ifges(7,6) * t148;
t33 = Ifges(6,1) * t63 + Ifges(6,4) * t62 + t145 * Ifges(6,5);
t32 = Ifges(6,4) * t63 + Ifges(6,2) * t62 + t145 * Ifges(6,6);
t14 = -mrSges(7,2) * t145 + mrSges(7,3) * t19;
t13 = mrSges(7,1) * t145 - mrSges(7,3) * t18;
t6 = Ifges(7,1) * t18 + Ifges(7,4) * t19 + t145 * Ifges(7,5);
t5 = Ifges(7,4) * t18 + Ifges(7,2) * t19 + t145 * Ifges(7,6);
t1 = [0.2e1 * (m(3) * qJ(2) + mrSges(3,3)) * (t170 ^ 2 + t172 ^ 2) * qJD(2) - (mrSges(4,3) * t226 - t174 * t86 + t177 * t87) * t144 + 0.2e1 * t197 * t192 + t69 * t226 + (t10 * t2 + t3 * t9 + t38 * t66) * t229 + (t11 * t30 + t12 * t31 + t67 * t91) * t230 + ((Ifges(5,5) * t177 + t194) * t150 + ((2 * Ifges(4,2)) + Ifges(7,3) + t235) * t148 - Ifges(6,5) * t98 - Ifges(6,6) * t97 + Ifges(7,5) * t65 + Ifges(7,6) * t64 + t117 * t228) * t145 + 0.2e1 * m(4) * (t117 * t89 - t212) + 0.2e1 * m(5) * (t39 * t79 + t40 * t78 - t212) + (-t194 * t144 + t89 * t228 + t179) * t148 + (t177 * t54 - t174 * t53 - 0.2e1 * Ifges(4,1) * t144 + (t148 * (-Ifges(5,5) * t174 - Ifges(5,6) * t177) - t174 * t87 - t177 * t86) * qJD(4) + 0.2e1 * (mrSges(4,3) + t189) * t90) * t150 + 0.2e1 * t9 * t13 + 0.2e1 * t10 * t14 + t19 * t34 + t18 * t35 + 0.2e1 * t38 * t37 + 0.2e1 * t31 * t44 + 0.2e1 * t30 * t45 + 0.2e1 * t2 * t48 + 0.2e1 * t3 * t49 + t62 * t51 + t63 * t52 + t64 * t5 + t65 * t6 + 0.2e1 * t66 * t8 + 0.2e1 * t67 * t68 + 0.2e1 * t78 * t82 + 0.2e1 * t79 * t83 + 0.2e1 * t12 * t84 + 0.2e1 * t11 * t85 + 0.2e1 * t91 * t36 - t97 * t32 - t98 * t33 + 0.2e1 * t39 * t110 + 0.2e1 * t40 * t111; t108 * t13 + t109 * t14 - t142 * t85 - t143 * t84 - t184 * t45 + t149 * t44 + t174 * t83 + t177 * t82 + t73 * t48 + t74 * t49 + (t110 * t177 - t111 * t174) * qJD(4) + m(7) * (t10 * t73 + t108 * t3 + t109 * t2 + t74 * t9) + m(6) * (-t11 * t184 + t12 * t149 - t142 * t30 - t143 * t31) + m(5) * (t174 * t39 + t177 * t40 + (-t174 * t78 + t177 * t79) * qJD(4)) + t192; 0.2e1 * m(6) * (t142 * t184 - t143 * t149) + 0.2e1 * m(7) * (t108 * t74 + t109 * t73); (t10 * t74 + t108 * t2 - t109 * t3 - t73 * t9) * mrSges(7,3) + (t155 * t221 - t144 * t220 - t40 * mrSges(5,3) + t90 * mrSges(5,2) + t54 / 0.2e1 + t211 / 0.2e1 + (-m(5) * t40 - t82) * pkin(8) + (-t86 / 0.2e1 - t209 / 0.2e1 - t79 * mrSges(5,3) + pkin(4) * t68 + t162 * t221 + t91 * t225 + (-m(5) * t79 - t110) * pkin(8)) * qJD(4)) * t174 + (t150 * t156 / 0.2e1 - t144 * t162 / 0.2e1 + t39 * mrSges(5,3) - t90 * mrSges(5,1) + t53 / 0.2e1 + t210 / 0.2e1 + (t87 / 0.2e1 - t78 * mrSges(5,3) + t150 * t220) * qJD(4) + (m(5) * (-t78 * qJD(4) + t39) + t83 - qJD(4) * t111) * pkin(8)) * t177 + m(6) * (t11 * t118 + t119 * t12 + t165 * t67 + t30 * t92 + t31 * t93) + (-m(5) * t90 - t69) * pkin(3) + m(7) * (t10 * t28 + t120 * t66 + t121 * t38 + t2 * t57 + t29 * t9 + t3 * t56) + (t216 + t234) * t148 / 0.2e1 + (-t11 * t149 - t12 * t184 - t142 * t31 + t143 * t30) * mrSges(6,3) + (Ifges(6,5) * t149 + Ifges(7,5) * t109 - Ifges(6,6) * t184 + Ifges(7,6) * t108) * t145 / 0.2e1 - t184 * t32 / 0.2e1 - t231 * t154 + t165 * t36 + t149 * t33 / 0.2e1 - t142 * t51 / 0.2e1 - t143 * t52 / 0.2e1 - Ifges(4,5) * t144 - Ifges(4,6) * t145 + t28 * t48 + t29 * t49 + t56 * t13 + t57 * t14 + t64 * t42 / 0.2e1 + t65 * t43 / 0.2e1 + t66 * t41 + t73 * t35 / 0.2e1 + t74 * t34 / 0.2e1 + t38 * t75 + t19 * t76 / 0.2e1 + t18 * t77 / 0.2e1 - t89 * mrSges(4,2) - t90 * mrSges(4,1) + t92 * t85 + t93 * t84 + t91 * t103 - t97 * t104 / 0.2e1 - t98 * t105 / 0.2e1 + t108 * t5 / 0.2e1 + t109 * t6 / 0.2e1 + t67 * t113 + t62 * t114 / 0.2e1 + t63 * t115 / 0.2e1 + t118 * t45 + t119 * t44 + t120 * t37 + t121 * t8; m(6) * (-t118 * t142 - t119 * t143 + t149 * t93 - t184 * t92) + m(7) * (t108 * t29 + t109 * t28 + t56 * t74 + t57 * t73); -0.2e1 * pkin(3) * t154 + 0.2e1 * t165 * t103 - t184 * t104 + t149 * t105 + t108 * t42 + t109 * t43 - t142 * t114 - t143 * t115 + 0.2e1 * t120 * t75 + 0.2e1 * t121 * t41 + t177 * t155 + t174 * t156 + t73 * t77 + t74 * t76 + (t177 * t162 + (0.2e1 * pkin(4) * t113 - t161) * t174) * qJD(4) + (t118 * t92 + t119 * t93 + t165 * t198) * t230 + (t120 * t121 + t28 * t57 + t29 * t56) * t229 + 0.2e1 * (t108 * t28 - t109 * t29 - t56 * t73 + t57 * t74) * mrSges(7,3) + 0.2e1 * (t118 * t143 - t119 * t142 - t149 * t92 - t184 * t93) * mrSges(6,3); t179 + m(7) * (t10 * t127 - t128 * t9 + t138 * t3 + t139 * t2) + (t169 * t44 + m(6) * (t11 * t171 + t12 * t169) + t171 * t45) * pkin(4) - t182 * Ifges(5,6) + t138 * t13 + t139 * t14 - Ifges(5,5) * t196 + t11 * mrSges(6,1) - t12 * mrSges(6,2) - t39 * mrSges(5,2) + t40 * mrSges(5,1) + t127 * t48 - t128 * t49 + t232; -t154 + m(7) * (-t108 * t128 + t109 * t127 + t138 * t74 + t139 * t73) + (-t142 * t171 - t143 * t169) * t225 + t180; m(7) * (t127 * t57 - t128 * t56 + t138 * t29 + t139 * t28) + t92 * mrSges(6,1) - t93 * mrSges(6,2) + (-t213 + (-mrSges(5,1) * t177 + mrSges(5,2) * t174) * pkin(8)) * qJD(4) + (m(6) * (t169 * t93 + t171 * t92) + (-t142 * t169 + t143 * t171) * mrSges(6,3)) * pkin(4) + (t108 * t127 + t109 * t128 - t138 * t73 + t139 * t74) * mrSges(7,3) + t185 + t234; 0.2e1 * m(7) * (t127 * t139 - t128 * t138) + 0.2e1 * t188; m(6) * t67 + m(7) * t38 + t36 + t8; 0; m(6) * t198 + m(7) * t120 - t180; 0; 0; t199 + t232; -t41; t185; t188; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
