% Calculate time derivative of joint inertia matrix for
% S6RRPPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
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
% Datum: 2019-03-09 09:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRR9_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR9_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR9_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR9_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR9_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:28:36
% EndTime: 2019-03-09 09:28:45
% DurationCPUTime: 4.04s
% Computational Cost: add. (3569->490), mult. (8941->702), div. (0->0), fcn. (7646->8), ass. (0->211)
t246 = Ifges(5,4) + Ifges(4,5);
t245 = Ifges(3,5) + Ifges(5,5);
t143 = sin(qJ(6));
t146 = cos(qJ(6));
t147 = cos(qJ(5));
t181 = qJD(6) * t147;
t144 = sin(qJ(5));
t187 = qJD(5) * t144;
t151 = t143 * t187 - t146 * t181;
t244 = t143 / 0.2e1;
t226 = t146 / 0.2e1;
t194 = t143 ^ 2 + t146 ^ 2;
t140 = cos(pkin(6));
t139 = sin(pkin(6));
t148 = cos(qJ(2));
t199 = t139 * t148;
t145 = sin(qJ(2));
t222 = pkin(1) * t140;
t125 = t145 * t222;
t88 = pkin(8) * t199 + t125;
t77 = -qJ(3) * t140 - t88;
t68 = pkin(3) * t199 - t77;
t52 = pkin(4) * t199 - pkin(9) * t140 + t68;
t141 = qJ(3) - pkin(9);
t142 = pkin(2) + qJ(4);
t165 = -t142 * t148 - pkin(1);
t55 = (-t141 * t145 + t165) * t139;
t215 = t144 * t52 + t147 * t55;
t243 = qJD(5) * t215;
t206 = t144 * mrSges(6,2);
t95 = (t147 * mrSges(6,1) - t206) * qJD(5);
t189 = qJD(3) * t145;
t193 = qJD(2) * t139;
t176 = t145 * t193;
t118 = pkin(2) * t176;
t196 = qJ(4) * t176 + t118;
t35 = (-t189 + (-qJD(2) * t141 - qJD(4)) * t148) * t139 + t196;
t231 = pkin(3) + pkin(8);
t180 = pkin(4) + t231;
t169 = t139 * t180;
t126 = t148 * t222;
t119 = qJD(2) * t126;
t130 = t140 * qJD(3);
t195 = t119 + t130;
t46 = -qJD(2) * t145 * t169 + t195;
t9 = -t144 * t35 + t147 * t46 - t243;
t242 = 2 * m(4);
t241 = 2 * m(5);
t240 = 2 * m(6);
t239 = 2 * m(7);
t238 = -0.2e1 * pkin(1);
t237 = 2 * mrSges(4,1);
t236 = -2 * mrSges(3,3);
t201 = qJ(3) * t145;
t78 = (-pkin(2) * t148 - pkin(1) - t201) * t139;
t235 = -0.2e1 * t78;
t200 = t139 * t145;
t86 = t140 * t147 + t144 * t200;
t65 = -t143 * t86 + t146 * t199;
t234 = t65 / 0.2e1;
t66 = t143 * t199 + t146 * t86;
t233 = t66 / 0.2e1;
t211 = Ifges(7,4) * t143;
t160 = Ifges(7,1) * t146 - t211;
t82 = Ifges(7,5) * t144 + t147 * t160;
t232 = t82 / 0.2e1;
t230 = Ifges(7,5) * t244 + Ifges(7,6) * t226;
t108 = Ifges(7,2) * t146 + t211;
t229 = t108 / 0.2e1;
t228 = -t143 / 0.2e1;
t227 = t144 / 0.2e1;
t225 = -t147 / 0.2e1;
t224 = t147 / 0.2e1;
t223 = m(7) * t144;
t221 = pkin(5) * t144;
t220 = pkin(10) * t147;
t83 = -pkin(8) * t176 + t119;
t219 = t83 * mrSges(3,2);
t218 = mrSges(4,1) + mrSges(5,1);
t153 = -t140 * t144 + t147 * t200;
t192 = qJD(2) * t148;
t175 = t139 * t192;
t63 = qJD(5) * t153 + t144 * t175;
t25 = -qJD(6) * t66 - t143 * t63 - t146 * t176;
t26 = qJD(6) * t65 - t143 * t176 + t63 * t146;
t12 = -mrSges(7,1) * t25 + mrSges(7,2) * t26;
t44 = -mrSges(6,1) * t176 - mrSges(6,3) * t63;
t217 = -t12 + t44;
t31 = -mrSges(7,1) * t65 + mrSges(7,2) * t66;
t71 = mrSges(6,1) * t199 - mrSges(6,3) * t86;
t216 = t31 - t71;
t64 = qJD(5) * t86 - t147 * t175;
t214 = Ifges(6,5) * t63 - Ifges(6,6) * t64;
t213 = Ifges(6,4) * t144;
t212 = Ifges(6,4) * t147;
t210 = Ifges(7,4) * t146;
t209 = Ifges(6,5) * t144;
t208 = Ifges(7,6) * t143;
t17 = pkin(10) * t199 + t215;
t87 = -pkin(8) * t200 + t126;
t79 = -pkin(2) * t140 - t87;
t166 = qJ(4) * t140 - t79;
t51 = (-pkin(3) - pkin(4)) * t200 + t166;
t27 = -pkin(5) * t153 - pkin(10) * t86 + t51;
t10 = -t143 * t17 + t146 * t27;
t207 = t10 * t143;
t84 = t88 * qJD(2);
t205 = t145 * t84;
t204 = t147 * mrSges(6,2);
t203 = t147 * mrSges(7,3);
t161 = mrSges(7,1) * t143 + mrSges(7,2) * t146;
t94 = t161 * qJD(6);
t202 = t147 * t94;
t198 = t143 * t144;
t197 = t144 * t146;
t136 = t144 ^ 2;
t191 = qJD(3) * t136;
t190 = qJD(3) * t144;
t188 = qJD(4) * t142;
t186 = qJD(5) * t146;
t185 = qJD(5) * t147;
t184 = qJD(6) * t143;
t183 = qJD(6) * t144;
t182 = qJD(6) * t146;
t138 = t147 ^ 2;
t133 = t138 * qJD(3);
t179 = 0.2e1 * t148;
t178 = Ifges(4,6) + Ifges(3,4) - Ifges(5,6);
t5 = Ifges(7,5) * t26 + Ifges(7,6) * t25 + Ifges(7,3) * t64;
t177 = (mrSges(4,2) - mrSges(3,1)) * t84;
t174 = t141 * t185;
t173 = t143 * t185;
t172 = t146 * t185;
t105 = -mrSges(7,1) * t146 + mrSges(7,2) * t143;
t168 = m(7) * pkin(5) + mrSges(6,1) - t105;
t167 = t175 * t245 + t176 * t246;
t30 = t64 * mrSges(6,1) + t63 * mrSges(6,2);
t164 = -t141 * t183 + qJD(4) + (pkin(5) * t147 + pkin(10) * t144) * qJD(5);
t129 = t140 * qJD(4);
t15 = t64 * pkin(5) - t63 * pkin(10) + t129 + (-t148 * t169 - t125) * qJD(2);
t8 = t144 * t46 + t147 * t35 + t185 * t52 - t187 * t55;
t3 = -pkin(10) * t176 + t8;
t1 = qJD(6) * t10 + t143 * t15 + t146 * t3;
t11 = t143 * t27 + t146 * t17;
t2 = -qJD(6) * t11 - t143 * t3 + t146 * t15;
t163 = t1 * t146 - t143 * t2;
t106 = t144 * mrSges(6,1) + t204;
t110 = Ifges(7,1) * t143 + t210;
t159 = -Ifges(7,2) * t143 + t210;
t13 = -mrSges(7,2) * t64 + mrSges(7,3) * t25;
t14 = mrSges(7,1) * t64 - mrSges(7,3) * t26;
t158 = t146 * t13 - t143 * t14;
t36 = mrSges(7,2) * t153 + mrSges(7,3) * t65;
t37 = -mrSges(7,1) * t153 - mrSges(7,3) * t66;
t157 = t143 * t37 - t146 * t36;
t18 = -t144 * t55 + t147 * t52;
t156 = t144 * t18 - t147 * t215;
t152 = t143 * t181 + t144 * t186;
t103 = t142 - t220 + t221;
t150 = qJD(6) * t103 + t174 + t190;
t149 = -t10 * t182 - t11 * t184 + t163;
t48 = -t152 * Ifges(7,5) + Ifges(7,6) * t151 + Ifges(7,3) * t185;
t128 = Ifges(7,5) * t182;
t121 = t141 * t133;
t111 = Ifges(6,1) * t147 - t213;
t109 = -Ifges(6,2) * t144 + t212;
t102 = mrSges(7,1) * t144 - t146 * t203;
t101 = -mrSges(7,2) * t144 - t143 * t203;
t100 = (-Ifges(6,1) * t144 - t212) * qJD(5);
t99 = t160 * qJD(6);
t98 = (-Ifges(6,2) * t147 - t213) * qJD(5);
t97 = t159 * qJD(6);
t96 = -Ifges(7,6) * t184 + t128;
t93 = mrSges(5,1) * t199 + mrSges(5,2) * t140;
t92 = -mrSges(4,1) * t199 - mrSges(4,3) * t140;
t91 = mrSges(5,1) * t200 - mrSges(5,3) * t140;
t89 = t161 * t147;
t81 = Ifges(7,6) * t144 + t147 * t159;
t80 = Ifges(7,3) * t144 + (Ifges(7,5) * t146 - t208) * t147;
t76 = -t130 - t83;
t75 = -mrSges(7,2) * t185 + mrSges(7,3) * t151;
t74 = mrSges(7,1) * t185 + mrSges(7,3) * t152;
t73 = t103 * t143 + t141 * t197;
t72 = t103 * t146 - t141 * t198;
t70 = -mrSges(6,2) * t199 + mrSges(6,3) * t153;
t69 = t118 + (-qJ(3) * t192 - t189) * t139;
t67 = (t165 - t201) * t139;
t59 = -t129 + (t199 * t231 + t125) * qJD(2);
t58 = -t176 * t231 + t195;
t57 = pkin(3) * t200 - t166;
t56 = -mrSges(7,1) * t151 - mrSges(7,2) * t152;
t54 = -mrSges(6,1) * t153 + mrSges(6,2) * t86;
t50 = -t110 * t181 + (Ifges(7,5) * t147 - t144 * t160) * qJD(5);
t49 = -t108 * t181 + (Ifges(7,6) * t147 - t144 * t159) * qJD(5);
t47 = -t129 + (t180 * t199 + t125) * qJD(2);
t45 = mrSges(6,2) * t176 - mrSges(6,3) * t64;
t42 = (-t189 + (-qJ(3) * qJD(2) - qJD(4)) * t148) * t139 + t196;
t40 = Ifges(6,1) * t86 + Ifges(6,4) * t153 + Ifges(6,5) * t199;
t39 = Ifges(6,4) * t86 + Ifges(6,2) * t153 + Ifges(6,6) * t199;
t33 = -t143 * t150 + t146 * t164;
t32 = t143 * t164 + t146 * t150;
t29 = Ifges(6,1) * t63 - Ifges(6,4) * t64 - Ifges(6,5) * t176;
t28 = Ifges(6,4) * t63 - Ifges(6,2) * t64 - Ifges(6,6) * t176;
t22 = Ifges(7,1) * t66 + Ifges(7,4) * t65 - Ifges(7,5) * t153;
t21 = Ifges(7,4) * t66 + Ifges(7,2) * t65 - Ifges(7,6) * t153;
t20 = Ifges(7,5) * t66 + Ifges(7,6) * t65 - Ifges(7,3) * t153;
t16 = -pkin(5) * t199 - t18;
t7 = Ifges(7,1) * t26 + Ifges(7,4) * t25 + Ifges(7,5) * t64;
t6 = Ifges(7,4) * t26 + Ifges(7,2) * t25 + Ifges(7,6) * t64;
t4 = pkin(5) * t176 - t9;
t19 = [(t167 + 0.2e1 * t177 - 0.2e1 * t219) * t140 + (t18 * t9 + t215 * t8 - t47 * t51) * t240 + 0.2e1 * t215 * t45 - (t5 - t28) * t153 + (t20 - t39) * t64 + 0.2e1 * m(3) * (t83 * t88 - t84 * t87) + (t148 * t214 + t205 * t237 + 0.2e1 * t42 * (-mrSges(5,2) * t145 - mrSges(5,3) * t148) + 0.2e1 * t69 * (mrSges(4,2) * t148 - mrSges(4,3) * t145) + 0.2e1 * (t148 * t83 + t205) * mrSges(3,3) + ((t79 * t237 + 0.2e1 * t57 * mrSges(5,1) - 0.2e1 * t67 * mrSges(5,2) + t87 * t236 + mrSges(4,3) * t235 + (mrSges(3,2) * t238 + t178 * t179) * t139 + (-(2 * Ifges(4,4)) + t245) * t140) * t148 + (t77 * t237 - 0.2e1 * t68 * mrSges(5,1) + mrSges(4,2) * t235 + t88 * t236 + 0.2e1 * t67 * mrSges(5,3) - Ifges(6,5) * t86 - Ifges(6,6) * t153 + (-(2 * Ifges(3,6)) + t246) * t140 + (-0.2e1 * t145 * t178 + mrSges(3,1) * t238 + (Ifges(3,1) - Ifges(3,2) + Ifges(4,2) - Ifges(5,2) - Ifges(4,3) + Ifges(5,3) - Ifges(6,3)) * t179) * t139) * t145) * qJD(2)) * t139 + (t1 * t11 + t10 * t2 + t16 * t4) * t239 + (t42 * t67 + t57 * t59 + t58 * t68) * t241 + (t69 * t78 + t76 * t77 + t79 * t84) * t242 + 0.2e1 * t11 * t13 + 0.2e1 * t10 * t14 + 0.2e1 * t16 * t12 + t25 * t21 + t26 * t22 + 0.2e1 * t4 * t31 + 0.2e1 * t1 * t36 + 0.2e1 * t2 * t37 + 0.2e1 * t18 * t44 + 0.2e1 * t51 * t30 - 0.2e1 * t47 * t54 + t63 * t40 + t65 * t6 + t66 * t7 + 0.2e1 * t8 * t70 + 0.2e1 * t9 * t71 + t86 * t29 + 0.2e1 * t59 * t91 + 0.2e1 * t76 * t92 + 0.2e1 * t58 * t93; m(6) * (t141 * t144 * t8 + qJD(4) * t51 - t142 * t47 + t190 * t215) - (t48 / 0.2e1 - t98 / 0.2e1) * t153 - t219 + t167 + (-t8 * mrSges(6,3) + qJD(3) * t70 + t141 * t45 + t5 / 0.2e1 - t28 / 0.2e1) * t144 + (t54 - t91) * qJD(4) + t177 + m(4) * (-pkin(2) * t84 - qJ(3) * t76 - qJD(3) * t77) + m(5) * (qJ(3) * t58 + qJD(3) * t68 - qJD(4) * t57 - t142 * t59) + m(7) * (t1 * t73 + t10 * t33 + t11 * t32 + t2 * t72) + (t80 / 0.2e1 - t109 / 0.2e1) * t64 + t50 * t233 + t49 * t234 + t26 * t232 + (-t22 * t197 / 0.2e1 + t21 * t198 / 0.2e1 - t144 * t40 / 0.2e1 + t20 * t224 + t39 * t225 + (-Ifges(6,6) * t147 - t209) * t199 / 0.2e1 + t156 * mrSges(6,3) + (-m(6) * t156 + t144 * t216 + t147 * t70 + t16 * t223) * t141) * qJD(5) + ((-mrSges(4,1) * pkin(2) - mrSges(5,1) * t142 - Ifges(4,4)) * t148 + (Ifges(6,5) * t225 + Ifges(6,6) * t227 - qJ(3) * t218 - Ifges(3,6)) * t145) * t193 + (-t9 * mrSges(6,3) + t6 * t228 + t7 * t226 + t29 / 0.2e1 + t217 * t141 + (-t146 * t21 / 0.2e1 + t22 * t228) * qJD(6) - t216 * qJD(3) + m(7) * (-qJD(3) * t16 - t141 * t4) + m(6) * (qJD(3) * t18 + t141 * t9)) * t147 + t32 * t36 + t33 * t37 + t16 * t56 + t58 * mrSges(5,2) - t59 * mrSges(5,3) + (-t92 + t93) * qJD(3) + t72 * t14 + t73 * t13 + t10 * t74 + t11 * t75 - t76 * mrSges(4,3) + t25 * t81 / 0.2e1 + t4 * t89 + t51 * t95 + t86 * t100 / 0.2e1 + t1 * t101 + t2 * t102 - t47 * t106 + t63 * t111 / 0.2e1 + t142 * t30; 0.2e1 * t32 * t101 + 0.2e1 * t33 * t102 + 0.2e1 * t142 * t95 + 0.2e1 * t72 * t74 + 0.2e1 * t73 * t75 + (t141 * t191 + t121 + t188) * t240 + t188 * t241 + (t32 * t73 + t33 * t72 + t121) * t239 + (t48 - t98 + (0.2e1 * t141 * t89 + t143 * t81 - t146 * t82 - t111) * qJD(5)) * t144 + (-0.2e1 * t141 * t56 - t143 * t49 + t146 * t50 + t100 + (-t143 * t82 - t146 * t81) * qJD(6) + (-0.2e1 * t141 ^ 2 * t223 - t109 + t80) * qJD(5)) * t147 + 0.2e1 * (mrSges(5,3) + t106) * qJD(4) + (0.2e1 * mrSges(5,2) + 0.2e1 * mrSges(4,3) - 0.2e1 * t89 * t147 + 0.2e1 * (-t136 - t138) * mrSges(6,3) + (t242 + t241) * qJ(3)) * qJD(3); -t143 * t13 - t146 * t14 + t157 * qJD(6) + t218 * t175 + m(7) * (-t1 * t143 - t146 * t2 + (-t11 * t146 + t207) * qJD(6)) + m(6) * t47 + m(5) * t59 + m(4) * t84 - t30; m(7) * (-t143 * t32 - t146 * t33 + (t143 * t72 - t146 * t73) * qJD(6)) - t101 * t182 - t143 * t75 + t102 * t184 - t146 * t74 - t95 + (-m(6) - m(5)) * qJD(4); 0; m(5) * t58 - mrSges(5,1) * t176 + ((-t157 + t70) * qJD(5) + m(7) * (-qJD(5) * t207 + t11 * t186 - t4) + m(6) * (t9 + t243) + t217) * t147 + (t45 + (-t143 * t36 - t146 * t37) * qJD(6) + t216 * qJD(5) + m(7) * (qJD(5) * t16 + t149) + m(6) * (-qJD(5) * t18 + t8) + t158) * t144; m(5) * qJD(3) + (-t56 + (t101 * t146 - t102 * t143) * qJD(5)) * t147 + m(7) * (t73 * t172 - t72 * t173 + t133) + m(6) * (t133 + t191) + (qJD(5) * t89 + m(7) * (-t143 * t33 + t146 * t32 - t182 * t72 - t184 * t73 - 0.2e1 * t174) - t101 * t184 + t146 * t75 - t102 * t182 - t143 * t74) * t144; 0; 0.2e1 * (-0.1e1 + t194) * t185 * t223; -Ifges(6,3) * t176 - t8 * mrSges(6,2) + t9 * mrSges(6,1) + t16 * t94 - t153 * t96 / 0.2e1 + t97 * t234 + t99 * t233 + t4 * t105 + t64 * t230 + t25 * t229 + t26 * t110 / 0.2e1 + t7 * t244 + t6 * t226 + (t21 * t228 + t22 * t226) * qJD(6) + (-m(7) * t4 - t12) * pkin(5) + ((-t10 * t146 - t11 * t143) * qJD(6) + t163) * mrSges(7,3) + (m(7) * t149 - t182 * t37 - t184 * t36 + t158) * pkin(10) + t214; -t141 * t202 - pkin(5) * t56 + t96 * t227 + (t147 * t168 - t206) * qJD(3) + (-t209 + (t230 - Ifges(6,6)) * t147 + (-t144 * t168 - t204) * t141) * qJD(5) + (-t110 * t187 / 0.2e1 + t32 * mrSges(7,3) + t99 * t224 + t49 / 0.2e1 + (-t72 * mrSges(7,3) + t108 * t225 + t232) * qJD(6) + (m(7) * (-qJD(6) * t72 + t32) + t75 - qJD(6) * t102) * pkin(10)) * t146 + (t187 * t229 - t33 * mrSges(7,3) + t97 * t225 + t50 / 0.2e1 + (t110 * t225 - t73 * mrSges(7,3) - t81 / 0.2e1) * qJD(6) + (m(7) * (-qJD(6) * t73 - t33) - qJD(6) * t101 - t74) * pkin(10)) * t143; 0; -t202 + (t144 * t105 + m(7) * (t194 * t220 - t221) + t194 * t203 - t106) * qJD(5); -0.2e1 * pkin(5) * t94 + t143 * t99 + t146 * t97 + (-t108 * t143 + t110 * t146) * qJD(6); mrSges(7,1) * t2 - mrSges(7,2) * t1 + t5; mrSges(7,1) * t33 - mrSges(7,2) * t32 + t48; t94; (t143 * t183 - t172) * mrSges(7,2) + (-t144 * t182 - t173) * mrSges(7,1); t128 + (pkin(10) * t105 - t208) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t19(1) t19(2) t19(4) t19(7) t19(11) t19(16); t19(2) t19(3) t19(5) t19(8) t19(12) t19(17); t19(4) t19(5) t19(6) t19(9) t19(13) t19(18); t19(7) t19(8) t19(9) t19(10) t19(14) t19(19); t19(11) t19(12) t19(13) t19(14) t19(15) t19(20); t19(16) t19(17) t19(18) t19(19) t19(20) t19(21);];
Mq  = res;
