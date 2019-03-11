% Calculate time derivative of joint inertia matrix for
% S6RRPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
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
% Datum: 2019-03-09 12:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRP7_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP7_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP7_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP7_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP7_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP7_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP7_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:16:48
% EndTime: 2019-03-09 12:16:57
% DurationCPUTime: 3.76s
% Computational Cost: add. (3793->403), mult. (8187->535), div. (0->0), fcn. (6775->6), ass. (0->170)
t118 = sin(qJ(5));
t121 = cos(qJ(5));
t175 = t118 ^ 2 + t121 ^ 2;
t215 = Ifges(7,4) + Ifges(6,5);
t251 = t118 * t215;
t119 = sin(qJ(4));
t122 = cos(qJ(4));
t120 = sin(qJ(2));
t220 = pkin(7) - pkin(8);
t95 = t220 * t120;
t123 = cos(qJ(2));
t96 = t220 * t123;
t243 = -t119 * t96 + t122 * t95;
t184 = qJD(4) * t243;
t250 = Ifges(7,2) + Ifges(6,3);
t249 = m(4) * pkin(7) + mrSges(4,2);
t241 = (mrSges(7,2) + mrSges(6,3)) * t175;
t129 = -mrSges(5,2) + t241;
t168 = qJD(4) * t122;
t248 = t175 * t168;
t166 = qJD(5) * t121;
t69 = t119 * t120 + t122 * t123;
t49 = (qJD(2) - qJD(4)) * t69;
t191 = t118 * t49;
t70 = -t119 * t123 + t120 * t122;
t135 = t70 * t166 + t191;
t167 = qJD(5) * t118;
t187 = t121 * t49;
t134 = t70 * t167 - t187;
t124 = -pkin(2) - pkin(3);
t83 = -t119 * qJ(3) + t122 * t124;
t62 = t122 * qJD(3) + qJD(4) * t83;
t247 = t175 * t62;
t195 = Ifges(7,5) * t121;
t141 = Ifges(7,3) * t118 + t195;
t74 = t141 * qJD(5);
t77 = Ifges(6,4) * t166 - Ifges(6,2) * t167;
t196 = Ifges(7,5) * t118;
t142 = Ifges(7,1) * t121 + t196;
t78 = t142 * qJD(5);
t79 = Ifges(6,1) * t166 - Ifges(6,4) * t167;
t89 = -Ifges(7,3) * t121 + t196;
t198 = Ifges(6,4) * t118;
t92 = Ifges(6,2) * t121 + t198;
t93 = Ifges(7,1) * t118 - t195;
t197 = Ifges(6,4) * t121;
t94 = Ifges(6,1) * t118 + t197;
t125 = (t118 * (t89 - t92) + t121 * (t93 + t94)) * qJD(5) + (t78 + t79) * t118 + (t77 - t74) * t121;
t88 = -t121 * mrSges(6,1) + t118 * mrSges(6,2);
t216 = mrSges(5,1) - t88;
t84 = t122 * qJ(3) + t119 * t124;
t63 = t119 * qJD(3) + qJD(4) * t84;
t244 = t216 * t63;
t242 = -t123 * pkin(2) - t120 * qJ(3);
t86 = -pkin(1) + t242;
t68 = t123 * pkin(3) - t86;
t32 = pkin(4) * t69 - pkin(9) * t70 + t68;
t53 = t119 * t95 + t122 * t96;
t206 = t118 * t32 + t121 * t53;
t177 = t206 * qJD(5);
t158 = t89 / 0.2e1 - t92 / 0.2e1;
t170 = qJD(4) * t119;
t173 = qJD(2) * t123;
t174 = qJD(2) * t120;
t48 = t119 * t173 + t120 * t168 - t122 * t174 - t123 * t170;
t176 = qJ(3) * t173 + t120 * qJD(3);
t58 = t124 * t174 + t176;
t15 = pkin(4) * t48 - pkin(9) * t49 + t58;
t147 = qJD(2) * t96;
t82 = t220 * t174;
t24 = t119 * t147 - t122 * t82 + t184;
t4 = -t118 * t24 + t121 * t15 - t177;
t2 = -pkin(5) * t48 - t4;
t8 = -Ifges(7,1) * t134 + Ifges(7,4) * t48 + Ifges(7,5) * t135;
t9 = -Ifges(6,1) * t134 - Ifges(6,4) * t135 + Ifges(6,5) * t48;
t240 = t158 * t49 + t2 * mrSges(7,2) + t8 / 0.2e1 + t9 / 0.2e1 - t4 * mrSges(6,3);
t3 = t118 * t15 + t121 * t24 + t32 * t166 - t167 * t53;
t1 = qJ(6) * t48 + qJD(6) * t69 + t3;
t157 = t93 / 0.2e1 + t94 / 0.2e1;
t6 = -Ifges(7,5) * t134 + Ifges(7,6) * t48 + Ifges(7,3) * t135;
t7 = -Ifges(6,4) * t134 - Ifges(6,2) * t135 + Ifges(6,6) * t48;
t239 = t157 * t49 + t1 * mrSges(7,2) + t3 * mrSges(6,3) - t6 / 0.2e1 + t7 / 0.2e1;
t139 = pkin(5) * t118 - qJ(6) * t121;
t65 = qJD(5) * t139 - t118 * qJD(6);
t238 = -t118 * (qJ(6) * mrSges(7,2) - Ifges(7,6)) - t121 * (pkin(5) * mrSges(7,2) - t215);
t103 = Ifges(6,6) * t167;
t194 = Ifges(6,6) * t121;
t23 = t139 * t70 - t243;
t182 = qJD(5) * t70;
t183 = qJD(4) * t53;
t25 = -t119 * t82 - t122 * t147 + t183;
t5 = (pkin(5) * t49 + qJ(6) * t182) * t118 + (-qJ(6) * t49 + (pkin(5) * qJD(5) - qJD(6)) * t70) * t121 + t25;
t72 = -mrSges(7,1) * t167 + mrSges(7,3) * t166;
t144 = t118 * mrSges(6,1) + t121 * mrSges(6,2);
t73 = t144 * qJD(5);
t87 = t121 * mrSges(7,1) + t118 * mrSges(7,3);
t235 = -t24 * mrSges(5,2) + Ifges(5,5) * t49 - (Ifges(5,6) - t194 / 0.2e1 + Ifges(7,6) * t121 / 0.2e1 - t251 / 0.2e1) * t48 + (Ifges(6,5) * t166 / 0.2e1 - t103 / 0.2e1 + (Ifges(7,4) * t121 + Ifges(7,6) * t118) * qJD(5) / 0.2e1) * t69 - t23 * t72 - t5 * t87 - t243 * t73;
t20 = -t118 * t53 + t121 * t32;
t13 = -pkin(5) * t69 - t20;
t186 = t121 * t70;
t43 = mrSges(6,1) * t69 - mrSges(6,3) * t186;
t44 = -mrSges(7,1) * t69 + mrSges(7,2) * t186;
t207 = -t43 + t44;
t18 = -mrSges(6,2) * t48 - mrSges(6,3) * t135;
t19 = -mrSges(7,2) * t135 + mrSges(7,3) * t48;
t212 = t18 + t19;
t234 = t207 * qJD(5) + m(7) * (t13 * qJD(5) + t1) + m(6) * (-t20 * qJD(5) + t3) + t212;
t12 = qJ(6) * t69 + t206;
t190 = t118 * t70;
t42 = -mrSges(6,2) * t69 - mrSges(6,3) * t190;
t45 = -mrSges(7,2) * t190 + mrSges(7,3) * t69;
t208 = t42 + t45;
t16 = mrSges(6,1) * t48 + mrSges(6,3) * t134;
t17 = -t48 * mrSges(7,1) - mrSges(7,2) * t134;
t213 = -t16 + t17;
t233 = -t208 * qJD(5) + m(7) * (-t12 * qJD(5) + t2) + m(6) * (-t4 - t177) + t213;
t232 = 2 * m(5);
t231 = 0.2e1 * m(6);
t230 = 0.2e1 * m(7);
t229 = -0.2e1 * pkin(1);
t228 = 0.2e1 * mrSges(5,2);
t227 = -2 * Ifges(5,4);
t226 = 0.2e1 * t25;
t225 = -0.2e1 * t243;
t224 = 0.2e1 * t58;
t223 = -0.2e1 * t73;
t222 = 0.2e1 * t86;
t219 = Ifges(6,6) * t69;
t218 = t25 * t243;
t217 = t243 * t63;
t214 = Ifges(4,5) - Ifges(3,4);
t26 = Ifges(7,6) * t69 + t141 * t70;
t27 = t219 + (-Ifges(6,2) * t118 + t197) * t70;
t211 = t26 - t27;
t28 = Ifges(7,4) * t69 + t142 * t70;
t29 = Ifges(6,5) * t69 + (Ifges(6,1) * t121 - t198) * t70;
t210 = t28 + t29;
t81 = -pkin(9) + t84;
t209 = t247 * t81;
t205 = t247 * pkin(9);
t204 = -t72 + t73;
t199 = t248 * pkin(9);
t189 = t119 * t62;
t185 = t122 * t63;
t171 = qJD(4) * t118;
t169 = qJD(4) * t121;
t165 = qJD(6) * t121;
t163 = 0.2e1 * t123;
t162 = t87 + t216;
t161 = t74 / 0.2e1 - t77 / 0.2e1;
t159 = t78 / 0.2e1 + t79 / 0.2e1;
t151 = t175 * t189 + t248 * t81;
t146 = (m(7) * pkin(9) + mrSges(7,2)) * t121;
t145 = -t123 * mrSges(4,1) - t120 * mrSges(4,3);
t143 = t118 * mrSges(7,1) - t121 * mrSges(7,3);
t140 = -t121 * pkin(5) - t118 * qJ(6);
t136 = t135 * Ifges(7,6) + t187 * t215 + t250 * t48;
t85 = -pkin(4) + t140;
t131 = -t26 / 0.2e1 + t27 / 0.2e1 + t12 * mrSges(7,2) + t206 * mrSges(6,3);
t130 = -t28 / 0.2e1 - t29 / 0.2e1 - t13 * mrSges(7,2) + t20 * mrSges(6,3);
t127 = m(7) * t140 - t87 + t88;
t126 = -m(7) * t139 - t143 - t144;
t80 = pkin(4) - t83;
t59 = -t83 - t85;
t35 = t144 * t70;
t34 = t143 * t70;
t33 = t63 - t65;
t11 = mrSges(6,1) * t135 - mrSges(6,2) * t134;
t10 = mrSges(7,1) * t135 + mrSges(7,3) * t134;
t14 = [0.2e1 * t13 * t17 + 0.2e1 * t12 * t19 + 0.2e1 * t20 * t16 + 0.2e1 * t206 * t18 + 0.2e1 * t23 * t10 + 0.2e1 * t5 * t34 + t35 * t226 + 0.2e1 * t3 * t42 + 0.2e1 * t4 * t43 + 0.2e1 * t2 * t44 + 0.2e1 * t1 * t45 + t11 * t225 + 0.2e1 * (mrSges(5,1) * t68 - mrSges(5,3) * t53) * t48 + (t1 * t12 + t13 * t2 + t23 * t5) * t230 + (t20 * t4 + t206 * t3 - t218) * t231 + (t24 * t53 + t58 * t68 - t218) * t232 + (mrSges(5,3) * t225 + t211 * t118 + t210 * t121 + t68 * t228) * t49 + (mrSges(5,1) * t224 - 0.2e1 * mrSges(5,3) * t24 + (-Ifges(6,6) * t118 + t227) * t49 + ((2 * Ifges(5,2)) + t250) * t48 + t136) * t69 + (mrSges(5,2) * t224 + mrSges(5,3) * t226 + 0.2e1 * Ifges(5,1) * t49 + (t8 + t9) * t121 + (t6 - t7) * t118 + (t227 + t215 * t121 + (-Ifges(6,6) + Ifges(7,6)) * t118) * t48 + ((t211 - t219) * t121 + (-t215 * t69 - t210) * t118) * qJD(5)) * t70 + ((mrSges(3,2) * t229 - 0.2e1 * t86 * mrSges(4,3) - t163 * t214) * t123 + (mrSges(3,1) * t229 + mrSges(4,1) * t222 + 0.2e1 * t214 * t120 + (Ifges(3,1) + Ifges(4,1) - Ifges(3,2) - Ifges(4,3)) * t163) * t120) * qJD(2) + (m(4) * t222 + 0.2e1 * t145) * (pkin(2) * t174 - t176); m(5) * (t24 * t84 - t25 * t83 + t53 * t62 - t217) + (-t84 * t48 - t83 * t49 - t62 * t69 + t63 * t70) * mrSges(5,3) + t216 * t25 + m(7) * (t23 * t33 + t5 * t59) + m(6) * (t25 * t80 - t217) + ((-pkin(2) * mrSges(4,2) + Ifges(4,4) + Ifges(3,5)) * t123 + (-qJ(3) * mrSges(4,2) - Ifges(3,6) + Ifges(4,6)) * t120 + (m(4) * t242 - t123 * mrSges(3,1) + t120 * mrSges(3,2) + t145) * pkin(7)) * qJD(2) + (t130 * qJD(5) + (-qJD(5) * t158 - t159) * t70 + (m(6) * t206 + m(7) * t12 + t208) * t62 + t234 * t81 - t239) * t121 + (t131 * qJD(5) + (qJD(5) * t157 - t161) * t70 + (-m(6) * t20 + m(7) * t13 + t207) * t62 + t233 * t81 - t240) * t118 + t33 * t34 + t59 * t10 + t63 * t35 + t80 * t11 - t235 + t249 * qJD(3) * t123; 0.2e1 * t33 * t87 + 0.2e1 * t59 * t72 + t80 * t223 + (t228 - 0.2e1 * t241) * t62 + (t33 * t59 + t209) * t230 + (t63 * t80 + t209) * t231 + (t62 * t84 - t63 * t83) * t232 + t125 + 0.2e1 * t244 + 0.2e1 * (m(4) * qJ(3) + mrSges(4,3)) * qJD(3); t249 * t173 + (-t49 * mrSges(5,3) - t10 - t11 + (-t69 * mrSges(5,3) + t118 * t207 + t121 * t208) * qJD(4) + m(5) * (-t25 + t183) + m(7) * (t12 * t169 + t13 * t171 - t5) + m(6) * (t169 * t206 - t171 * t20 - t25)) * t122 + (-t48 * mrSges(5,3) + t212 * t121 + t213 * t118 + (t70 * mrSges(5,3) + t34 + t35) * qJD(4) + (-t118 * t208 + t121 * t207) * qJD(5) + m(5) * (t24 - t184) + m(7) * (qJD(4) * t23 + t1 * t121 + t118 * t2 - t12 * t167 + t13 * t166) + m(6) * (-t118 * t4 + t121 * t3 - t166 * t20 - t167 * t206 - t184)) * t119; t204 * t122 + (t162 * t119 - t129 * t122) * qJD(4) + m(7) * (-t122 * t33 + t170 * t59 + t151) + m(6) * (t170 * t80 + t151 - t185) + m(5) * (t189 - t185 + (-t119 * t83 + t122 * t84) * qJD(4)); 0.4e1 * (m(7) / 0.2e1 + m(6) / 0.2e1) * (-0.1e1 + t175) * t119 * t168; m(7) * (t23 * t65 + t5 * t85) - pkin(4) * t11 + t65 * t34 + t85 * t10 + (-m(6) * pkin(4) - t216) * t25 + (t159 * t70 + (t158 * t70 - t130) * qJD(5) + t234 * pkin(9) + t239) * t121 + (t161 * t70 + (-t157 * t70 - t131) * qJD(5) + t233 * pkin(9) + t240) * t118 + t235; (t65 - t33) * t87 + (pkin(4) + t80) * t73 + (-t59 + t85) * t72 - t244 + m(6) * (-pkin(4) * t63 + t205) + m(7) * (t33 * t85 + t59 * t65 + t205) + t129 * t62 - t125; -t162 * t170 + m(7) * (t170 * t85 + t199) + m(6) * (-pkin(4) * t170 + t199) + (-m(7) * t65 + qJD(4) * t129 - t204) * t122; pkin(4) * t223 - 0.2e1 * t85 * t72 + 0.2e1 * (m(7) * t85 - t87) * t65 + t125; -Ifges(6,6) * t191 - pkin(5) * t17 + m(7) * (-pkin(5) * t2 + qJ(6) * t1 + qJD(6) * t12) + qJD(6) * t45 + qJ(6) * t19 + t1 * mrSges(7,3) + t4 * mrSges(6,1) - t2 * mrSges(7,1) - t3 * mrSges(6,2) + (-t194 - t251) * t182 + t136; t103 + (m(7) * t81 - mrSges(7,2)) * t165 + t126 * t62 + (t127 * t81 - t238) * qJD(5); t126 * t168 + (m(7) * t165 + qJD(5) * t127) * t119; -t103 + qJD(6) * t146 + (pkin(9) * t127 + t238) * qJD(5); 0.2e1 * (m(7) * qJ(6) + mrSges(7,3)) * qJD(6); m(7) * t2 + t17; -mrSges(7,2) * t166 + (t118 * t62 + t166 * t81) * m(7); (t118 * t168 + t119 * t166) * m(7); qJD(5) * t146; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t14(1) t14(2) t14(4) t14(7) t14(11) t14(16); t14(2) t14(3) t14(5) t14(8) t14(12) t14(17); t14(4) t14(5) t14(6) t14(9) t14(13) t14(18); t14(7) t14(8) t14(9) t14(10) t14(14) t14(19); t14(11) t14(12) t14(13) t14(14) t14(15) t14(20); t14(16) t14(17) t14(18) t14(19) t14(20) t14(21);];
Mq  = res;
