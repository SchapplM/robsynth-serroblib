% Calculate time derivative of joint inertia matrix for
% S5RRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR10_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR10_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR10_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR10_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR10_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR10_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR10_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:23:50
% EndTime: 2019-12-31 20:23:59
% DurationCPUTime: 3.26s
% Computational Cost: add. (4260->391), mult. (11956->597), div. (0->0), fcn. (11415->10), ass. (0->184)
t132 = sin(pkin(5));
t226 = 0.2e1 * t132;
t135 = sin(qJ(5));
t138 = cos(qJ(5));
t136 = sin(qJ(4));
t172 = qJD(5) * t136;
t139 = cos(qJ(4));
t174 = qJD(4) * t139;
t225 = -t135 * t172 + t138 * t174;
t134 = cos(pkin(5));
t131 = sin(pkin(10));
t133 = cos(pkin(10));
t140 = cos(qJ(2));
t203 = pkin(1) * t134;
t122 = t140 * t203;
t137 = sin(qJ(2));
t193 = -pkin(7) - qJ(3);
t159 = t193 * t137;
t73 = pkin(2) * t134 + t132 * t159 + t122;
t181 = t132 * t140;
t121 = t137 * t203;
t93 = pkin(7) * t181 + t121;
t82 = qJ(3) * t181 + t93;
t54 = t131 * t73 + t133 * t82;
t43 = pkin(8) * t134 + t54;
t180 = t133 * t140;
t182 = t132 * t137;
t85 = t131 * t182 - t132 * t180;
t86 = (t131 * t140 + t133 * t137) * t132;
t96 = (-pkin(2) * t140 - pkin(1)) * t132;
t59 = t85 * pkin(3) - t86 * pkin(8) + t96;
t192 = t136 * t59 + t139 * t43;
t219 = qJD(4) * t192;
t221 = (t193 * t181 - t121) * qJD(2) - qJD(3) * t182;
t115 = qJD(2) * t122;
t66 = t115 + (qJD(2) * t159 + qJD(3) * t140) * t132;
t40 = t221 * t131 + t133 * t66;
t176 = qJD(2) * t132;
t166 = t137 * t176;
t155 = pkin(2) * t166;
t83 = qJD(2) * t86;
t84 = (-t131 * t137 + t180) * t176;
t56 = pkin(3) * t83 - pkin(8) * t84 + t155;
t12 = -t136 * t40 + t139 * t56 - t219;
t146 = t134 * t139 - t136 * t86;
t52 = qJD(4) * t146 + t139 * t84;
t72 = t134 * t136 + t139 * t86;
t58 = t135 * t85 + t138 * t72;
t21 = -qJD(5) * t58 - t135 * t52 + t138 * t83;
t57 = -t135 * t72 + t138 * t85;
t22 = qJD(5) * t57 + t135 * t83 + t138 * t52;
t6 = -mrSges(6,1) * t21 + mrSges(6,2) * t22;
t8 = -pkin(4) * t83 - t12;
t167 = -m(6) * t8 - t6;
t37 = mrSges(5,1) * t83 - mrSges(5,3) * t52;
t224 = t167 + m(5) * (t12 + t219) + t37;
t171 = qJD(5) * t138;
t127 = Ifges(6,5) * t171;
t173 = qJD(5) * t135;
t223 = Ifges(6,6) * t173 / 0.2e1 - t127 / 0.2e1;
t222 = t135 / 0.2e1;
t204 = t138 / 0.2e1;
t177 = t135 ^ 2 + t138 ^ 2;
t144 = -t135 * t174 - t136 * t171;
t217 = 0.2e1 * m(6);
t216 = -2 * mrSges(3,3);
t215 = -2 * mrSges(4,3);
t214 = -2 * Ifges(4,4);
t213 = 0.2e1 * t96;
t212 = t21 / 0.2e1;
t188 = Ifges(6,4) * t135;
t109 = Ifges(6,2) * t138 + t188;
t187 = Ifges(6,4) * t138;
t150 = -Ifges(6,2) * t135 + t187;
t63 = -t109 * t172 + (Ifges(6,6) * t136 + t150 * t139) * qJD(4);
t211 = t63 / 0.2e1;
t111 = Ifges(6,1) * t135 + t187;
t151 = Ifges(6,1) * t138 - t188;
t64 = -t111 * t172 + (Ifges(6,5) * t136 + t151 * t139) * qJD(4);
t210 = t64 / 0.2e1;
t89 = -Ifges(6,5) * t139 + t136 * t151;
t209 = t89 / 0.2e1;
t208 = Ifges(6,5) * t222 + Ifges(6,6) * t204;
t207 = t111 / 0.2e1;
t206 = -t135 / 0.2e1;
t205 = -t138 / 0.2e1;
t202 = pkin(4) * t136;
t201 = pkin(9) * t139;
t39 = t131 * t66 - t133 * t221;
t200 = t39 * mrSges(4,1);
t199 = t40 * mrSges(4,2);
t198 = t83 * Ifges(5,5);
t197 = t83 * Ifges(5,6);
t196 = t85 * Ifges(5,6);
t90 = -pkin(7) * t166 + t115;
t195 = t90 * mrSges(3,2);
t91 = t93 * qJD(2);
t194 = t91 * mrSges(3,1);
t191 = mrSges(6,3) * t136;
t190 = Ifges(5,4) * t136;
t189 = Ifges(5,4) * t139;
t186 = Ifges(6,6) * t135;
t28 = -t136 * t43 + t139 * t59;
t23 = -pkin(4) * t85 - t28;
t185 = qJD(4) * t23;
t124 = pkin(2) * t131 + pkin(8);
t179 = t135 * t139;
t125 = -pkin(2) * t133 - pkin(3);
t95 = -pkin(4) * t139 - pkin(9) * t136 + t125;
t69 = -t124 * t179 + t138 * t95;
t184 = qJD(5) * t69;
t178 = t138 * t139;
t70 = t124 * t178 + t135 * t95;
t183 = qJD(5) * t70;
t175 = qJD(4) * t136;
t170 = t139 * t217;
t51 = qJD(4) * t72 + t136 * t84;
t3 = Ifges(6,5) * t22 + Ifges(6,6) * t21 + Ifges(6,3) * t51;
t169 = Ifges(5,5) * t52 - Ifges(5,6) * t51 + Ifges(5,3) * t83;
t168 = Ifges(3,5) * t140 * t176 + Ifges(4,5) * t84 - Ifges(4,6) * t83;
t165 = t124 * t175;
t53 = -t131 * t82 + t133 * t73;
t106 = (-t201 + t202) * qJD(4);
t44 = t106 * t135 - t138 * t165 + t184;
t157 = t44 - t184;
t45 = t106 * t138 + t135 * t165 - t183;
t156 = -t45 - t183;
t15 = pkin(4) * t51 - pkin(9) * t52 + t39;
t11 = t136 * t56 + t139 * t40 + t59 * t174 - t175 * t43;
t7 = pkin(9) * t83 + t11;
t24 = pkin(9) * t85 + t192;
t42 = -pkin(3) * t134 - t53;
t27 = -pkin(4) * t146 - pkin(9) * t72 + t42;
t9 = -t135 * t24 + t138 * t27;
t1 = qJD(5) * t9 + t135 * t15 + t138 * t7;
t10 = t135 * t27 + t138 * t24;
t2 = -qJD(5) * t10 - t135 * t7 + t138 * t15;
t154 = t1 * t138 - t135 * t2;
t153 = t136 * mrSges(5,1) + t139 * mrSges(5,2);
t107 = -mrSges(6,1) * t138 + mrSges(6,2) * t135;
t152 = mrSges(6,1) * t135 + mrSges(6,2) * t138;
t13 = -mrSges(6,2) * t51 + mrSges(6,3) * t21;
t14 = mrSges(6,1) * t51 - mrSges(6,3) * t22;
t149 = t138 * t13 - t135 * t14;
t17 = Ifges(6,4) * t58 + Ifges(6,2) * t57 - Ifges(6,6) * t146;
t18 = Ifges(6,1) * t58 + Ifges(6,4) * t57 - Ifges(6,5) * t146;
t145 = t17 * t206 + t18 * t204;
t31 = -mrSges(6,1) * t57 + mrSges(6,2) * t58;
t36 = -mrSges(5,2) * t83 - mrSges(5,3) * t51;
t61 = mrSges(5,1) * t85 - mrSges(5,3) * t72;
t142 = t36 + m(5) * (-t28 * qJD(4) + t11) + (t31 - t61) * qJD(4);
t141 = -t10 * t173 - t171 * t9 + t154;
t62 = Ifges(6,5) * t225 + Ifges(6,6) * t144 + Ifges(6,3) * t175;
t128 = Ifges(5,5) * t174;
t112 = Ifges(5,1) * t136 + t189;
t110 = Ifges(5,2) * t139 + t190;
t105 = -mrSges(6,1) * t139 - t138 * t191;
t104 = mrSges(6,2) * t139 - t135 * t191;
t103 = (Ifges(5,1) * t139 - t190) * qJD(4);
t102 = t151 * qJD(5);
t101 = (-Ifges(5,2) * t136 + t189) * qJD(4);
t100 = t150 * qJD(5);
t98 = t153 * qJD(4);
t97 = t152 * qJD(5);
t94 = t152 * t136;
t92 = -pkin(7) * t182 + t122;
t88 = -Ifges(6,6) * t139 + t136 * t150;
t87 = -Ifges(6,3) * t139 + (Ifges(6,5) * t138 - t186) * t136;
t81 = -mrSges(6,2) * t175 + mrSges(6,3) * t144;
t80 = mrSges(6,1) * t175 - mrSges(6,3) * t225;
t76 = t84 * mrSges(4,2);
t67 = t144 * mrSges(6,1) - mrSges(6,2) * t225;
t60 = -mrSges(5,2) * t85 + mrSges(5,3) * t146;
t35 = Ifges(5,1) * t72 + Ifges(5,4) * t146 + Ifges(5,5) * t85;
t34 = Ifges(5,4) * t72 + Ifges(5,2) * t146 + t196;
t33 = -mrSges(6,1) * t146 - mrSges(6,3) * t58;
t32 = mrSges(6,2) * t146 + mrSges(6,3) * t57;
t30 = mrSges(5,1) * t51 + mrSges(5,2) * t52;
t26 = Ifges(5,1) * t52 - Ifges(5,4) * t51 + t198;
t25 = Ifges(5,4) * t52 - Ifges(5,2) * t51 + t197;
t16 = Ifges(6,5) * t58 + Ifges(6,6) * t57 - Ifges(6,3) * t146;
t5 = Ifges(6,1) * t22 + Ifges(6,4) * t21 + Ifges(6,5) * t51;
t4 = Ifges(6,4) * t22 + Ifges(6,2) * t21 + Ifges(6,6) * t51;
t19 = [(0.2e1 * (t137 * t91 + t140 * t90) * mrSges(3,3) + ((t92 * t216 + Ifges(3,5) * t134 + (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t140) * t226) * t140 + (-0.2e1 * Ifges(3,6) * t134 + 0.2e1 * pkin(2) * (mrSges(4,1) * t85 + mrSges(4,2) * t86) + t93 * t216 + m(4) * pkin(2) * t213 + (-pkin(1) * mrSges(3,1) - Ifges(3,4) * t137 + (Ifges(3,1) - Ifges(3,2)) * t140) * t226) * t137) * qJD(2)) * t132 + 0.2e1 * m(4) * (-t39 * t53 + t40 * t54) + (0.2e1 * Ifges(4,1) * t86 + Ifges(4,5) * t134 + t85 * t214 + t53 * t215) * t84 + t85 * t169 + 0.2e1 * (t39 * t86 - t40 * t85) * mrSges(4,3) + 0.2e1 * m(5) * (t11 * t192 + t12 * t28 + t39 * t42) + 0.2e1 * t192 * t36 + 0.2e1 * m(3) * (t90 * t93 - t91 * t92) + (t168 - 0.2e1 * t194 - 0.2e1 * t195 - 0.2e1 * t199 - 0.2e1 * t200) * t134 + t76 * t213 + (t1 * t10 + t2 * t9 + t23 * t8) * t217 + (mrSges(4,1) * t213 + t54 * t215 + t86 * t214 + Ifges(5,5) * t72 - Ifges(4,6) * t134 + Ifges(5,6) * t146 + ((2 * Ifges(4,2)) + Ifges(5,3)) * t85) * t83 + 0.2e1 * t39 * (-mrSges(5,1) * t146 + mrSges(5,2) * t72) - t146 * t3 + t146 * t25 + (t16 - t34) * t51 + 0.2e1 * t10 * t13 + 0.2e1 * t9 * t14 + t21 * t17 + t22 * t18 + 0.2e1 * t23 * t6 + 0.2e1 * t8 * t31 + 0.2e1 * t1 * t32 + 0.2e1 * t2 * t33 + 0.2e1 * t28 * t37 + 0.2e1 * t42 * t30 + t52 * t35 + t57 * t4 + t58 * t5 + 0.2e1 * t11 * t60 + 0.2e1 * t12 * t61 + t72 * t26; (-t12 * mrSges(5,3) + t4 * t206 + t5 * t204 + t26 / 0.2e1 + t39 * mrSges(5,2) + t198 / 0.2e1 + (t17 * t205 + t18 * t206) * qJD(5) + (-t192 * mrSges(5,3) + t16 / 0.2e1 - t34 / 0.2e1 - t196 / 0.2e1) * qJD(4) + (-qJD(4) * t60 - t224) * t124) * t136 - Ifges(3,6) * t166 + t168 + t85 * t128 / 0.2e1 + (t11 * mrSges(5,3) - t3 / 0.2e1 + t25 / 0.2e1 - t39 * mrSges(5,1) + t197 / 0.2e1 + (-t28 * mrSges(5,3) + t35 / 0.2e1 + t145) * qJD(4) + (m(6) * t185 + t142) * t124) * t139 + (m(4) * (t131 * t40 - t133 * t39) + (-t131 * t83 - t133 * t84) * mrSges(4,3)) * pkin(2) + t88 * t212 + t22 * t209 + t58 * t210 + t57 * t211 + (m(5) * t39 + t30) * t125 - (t62 / 0.2e1 - t101 / 0.2e1) * t146 + (t87 / 0.2e1 - t110 / 0.2e1) * t51 - t194 - t195 - t199 - t200 + m(6) * (t1 * t70 + t44 * t10 + t2 * t69 + t45 * t9) + t44 * t32 + t45 * t33 - t23 * t67 + t69 * t14 + t70 * t13 + t9 * t80 + t10 * t81 + t8 * t94 + t42 * t98 + t72 * t103 / 0.2e1 + t1 * t104 + t2 * t105 + t52 * t112 / 0.2e1; 0.2e1 * t45 * t105 + 0.2e1 * t69 * t80 + (t44 * t70 + t45 * t69) * t217 + 0.2e1 * t44 * t104 + 0.2e1 * t70 * t81 + 0.2e1 * t125 * t98 + (t101 - t62 + (0.2e1 * t124 * t94 - t135 * t88 + t138 * t89 + t112) * qJD(4)) * t139 + (-0.2e1 * t124 * t67 - t135 * t63 + t138 * t64 + t103 + (-t135 * t89 - t138 * t88) * qJD(5) + (t124 ^ 2 * t170 - t110 + t87) * qJD(4)) * t136; m(4) * t155 + t83 * mrSges(4,1) + t76 + ((-t135 * t33 + t138 * t32 + t60 + m(6) * (t10 * t138 - t135 * t9)) * qJD(4) + t224) * t139 + ((-t135 * t32 - t138 * t33) * qJD(5) + m(6) * (t141 + t185) + t142 + t149) * t136; t139 * t67 + (-t105 * t171 - t135 * t80 + m(6) * (-t135 * t45 + t138 * t44 - t171 * t69 - t173 * t70) - t104 * t173 + t138 * t81) * t136 + (t136 * t94 - t105 * t179 + m(6) * (t178 * t70 - t179 * t69 + (t136 ^ 2 - t139 ^ 2) * t124) + t104 * t178) * qJD(4); (-0.1e1 + t177) * t170 * t175; -t11 * mrSges(5,2) + t12 * mrSges(5,1) + t23 * t97 + t146 * t223 + t57 * t100 / 0.2e1 + t58 * t102 / 0.2e1 + t8 * t107 + t51 * t208 + t109 * t212 + t22 * t207 + t5 * t222 + t4 * t204 + t145 * qJD(5) + t167 * pkin(4) + ((-t10 * t135 - t138 * t9) * qJD(5) + t154) * mrSges(6,3) + (m(6) * t141 - t171 * t33 - t173 * t32 + t149) * pkin(9) + t169; pkin(4) * t67 + t128 + (t223 + (-m(6) * pkin(4) - mrSges(5,1) + t107) * t124 * qJD(4)) * t139 + (t174 * t207 + qJD(5) * t209 + t211 + t157 * mrSges(6,3) + (m(6) * t157 - qJD(5) * t105 + t81) * pkin(9)) * t138 + (-t109 * t174 / 0.2e1 - qJD(5) * t88 / 0.2e1 + t210 + t156 * mrSges(6,3) + (m(6) * t156 - qJD(5) * t104 - t80) * pkin(9)) * t135 + (t124 * t97 + t100 * t206 + t102 * t204 + (t109 * t205 + t111 * t206) * qJD(5) + (t124 * mrSges(5,2) - Ifges(5,6) + t208) * qJD(4)) * t136; -t139 * t97 + (t136 * t107 + m(6) * (t177 * t201 - t202) + t177 * t139 * mrSges(6,3) - t153) * qJD(4); -0.2e1 * pkin(4) * t97 + t100 * t138 + t102 * t135 + (-t109 * t135 + t111 * t138) * qJD(5); mrSges(6,1) * t2 - mrSges(6,2) * t1 + t3; mrSges(6,1) * t45 - mrSges(6,2) * t44 + t62; t67; t127 + (t107 * pkin(9) - t186) * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t19(1), t19(2), t19(4), t19(7), t19(11); t19(2), t19(3), t19(5), t19(8), t19(12); t19(4), t19(5), t19(6), t19(9), t19(13); t19(7), t19(8), t19(9), t19(10), t19(14); t19(11), t19(12), t19(13), t19(14), t19(15);];
Mq = res;
