% Calculate time derivative of joint inertia matrix for
% S5RPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR14_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR14_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR14_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RPRRR14_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR14_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR14_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR14_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:16:56
% EndTime: 2019-12-31 19:17:07
% DurationCPUTime: 3.86s
% Computational Cost: add. (6828->419), mult. (20093->648), div. (0->0), fcn. (20868->12), ass. (0->187)
t140 = sin(qJ(5));
t143 = cos(qJ(5));
t141 = sin(qJ(4));
t176 = qJD(5) * t141;
t144 = cos(qJ(4));
t178 = qJD(4) * t144;
t148 = -t140 * t176 + t143 * t178;
t134 = sin(pkin(11));
t136 = sin(pkin(5));
t137 = cos(pkin(11));
t139 = cos(pkin(5));
t142 = sin(qJ(3));
t138 = cos(pkin(6));
t145 = cos(qJ(3));
t184 = t138 * t145;
t135 = sin(pkin(6));
t187 = t135 * t145;
t229 = (-t134 * t142 + t137 * t184) * t136 + t139 * t187;
t228 = -2 * Ifges(4,4);
t175 = qJD(5) * t143;
t132 = Ifges(6,5) * t175;
t177 = qJD(5) * t140;
t227 = Ifges(6,6) * t177 / 0.2e1 - t132 / 0.2e1;
t226 = t140 / 0.2e1;
t210 = t143 / 0.2e1;
t209 = pkin(1) * t139;
t129 = t137 * t209;
t190 = t134 * t136;
t86 = pkin(2) * t139 + t129 + (-pkin(8) * t138 - qJ(2)) * t190;
t95 = (-pkin(8) * t134 * t135 - pkin(2) * t137 - pkin(1)) * t136;
t156 = t135 * t95 + t138 * t86;
t186 = t136 * t137;
t183 = qJ(2) * t186 + t134 * t209;
t83 = (t135 * t139 + t138 * t186) * pkin(8) + t183;
t48 = -t142 * t83 + t156 * t145;
t64 = -t135 * t86 + t138 * t95;
t185 = t138 * t142;
t188 = t135 * t142;
t85 = t139 * t188 + (t134 * t145 + t137 * t185) * t136;
t43 = -pkin(3) * t229 - pkin(9) * t85 + t64;
t102 = -t135 * t186 + t138 * t139;
t76 = t145 * t83;
t49 = t86 * t185 + t95 * t188 + t76;
t46 = pkin(9) * t102 + t49;
t203 = t141 * t43 + t144 * t46;
t224 = qJD(4) * t203;
t182 = qJD(2) * t136;
t39 = (-t134 * t185 + t137 * t145) * t182 + t48 * qJD(3);
t168 = t134 * t182;
t161 = t135 * t168;
t80 = t229 * qJD(3);
t81 = t85 * qJD(3);
t60 = pkin(3) * t81 - pkin(9) * t80 + t161;
t11 = -t141 * t39 + t144 * t60 - t224;
t118 = -mrSges(6,1) * t143 + mrSges(6,2) * t140;
t223 = -m(6) * pkin(4) - mrSges(5,1) + t118;
t222 = 0.2e1 * m(6);
t221 = 0.2e1 * pkin(9);
t220 = -2 * mrSges(4,3);
t219 = m(5) * pkin(3);
t153 = t102 * t144 - t141 * t85;
t54 = qJD(4) * t153 + t144 * t80;
t66 = t102 * t141 + t144 * t85;
t56 = -t140 * t229 + t143 * t66;
t26 = -qJD(5) * t56 - t140 * t54 + t143 * t81;
t218 = t26 / 0.2e1;
t198 = Ifges(6,4) * t140;
t121 = Ifges(6,2) * t143 + t198;
t197 = Ifges(6,4) * t143;
t157 = -Ifges(6,2) * t140 + t197;
t72 = -t121 * t176 + (Ifges(6,6) * t141 + t144 * t157) * qJD(4);
t217 = t72 / 0.2e1;
t123 = Ifges(6,1) * t140 + t197;
t158 = Ifges(6,1) * t143 - t198;
t73 = -t123 * t176 + (Ifges(6,5) * t141 + t144 * t158) * qJD(4);
t216 = t73 / 0.2e1;
t101 = -Ifges(6,5) * t144 + t141 * t158;
t215 = t101 / 0.2e1;
t214 = Ifges(6,5) * t226 + Ifges(6,6) * t210;
t213 = t123 / 0.2e1;
t212 = -t140 / 0.2e1;
t211 = -t143 / 0.2e1;
t208 = pkin(9) * t144;
t207 = t81 * Ifges(5,5);
t206 = t81 * Ifges(5,6);
t205 = t229 * Ifges(5,6);
t55 = -t140 * t66 - t143 * t229;
t31 = -mrSges(6,1) * t55 + mrSges(6,2) * t56;
t58 = -mrSges(5,1) * t229 - mrSges(5,3) * t66;
t204 = t31 - t58;
t202 = Ifges(4,5) * t80 - Ifges(4,6) * t81;
t201 = mrSges(6,3) * t141;
t200 = Ifges(5,4) * t141;
t199 = Ifges(5,4) * t144;
t196 = Ifges(6,6) * t140;
t103 = -t144 * t138 + t141 * t188;
t104 = t138 * t141 + t144 * t188;
t180 = qJD(3) * t145;
t166 = t135 * t180;
t89 = qJD(4) * t104 + t141 * t166;
t195 = t103 * t89;
t40 = (t134 * t184 + t137 * t142) * t182 + (t142 * t156 + t76) * qJD(3);
t194 = t145 * t40;
t88 = -qJD(4) * t103 + t144 * t166;
t193 = t88 * t144;
t192 = t89 * t141;
t191 = -mrSges(5,1) * t144 + mrSges(5,2) * t141 - mrSges(4,1);
t181 = qJD(3) * t142;
t179 = qJD(4) * t141;
t174 = qJD(5) * t144;
t27 = qJD(5) * t55 + t140 * t81 + t143 * t54;
t53 = qJD(4) * t66 + t141 * t80;
t3 = Ifges(6,5) * t27 + Ifges(6,6) * t26 + Ifges(6,3) * t53;
t173 = Ifges(5,5) * t54 - Ifges(5,6) * t53 + Ifges(5,3) * t81;
t12 = -mrSges(6,1) * t26 + mrSges(6,2) * t27;
t9 = -pkin(4) * t81 - t11;
t171 = -m(6) * t9 - t12;
t167 = t135 * t181;
t116 = (pkin(4) * t141 - pkin(10) * t144) * qJD(4);
t117 = -pkin(4) * t144 - pkin(10) * t141 - pkin(3);
t69 = t117 * t175 + t116 * t140 + (-t140 * t174 - t143 * t179) * pkin(9);
t96 = t117 * t143 - t140 * t208;
t163 = -qJD(5) * t96 + t69;
t70 = -t117 * t177 + t116 * t143 + (t140 * t179 - t143 * t174) * pkin(9);
t97 = t117 * t140 + t143 * t208;
t162 = -qJD(5) * t97 - t70;
t15 = t53 * pkin(4) - t54 * pkin(10) + t40;
t17 = -pkin(10) * t229 + t203;
t45 = -t102 * pkin(3) - t48;
t25 = -pkin(4) * t153 - t66 * pkin(10) + t45;
t6 = -t140 * t17 + t143 * t25;
t10 = t141 * t60 + t144 * t39 + t43 * t178 - t179 * t46;
t8 = pkin(10) * t81 + t10;
t1 = qJD(5) * t6 + t140 * t15 + t143 * t8;
t7 = t140 * t25 + t143 * t17;
t2 = -qJD(5) * t7 - t140 * t8 + t143 * t15;
t160 = t1 * t143 - t2 * t140;
t159 = mrSges(6,1) * t140 + mrSges(6,2) * t143;
t18 = -t141 * t46 + t144 * t43;
t21 = Ifges(6,4) * t56 + Ifges(6,2) * t55 - Ifges(6,6) * t153;
t22 = Ifges(6,1) * t56 + Ifges(6,4) * t55 - Ifges(6,5) * t153;
t152 = t21 * t212 + t210 * t22;
t151 = t103 * t178 + t192;
t90 = -t140 * t104 - t143 * t187;
t150 = -t143 * t104 + t140 * t187;
t147 = t140 * t178 + t141 * t175;
t71 = Ifges(6,5) * t148 - Ifges(6,6) * t147 + Ifges(6,3) * t179;
t133 = Ifges(5,5) * t178;
t124 = Ifges(5,1) * t141 + t199;
t122 = Ifges(5,2) * t144 + t200;
t115 = -mrSges(6,1) * t144 - t143 * t201;
t114 = mrSges(6,2) * t144 - t140 * t201;
t112 = (Ifges(5,1) * t144 - t200) * qJD(4);
t111 = t158 * qJD(5);
t110 = (-Ifges(5,2) * t141 + t199) * qJD(4);
t109 = t157 * qJD(5);
t107 = (mrSges(5,1) * t141 + mrSges(5,2) * t144) * qJD(4);
t106 = t159 * qJD(5);
t105 = t159 * t141;
t100 = -Ifges(6,6) * t144 + t141 * t157;
t99 = -Ifges(6,3) * t144 + (Ifges(6,5) * t143 - t196) * t141;
t94 = -mrSges(6,2) * t179 - mrSges(6,3) * t147;
t93 = mrSges(6,1) * t179 - mrSges(6,3) * t148;
t82 = mrSges(6,1) * t147 + mrSges(6,2) * t148;
t68 = mrSges(4,1) * t102 - mrSges(4,3) * t85;
t67 = -mrSges(4,2) * t102 + mrSges(4,3) * t229;
t63 = qJD(5) * t150 - t140 * t88 + t143 * t167;
t62 = qJD(5) * t90 + t140 * t167 + t143 * t88;
t61 = mrSges(4,1) * t81 + mrSges(4,2) * t80;
t57 = mrSges(5,2) * t229 + mrSges(5,3) * t153;
t47 = -mrSges(5,1) * t153 + mrSges(5,2) * t66;
t37 = mrSges(5,1) * t81 - mrSges(5,3) * t54;
t36 = -mrSges(5,2) * t81 - mrSges(5,3) * t53;
t35 = Ifges(5,1) * t66 + Ifges(5,4) * t153 - Ifges(5,5) * t229;
t34 = Ifges(5,4) * t66 + Ifges(5,2) * t153 - t205;
t33 = -mrSges(6,1) * t153 - mrSges(6,3) * t56;
t32 = mrSges(6,2) * t153 + mrSges(6,3) * t55;
t30 = mrSges(5,1) * t53 + mrSges(5,2) * t54;
t29 = Ifges(5,1) * t54 - Ifges(5,4) * t53 + t207;
t28 = Ifges(5,4) * t54 - Ifges(5,2) * t53 + t206;
t20 = Ifges(6,5) * t56 + Ifges(6,6) * t55 - Ifges(6,3) * t153;
t16 = pkin(4) * t229 - t18;
t14 = mrSges(6,1) * t53 - mrSges(6,3) * t27;
t13 = -mrSges(6,2) * t53 + mrSges(6,3) * t26;
t5 = Ifges(6,1) * t27 + Ifges(6,4) * t26 + Ifges(6,5) * t53;
t4 = Ifges(6,4) * t27 + Ifges(6,2) * t26 + Ifges(6,6) * t53;
t19 = [0.2e1 * m(4) * (t161 * t64 + t39 * t49 - t40 * t48) + t102 * t202 + (t1 * t7 + t16 * t9 + t2 * t6) * t222 - 0.2e1 * (mrSges(3,1) * t139 - mrSges(3,3) * t190) * t168 - (t3 - t28) * t153 + (0.2e1 * Ifges(4,1) * t85 + Ifges(4,5) * t102 + t220 * t48 - t228 * t229) * t80 + 0.2e1 * (-mrSges(4,1) * t229 + mrSges(4,2) * t85) * t161 + (t85 * t228 + Ifges(5,5) * t66 - Ifges(4,6) * t102 + Ifges(5,6) * t153 + t220 * t49 - ((2 * Ifges(4,2)) + Ifges(5,3)) * t229) * t81 - t229 * t173 + (t20 - t34) * t53 + 0.2e1 * t7 * t13 + 0.2e1 * t6 * t14 + 0.2e1 * t16 * t12 + t26 * t21 + t27 * t22 + 0.2e1 * t9 * t31 + 0.2e1 * t1 * t32 + 0.2e1 * t2 * t33 + 0.2e1 * t18 * t37 + 0.2e1 * t45 * t30 + 0.2e1 * t40 * t47 + t54 * t35 + t55 * t4 + t56 * t5 + 0.2e1 * t10 * t57 + 0.2e1 * t11 * t58 + 0.2e1 * t64 * t61 + t66 * t29 + 0.2e1 * t39 * t67 - 0.2e1 * t40 * t68 + 0.2e1 * (t137 * (-mrSges(3,2) * t139 + mrSges(3,3) * t186) + m(3) * (t183 * t137 + (qJ(2) * t190 - t129) * t134)) * t182 + 0.2e1 * m(5) * (t10 * t203 + t11 * t18 + t40 * t45) + 0.2e1 * t203 * t36; t104 * t36 - t150 * t13 + t138 * t61 + t90 * t14 + t62 * t32 + t63 * t33 + t88 * t57 + t204 * t89 + (t12 - t37) * t103 + m(6) * (-t1 * t150 + t103 * t9 + t16 * t89 + t2 * t90 + t6 * t63 + t62 * t7) + m(5) * (t104 * t10 - t103 * t11 - t89 * t18 + t203 * t88) + (-t145 * t30 + (-t142 * t81 - t145 * t80) * mrSges(4,3) + (t145 * t67 + (t47 - t68) * t142) * qJD(3) + m(5) * (t181 * t45 - t194) + m(4) * (t138 * t168 + t142 * t39 + t180 * t49 - t181 * t48 - t194)) * t135; 0.2e1 * m(5) * (-t135 ^ 2 * t142 * t180 + t104 * t88 + t195) + 0.2e1 * m(6) * (-t150 * t62 + t63 * t90 + t195); (t5 * t210 + t4 * t212 - t11 * mrSges(5,3) + t29 / 0.2e1 + t207 / 0.2e1 + (t21 * t211 + t212 * t22) * qJD(5) + (t20 / 0.2e1 - t34 / 0.2e1 - t203 * mrSges(5,3) + t205 / 0.2e1) * qJD(4) + (-qJD(4) * t57 - t37 + m(5) * (-t11 - t224) - t171) * pkin(9)) * t141 + (t99 / 0.2e1 - t122 / 0.2e1) * t53 + (t191 - t219) * t40 + t27 * t215 + t56 * t216 + t55 * t217 + t100 * t218 - (t71 / 0.2e1 - t110 / 0.2e1) * t153 - t229 * t133 / 0.2e1 + t202 + (t10 * mrSges(5,3) - t3 / 0.2e1 + t28 / 0.2e1 + t206 / 0.2e1 + (t35 / 0.2e1 - t18 * mrSges(5,3) + t152) * qJD(4) + (m(5) * t10 + t36 + (-m(5) * t18 + m(6) * t16 + t204) * qJD(4)) * pkin(9)) * t144 + m(6) * (t1 * t97 + t2 * t96 + t6 * t70 + t69 * t7) - pkin(3) * t30 - t39 * mrSges(4,2) + t69 * t32 + t70 * t33 + t16 * t82 + t6 * t93 + t7 * t94 + t96 * t14 + t97 * t13 + t9 * t105 + t45 * t107 + t66 * t112 / 0.2e1 + t1 * t114 + t2 * t115 + t54 * t124 / 0.2e1; t103 * t82 + t89 * t105 + t62 * t114 + t63 * t115 + t90 * t93 - t150 * t94 + (-t145 * t107 + (-mrSges(4,2) * t145 + t142 * t191) * qJD(3)) * t135 - t167 * t219 + m(6) * (-t150 * t69 + t62 * t97 + t63 * t96 + t70 * t90) + (m(5) * (-t104 * t179 + t151 + t193) / 0.2e1 + m(6) * t151 / 0.2e1) * t221 + (t192 + t193 + (t103 * t144 - t104 * t141) * qJD(4)) * mrSges(5,3); 0.2e1 * t69 * t114 + 0.2e1 * t97 * t94 + 0.2e1 * t70 * t115 + 0.2e1 * t96 * t93 - 0.2e1 * pkin(3) * t107 + (t69 * t97 + t70 * t96) * t222 + (t110 - t71 + (-t100 * t140 + t101 * t143 + t105 * t221 + t124) * qJD(4)) * t144 + (t82 * t221 - t140 * t72 + t143 * t73 + t112 + (-t100 * t143 - t101 * t140) * qJD(5) + (pkin(9) ^ 2 * t144 * t222 - t122 + t99) * qJD(4)) * t141; -t10 * mrSges(5,2) + t11 * mrSges(5,1) + t16 * t106 + t153 * t227 + t55 * t109 / 0.2e1 + t56 * t111 / 0.2e1 + t9 * t118 + t53 * t214 + t121 * t218 + t27 * t213 + t5 * t226 + t4 * t210 + t152 * qJD(5) + t171 * pkin(4) + ((-t140 * t7 - t143 * t6) * qJD(5) + t160) * mrSges(6,3) + (m(6) * (-t175 * t6 - t177 * t7 + t160) + t143 * t13 - t140 * t14 - t33 * t175 - t32 * t177) * pkin(10) + t173; -t88 * mrSges(5,2) + t103 * t106 + (m(6) * pkin(10) + mrSges(6,3)) * (-t63 * t140 + t62 * t143 + (t140 * t150 - t143 * t90) * qJD(5)) + t223 * t89; -pkin(4) * t82 + t133 + (t223 * qJD(4) * pkin(9) + t227) * t144 + (qJD(5) * t215 + t178 * t213 + t217 + t163 * mrSges(6,3) + (m(6) * t163 - qJD(5) * t115 + t94) * pkin(10)) * t143 + (-qJD(5) * t100 / 0.2e1 - t121 * t178 / 0.2e1 + t216 + t162 * mrSges(6,3) + (m(6) * t162 - qJD(5) * t114 - t93) * pkin(10)) * t140 + (t111 * t210 + t109 * t212 + pkin(9) * t106 + (t121 * t211 + t123 * t212) * qJD(5) + (pkin(9) * mrSges(5,2) - Ifges(5,6) + t214) * qJD(4)) * t141; -0.2e1 * pkin(4) * t106 + t109 * t143 + t111 * t140 + (-t121 * t140 + t123 * t143) * qJD(5); mrSges(6,1) * t2 - mrSges(6,2) * t1 + t3; mrSges(6,1) * t63 - mrSges(6,2) * t62; mrSges(6,1) * t70 - mrSges(6,2) * t69 + t71; t132 + (pkin(10) * t118 - t196) * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t19(1), t19(2), t19(4), t19(7), t19(11); t19(2), t19(3), t19(5), t19(8), t19(12); t19(4), t19(5), t19(6), t19(9), t19(13); t19(7), t19(8), t19(9), t19(10), t19(14); t19(11), t19(12), t19(13), t19(14), t19(15);];
Mq = res;
