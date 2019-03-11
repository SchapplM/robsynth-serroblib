% Calculate time derivative of joint inertia matrix for
% S6RRPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3,theta4]';
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
% Datum: 2019-03-09 08:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRR2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR2_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR2_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR2_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR2_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR2_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR2_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:49:53
% EndTime: 2019-03-09 08:49:58
% DurationCPUTime: 2.88s
% Computational Cost: add. (8326->400), mult. (18269->607), div. (0->0), fcn. (18556->10), ass. (0->167)
t155 = sin(pkin(11));
t156 = sin(pkin(10));
t199 = pkin(2) * t156;
t149 = qJ(4) + t199;
t196 = pkin(8) + t149;
t137 = t196 * t155;
t157 = cos(pkin(11));
t138 = t196 * t157;
t160 = sin(qJ(5));
t163 = cos(qJ(5));
t102 = -t160 * t137 + t163 * t138;
t181 = t157 * t163;
t167 = t155 * t160 - t181;
t135 = t167 * qJD(5);
t158 = cos(pkin(10));
t161 = sin(qJ(2));
t164 = cos(qJ(2));
t142 = t156 * t164 + t158 * t161;
t133 = t142 * qJD(2);
t140 = t156 * t161 - t158 * t164;
t134 = t140 * qJD(2);
t143 = t155 * t163 + t157 * t160;
t136 = t143 * qJD(5);
t61 = t134 * t167 - t136 * t142;
t62 = t134 * t143 + t135 * t142;
t212 = Ifges(6,5) * t61 + Ifges(6,6) * t62 + Ifges(6,3) * t133;
t211 = 2 * m(5);
t210 = 2 * m(6);
t209 = 2 * m(7);
t154 = t157 ^ 2;
t152 = -pkin(2) * t164 - pkin(1);
t207 = 0.2e1 * t152;
t206 = m(4) * pkin(2);
t205 = m(7) * pkin(5);
t159 = sin(qJ(6));
t162 = cos(qJ(6));
t107 = -t143 * t159 - t162 * t167;
t204 = t107 / 0.2e1;
t108 = t143 * t162 - t159 * t167;
t203 = t108 / 0.2e1;
t201 = -t167 / 0.2e1;
t200 = t143 / 0.2e1;
t198 = pkin(2) * t158;
t197 = pkin(5) * t136;
t195 = -qJ(3) - pkin(7);
t182 = t142 * t157;
t106 = t140 * pkin(3) - t142 * qJ(4) + t152;
t146 = t195 * t161;
t147 = t195 * t164;
t115 = t146 * t156 - t147 * t158;
t76 = t157 * t106 - t115 * t155;
t53 = pkin(4) * t140 - pkin(8) * t182 + t76;
t183 = t142 * t155;
t77 = t155 * t106 + t157 * t115;
t60 = -pkin(8) * t183 + t77;
t26 = t160 * t53 + t163 * t60;
t69 = qJD(6) * t107 - t135 * t162 - t136 * t159;
t70 = -qJD(6) * t108 + t135 * t159 - t136 * t162;
t194 = Ifges(7,5) * t69 + Ifges(7,6) * t70;
t173 = pkin(2) * qJD(2) * t161;
t83 = pkin(3) * t133 + qJ(4) * t134 - qJD(4) * t142 + t173;
t172 = qJD(2) * t195;
t132 = qJD(3) * t164 + t161 * t172;
t165 = -t161 * qJD(3) + t164 * t172;
t94 = t158 * t132 + t156 * t165;
t48 = t155 * t83 + t157 * t94;
t193 = Ifges(5,4) * t155;
t192 = Ifges(5,4) * t157;
t114 = -t158 * t146 - t147 * t156;
t93 = t132 * t156 - t158 * t165;
t191 = t114 * t93;
t190 = t133 * Ifges(5,5);
t189 = t133 * Ifges(5,6);
t188 = t136 * mrSges(6,1);
t187 = t155 * Ifges(5,2);
t186 = t160 * t60;
t185 = t134 * t155;
t184 = t134 * t157;
t92 = -mrSges(5,1) * t185 - mrSges(5,2) * t184;
t180 = -Ifges(6,5) * t135 - Ifges(6,6) * t136;
t178 = qJD(5) * t163;
t177 = qJD(6) * t159;
t176 = qJD(6) * t162;
t175 = 0.2e1 * t164;
t95 = t143 * t142;
t96 = t167 * t142;
t63 = t159 * t96 - t162 * t95;
t19 = qJD(6) * t63 + t159 * t62 + t162 * t61;
t64 = -t159 * t95 - t162 * t96;
t20 = -qJD(6) * t64 - t159 * t61 + t162 * t62;
t174 = Ifges(7,5) * t19 + Ifges(7,6) * t20 + Ifges(7,3) * t133;
t151 = -pkin(3) - t198;
t31 = -t62 * mrSges(6,1) + t61 * mrSges(6,2);
t8 = -t20 * mrSges(7,1) + t19 * mrSges(7,2);
t33 = -t70 * mrSges(7,1) + t69 * mrSges(7,2);
t47 = -t155 * t94 + t157 * t83;
t25 = t163 * t53 - t186;
t171 = t133 * mrSges(4,1) - t134 * mrSges(4,2);
t101 = -t163 * t137 - t138 * t160;
t127 = t135 * mrSges(6,2);
t170 = -t127 + t33;
t78 = -pkin(4) * t185 + t93;
t91 = pkin(4) * t183 + t114;
t169 = Ifges(5,5) * t157 - Ifges(5,6) * t155;
t23 = pkin(5) * t140 + pkin(9) * t96 + t25;
t24 = -pkin(9) * t95 + t26;
t9 = -t159 * t24 + t162 * t23;
t10 = t159 * t23 + t162 * t24;
t88 = -pkin(9) * t143 + t101;
t89 = -pkin(9) * t167 + t102;
t42 = -t159 * t89 + t162 * t88;
t43 = t159 * t88 + t162 * t89;
t86 = -t137 * t178 + qJD(4) * t181 + (-qJD(4) * t155 - qJD(5) * t138) * t160;
t71 = -pkin(9) * t136 + t86;
t87 = -t143 * qJD(4) - qJD(5) * t102;
t72 = pkin(9) * t135 + t87;
t21 = qJD(6) * t42 + t159 * t72 + t162 * t71;
t22 = -qJD(6) * t43 - t159 * t71 + t162 * t72;
t168 = t22 * mrSges(7,1) - t21 * mrSges(7,2) + t194;
t145 = -pkin(4) * t157 + t151;
t39 = pkin(4) * t133 + pkin(8) * t184 + t47;
t41 = pkin(8) * t185 + t48;
t12 = -qJD(5) * t26 - t160 * t41 + t163 * t39;
t4 = pkin(5) * t133 - pkin(9) * t61 + t12;
t11 = -qJD(5) * t186 + t160 * t39 + t163 * t41 + t53 * t178;
t5 = pkin(9) * t62 + t11;
t2 = qJD(6) * t9 + t159 * t4 + t162 * t5;
t3 = -qJD(6) * t10 - t159 * t5 + t162 * t4;
t166 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + t174;
t139 = (-mrSges(7,1) * t159 - mrSges(7,2) * t162) * qJD(6) * pkin(5);
t116 = pkin(5) * t167 + t145;
t112 = Ifges(6,1) * t143 - Ifges(6,4) * t167;
t111 = Ifges(6,4) * t143 - Ifges(6,2) * t167;
t110 = mrSges(5,1) * t140 - mrSges(5,3) * t182;
t109 = -mrSges(5,2) * t140 - mrSges(5,3) * t183;
t105 = -Ifges(6,1) * t135 - Ifges(6,4) * t136;
t104 = -Ifges(6,4) * t135 - Ifges(6,2) * t136;
t103 = -t127 + t188;
t100 = mrSges(5,1) * t133 + mrSges(5,3) * t184;
t99 = -mrSges(5,2) * t133 + mrSges(5,3) * t185;
t85 = mrSges(6,1) * t140 + mrSges(6,3) * t96;
t84 = -mrSges(6,2) * t140 - mrSges(6,3) * t95;
t82 = t190 - (t157 * Ifges(5,1) - t193) * t134;
t81 = t189 - (-t187 + t192) * t134;
t75 = Ifges(7,1) * t108 + Ifges(7,4) * t107;
t74 = Ifges(7,4) * t108 + Ifges(7,2) * t107;
t73 = -mrSges(7,1) * t107 + mrSges(7,2) * t108;
t65 = pkin(5) * t95 + t91;
t55 = -Ifges(6,1) * t96 - Ifges(6,4) * t95 + Ifges(6,5) * t140;
t54 = -Ifges(6,4) * t96 - Ifges(6,2) * t95 + Ifges(6,6) * t140;
t52 = mrSges(7,1) * t140 - mrSges(7,3) * t64;
t51 = -mrSges(7,2) * t140 + mrSges(7,3) * t63;
t46 = -mrSges(6,2) * t133 + mrSges(6,3) * t62;
t45 = mrSges(6,1) * t133 - mrSges(6,3) * t61;
t36 = -pkin(5) * t62 + t78;
t35 = Ifges(7,1) * t69 + Ifges(7,4) * t70;
t34 = Ifges(7,4) * t69 + Ifges(7,2) * t70;
t32 = -mrSges(7,1) * t63 + mrSges(7,2) * t64;
t30 = Ifges(7,1) * t64 + Ifges(7,4) * t63 + Ifges(7,5) * t140;
t29 = Ifges(7,4) * t64 + Ifges(7,2) * t63 + Ifges(7,6) * t140;
t28 = Ifges(6,1) * t61 + Ifges(6,4) * t62 + t133 * Ifges(6,5);
t27 = Ifges(6,4) * t61 + Ifges(6,2) * t62 + t133 * Ifges(6,6);
t14 = -mrSges(7,2) * t133 + mrSges(7,3) * t20;
t13 = mrSges(7,1) * t133 - mrSges(7,3) * t19;
t7 = Ifges(7,1) * t19 + Ifges(7,4) * t20 + t133 * Ifges(7,5);
t6 = Ifges(7,4) * t19 + Ifges(7,2) * t20 + t133 * Ifges(7,6);
t1 = [(-0.2e1 * t94 * mrSges(4,3) - 0.2e1 * (-Ifges(4,4) + t169) * t134 + ((2 * Ifges(4,2)) + (2 * Ifges(5,3)) + Ifges(6,3) + Ifges(7,3)) * t133 + t174 + t212) * t140 + 0.2e1 * (-t114 * t134 - t115 * t133) * mrSges(4,3) + (t157 * t82 - t155 * t81 + (-0.2e1 * Ifges(4,4) + t169) * t133 - (Ifges(5,1) * t154 + (2 * Ifges(4,1)) + (t187 - 0.2e1 * t192) * t155) * t134 + 0.2e1 * (mrSges(5,1) * t155 + mrSges(5,2) * t157 + mrSges(4,3)) * t93) * t142 + t133 * (-Ifges(6,5) * t96 - Ifges(6,6) * t95) + 0.2e1 * t78 * (mrSges(6,1) * t95 - mrSges(6,2) * t96) + t133 * (Ifges(7,5) * t64 + Ifges(7,6) * t63) + ((-mrSges(3,2) * pkin(1) + Ifges(3,4) * t164) * t175 + (t206 * t207 - 0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * pkin(2) * (mrSges(4,1) * t140 + mrSges(4,2) * t142) - 0.2e1 * Ifges(3,4) * t161 + (Ifges(3,1) - Ifges(3,2)) * t175) * t161) * qJD(2) + 0.2e1 * t48 * t109 + 0.2e1 * t47 * t110 + 0.2e1 * t114 * t92 + t171 * t207 + (t10 * t2 + t3 * t9 + t36 * t65) * t209 + (t11 * t26 + t12 * t25 + t78 * t91) * t210 + (t47 * t76 + t48 * t77 + t191) * t211 + 0.2e1 * t9 * t13 + 0.2e1 * t10 * t14 + t20 * t29 + t19 * t30 + 0.2e1 * t36 * t32 + 0.2e1 * t25 * t45 + 0.2e1 * t26 * t46 + 0.2e1 * t2 * t51 + 0.2e1 * t3 * t52 + t61 * t55 + t62 * t54 + t63 * t6 + t64 * t7 + 0.2e1 * t65 * t8 + 0.2e1 * t11 * t84 + 0.2e1 * t12 * t85 + 0.2e1 * m(4) * (t115 * t94 + t191) + 0.2e1 * t91 * t31 - t95 * t27 - t96 * t28 + 0.2e1 * t77 * t99 + 0.2e1 * t76 * t100; (Ifges(3,5) * t164 - Ifges(3,6) * t161 + (-mrSges(3,1) * t164 + mrSges(3,2) * t161) * pkin(7)) * qJD(2) - (-pkin(5) * t32 + t54 / 0.2e1) * t136 + m(7) * (t10 * t21 + t116 * t36 + t197 * t65 + t2 * t43 + t22 * t9 + t3 * t42) + (t194 + t180) * t140 / 0.2e1 + (-t11 * t167 - t12 * t143 + t135 * t25 - t136 * t26) * mrSges(6,3) + t78 * (mrSges(6,1) * t167 + mrSges(6,2) * t143) + m(6) * (t101 * t12 + t102 * t11 + t145 * t78 + t25 * t87 + t26 * t86) + m(5) * (t151 * t93 + (-t155 * t47 + t157 * t48) * t149 + (-t155 * t76 + t157 * t77) * qJD(4)) + t151 * t92 + t145 * t31 + (-mrSges(4,3) * t199 + Ifges(6,5) * t200 + Ifges(7,5) * t203 + Ifges(6,6) * t201 + Ifges(7,6) * t204 - Ifges(4,6)) * t133 - t135 * t55 / 0.2e1 - (t157 * (Ifges(5,1) * t155 + t192) / 0.2e1 - t155 * (Ifges(5,2) * t157 + t193) / 0.2e1 + Ifges(4,5) - mrSges(4,3) * t198) * t134 + t62 * t111 / 0.2e1 + t61 * t112 / 0.2e1 + t116 * t8 + t28 * t200 + t27 * t201 + t7 * t203 + t6 * t204 + (t156 * t94 - t158 * t93) * t206 + (t10 * t70 + t107 * t2 - t108 * t3 - t69 * t9) * mrSges(7,3) + t42 * t13 + t43 * t14 + t21 * t51 + t22 * t52 + t63 * t34 / 0.2e1 + t64 * t35 / 0.2e1 + t65 * t33 + t69 * t30 / 0.2e1 + t70 * t29 / 0.2e1 + t36 * t73 + t20 * t74 / 0.2e1 + t19 * t75 / 0.2e1 + t86 * t84 + t87 * t85 + (t149 * t99 + qJD(4) * t109 + t48 * mrSges(5,3) + t189 / 0.2e1 - t93 * mrSges(5,1) + t81 / 0.2e1) * t157 + (-t149 * t100 - qJD(4) * t110 - t47 * mrSges(5,3) + t82 / 0.2e1 + t190 / 0.2e1 + t93 * mrSges(5,2)) * t155 - t93 * mrSges(4,1) - t94 * mrSges(4,2) + t101 * t45 + t102 * t46 + t91 * t103 - t95 * t104 / 0.2e1 - t96 * t105 / 0.2e1; 0.2e1 * t145 * t103 - t167 * t104 + t143 * t105 + t107 * t34 + t108 * t35 - t135 * t112 + 0.2e1 * t116 * t33 + t69 * t75 + t70 * t74 - (-0.2e1 * pkin(5) * t73 + t111) * t136 + (t116 * t197 + t21 * t43 + t22 * t42) * t209 + (t101 * t87 + t102 * t86) * t210 + 0.2e1 * (t107 * t21 - t108 * t22 - t42 * t69 + t43 * t70) * mrSges(7,3) + 0.2e1 * (t101 * t135 - t102 * t136 - t143 * t87 - t167 * t86) * mrSges(6,3) + (t149 * t211 + 0.2e1 * mrSges(5,3)) * qJD(4) * (t155 ^ 2 + t154); m(4) * t173 + t157 * t100 + t107 * t13 + t108 * t14 - t135 * t84 - t136 * t85 - t167 * t45 + t143 * t46 + t155 * t99 + t69 * t51 + t70 * t52 + m(7) * (t10 * t69 + t107 * t3 + t108 * t2 + t70 * t9) + m(6) * (t11 * t143 - t12 * t167 - t135 * t26 - t136 * t25) + m(5) * (t155 * t48 + t157 * t47) + t171; m(7) * (t107 * t22 + t108 * t21 + t42 * t70 + t43 * t69) + m(6) * (-t101 * t136 - t102 * t135 + t143 * t86 - t167 * t87); 0.2e1 * m(6) * (-t135 * t143 + t136 * t167) + 0.2e1 * m(7) * (t107 * t70 + t108 * t69); m(5) * t93 + m(6) * t78 + m(7) * t36 + t31 + t8 + t92; -(-mrSges(6,1) - t205) * t136 + t170; 0; 0; t12 * mrSges(6,1) - t11 * mrSges(6,2) + (-t52 * t177 + t162 * t13 + m(7) * (t10 * t176 + t159 * t2 + t162 * t3 - t177 * t9) + t51 * t176 + t159 * t14) * pkin(5) + t166 + t212; t87 * mrSges(6,1) - t86 * mrSges(6,2) + (m(7) * (t159 * t21 + t162 * t22 + (-t159 * t42 + t162 * t43) * qJD(6)) + (t159 * t70 - t162 * t69 + (t107 * t162 + t108 * t159) * qJD(6)) * mrSges(7,3)) * pkin(5) + t168 + t180; -t188 + (t159 * t69 + t162 * t70 + (-t107 * t159 + t108 * t162) * qJD(6)) * t205 - t170; 0; 0.2e1 * t139; t166; t168; -t33; 0; t139; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
