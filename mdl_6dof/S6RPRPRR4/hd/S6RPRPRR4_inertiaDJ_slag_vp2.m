% Calculate time derivative of joint inertia matrix for
% S6RPRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
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
% Datum: 2019-03-09 03:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRR4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR4_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR4_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR4_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR4_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR4_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR4_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:44:29
% EndTime: 2019-03-09 03:44:34
% DurationCPUTime: 2.17s
% Computational Cost: add. (2811->339), mult. (5862->500), div. (0->0), fcn. (4706->8), ass. (0->157)
t95 = sin(pkin(10)) * pkin(1) + pkin(7);
t172 = pkin(4) + t95;
t107 = sin(qJ(5));
t108 = sin(qJ(3));
t155 = qJD(3) * t108;
t140 = t107 * t155;
t110 = cos(qJ(5));
t111 = cos(qJ(3));
t150 = qJD(5) * t111;
t141 = t110 * t150;
t116 = t140 - t141;
t142 = t107 * t150;
t154 = qJD(3) * t110;
t117 = t108 * t154 + t142;
t35 = -t117 * mrSges(6,1) + t116 * mrSges(6,2);
t106 = sin(qJ(6));
t109 = cos(qJ(6));
t158 = t109 * t110;
t120 = t106 * t107 - t158;
t121 = t106 * t110 + t109 * t107;
t187 = qJD(5) + qJD(6);
t184 = t187 * t111;
t24 = t120 * t184 + t121 * t155;
t25 = -t120 * t155 + t121 * t184;
t6 = -t25 * mrSges(7,1) + t24 * mrSges(7,2);
t194 = -t35 - t6;
t193 = 2 * qJ(4);
t143 = -cos(pkin(10)) * pkin(1) - pkin(2);
t192 = 0.2e1 * t143;
t191 = m(5) * t95;
t190 = -mrSges(4,1) + mrSges(5,2);
t189 = mrSges(4,2) - mrSges(5,3);
t151 = qJD(5) * t110;
t152 = qJD(5) * t107;
t147 = t108 * qJD(4);
t134 = pkin(3) * t155 - t147;
t52 = (pkin(8) * t108 - qJ(4) * t111) * qJD(3) + t134;
t112 = -pkin(3) - pkin(8);
t159 = qJ(4) * t108;
t119 = t143 - t159;
t57 = t112 * t111 + t119;
t153 = qJD(3) * t111;
t66 = t172 * t153;
t70 = t172 * t108;
t11 = t107 * t66 + t110 * t52 + t70 * t151 - t57 * t152;
t10 = t117 * pkin(9) + t11;
t135 = -t107 * t52 + t110 * t66;
t137 = pkin(9) * t111 - t57;
t63 = t107 * t70;
t7 = (-pkin(9) * t107 * t108 + pkin(5) * t111) * qJD(3) + (t137 * t110 - t63) * qJD(5) + t135;
t64 = t110 * t70;
t26 = pkin(5) * t108 + t137 * t107 + t64;
t157 = t110 * t111;
t31 = t110 * t57 + t63;
t27 = -pkin(9) * t157 + t31;
t8 = -t106 * t27 + t109 * t26;
t2 = t8 * qJD(6) + t10 * t109 + t106 * t7;
t9 = t106 * t26 + t109 * t27;
t3 = -t9 * qJD(6) - t10 * t106 + t109 * t7;
t188 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t12 = -t31 * qJD(5) + t135;
t126 = t11 * t107 + t12 * t110;
t80 = mrSges(6,3) * t107 * t111 + mrSges(6,1) * t108;
t81 = -mrSges(6,2) * t108 - mrSges(6,3) * t157;
t123 = t107 * t80 - t110 * t81;
t50 = -mrSges(6,2) * t153 + t117 * mrSges(6,3);
t51 = mrSges(6,1) * t153 - t116 * mrSges(6,3);
t186 = t123 * qJD(5) - t107 * t50 - t110 * t51;
t149 = qJD(6) * t106;
t40 = -t106 * t152 - t107 * t149 + t158 * t187;
t41 = t187 * t121;
t185 = (-t106 * t120 - t109 * t121) * qJD(6) - t106 * t40 + t109 * t41;
t183 = 2 * m(7);
t182 = 2 * mrSges(5,3);
t171 = pkin(3) * t111;
t181 = -0.2e1 * t119 + 0.2e1 * t171;
t180 = m(7) * pkin(5);
t179 = -t121 / 0.2e1;
t178 = -t120 / 0.2e1;
t67 = qJ(4) * t153 - t134;
t177 = m(5) * t67;
t175 = -t107 / 0.2e1;
t174 = t108 / 0.2e1;
t173 = -t110 / 0.2e1;
t170 = pkin(9) - t112;
t169 = -Ifges(7,5) * t41 - Ifges(7,6) * t40;
t168 = Ifges(6,4) * t107;
t167 = Ifges(6,4) * t110;
t131 = Ifges(6,1) * t107 + t167;
t164 = t108 * Ifges(6,5);
t59 = -t131 * t111 + t164;
t166 = t107 * t59;
t86 = Ifges(6,1) * t110 - t168;
t165 = t107 * t86;
t130 = Ifges(6,2) * t110 + t168;
t58 = t108 * Ifges(6,6) - t130 * t111;
t162 = t110 * t58;
t85 = -Ifges(6,2) * t107 + t167;
t161 = t110 * t85;
t71 = t172 * t111;
t156 = t107 ^ 2 + t110 ^ 2;
t148 = qJD(6) * t109;
t146 = 2 * mrSges(7,3);
t145 = Ifges(7,5) * t24 + Ifges(7,6) * t25 + Ifges(7,3) * t153;
t144 = mrSges(5,1) + t191;
t93 = t108 * t153;
t138 = -t41 * mrSges(7,1) - t40 * mrSges(7,2);
t136 = m(6) * t156;
t83 = t170 * t110;
t133 = -t120 * t41 - t121 * t40;
t132 = mrSges(6,1) * t110 - mrSges(6,2) * t107;
t84 = mrSges(6,1) * t107 + mrSges(6,2) * t110;
t129 = -Ifges(6,5) * t107 - Ifges(6,6) * t110;
t82 = t170 * t107;
t46 = -t106 * t83 - t109 * t82;
t45 = t106 * t82 - t109 * t83;
t30 = -t107 * t57 + t64;
t125 = t30 * t107 - t31 * t110;
t73 = t170 * t152;
t74 = qJD(5) * t83;
t16 = t45 * qJD(6) + t106 * t73 - t109 * t74;
t17 = -t46 * qJD(6) + t106 * t74 + t109 * t73;
t122 = t17 * mrSges(7,1) - t16 * mrSges(7,2) + t169;
t118 = Ifges(6,5) * t140 + t117 * Ifges(6,6) + Ifges(6,3) * t153 + t145;
t115 = t120 * t3 - t121 * t2 - t9 * t40 + t8 * t41;
t114 = t120 * t17 - t121 * t16 - t40 * t46 + t41 * t45;
t61 = t120 * t111;
t62 = t121 * t111;
t113 = t120 * t25 - t121 * t24 + t40 * t62 + t41 * t61;
t96 = pkin(5) * t107 + qJ(4);
t91 = pkin(5) * t151 + qJD(4);
t79 = t131 * qJD(5);
t78 = t130 * qJD(5);
t77 = t132 * qJD(5);
t72 = (-mrSges(7,1) * t106 - mrSges(7,2) * t109) * qJD(6) * pkin(5);
t69 = t132 * t111;
t65 = t172 * t155;
t54 = pkin(5) * t157 + t71;
t49 = mrSges(7,1) * t108 + mrSges(7,3) * t62;
t48 = -mrSges(7,2) * t108 + mrSges(7,3) * t61;
t44 = -Ifges(7,1) * t120 - Ifges(7,4) * t121;
t43 = -Ifges(7,4) * t120 - Ifges(7,2) * t121;
t42 = mrSges(7,1) * t121 - mrSges(7,2) * t120;
t36 = -pkin(5) * t142 + (-pkin(5) * t110 - t172) * t155;
t34 = -mrSges(7,1) * t61 - mrSges(7,2) * t62;
t33 = -t86 * t150 + (t111 * Ifges(6,5) + t131 * t108) * qJD(3);
t32 = -t85 * t150 + (t111 * Ifges(6,6) + t130 * t108) * qJD(3);
t29 = -Ifges(7,1) * t62 + Ifges(7,4) * t61 + Ifges(7,5) * t108;
t28 = -Ifges(7,4) * t62 + Ifges(7,2) * t61 + Ifges(7,6) * t108;
t20 = -Ifges(7,1) * t41 - Ifges(7,4) * t40;
t19 = -Ifges(7,4) * t41 - Ifges(7,2) * t40;
t18 = mrSges(7,1) * t40 - mrSges(7,2) * t41;
t14 = -mrSges(7,2) * t153 + mrSges(7,3) * t25;
t13 = mrSges(7,1) * t153 - mrSges(7,3) * t24;
t5 = Ifges(7,1) * t24 + Ifges(7,4) * t25 + Ifges(7,5) * t153;
t4 = Ifges(7,4) * t24 + Ifges(7,2) * t25 + Ifges(7,6) * t153;
t1 = [0.2e1 * t8 * t13 + 0.2e1 * t9 * t14 + t177 * t181 + t25 * t28 + t24 * t29 + 0.2e1 * t36 * t34 + 0.2e1 * t2 * t48 + 0.2e1 * t3 * t49 + 0.2e1 * t31 * t50 + 0.2e1 * t30 * t51 + 0.2e1 * t54 * t6 + t61 * t4 - t62 * t5 - 0.2e1 * t65 * t69 + 0.2e1 * t71 * t35 + 0.2e1 * t12 * t80 + 0.2e1 * t11 * t81 + (t2 * t9 + t3 * t8 + t36 * t54) * t183 + 0.2e1 * m(6) * (t11 * t31 + t12 * t30 - t65 * t71) + (t67 * t182 + t118) * t108 + (-0.2e1 * t67 * mrSges(5,2) - t107 * t33 - t110 * t32 + (t107 * t58 + (-t59 - t164) * t110) * qJD(5)) * t111 + ((t166 + mrSges(5,2) * t181 + mrSges(4,1) * t192 + t162 + 0.2e1 * (-Ifges(5,6) - Ifges(4,4)) * t108) * t108 + (-Ifges(7,5) * t62 + Ifges(7,6) * t61 + mrSges(5,3) * t181 + mrSges(4,2) * t192 + (0.2e1 * Ifges(4,4) + 0.2e1 * Ifges(5,6) + t129) * t111 + (-(2 * Ifges(5,3)) - (2 * Ifges(4,2)) + (2 * Ifges(5,2)) + (2 * Ifges(4,1)) + Ifges(7,3) + Ifges(6,3)) * t108) * t111) * qJD(3); t24 * t48 - t62 * t14 + t25 * t49 + t61 * t13 + m(7) * (-t2 * t62 + t24 * t9 + t25 * t8 + t3 * t61) + ((t107 * t81 + t110 * t80) * qJD(3) + m(7) * t36 + m(6) * (qJD(3) * t107 * t31 + t30 * t154 - t65) - t194) * t108 + (m(6) * (-t31 * t151 + t30 * t152 - t126) + (m(6) * t71 + m(7) * t54 + t34 + t69) * qJD(3) + t186) * t111; 0.2e1 * m(7) * (-t24 * t62 + t25 * t61 + t93) + 0.2e1 * m(6) * (-t156 + 0.1e1) * t93; t115 * mrSges(7,3) + (t144 * qJD(4) - t78 * t173 - t79 * t175) * t111 + (-t12 * mrSges(6,3) + t112 * t51 + t33 / 0.2e1) * t110 + (-t11 * mrSges(6,3) + t112 * t50 - t32 / 0.2e1) * t107 + ((-t159 - t171) * t191 + (Ifges(7,5) * t178 + Ifges(7,6) * t179 + Ifges(6,5) * t110 / 0.2e1 + Ifges(6,6) * t175 - Ifges(5,4) + Ifges(4,5) - pkin(3) * mrSges(5,1) + t190 * t95) * t111 + (Ifges(5,5) - Ifges(4,6) - qJ(4) * mrSges(5,1) + t165 / 0.2e1 + t161 / 0.2e1 + t189 * t95) * t108) * qJD(3) + m(6) * (-qJ(4) * t65 + qJD(4) * t71 + t112 * t126) + qJ(4) * t35 - t40 * t28 / 0.2e1 - t41 * t29 / 0.2e1 + t36 * t42 + t25 * t43 / 0.2e1 + t24 * t44 / 0.2e1 + t45 * t13 + t46 * t14 + t16 * t48 + t17 * t49 + t54 * t18 + t61 * t19 / 0.2e1 - t62 * t20 / 0.2e1 + qJD(4) * t69 + t71 * t77 - t65 * t84 + t91 * t34 + t96 * t6 + t169 * t174 + t5 * t178 + t4 * t179 + (-t166 / 0.2e1 - t162 / 0.2e1 + t129 * t174 + (t107 * t85 / 0.2e1 + t86 * t173) * t111 + t125 * mrSges(6,3) + (-m(6) * t125 - t123) * t112) * qJD(5) + m(7) * (t16 * t9 + t17 * t8 + t2 * t46 + t3 * t45 + t36 * t96 + t54 * t91); (t18 + t77) * t108 + m(7) * (t108 * t91 - t16 * t62 + t17 * t61 + t24 * t46 + t25 * t45) + m(6) * t147 + t177 + t113 * mrSges(7,3) + ((m(6) * qJ(4) + m(7) * t96 - t189 + t42 + t84) * t111 + (-t156 * mrSges(6,3) + t112 * t136 + t190) * t108) * qJD(3); -t40 * t43 - t121 * t19 + 0.2e1 * t91 * t42 + 0.2e1 * t96 * t18 - t41 * t44 - t120 * t20 + (t16 * t46 + t17 * t45 + t91 * t96) * t183 + t77 * t193 - t110 * t79 + t107 * t78 + (-t161 - t165) * qJD(5) + (t182 + 0.2e1 * t84 + (m(5) + m(6)) * t193) * qJD(4) + t114 * t146; -t120 * t13 + t121 * t14 + t40 * t48 - t41 * t49 + t144 * t153 - m(7) * t115 + m(6) * (-t125 * qJD(5) + t126) - t186; -m(7) * t113 + (t136 + m(5)) * t155; -m(7) * t114 + t133 * t146; -0.2e1 * m(7) * t133; -Ifges(6,5) * t141 + t12 * mrSges(6,1) - t11 * mrSges(6,2) + (-t49 * t149 + t109 * t13 + m(7) * (t106 * t2 + t109 * t3 + t9 * t148 - t8 * t149) + t48 * t148 + t106 * t14) * pkin(5) + t118 + t188; (t106 * t24 + t109 * t25 + (-t106 * t61 - t109 * t62) * qJD(6)) * t180 + t194; ((-mrSges(6,2) * t112 - Ifges(6,6)) * t110 + (-mrSges(6,1) * t112 - Ifges(6,5)) * t107) * qJD(5) + (m(7) * (t106 * t16 + t109 * t17 + (-t106 * t45 + t109 * t46) * qJD(6)) + t185 * mrSges(7,3)) * pkin(5) + t122; -t84 * qJD(5) - t180 * t185 + t138; 0.2e1 * t72; t145 + t188; -t6; t122; t138; t72; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
