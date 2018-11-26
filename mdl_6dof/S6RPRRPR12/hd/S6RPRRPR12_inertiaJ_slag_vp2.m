% Calculate joint inertia matrix for
% S6RPRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2]';
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
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 16:23
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPR12_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR12_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRPR12_inertiaJ_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR12_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR12_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR12_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:23:11
% EndTime: 2018-11-23 16:23:12
% DurationCPUTime: 1.68s
% Computational Cost: add. (3727->359), mult. (9804->504), div. (0->0), fcn. (10981->12), ass. (0->138)
t183 = Ifges(6,1) + Ifges(5,3);
t126 = sin(qJ(4));
t129 = cos(qJ(4));
t182 = t126 ^ 2 + t129 ^ 2;
t125 = sin(qJ(6));
t172 = -t125 / 0.2e1;
t167 = mrSges(6,1) + mrSges(5,3);
t180 = -m(6) * pkin(4) + mrSges(6,2);
t120 = sin(pkin(7));
t123 = cos(pkin(7));
t127 = sin(qJ(3));
t130 = cos(qJ(3));
t124 = cos(pkin(6));
t121 = sin(pkin(6));
t122 = cos(pkin(12));
t154 = t121 * t122;
t146 = t123 * t154;
t119 = sin(pkin(12));
t168 = pkin(1) * t124;
t77 = qJ(2) * t154 + t119 * t168;
t53 = (t120 * t124 + t146) * pkin(9) + t77;
t104 = t122 * t168;
t157 = t119 * t121;
t61 = pkin(2) * t124 + t104 + (-pkin(9) * t123 - qJ(2)) * t157;
t68 = (-pkin(9) * t119 * t120 - pkin(2) * t122 - pkin(1)) * t121;
t29 = -t127 * t53 + (t120 * t68 + t123 * t61) * t130;
t156 = t120 * t127;
t80 = t123 * t126 + t129 * t156;
t75 = t80 ^ 2;
t179 = 2 * mrSges(3,1);
t178 = 0.2e1 * t124;
t128 = cos(qJ(6));
t153 = t123 * t127;
t60 = t124 * t156 + (t119 * t130 + t122 * t153) * t121;
t74 = -t120 * t154 + t123 * t124;
t46 = t126 * t60 - t129 * t74;
t155 = t120 * t130;
t59 = -t124 * t155 + t127 * t157 - t130 * t146;
t33 = -t125 * t59 + t128 * t46;
t177 = t33 / 0.2e1;
t161 = Ifges(7,4) * t128;
t72 = Ifges(7,5) * t126 + (-Ifges(7,1) * t125 - t161) * t129;
t176 = t72 / 0.2e1;
t109 = Ifges(7,5) * t128;
t175 = Ifges(7,6) * t172 + t109 / 0.2e1;
t174 = pkin(4) + pkin(11);
t173 = pkin(5) + pkin(10);
t171 = -t128 / 0.2e1;
t170 = Ifges(6,4) * t59;
t169 = Ifges(6,5) * t59;
t90 = mrSges(7,1) * t125 + mrSges(7,2) * t128;
t166 = mrSges(6,3) + t90;
t40 = -t120 * t61 + t123 * t68;
t22 = pkin(3) * t59 - pkin(10) * t60 + t40;
t30 = t130 * t53 + t61 * t153 + t68 * t156;
t26 = pkin(10) * t74 + t30;
t9 = t126 * t22 + t129 * t26;
t35 = mrSges(6,1) * t46 - mrSges(6,3) * t59;
t37 = -mrSges(5,2) * t59 - mrSges(5,3) * t46;
t165 = -t35 + t37;
t47 = t126 * t74 + t129 * t60;
t36 = t47 * mrSges(6,1) + t59 * mrSges(6,2);
t38 = mrSges(5,1) * t59 - mrSges(5,3) * t47;
t164 = t36 - t38;
t163 = mrSges(7,3) * t129;
t162 = Ifges(7,4) * t125;
t160 = qJ(5) * t80;
t86 = -mrSges(7,2) * t126 - t128 * t163;
t159 = t125 * t86;
t78 = -t129 * t123 + t126 * t156;
t158 = t126 * t78;
t152 = t128 * t174;
t151 = Ifges(5,5) * t126 + Ifges(5,6) * t129;
t150 = t182 * pkin(10) ^ 2;
t149 = -t125 ^ 2 - t128 ^ 2;
t34 = t125 * t46 + t128 * t59;
t10 = Ifges(7,5) * t34 + Ifges(7,6) * t33 + Ifges(7,3) * t47;
t147 = Ifges(4,5) * t60 - Ifges(4,6) * t59 + Ifges(4,3) * t74;
t145 = m(7) * t149;
t144 = -qJ(5) * t126 - pkin(3);
t8 = -t126 * t26 + t129 * t22;
t143 = t149 * mrSges(7,3);
t6 = -qJ(5) * t59 - t9;
t3 = pkin(5) * t47 - t174 * t59 - t8;
t25 = -pkin(3) * t74 - t29;
t134 = -qJ(5) * t47 + t25;
t5 = t174 * t46 + t134;
t1 = -t125 * t5 + t128 * t3;
t2 = t125 * t3 + t128 * t5;
t141 = t1 * t128 + t125 * t2;
t140 = mrSges(7,1) * t128 - t125 * mrSges(7,2);
t83 = -t174 * t129 + t144;
t99 = t173 * t126;
t62 = -t125 * t83 + t128 * t99;
t63 = t125 * t99 + t128 * t83;
t138 = t125 * t63 + t128 * t62;
t65 = t125 * t155 + t128 * t78;
t66 = t125 * t78 - t128 * t155;
t137 = t125 * t66 + t128 * t65;
t136 = t183 * t59 + (-Ifges(6,4) + Ifges(5,5)) * t47 + (Ifges(6,5) - Ifges(5,6)) * t46;
t135 = (t129 * t80 + t158) * pkin(10);
t70 = Ifges(7,3) * t126 + (-Ifges(7,5) * t125 - Ifges(7,6) * t128) * t129;
t132 = qJ(5) ^ 2;
t113 = t120 ^ 2;
t105 = t113 * t130 ^ 2;
t102 = mrSges(3,2) * t157;
t100 = t173 * t129;
t97 = Ifges(5,1) * t126 + Ifges(5,4) * t129;
t96 = Ifges(7,1) * t128 - t162;
t95 = Ifges(5,4) * t126 + Ifges(5,2) * t129;
t94 = -Ifges(7,2) * t125 + t161;
t92 = -Ifges(6,2) * t126 - Ifges(6,6) * t129;
t91 = -Ifges(6,6) * t126 - Ifges(6,3) * t129;
t89 = -mrSges(5,1) * t129 + mrSges(5,2) * t126;
t88 = mrSges(6,2) * t129 - mrSges(6,3) * t126;
t87 = -pkin(4) * t129 + t144;
t85 = mrSges(7,1) * t126 + t125 * t163;
t82 = t140 * t129;
t76 = -qJ(2) * t157 + t104;
t71 = Ifges(7,6) * t126 + (-Ifges(7,2) * t128 - t162) * t129;
t49 = mrSges(4,1) * t74 - mrSges(4,3) * t60;
t48 = -mrSges(4,2) * t74 - mrSges(4,3) * t59;
t39 = mrSges(4,1) * t59 + mrSges(4,2) * t60;
t28 = -mrSges(6,2) * t46 - mrSges(6,3) * t47;
t27 = mrSges(5,1) * t46 + mrSges(5,2) * t47;
t20 = Ifges(5,1) * t47 - Ifges(5,4) * t46 + Ifges(5,5) * t59;
t19 = Ifges(5,4) * t47 - Ifges(5,2) * t46 + Ifges(5,6) * t59;
t18 = -Ifges(6,2) * t47 + Ifges(6,6) * t46 + t170;
t17 = -Ifges(6,6) * t47 + Ifges(6,3) * t46 + t169;
t16 = mrSges(7,1) * t47 - mrSges(7,3) * t34;
t15 = -mrSges(7,2) * t47 + mrSges(7,3) * t33;
t14 = -mrSges(7,1) * t33 + mrSges(7,2) * t34;
t13 = pkin(4) * t46 + t134;
t12 = Ifges(7,1) * t34 + Ifges(7,4) * t33 + Ifges(7,5) * t47;
t11 = Ifges(7,4) * t34 + Ifges(7,2) * t33 + Ifges(7,6) * t47;
t7 = -pkin(4) * t59 - t8;
t4 = -pkin(5) * t46 - t6;
t21 = [t60 * (Ifges(4,1) * t60 + Ifges(4,5) * t74) + (-0.2e1 * Ifges(4,4) * t60 + Ifges(4,2) * t59 - Ifges(4,6) * t74 + t136) * t59 + m(3) * (pkin(1) ^ 2 * t121 ^ 2 + t76 ^ 2 + t77 ^ 2) + m(4) * (t29 ^ 2 + t30 ^ 2 + t40 ^ 2) + m(5) * (t25 ^ 2 + t8 ^ 2 + t9 ^ 2) + m(6) * (t13 ^ 2 + t6 ^ 2 + t7 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t4 ^ 2) + (t20 + t10 - t18) * t47 + (-t19 + t17) * t46 + t74 * t147 + 0.2e1 * t30 * t48 + 0.2e1 * t29 * t49 + 0.2e1 * t40 * t39 + t33 * t11 + t34 * t12 + 0.2e1 * t6 * t35 + 0.2e1 * t7 * t36 + 0.2e1 * t9 * t37 + 0.2e1 * t8 * t38 + 0.2e1 * t25 * t27 + 0.2e1 * t13 * t28 + 0.2e1 * t4 * t14 + 0.2e1 * t2 * t15 + 0.2e1 * t1 * t16 + Ifges(2,3) + (-0.2e1 * mrSges(3,2) * t77 + Ifges(3,3) * t124 + t76 * t179) * t124 + (-0.2e1 * pkin(1) * t102 + (-0.2e1 * mrSges(3,3) * t76 + Ifges(3,1) * t157 + Ifges(3,5) * t178) * t119 + (0.2e1 * t77 * mrSges(3,3) + Ifges(3,6) * t178 + (0.2e1 * Ifges(3,4) * t119 + Ifges(3,2) * t122 + pkin(1) * t179) * t121) * t122) * t121; t123 * t39 + t66 * t15 + t65 * t16 + t102 + t164 * t78 + (-m(3) * pkin(1) - mrSges(3,1) * t122) * t121 + (t14 + t165) * t80 + (t127 * t48 + (-t27 - t28 + t49) * t130) * t120 + m(7) * (t1 * t65 + t2 * t66 + t4 * t80) + m(6) * (-t13 * t155 - t6 * t80 + t7 * t78) + m(5) * (-t25 * t155 - t78 * t8 + t80 * t9) + m(4) * (t123 * t40 + (t127 * t30 + t130 * t29) * t120); m(3) + m(7) * (t65 ^ 2 + t66 ^ 2 + t75) + m(4) * (t113 * t127 ^ 2 + t123 ^ 2 + t105) + (m(6) + m(5)) * (t78 ^ 2 + t105 + t75); m(7) * (t1 * t62 + t100 * t4 + t2 * t63) + (t70 / 0.2e1 - t92 / 0.2e1 + t97 / 0.2e1) * t47 + (t91 / 0.2e1 - t95 / 0.2e1) * t46 + t147 + t62 * t16 + t63 * t15 + t71 * t177 + t34 * t176 + t4 * t82 + t1 * t85 + t2 * t86 + t87 * t28 + t13 * t88 + t25 * t89 + t100 * t14 + t59 * t151 / 0.2e1 - pkin(3) * t27 + t29 * mrSges(4,1) - t30 * mrSges(4,2) + (t7 * mrSges(6,1) - t8 * mrSges(5,3) - t18 / 0.2e1 + t20 / 0.2e1 + t10 / 0.2e1 - t170 / 0.2e1 + t164 * pkin(10)) * t126 + (t11 * t171 + t12 * t172 + t9 * mrSges(5,3) - t6 * mrSges(6,1) - t17 / 0.2e1 + t19 / 0.2e1 - t169 / 0.2e1 + t165 * pkin(10)) * t129 + m(5) * (-pkin(3) * t25 + (-t126 * t8 + t129 * t9) * pkin(10)) + m(6) * (t13 * t87 + (t126 * t7 - t129 * t6) * pkin(10)); t65 * t85 + t66 * t86 + t167 * t158 + (t167 * t129 + t82) * t80 + (-t127 * mrSges(4,2) + (mrSges(4,1) - t88 - t89) * t130) * t120 + m(7) * (t100 * t80 + t62 * t65 + t63 * t66) + m(6) * (-t87 * t155 + t135) + m(5) * (pkin(3) * t155 + t135); -0.2e1 * pkin(3) * t89 + 0.2e1 * t100 * t82 + 0.2e1 * t62 * t85 + 0.2e1 * t63 * t86 + 0.2e1 * t87 * t88 + Ifges(4,3) + (t70 + t97 - t92) * t126 + m(7) * (t100 ^ 2 + t62 ^ 2 + t63 ^ 2) + m(6) * (t87 ^ 2 + t150) + m(5) * (pkin(3) ^ 2 + t150) + (-t125 * t72 - t128 * t71 - t91 + t95) * t129 + 0.2e1 * t167 * pkin(10) * t182; (-t1 * mrSges(7,3) - t174 * t16 + t12 / 0.2e1) * t128 + (-t2 * mrSges(7,3) - t174 * t15 - t11 / 0.2e1) * t125 + (-t35 + t14) * qJ(5) + m(6) * (-pkin(4) * t7 - qJ(5) * t6) + t4 * t90 + t47 * t175 + t94 * t177 + t34 * t96 / 0.2e1 - pkin(4) * t36 + t8 * mrSges(5,1) - t9 * mrSges(5,2) - t6 * mrSges(6,3) + t7 * mrSges(6,2) + t136 + m(7) * (qJ(5) * t4 - t141 * t174); (-mrSges(5,1) + mrSges(6,2)) * t78 - t137 * mrSges(7,3) + (-mrSges(5,2) + t166) * t80 + m(7) * (-t137 * t174 + t160) + m(6) * (-pkin(4) * t78 + t160); -t174 * t159 - t85 * t152 + qJ(5) * t82 + m(7) * (qJ(5) * t100 - t138 * t174) + t128 * t176 + t71 * t172 + t100 * t90 - t138 * mrSges(7,3) + (-pkin(4) * mrSges(6,1) - Ifges(6,4) + t175) * t126 + (qJ(5) * mrSges(6,1) + t94 * t171 + t96 * t172 - Ifges(6,5)) * t129 + ((m(6) * qJ(5) - mrSges(5,2) + mrSges(6,3)) * t129 + (-mrSges(5,1) + t180) * t126) * pkin(10) + t151; -0.2e1 * pkin(4) * mrSges(6,2) - t125 * t94 + t128 * t96 + m(6) * (pkin(4) ^ 2 + t132) + m(7) * (-t149 * t174 ^ 2 + t132) + 0.2e1 * t166 * qJ(5) - 0.2e1 * t174 * t143 + t183; m(6) * t7 + m(7) * t141 + t125 * t15 + t128 * t16 + t36; m(6) * t78 + m(7) * t137; m(7) * t138 + t159 + t128 * t85 + (m(6) * pkin(10) + mrSges(6,1)) * t126; t145 * t174 + t143 + t180; m(6) - t145; mrSges(7,1) * t1 - mrSges(7,2) * t2 + t10; mrSges(7,1) * t65 - mrSges(7,2) * t66; mrSges(7,1) * t62 - mrSges(7,2) * t63 + t70; -mrSges(7,1) * t152 + t109 + (mrSges(7,2) * t174 - Ifges(7,6)) * t125; t140; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t21(1) t21(2) t21(4) t21(7) t21(11) t21(16); t21(2) t21(3) t21(5) t21(8) t21(12) t21(17); t21(4) t21(5) t21(6) t21(9) t21(13) t21(18); t21(7) t21(8) t21(9) t21(10) t21(14) t21(19); t21(11) t21(12) t21(13) t21(14) t21(15) t21(20); t21(16) t21(17) t21(18) t21(19) t21(20) t21(21);];
Mq  = res;
