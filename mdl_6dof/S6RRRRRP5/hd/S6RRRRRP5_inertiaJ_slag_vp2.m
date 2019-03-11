% Calculate joint inertia matrix for
% S6RRRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRP5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP5_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP5_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP5_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP5_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP5_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:18:12
% EndTime: 2019-03-10 01:18:16
% DurationCPUTime: 1.57s
% Computational Cost: add. (2819->342), mult. (5652->473), div. (0->0), fcn. (5978->8), ass. (0->141)
t176 = 2 * pkin(7);
t175 = mrSges(6,2) + mrSges(7,2);
t131 = sin(qJ(2));
t129 = sin(qJ(4));
t130 = sin(qJ(3));
t133 = cos(qJ(4));
t134 = cos(qJ(3));
t97 = t129 * t134 + t130 * t133;
t85 = t97 * t131;
t96 = -t129 * t130 + t133 * t134;
t86 = t96 * t131;
t174 = -Ifges(5,5) * t86 + Ifges(5,6) * t85;
t116 = pkin(3) * t133 + pkin(4);
t128 = sin(qJ(5));
t132 = cos(qJ(5));
t169 = pkin(3) * t129;
t89 = t132 * t116 - t128 * t169;
t87 = pkin(5) + t89;
t81 = t87 * mrSges(7,1);
t84 = t89 * mrSges(6,1);
t173 = t81 + t84;
t172 = m(7) * pkin(5);
t171 = -pkin(9) - pkin(8);
t170 = m(7) * t128;
t168 = pkin(4) * t128;
t167 = pkin(4) * t132;
t135 = cos(qJ(2));
t166 = pkin(7) * t135;
t123 = t131 * pkin(7);
t61 = -t128 * t97 + t132 * t96;
t90 = t116 * t128 + t132 * t169;
t165 = t61 * t90;
t163 = Ifges(6,3) + Ifges(7,3);
t157 = t131 * t134;
t104 = -pkin(2) * t135 - pkin(8) * t131 - pkin(1);
t95 = t134 * t104;
t63 = -pkin(9) * t157 + t95 + (-pkin(7) * t130 - pkin(3)) * t135;
t158 = t130 * t131;
t76 = t130 * t104 + t134 * t166;
t69 = -pkin(9) * t158 + t76;
t32 = -t129 * t69 + t133 * t63;
t16 = -pkin(4) * t135 - pkin(10) * t86 + t32;
t33 = t129 * t63 + t133 * t69;
t23 = -pkin(10) * t85 + t33;
t6 = t128 * t16 + t132 * t23;
t46 = -t128 * t86 - t132 * t85;
t35 = mrSges(7,2) * t135 + mrSges(7,3) * t46;
t36 = mrSges(6,2) * t135 + mrSges(6,3) * t46;
t162 = t35 + t36;
t108 = t171 * t130;
t109 = t171 * t134;
t71 = t133 * t108 + t109 * t129;
t52 = -pkin(10) * t97 + t71;
t72 = t129 * t108 - t133 * t109;
t53 = pkin(10) * t96 + t72;
t22 = t128 * t52 + t132 * t53;
t161 = Ifges(4,4) * t130;
t160 = Ifges(4,4) * t134;
t159 = t128 * t61;
t103 = pkin(3) * t158 + t123;
t156 = t130 ^ 2 + t134 ^ 2;
t155 = -Ifges(5,3) - t163;
t47 = -t128 * t85 + t132 * t86;
t5 = -t128 * t23 + t132 * t16;
t2 = -pkin(5) * t135 - qJ(6) * t47 + t5;
t37 = -mrSges(7,1) * t135 - mrSges(7,3) * t47;
t154 = m(7) * t2 + t37;
t117 = -pkin(3) * t134 - pkin(2);
t153 = t175 * t90;
t17 = -t46 * mrSges(7,1) + t47 * mrSges(7,2);
t62 = t128 * t96 + t132 * t97;
t25 = -t61 * mrSges(7,1) + t62 * mrSges(7,2);
t21 = -t128 * t53 + t132 * t52;
t151 = Ifges(4,3) - t155;
t150 = (-Ifges(6,5) - Ifges(7,5)) * t47 + (-Ifges(6,6) - Ifges(7,6)) * t46;
t115 = pkin(5) + t167;
t111 = t115 * mrSges(7,1);
t118 = mrSges(6,1) * t167;
t149 = t111 + t118 + t163;
t148 = t175 * t168;
t64 = pkin(4) * t85 + t103;
t8 = -qJ(6) * t62 + t21;
t147 = m(7) * t8 - t62 * mrSges(7,3);
t146 = mrSges(4,1) * t130 + mrSges(4,2) * t134;
t78 = -pkin(4) * t96 + t117;
t145 = (mrSges(5,1) * t133 - mrSges(5,2) * t129) * pkin(3);
t3 = qJ(6) * t46 + t6;
t144 = t5 * mrSges(6,1) + t2 * mrSges(7,1) - t6 * mrSges(6,2) - t3 * mrSges(7,2) - t150;
t57 = Ifges(7,6) * t61;
t58 = Ifges(6,6) * t61;
t59 = Ifges(7,5) * t62;
t60 = Ifges(6,5) * t62;
t9 = qJ(6) * t61 + t22;
t143 = t21 * mrSges(6,1) + t8 * mrSges(7,1) - t22 * mrSges(6,2) - t9 * mrSges(7,2) + t57 + t58 + t59 + t60;
t142 = t32 * mrSges(5,1) - t33 * mrSges(5,2) + t144 - t174;
t92 = Ifges(5,6) * t96;
t93 = Ifges(5,5) * t97;
t141 = t71 * mrSges(5,1) - t72 * mrSges(5,2) + t143 + t92 + t93;
t139 = pkin(4) ^ 2;
t138 = pkin(7) ^ 2;
t136 = pkin(5) * mrSges(7,1);
t127 = t135 ^ 2;
t125 = t131 ^ 2;
t122 = t125 * t138;
t121 = t128 ^ 2 * t139;
t120 = Ifges(4,5) * t130;
t119 = Ifges(4,6) * t134;
t112 = Ifges(4,5) * t157;
t107 = Ifges(4,1) * t130 + t160;
t106 = Ifges(4,2) * t134 + t161;
t105 = -mrSges(4,1) * t134 + mrSges(4,2) * t130;
t102 = -mrSges(4,1) * t135 - mrSges(4,3) * t157;
t101 = mrSges(4,2) * t135 - mrSges(4,3) * t158;
t91 = t146 * t131;
t88 = t90 ^ 2;
t83 = -Ifges(4,5) * t135 + (Ifges(4,1) * t134 - t161) * t131;
t82 = -Ifges(4,6) * t135 + (-Ifges(4,2) * t130 + t160) * t131;
t77 = t90 * t168;
t75 = -t130 * t166 + t95;
t74 = -mrSges(5,1) * t135 - mrSges(5,3) * t86;
t73 = mrSges(5,2) * t135 - mrSges(5,3) * t85;
t68 = Ifges(5,1) * t97 + Ifges(5,4) * t96;
t67 = Ifges(5,4) * t97 + Ifges(5,2) * t96;
t66 = -mrSges(5,1) * t96 + mrSges(5,2) * t97;
t51 = mrSges(5,1) * t85 + mrSges(5,2) * t86;
t45 = Ifges(5,1) * t86 - Ifges(5,4) * t85 - Ifges(5,5) * t135;
t44 = Ifges(5,4) * t86 - Ifges(5,2) * t85 - Ifges(5,6) * t135;
t38 = -mrSges(6,1) * t135 - mrSges(6,3) * t47;
t34 = -pkin(5) * t61 + t78;
t30 = Ifges(6,1) * t62 + Ifges(6,4) * t61;
t29 = Ifges(7,1) * t62 + Ifges(7,4) * t61;
t28 = Ifges(6,4) * t62 + Ifges(6,2) * t61;
t27 = Ifges(7,4) * t62 + Ifges(7,2) * t61;
t26 = -mrSges(6,1) * t61 + mrSges(6,2) * t62;
t24 = -pkin(5) * t46 + t64;
t18 = -mrSges(6,1) * t46 + mrSges(6,2) * t47;
t13 = Ifges(6,1) * t47 + Ifges(6,4) * t46 - Ifges(6,5) * t135;
t12 = Ifges(7,1) * t47 + Ifges(7,4) * t46 - Ifges(7,5) * t135;
t11 = Ifges(6,4) * t47 + Ifges(6,2) * t46 - Ifges(6,6) * t135;
t10 = Ifges(7,4) * t47 + Ifges(7,2) * t46 - Ifges(7,6) * t135;
t1 = [0.2e1 * t76 * t101 + 0.2e1 * t75 * t102 + 0.2e1 * t103 * t51 + 0.2e1 * t24 * t17 + 0.2e1 * t64 * t18 + 0.2e1 * t2 * t37 + 0.2e1 * t3 * t35 + 0.2e1 * t32 * t74 + 0.2e1 * t33 * t73 + 0.2e1 * t6 * t36 + 0.2e1 * t5 * t38 - t85 * t44 + t86 * t45 + Ifges(2,3) + (t12 + t13) * t47 + (t10 + t11) * t46 + (t125 + t127) * mrSges(3,3) * t176 + (-0.2e1 * mrSges(3,2) * pkin(1) + Ifges(3,1) * t131 - t130 * t82 + t134 * t83 + t91 * t176) * t131 + m(3) * (pkin(1) ^ 2 + t127 * t138 + t122) + m(4) * (t75 ^ 2 + t76 ^ 2 + t122) + m(5) * (t103 ^ 2 + t32 ^ 2 + t33 ^ 2) + m(6) * (t5 ^ 2 + t6 ^ 2 + t64 ^ 2) + m(7) * (t2 ^ 2 + t24 ^ 2 + t3 ^ 2) + (0.2e1 * pkin(1) * mrSges(3,1) - t112 + (Ifges(3,2) + t151) * t135 + (Ifges(4,6) * t130 + (2 * Ifges(3,4))) * t131 + t150 + t174) * t135; m(5) * (t103 * t117 + t32 * t71 + t33 * t72) + m(6) * (t21 * t5 + t22 * t6 + t64 * t78) + m(7) * (t2 * t8 + t24 * t34 + t3 * t9) + (t12 / 0.2e1 + t13 / 0.2e1 - t5 * mrSges(6,3) - t2 * mrSges(7,3)) * t62 + (t10 / 0.2e1 + t11 / 0.2e1 + t6 * mrSges(6,3) + t3 * mrSges(7,3)) * t61 + t117 * t51 + t103 * t66 - pkin(2) * t91 + t96 * t44 / 0.2e1 + t97 * t45 / 0.2e1 - t85 * t67 / 0.2e1 + t86 * t68 / 0.2e1 + t64 * t26 + t72 * t73 + t71 * t74 + t78 * t18 + (Ifges(3,5) + t134 * t107 / 0.2e1 - t130 * t106 / 0.2e1 + (t105 - mrSges(3,1)) * pkin(7)) * t131 + m(4) * (-pkin(2) * t123 + (-t130 * t75 + t134 * t76) * pkin(8)) + (t29 / 0.2e1 + t30 / 0.2e1) * t47 + (t27 / 0.2e1 + t28 / 0.2e1) * t46 + t34 * t17 + t9 * t35 + t22 * t36 + t8 * t37 + t21 * t38 + t24 * t25 + (-t120 / 0.2e1 - t119 / 0.2e1 - t93 / 0.2e1 - t92 / 0.2e1 - t59 / 0.2e1 - t57 / 0.2e1 - t60 / 0.2e1 - t58 / 0.2e1 + Ifges(3,6) - pkin(7) * mrSges(3,2)) * t135 + (-t32 * t97 + t33 * t96) * mrSges(5,3) + (t82 / 0.2e1 + pkin(8) * t101 + t76 * mrSges(4,3)) * t134 + (t83 / 0.2e1 - pkin(8) * t102 - t75 * mrSges(4,3)) * t130; -0.2e1 * pkin(2) * t105 + t134 * t106 + t130 * t107 + 0.2e1 * t117 * t66 + 0.2e1 * t34 * t25 + 0.2e1 * t78 * t26 + t96 * t67 + t97 * t68 + Ifges(3,3) + (-0.2e1 * mrSges(6,3) * t21 - 0.2e1 * mrSges(7,3) * t8 + t29 + t30) * t62 + (0.2e1 * mrSges(6,3) * t22 + 0.2e1 * mrSges(7,3) * t9 + t27 + t28) * t61 + m(4) * (pkin(8) ^ 2 * t156 + pkin(2) ^ 2) + m(5) * (t117 ^ 2 + t71 ^ 2 + t72 ^ 2) + m(6) * (t21 ^ 2 + t22 ^ 2 + t78 ^ 2) + m(7) * (t34 ^ 2 + t8 ^ 2 + t9 ^ 2) + 0.2e1 * (-t71 * t97 + t72 * t96) * mrSges(5,3) + 0.2e1 * t156 * pkin(8) * mrSges(4,3); t162 * t90 + t142 + m(7) * (t2 * t87 + t3 * t90) + m(6) * (t5 * t89 + t6 * t90) - t151 * t135 + t112 + t89 * t38 + t87 * t37 + t75 * mrSges(4,1) - t76 * mrSges(4,2) - Ifges(4,6) * t158 + (m(5) * (t129 * t33 + t133 * t32) + t133 * t74 + t129 * t73) * pkin(3); -t146 * pkin(8) + (-t62 * t87 + t165) * mrSges(7,3) + (-t62 * t89 + t165) * mrSges(6,3) + t141 + t119 + t120 + m(7) * (t8 * t87 + t9 * t90) + m(6) * (t21 * t89 + t22 * t90) + (m(5) * (t129 * t72 + t133 * t71) + (t129 * t96 - t133 * t97) * mrSges(5,3)) * pkin(3); 0.2e1 * t81 + 0.2e1 * t84 + m(7) * (t87 ^ 2 + t88) + m(6) * (t89 ^ 2 + t88) + m(5) * (t129 ^ 2 + t133 ^ 2) * pkin(3) ^ 2 + t151 - 0.2e1 * t153 + 0.2e1 * t145; t155 * t135 + t142 + t154 * t115 + (t132 * t38 + t162 * t128 + t3 * t170 + m(6) * (t128 * t6 + t132 * t5)) * pkin(4); t141 + t147 * t115 + (mrSges(7,3) * t159 + (-t132 * t62 + t159) * mrSges(6,3) + t9 * t170 + m(6) * (t128 * t22 + t132 * t21)) * pkin(4); Ifges(5,3) + t145 + m(7) * (t115 * t87 + t77) + m(6) * (t89 * t167 + t77) + t149 + t175 * (-t90 - t168) + t173; 0.2e1 * t111 + 0.2e1 * t118 - 0.2e1 * t148 + m(6) * (t132 ^ 2 * t139 + t121) + m(7) * (t115 ^ 2 + t121) - t155; t154 * pkin(5) - t163 * t135 + t144; pkin(5) * t147 + t143; t87 * t172 + t136 - t153 + t163 + t173; t115 * t172 + t136 - t148 + t149; m(7) * pkin(5) ^ 2 + 0.2e1 * t136 + t163; m(7) * t24 + t17; m(7) * t34 + t25; 0; 0; 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
