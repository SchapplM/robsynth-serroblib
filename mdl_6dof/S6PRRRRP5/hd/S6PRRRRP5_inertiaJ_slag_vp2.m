% Calculate joint inertia matrix for
% S6PRRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-03-09 00:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRRP5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP5_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP5_inertiaJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP5_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP5_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP5_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:20:17
% EndTime: 2019-03-09 00:20:22
% DurationCPUTime: 1.78s
% Computational Cost: add. (1933->380), mult. (4729->518), div. (0->0), fcn. (5065->12), ass. (0->136)
t132 = sin(qJ(5));
t136 = cos(qJ(5));
t89 = -mrSges(6,1) * t136 + mrSges(6,2) * t132;
t181 = -m(6) * pkin(4) + t89;
t180 = -mrSges(5,1) + t181;
t179 = 2 * pkin(10);
t178 = -2 * mrSges(7,3);
t128 = sin(pkin(7));
t138 = cos(qJ(3));
t157 = t128 * t138;
t130 = cos(pkin(7));
t133 = sin(qJ(4));
t137 = cos(qJ(4));
t134 = sin(qJ(3));
t158 = t128 * t134;
t70 = t130 * t133 + t137 * t158;
t44 = -t132 * t70 - t136 * t157;
t45 = -t132 * t157 + t136 * t70;
t69 = -t137 * t130 + t133 * t158;
t10 = Ifges(6,5) * t45 + Ifges(6,6) * t44 + Ifges(6,3) * t69;
t9 = Ifges(7,5) * t45 + Ifges(7,6) * t44 + Ifges(7,3) * t69;
t177 = t10 + t9;
t113 = -pkin(5) * t136 - pkin(4);
t88 = -t136 * mrSges(7,1) + t132 * mrSges(7,2);
t176 = m(7) * t113 + t88;
t129 = sin(pkin(6));
t131 = cos(pkin(6));
t135 = sin(qJ(2));
t139 = cos(qJ(2));
t156 = t130 * t139;
t38 = t131 * t158 + (t134 * t156 + t135 * t138) * t129;
t68 = -t128 * t129 * t139 + t130 * t131;
t20 = t133 * t38 - t68 * t137;
t175 = t20 ^ 2;
t36 = -t131 * t157 + (t134 * t135 - t138 * t156) * t129;
t174 = t36 ^ 2;
t172 = m(7) * pkin(5);
t171 = pkin(2) * t138;
t170 = pkin(10) * t137;
t169 = -Ifges(6,3) - Ifges(7,3);
t168 = -qJ(6) - pkin(11);
t103 = pkin(9) * t158;
t55 = t103 + (-pkin(3) - t171) * t130;
t23 = pkin(4) * t69 - pkin(11) * t70 + t55;
t73 = t130 * t134 * pkin(2) + pkin(9) * t157;
t56 = pkin(10) * t130 + t73;
t57 = (-pkin(3) * t138 - pkin(10) * t134 - pkin(2)) * t128;
t31 = t133 * t57 + t137 * t56;
t25 = -pkin(11) * t157 + t31;
t5 = t132 * t23 + t136 * t25;
t16 = -mrSges(6,1) * t44 + mrSges(6,2) * t45;
t48 = -mrSges(5,1) * t157 - mrSges(5,3) * t70;
t167 = -t48 + t16;
t166 = -Ifges(5,5) * t70 + Ifges(5,6) * t69;
t164 = Ifges(6,4) * t132;
t163 = Ifges(6,4) * t136;
t162 = Ifges(7,4) * t132;
t161 = Ifges(7,4) * t136;
t160 = t133 * t20;
t22 = t133 * t68 + t137 * t38;
t159 = t137 * t22;
t86 = -pkin(4) * t137 - pkin(11) * t133 - pkin(3);
t53 = t132 * t86 + t136 * t170;
t155 = t132 * t133;
t154 = t133 * t136;
t74 = mrSges(7,1) * t155 + mrSges(7,2) * t154;
t92 = Ifges(7,5) * t132 + Ifges(7,6) * t136;
t93 = Ifges(6,5) * t132 + Ifges(6,6) * t136;
t153 = Ifges(5,5) * t133 + Ifges(5,6) * t137;
t152 = t132 ^ 2 + t136 ^ 2;
t11 = Ifges(7,4) * t45 + Ifges(7,2) * t44 + Ifges(7,6) * t69;
t12 = Ifges(6,4) * t45 + Ifges(6,2) * t44 + Ifges(6,6) * t69;
t151 = -t11 / 0.2e1 - t12 / 0.2e1;
t13 = Ifges(7,1) * t45 + Ifges(7,4) * t44 + Ifges(7,5) * t69;
t14 = Ifges(6,1) * t45 + Ifges(6,4) * t44 + Ifges(6,5) * t69;
t150 = t13 / 0.2e1 + t14 / 0.2e1;
t60 = -Ifges(7,6) * t137 + (-Ifges(7,2) * t132 + t161) * t133;
t61 = -Ifges(6,6) * t137 + (-Ifges(6,2) * t132 + t163) * t133;
t149 = t60 / 0.2e1 + t61 / 0.2e1;
t62 = -Ifges(7,5) * t137 + (Ifges(7,1) * t136 - t162) * t133;
t63 = -Ifges(6,5) * t137 + (Ifges(6,1) * t136 - t164) * t133;
t148 = t62 / 0.2e1 + t63 / 0.2e1;
t147 = t92 / 0.2e1 + t93 / 0.2e1;
t94 = Ifges(7,2) * t136 + t162;
t95 = Ifges(6,2) * t136 + t164;
t146 = t95 / 0.2e1 + t94 / 0.2e1;
t97 = Ifges(7,1) * t132 + t161;
t98 = Ifges(6,1) * t132 + t163;
t145 = t97 / 0.2e1 + t98 / 0.2e1;
t144 = Ifges(4,5) * t158 + Ifges(4,6) * t157 + Ifges(4,3) * t130;
t15 = -t44 * mrSges(7,1) + t45 * mrSges(7,2);
t4 = -t132 * t25 + t136 * t23;
t30 = -t133 * t56 + t137 * t57;
t24 = pkin(4) * t157 - t30;
t142 = mrSges(6,1) * t132 + mrSges(6,2) * t136;
t141 = pkin(10) ^ 2;
t127 = t137 ^ 2;
t125 = t133 ^ 2;
t122 = t125 * t141;
t110 = Ifges(6,5) * t154;
t109 = Ifges(7,5) * t154;
t99 = Ifges(5,1) * t133 + Ifges(5,4) * t137;
t96 = Ifges(5,4) * t133 + Ifges(5,2) * t137;
t91 = t168 * t136;
t90 = -mrSges(5,1) * t137 + mrSges(5,2) * t133;
t87 = t168 * t132;
t85 = (pkin(5) * t132 + pkin(10)) * t133;
t83 = -mrSges(6,1) * t137 - mrSges(6,3) * t154;
t82 = -mrSges(7,1) * t137 - mrSges(7,3) * t154;
t81 = mrSges(6,2) * t137 - mrSges(6,3) * t155;
t80 = mrSges(7,2) * t137 - mrSges(7,3) * t155;
t79 = -mrSges(4,2) * t130 + mrSges(4,3) * t157;
t78 = mrSges(4,1) * t130 - mrSges(4,3) * t158;
t77 = t136 * t86;
t75 = t142 * t133;
t72 = t130 * t171 - t103;
t71 = (-mrSges(4,1) * t138 + mrSges(4,2) * t134) * t128;
t59 = -Ifges(6,6) * t155 - Ifges(6,3) * t137 + t110;
t58 = -Ifges(7,6) * t155 - Ifges(7,3) * t137 + t109;
t52 = -t132 * t170 + t77;
t47 = mrSges(5,2) * t157 - mrSges(5,3) * t69;
t46 = -qJ(6) * t155 + t53;
t35 = -qJ(6) * t154 + t77 + (-pkin(10) * t132 - pkin(5)) * t137;
t34 = mrSges(5,1) * t69 + mrSges(5,2) * t70;
t33 = Ifges(5,1) * t70 - Ifges(5,4) * t69 - Ifges(5,5) * t157;
t32 = Ifges(5,4) * t70 - Ifges(5,2) * t69 - Ifges(5,6) * t157;
t29 = mrSges(6,1) * t69 - mrSges(6,3) * t45;
t28 = mrSges(7,1) * t69 - mrSges(7,3) * t45;
t27 = -mrSges(6,2) * t69 + mrSges(6,3) * t44;
t26 = -mrSges(7,2) * t69 + mrSges(7,3) * t44;
t8 = -pkin(5) * t44 + t24;
t7 = t132 * t36 + t136 * t22;
t6 = -t132 * t22 + t136 * t36;
t3 = qJ(6) * t44 + t5;
t2 = pkin(5) * t69 - qJ(6) * t45 + t4;
t1 = [m(2) + m(5) * (t22 ^ 2 + t174 + t175) + m(4) * (t38 ^ 2 + t68 ^ 2 + t174) + m(3) * (t131 ^ 2 + (t135 ^ 2 + t139 ^ 2) * t129 ^ 2) + 0.2e1 * (m(7) / 0.2e1 + m(6) / 0.2e1) * (t6 ^ 2 + t7 ^ 2 + t175); t22 * t47 + t38 * t79 + t68 * t71 + (t26 + t27) * t7 + (t28 + t29) * t6 + (t34 - t78) * t36 + (mrSges(3,1) * t139 - mrSges(3,2) * t135) * t129 + (t15 + t167) * t20 + m(7) * (t2 * t6 + t20 * t8 + t3 * t7) + m(6) * (t20 * t24 + t4 * t6 + t5 * t7) + m(5) * (-t20 * t30 + t22 * t31 + t36 * t55) + m(4) * (-pkin(2) * t128 * t68 - t36 * t72 + t38 * t73); t130 * t144 + t70 * t33 + 0.2e1 * t72 * t78 + 0.2e1 * t73 * t79 + 0.2e1 * t55 * t34 + 0.2e1 * t31 * t47 + 0.2e1 * t30 * t48 + 0.2e1 * t3 * t26 + 0.2e1 * t5 * t27 + 0.2e1 * t2 * t28 + 0.2e1 * t4 * t29 + 0.2e1 * t24 * t16 + 0.2e1 * t8 * t15 + Ifges(3,3) + (t13 + t14) * t45 + (t11 + t12) * t44 + (-t32 + t177) * t69 + (-0.2e1 * pkin(2) * t71 + (Ifges(4,1) * t158 + Ifges(4,5) * t130) * t134 + (Ifges(4,6) * t130 + (0.2e1 * Ifges(4,4) * t134 + (Ifges(4,2) + Ifges(5,3)) * t138) * t128 + t166) * t138) * t128 + m(4) * (pkin(2) ^ 2 * t128 ^ 2 + t72 ^ 2 + t73 ^ 2) + m(5) * (t30 ^ 2 + t31 ^ 2 + t55 ^ 2) + m(6) * (t24 ^ 2 + t4 ^ 2 + t5 ^ 2) + m(7) * (t2 ^ 2 + t3 ^ 2 + t8 ^ 2); mrSges(5,3) * t159 - t38 * mrSges(4,2) + (t80 + t81) * t7 + (t82 + t83) * t6 + (-mrSges(4,1) + t90) * t36 + (t133 * mrSges(5,3) + t74 + t75) * t20 + m(7) * (t20 * t85 + t35 * t6 + t46 * t7) + m(6) * (pkin(10) * t160 + t52 * t6 + t53 * t7) + m(5) * (-pkin(3) * t36 + (t159 + t160) * pkin(10)); t144 + (t33 / 0.2e1 - t30 * mrSges(5,3) + t150 * t136 + t151 * t132 + t167 * pkin(10)) * t133 - t153 * t157 / 0.2e1 + t148 * t45 + t149 * t44 + m(5) * (-pkin(3) * t55 + (-t30 * t133 + t31 * t137) * pkin(10)) + t70 * t99 / 0.2e1 + t3 * t80 + t5 * t81 + t2 * t82 + t4 * t83 + t85 * t15 + t55 * t90 + t72 * mrSges(4,1) - t73 * mrSges(4,2) + t8 * t74 + t24 * t75 + t46 * t26 + t52 * t29 + t53 * t27 + t35 * t28 - pkin(3) * t34 + (t32 / 0.2e1 - t9 / 0.2e1 - t10 / 0.2e1 + pkin(10) * t47 + t31 * mrSges(5,3)) * t137 + m(7) * (t2 * t35 + t3 * t46 + t8 * t85) + m(6) * (pkin(10) * t133 * t24 + t4 * t52 + t5 * t53) + (-t96 / 0.2e1 + t58 / 0.2e1 + t59 / 0.2e1) * t69; -0.2e1 * pkin(3) * t90 + 0.2e1 * t35 * t82 + 0.2e1 * t46 * t80 + 0.2e1 * t52 * t83 + 0.2e1 * t53 * t81 + 0.2e1 * t85 * t74 + Ifges(4,3) + (t125 + t127) * mrSges(5,3) * t179 + (-t59 - t58 + t96) * t137 + m(6) * (t52 ^ 2 + t53 ^ 2 + t122) + m(7) * (t35 ^ 2 + t46 ^ 2 + t85 ^ 2) + m(5) * (pkin(3) ^ 2 + t127 * t141 + t122) + (t75 * t179 + t99 + (t62 + t63) * t136 + (-t60 - t61) * t132) * t133; -t22 * mrSges(5,2) + m(7) * (t6 * t87 - t7 * t91) + (t176 + t180) * t20 + (m(6) * pkin(11) + mrSges(6,3) + mrSges(7,3)) * (-t132 * t6 + t136 * t7); -Ifges(5,3) * t157 + t30 * mrSges(5,1) - t31 * mrSges(5,2) - pkin(4) * t16 + t113 * t15 - t91 * t26 + t87 * t28 + t8 * t88 + t147 * t69 + t145 * t45 + t146 * t44 + m(7) * (t113 * t8 + t2 * t87 - t3 * t91) + (t3 * mrSges(7,3) + t5 * mrSges(6,3) + (m(6) * t5 + t27) * pkin(11) - t151) * t136 + (-t4 * mrSges(6,3) - t2 * mrSges(7,3) + (-m(6) * t4 - t29) * pkin(11) + t150) * t132 - t166 + t181 * t24; t113 * t74 + t87 * t82 + t85 * t88 - t91 * t80 - pkin(4) * t75 + m(7) * (t113 * t85 + t35 * t87 - t46 * t91) - t147 * t137 + (-t137 * mrSges(5,2) + t133 * t180) * pkin(10) + (t46 * mrSges(7,3) + t53 * mrSges(6,3) + t145 * t133 + (m(6) * t53 + t81) * pkin(11) + t149) * t136 + (-t52 * mrSges(6,3) - t35 * mrSges(7,3) - t146 * t133 + (-m(6) * t52 - t83) * pkin(11) + t148) * t132 + t153; -0.2e1 * pkin(4) * t89 + 0.2e1 * t113 * t88 + Ifges(5,3) + 0.2e1 * t152 * pkin(11) * mrSges(6,3) + m(7) * (t113 ^ 2 + t87 ^ 2 + t91 ^ 2) + m(6) * (pkin(11) ^ 2 * t152 + pkin(4) ^ 2) + (t178 * t91 + t94 + t95) * t136 + (t178 * t87 + t97 + t98) * t132; (-mrSges(6,2) - mrSges(7,2)) * t7 + (mrSges(6,1) + mrSges(7,1) + t172) * t6; mrSges(6,1) * t4 + mrSges(7,1) * t2 - mrSges(6,2) * t5 - mrSges(7,2) * t3 + (m(7) * t2 + t28) * pkin(5) + t177; mrSges(6,1) * t52 + mrSges(7,1) * t35 - mrSges(6,2) * t53 - mrSges(7,2) * t46 + t109 + t110 + t169 * t137 + (-Ifges(6,6) - Ifges(7,6)) * t155 + (m(7) * t35 + t82) * pkin(5); mrSges(7,1) * t87 + mrSges(7,2) * t91 - t142 * pkin(11) + (m(7) * t87 - t132 * mrSges(7,3)) * pkin(5) + t93 + t92; (0.2e1 * mrSges(7,1) + t172) * pkin(5) - t169; m(7) * t20; m(7) * t8 + t15; m(7) * t85 + t74; t176; 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
