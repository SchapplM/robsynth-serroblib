% Calculate joint inertia matrix for
% S6RRRPRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
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
% Datum: 2019-03-09 18:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRP12_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP12_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP12_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP12_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP12_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP12_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:51:54
% EndTime: 2019-03-09 17:51:59
% DurationCPUTime: 1.82s
% Computational Cost: add. (2010->388), mult. (4369->511), div. (0->0), fcn. (4354->8), ass. (0->141)
t189 = pkin(4) + pkin(9);
t131 = sin(qJ(3));
t134 = cos(qJ(3));
t188 = t131 ^ 2 + t134 ^ 2;
t130 = sin(qJ(5));
t133 = cos(qJ(5));
t162 = -t130 ^ 2 - t133 ^ 2;
t128 = sin(pkin(6));
t135 = cos(qJ(2));
t169 = t128 * t135;
t129 = cos(pkin(6));
t132 = sin(qJ(2));
t170 = t128 * t132;
t69 = -t129 * t134 + t131 * t170;
t43 = t130 * t169 + t133 * t69;
t44 = t130 * t69 - t133 * t169;
t70 = t129 * t131 + t134 * t170;
t8 = Ifges(6,5) * t44 + Ifges(6,6) * t43 + Ifges(6,3) * t70;
t9 = Ifges(7,4) * t44 + Ifges(7,2) * t70 - Ifges(7,6) * t43;
t187 = t8 + t9;
t186 = -m(5) * pkin(3) + mrSges(5,2);
t185 = (mrSges(7,2) + mrSges(6,3)) * t162;
t184 = pkin(3) + pkin(10);
t72 = t129 * t132 * pkin(1) + pkin(8) * t169;
t54 = pkin(9) * t129 + t72;
t55 = (-pkin(2) * t135 - pkin(9) * t132 - pkin(1)) * t128;
t26 = -t131 * t54 + t134 * t55;
t21 = pkin(3) * t169 - t26;
t13 = pkin(4) * t70 + pkin(10) * t169 + t21;
t103 = pkin(8) * t170;
t183 = pkin(1) * t135;
t53 = t103 + (-pkin(2) - t183) * t129;
t140 = -qJ(4) * t70 + t53;
t15 = t184 * t69 + t140;
t4 = t130 * t13 + t133 * t15;
t71 = t129 * t183 - t103;
t182 = t71 * mrSges(3,1);
t181 = t72 * mrSges(3,2);
t180 = Ifges(5,1) + Ifges(4,3);
t22 = mrSges(7,2) * t43 + mrSges(7,3) * t70;
t23 = -mrSges(6,2) * t70 + mrSges(6,3) * t43;
t179 = t22 + t23;
t24 = mrSges(6,1) * t70 - mrSges(6,3) * t44;
t25 = -t70 * mrSges(7,1) + t44 * mrSges(7,2);
t178 = t24 - t25;
t27 = t131 * t55 + t134 * t54;
t152 = -qJ(4) * t131 - pkin(2);
t75 = -t184 * t134 + t152;
t98 = t189 * t131;
t37 = t130 * t98 + t133 * t75;
t168 = t130 * t134;
t78 = mrSges(6,1) * t131 + mrSges(6,3) * t168;
t79 = -t131 * mrSges(7,1) - mrSges(7,2) * t168;
t177 = t78 - t79;
t167 = t133 * t134;
t80 = -mrSges(6,2) * t131 - mrSges(6,3) * t167;
t81 = -mrSges(7,2) * t167 + mrSges(7,3) * t131;
t176 = t80 + t81;
t175 = Ifges(6,4) * t130;
t174 = Ifges(6,4) * t133;
t173 = Ifges(7,5) * t130;
t172 = Ifges(7,5) * t133;
t171 = Ifges(6,6) * t133;
t166 = Ifges(7,2) * t131 + Ifges(7,6) * t167;
t165 = t162 * t184 ^ 2;
t164 = Ifges(4,5) * t131 + Ifges(4,6) * t134;
t92 = Ifges(7,4) * t133 + Ifges(7,6) * t130;
t163 = t188 * pkin(9) ^ 2;
t99 = t189 * t134;
t10 = Ifges(6,4) * t44 + Ifges(6,2) * t43 + Ifges(6,6) * t70;
t7 = Ifges(7,5) * t44 + Ifges(7,6) * t70 - Ifges(7,3) * t43;
t161 = t7 / 0.2e1 - t10 / 0.2e1;
t11 = Ifges(7,1) * t44 + Ifges(7,4) * t70 - Ifges(7,5) * t43;
t12 = Ifges(6,1) * t44 + Ifges(6,4) * t43 + Ifges(6,5) * t70;
t160 = -t11 / 0.2e1 - t12 / 0.2e1;
t56 = Ifges(7,6) * t131 + (Ifges(7,3) * t133 - t173) * t134;
t59 = Ifges(6,6) * t131 + (-Ifges(6,2) * t133 - t175) * t134;
t159 = -t56 / 0.2e1 + t59 / 0.2e1;
t60 = Ifges(7,4) * t131 + (-Ifges(7,1) * t130 + t172) * t134;
t61 = Ifges(6,5) * t131 + (-Ifges(6,1) * t130 - t174) * t134;
t158 = t60 / 0.2e1 + t61 / 0.2e1;
t91 = Ifges(6,5) * t133 - Ifges(6,6) * t130;
t157 = t91 / 0.2e1 + t92 / 0.2e1;
t90 = Ifges(7,3) * t130 + t172;
t93 = -Ifges(6,2) * t130 + t174;
t156 = -t93 / 0.2e1 + t90 / 0.2e1;
t95 = Ifges(7,1) * t133 + t173;
t96 = Ifges(6,1) * t133 - t175;
t155 = t95 / 0.2e1 + t96 / 0.2e1;
t154 = Ifges(3,5) * t170 + Ifges(3,6) * t169 + Ifges(3,3) * t129;
t153 = t169 / 0.2e1;
t151 = (Ifges(5,4) - Ifges(4,5)) * t70 + (-Ifges(5,5) + Ifges(4,6)) * t69;
t1 = qJ(6) * t70 + t4;
t3 = t13 * t133 - t130 * t15;
t2 = -pkin(5) * t70 - t3;
t148 = t1 * t130 - t133 * t2;
t147 = t130 * t4 + t133 * t3;
t146 = t133 * mrSges(6,1) - t130 * mrSges(6,2);
t145 = t133 * mrSges(7,1) + t130 * mrSges(7,3);
t144 = pkin(5) * t133 + qJ(6) * t130;
t34 = qJ(6) * t131 + t37;
t36 = -t130 * t75 + t133 * t98;
t35 = -pkin(5) * t131 - t36;
t143 = t130 * t34 - t133 * t35;
t142 = t130 * t37 + t133 * t36;
t20 = qJ(4) * t169 - t27;
t46 = t70 * mrSges(5,1) - mrSges(5,2) * t169;
t141 = -0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * t162;
t16 = -pkin(4) * t69 - t20;
t139 = m(7) * t144 + t145 + t146;
t137 = qJ(4) ^ 2;
t112 = Ifges(6,3) * t131;
t97 = Ifges(4,1) * t131 + Ifges(4,4) * t134;
t94 = Ifges(4,4) * t131 + Ifges(4,2) * t134;
t89 = -Ifges(5,2) * t131 - Ifges(5,6) * t134;
t88 = -Ifges(5,6) * t131 - Ifges(5,3) * t134;
t87 = mrSges(6,1) * t130 + mrSges(6,2) * t133;
t86 = mrSges(7,1) * t130 - mrSges(7,3) * t133;
t85 = -mrSges(4,1) * t134 + mrSges(4,2) * t131;
t84 = mrSges(5,2) * t134 - mrSges(5,3) * t131;
t83 = -pkin(3) * t134 + t152;
t82 = pkin(5) * t130 - qJ(6) * t133 + qJ(4);
t74 = t146 * t134;
t73 = t145 * t134;
t58 = -Ifges(7,4) * t168 + t166;
t57 = t112 + (-Ifges(6,5) * t130 - t171) * t134;
t49 = t144 * t134 + t99;
t48 = -mrSges(4,1) * t169 - mrSges(4,3) * t70;
t47 = mrSges(4,2) * t169 - mrSges(4,3) * t69;
t45 = mrSges(5,1) * t69 + mrSges(5,3) * t169;
t33 = -mrSges(5,2) * t69 - mrSges(5,3) * t70;
t32 = mrSges(4,1) * t69 + mrSges(4,2) * t70;
t31 = Ifges(4,1) * t70 - Ifges(4,4) * t69 - Ifges(4,5) * t169;
t30 = Ifges(4,4) * t70 - Ifges(4,2) * t69 - Ifges(4,6) * t169;
t29 = -Ifges(5,4) * t169 - Ifges(5,2) * t70 + Ifges(5,6) * t69;
t28 = -Ifges(5,5) * t169 - Ifges(5,6) * t70 + Ifges(5,3) * t69;
t19 = pkin(3) * t69 + t140;
t18 = -mrSges(6,1) * t43 + mrSges(6,2) * t44;
t17 = -mrSges(7,1) * t43 - mrSges(7,3) * t44;
t5 = -pkin(5) * t43 - qJ(6) * t44 + t16;
t6 = [0.2e1 * t1 * t22 + 0.2e1 * t16 * t18 + 0.2e1 * t5 * t17 + 0.2e1 * t19 * t33 + 0.2e1 * t2 * t25 + 0.2e1 * t20 * t45 + 0.2e1 * t21 * t46 + 0.2e1 * t4 * t23 + 0.2e1 * t3 * t24 + 0.2e1 * t26 * t48 + 0.2e1 * t27 * t47 + 0.2e1 * t53 * t32 + Ifges(2,3) + (t28 - t30) * t69 + (t11 + t12) * t44 + (-t7 + t10) * t43 + (t154 - 0.2e1 * t181 + 0.2e1 * t182) * t129 + (-t29 + t31 + t187) * t70 + ((-0.2e1 * t71 * mrSges(3,3) + Ifges(3,5) * t129 + (-0.2e1 * mrSges(3,2) * pkin(1) + Ifges(3,1) * t132) * t128) * t132 + (0.2e1 * t72 * mrSges(3,3) + Ifges(3,6) * t129 + (0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * Ifges(3,4) * t132 + (Ifges(3,2) + t180) * t135) * t128 + t151) * t135) * t128 + m(3) * (pkin(1) ^ 2 * t128 ^ 2 + t71 ^ 2 + t72 ^ 2) + m(4) * (t26 ^ 2 + t27 ^ 2 + t53 ^ 2) + m(6) * (t16 ^ 2 + t3 ^ 2 + t4 ^ 2) + m(5) * (t19 ^ 2 + t20 ^ 2 + t21 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2); -t164 * t169 / 0.2e1 + (-t28 / 0.2e1 + t30 / 0.2e1 + Ifges(5,5) * t153 + t27 * mrSges(4,3) - t20 * mrSges(5,1) + t161 * t133 + t160 * t130 + (-t45 + t47) * pkin(9)) * t134 + t158 * t44 + t159 * t43 + (Ifges(5,4) * t153 - t26 * mrSges(4,3) + t21 * mrSges(5,1) - t29 / 0.2e1 + t31 / 0.2e1 + t8 / 0.2e1 + t9 / 0.2e1 + (t46 - t48) * pkin(9)) * t131 + t154 + m(5) * (t19 * t83 + (t131 * t21 - t134 * t20) * pkin(9)) + m(4) * (-pkin(2) * t53 + (-t131 * t26 + t134 * t27) * pkin(9)) + t19 * t84 + t53 * t85 + t99 * t18 + t5 * t73 + t16 * t74 + t3 * t78 + t2 * t79 + t4 * t80 + t1 * t81 + t83 * t33 + t49 * t17 + t34 * t22 + t35 * t25 + t36 * t24 + t37 * t23 - pkin(2) * t32 + m(6) * (t16 * t99 + t3 * t36 + t37 * t4) + m(7) * (t1 * t34 + t2 * t35 + t49 * t5) + (-t89 / 0.2e1 + t97 / 0.2e1 + t57 / 0.2e1 + t58 / 0.2e1) * t70 + (t88 / 0.2e1 - t94 / 0.2e1) * t69 - t181 + t182; -0.2e1 * pkin(2) * t85 + 0.2e1 * t34 * t81 + 0.2e1 * t35 * t79 + 0.2e1 * t36 * t78 + 0.2e1 * t37 * t80 + 0.2e1 * t49 * t73 + 0.2e1 * t99 * t74 + 0.2e1 * t83 * t84 + Ifges(3,3) + (-t89 + t97 + t57 + t58) * t131 + m(5) * (t83 ^ 2 + t163) + m(4) * (pkin(2) ^ 2 + t163) + m(6) * (t36 ^ 2 + t37 ^ 2 + t99 ^ 2) + m(7) * (t34 ^ 2 + t35 ^ 2 + t49 ^ 2) + (-t88 + t94 + (t56 - t59) * t133 + (-t60 - t61) * t130) * t134 + 0.2e1 * (mrSges(5,1) + mrSges(4,3)) * pkin(9) * t188; (-t1 * mrSges(7,2) - t4 * mrSges(6,3) - t179 * t184 + t161) * t130 + (t2 * mrSges(7,2) - t3 * mrSges(6,3) - t178 * t184 - t160) * t133 - t151 - t180 * t169 + t5 * t86 + t16 * t87 + t82 * t17 - pkin(3) * t46 + t26 * mrSges(4,1) - t27 * mrSges(4,2) - t20 * mrSges(5,3) + t21 * mrSges(5,2) + (-t45 + t18) * qJ(4) + m(5) * (-pkin(3) * t21 - qJ(4) * t20) + t157 * t70 + t155 * t44 - t156 * t43 + m(6) * (qJ(4) * t16 - t147 * t184) + m(7) * (-t148 * t184 + t5 * t82); qJ(4) * t74 + t49 * t86 + t82 * t73 + t99 * t87 + (t35 * mrSges(7,2) - t36 * mrSges(6,3) - t177 * t184 + t158) * t133 + (-t34 * mrSges(7,2) - t37 * mrSges(6,3) - t176 * t184 - t159) * t130 + m(6) * (qJ(4) * t99 - t142 * t184) + m(7) * (-t143 * t184 + t49 * t82) + (-pkin(3) * mrSges(5,1) - Ifges(5,4) + t157) * t131 + (qJ(4) * mrSges(5,1) - t155 * t130 + t156 * t133 - Ifges(5,5)) * t134 + ((m(5) * qJ(4) - mrSges(4,2) + mrSges(5,3)) * t134 + (-mrSges(4,1) + t186) * t131) * pkin(9) + t164; -0.2e1 * pkin(3) * mrSges(5,2) + 0.2e1 * t82 * t86 + (t96 + t95) * t133 + (-t93 + t90) * t130 + 0.2e1 * (t87 + mrSges(5,3)) * qJ(4) + m(7) * (t82 ^ 2 - t165) + m(6) * (t137 - t165) + m(5) * (pkin(3) ^ 2 + t137) + t180 - 0.2e1 * t184 * t185; m(5) * t21 + m(6) * t147 + m(7) * t148 + t179 * t130 + t178 * t133 + t46; t177 * t133 + (m(5) * pkin(9) + mrSges(5,1)) * t131 + t176 * t130 + m(7) * t143 + m(6) * t142; -t141 * t184 + t185 + t186; m(5) + t141; -pkin(5) * t25 + m(7) * (-pkin(5) * t2 + qJ(6) * t1) + qJ(6) * t22 + t1 * mrSges(7,3) - t4 * mrSges(6,2) + t3 * mrSges(6,1) - t2 * mrSges(7,1) + t187; -pkin(5) * t79 + m(7) * (-pkin(5) * t35 + qJ(6) * t34) + qJ(6) * t81 + t34 * mrSges(7,3) - t35 * mrSges(7,1) - t37 * mrSges(6,2) + t36 * mrSges(6,1) + t112 + (-t171 + (-Ifges(7,4) - Ifges(6,5)) * t130) * t134 + t166; -t144 * mrSges(7,2) - t139 * t184 + t91 + t92; t139; Ifges(7,2) + Ifges(6,3) + 0.2e1 * pkin(5) * mrSges(7,1) + 0.2e1 * qJ(6) * mrSges(7,3) + m(7) * (pkin(5) ^ 2 + qJ(6) ^ 2); m(7) * t2 + t25; m(7) * t35 + t79; (m(7) * t184 + mrSges(7,2)) * t133; -m(7) * t133; -m(7) * pkin(5) - mrSges(7,1); m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t6(1) t6(2) t6(4) t6(7) t6(11) t6(16); t6(2) t6(3) t6(5) t6(8) t6(12) t6(17); t6(4) t6(5) t6(6) t6(9) t6(13) t6(18); t6(7) t6(8) t6(9) t6(10) t6(14) t6(19); t6(11) t6(12) t6(13) t6(14) t6(15) t6(20); t6(16) t6(17) t6(18) t6(19) t6(20) t6(21);];
Mq  = res;
