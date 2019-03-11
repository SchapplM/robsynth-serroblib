% Calculate joint inertia matrix for
% S6PRRRRP6
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
% Datum: 2019-03-09 00:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRRP6_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP6_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP6_inertiaJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP6_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP6_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP6_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:28:42
% EndTime: 2019-03-09 00:28:45
% DurationCPUTime: 1.63s
% Computational Cost: add. (1957->379), mult. (4789->521), div. (0->0), fcn. (5103->12), ass. (0->140)
t133 = sin(qJ(5));
t137 = cos(qJ(5));
t192 = t133 ^ 2 + t137 ^ 2;
t134 = sin(qJ(4));
t138 = cos(qJ(4));
t130 = sin(pkin(6));
t132 = cos(pkin(6));
t135 = sin(qJ(3));
t136 = sin(qJ(2));
t139 = cos(qJ(3));
t131 = cos(pkin(7));
t140 = cos(qJ(2));
t164 = t131 * t140;
t129 = sin(pkin(7));
t166 = t129 * t135;
t40 = t132 * t166 + (t135 * t164 + t136 * t139) * t130;
t72 = -t129 * t130 * t140 + t131 * t132;
t24 = t134 * t72 + t138 * t40;
t165 = t129 * t139;
t38 = -t132 * t165 + (t135 * t136 - t139 * t164) * t130;
t10 = t133 * t38 + t137 * t24;
t8 = t133 * t24 - t38 * t137;
t191 = t10 * t137 + t133 * t8;
t190 = 2 * pkin(10);
t189 = mrSges(6,3) + mrSges(7,2);
t74 = t131 * t134 + t138 * t166;
t46 = t133 * t74 + t137 * t165;
t47 = -t133 * t165 + t137 * t74;
t73 = -t138 * t131 + t134 * t166;
t12 = Ifges(6,5) * t47 - Ifges(6,6) * t46 + Ifges(6,3) * t73;
t13 = Ifges(7,4) * t47 + Ifges(7,2) * t73 + Ifges(7,6) * t46;
t188 = t12 + t13;
t187 = -m(7) * pkin(5) - mrSges(7,1);
t22 = t134 * t40 - t72 * t138;
t21 = t22 ^ 2;
t186 = t38 ^ 2;
t185 = pkin(2) * t139;
t184 = pkin(10) * t138;
t182 = -Ifges(7,2) - Ifges(6,3);
t104 = pkin(9) * t166;
t58 = t104 + (-pkin(3) - t185) * t131;
t25 = pkin(4) * t73 - pkin(11) * t74 + t58;
t77 = t131 * t135 * pkin(2) + pkin(9) * t165;
t59 = pkin(10) * t131 + t77;
t60 = (-pkin(3) * t139 - pkin(10) * t135 - pkin(2)) * t129;
t33 = t134 * t60 + t138 * t59;
t27 = -pkin(11) * t165 + t33;
t4 = t133 * t25 + t137 * t27;
t28 = -mrSges(7,2) * t46 + mrSges(7,3) * t73;
t29 = -mrSges(6,2) * t73 - mrSges(6,3) * t46;
t181 = t28 + t29;
t30 = mrSges(6,1) * t73 - mrSges(6,3) * t47;
t31 = -t73 * mrSges(7,1) + t47 * mrSges(7,2);
t180 = -t30 + t31;
t18 = mrSges(6,1) * t46 + mrSges(6,2) * t47;
t49 = -mrSges(5,1) * t165 - mrSges(5,3) * t74;
t179 = -t49 + t18;
t178 = -Ifges(5,5) * t74 + Ifges(5,6) * t73;
t163 = t133 * t134;
t83 = mrSges(6,2) * t138 - mrSges(6,3) * t163;
t86 = -mrSges(7,2) * t163 - mrSges(7,3) * t138;
t177 = t83 + t86;
t162 = t134 * t137;
t84 = -mrSges(6,1) * t138 - mrSges(6,3) * t162;
t85 = t138 * mrSges(7,1) + mrSges(7,2) * t162;
t176 = -t84 + t85;
t91 = -mrSges(6,1) * t137 + mrSges(6,2) * t133;
t175 = t91 - mrSges(5,1);
t174 = Ifges(6,4) * t133;
t173 = Ifges(6,4) * t137;
t172 = Ifges(7,5) * t133;
t171 = Ifges(7,5) * t137;
t169 = t134 * t22;
t89 = -pkin(4) * t138 - pkin(11) * t134 - pkin(3);
t168 = t137 * t89;
t167 = t138 * t24;
t56 = t133 * t89 + t137 * t184;
t161 = Ifges(7,4) * t162 + Ifges(7,6) * t163;
t94 = Ifges(6,5) * t133 + Ifges(6,6) * t137;
t160 = Ifges(5,5) * t134 + Ifges(5,6) * t138;
t159 = t192 * pkin(11) ^ 2;
t11 = Ifges(7,5) * t47 + Ifges(7,6) * t73 + Ifges(7,3) * t46;
t14 = Ifges(6,4) * t47 - Ifges(6,2) * t46 + Ifges(6,6) * t73;
t158 = t11 / 0.2e1 - t14 / 0.2e1;
t15 = Ifges(7,1) * t47 + Ifges(7,4) * t73 + Ifges(7,5) * t46;
t16 = Ifges(6,1) * t47 - Ifges(6,4) * t46 + Ifges(6,5) * t73;
t157 = t15 / 0.2e1 + t16 / 0.2e1;
t61 = -Ifges(7,6) * t138 + (Ifges(7,3) * t133 + t171) * t134;
t64 = -Ifges(6,6) * t138 + (-Ifges(6,2) * t133 + t173) * t134;
t156 = t61 / 0.2e1 - t64 / 0.2e1;
t65 = -Ifges(7,4) * t138 + (Ifges(7,1) * t137 + t172) * t134;
t66 = -Ifges(6,5) * t138 + (Ifges(6,1) * t137 - t174) * t134;
t155 = t65 / 0.2e1 + t66 / 0.2e1;
t93 = -Ifges(7,3) * t137 + t172;
t96 = Ifges(6,2) * t137 + t174;
t154 = t93 / 0.2e1 - t96 / 0.2e1;
t95 = Ifges(7,4) * t133 - Ifges(7,6) * t137;
t153 = t94 / 0.2e1 + t95 / 0.2e1;
t98 = Ifges(7,1) * t133 - t171;
t99 = Ifges(6,1) * t133 + t173;
t152 = t98 / 0.2e1 + t99 / 0.2e1;
t150 = Ifges(4,5) * t166 + Ifges(4,6) * t165 + Ifges(4,3) * t131;
t32 = -t134 * t59 + t138 * t60;
t149 = t191 * pkin(11);
t26 = pkin(4) * t165 - t32;
t147 = Ifges(6,5) * t162 - Ifges(6,6) * t163;
t145 = t133 * mrSges(6,1) + t137 * mrSges(6,2);
t144 = t133 * mrSges(7,1) - t137 * mrSges(7,3);
t143 = -pkin(5) * t133 + qJ(6) * t137;
t3 = -t133 * t27 + t137 * t25;
t142 = pkin(10) ^ 2;
t128 = t138 ^ 2;
t126 = t134 ^ 2;
t122 = t126 * t142;
t100 = Ifges(5,1) * t134 + Ifges(5,4) * t138;
t97 = Ifges(5,4) * t134 + Ifges(5,2) * t138;
t92 = -mrSges(5,1) * t138 + mrSges(5,2) * t134;
t90 = -mrSges(7,1) * t137 - mrSges(7,3) * t133;
t88 = -pkin(5) * t137 - qJ(6) * t133 - pkin(4);
t82 = -mrSges(4,2) * t131 + mrSges(4,3) * t165;
t81 = mrSges(4,1) * t131 - mrSges(4,3) * t166;
t79 = t145 * t134;
t78 = t144 * t134;
t76 = t131 * t185 - t104;
t75 = (-mrSges(4,1) * t139 + mrSges(4,2) * t135) * t129;
t67 = (pkin(10) - t143) * t134;
t63 = -Ifges(7,2) * t138 + t161;
t62 = -Ifges(6,3) * t138 + t147;
t55 = -t133 * t184 + t168;
t51 = -t168 + (pkin(10) * t133 + pkin(5)) * t138;
t50 = -qJ(6) * t138 + t56;
t48 = mrSges(5,2) * t165 - mrSges(5,3) * t73;
t36 = mrSges(5,1) * t73 + mrSges(5,2) * t74;
t35 = Ifges(5,1) * t74 - Ifges(5,4) * t73 - Ifges(5,5) * t165;
t34 = Ifges(5,4) * t74 - Ifges(5,2) * t73 - Ifges(5,6) * t165;
t17 = mrSges(7,1) * t46 - mrSges(7,3) * t47;
t5 = pkin(5) * t46 - qJ(6) * t47 + t26;
t2 = -pkin(5) * t73 - t3;
t1 = qJ(6) * t73 + t4;
t6 = [m(2) + m(5) * (t24 ^ 2 + t186 + t21) + m(4) * (t40 ^ 2 + t72 ^ 2 + t186) + m(3) * (t132 ^ 2 + (t136 ^ 2 + t140 ^ 2) * t130 ^ 2) + (m(6) + m(7)) * (t10 ^ 2 + t8 ^ 2 + t21); t24 * t48 + t40 * t82 + t72 * t75 + t180 * t8 + (t36 - t81) * t38 + (mrSges(3,1) * t140 - mrSges(3,2) * t136) * t130 + t181 * t10 + (t17 + t179) * t22 + m(6) * (t10 * t4 + t22 * t26 - t3 * t8) + m(7) * (t1 * t10 + t2 * t8 + t22 * t5) + m(5) * (-t22 * t32 + t24 * t33 + t38 * t58) + m(4) * (-pkin(2) * t129 * t72 - t38 * t76 + t40 * t77); t131 * t150 + 0.2e1 * t76 * t81 + 0.2e1 * t77 * t82 + t74 * t35 + 0.2e1 * t58 * t36 + 0.2e1 * t33 * t48 + 0.2e1 * t32 * t49 + 0.2e1 * t1 * t28 + 0.2e1 * t4 * t29 + 0.2e1 * t3 * t30 + 0.2e1 * t2 * t31 + 0.2e1 * t5 * t17 + 0.2e1 * t26 * t18 + Ifges(3,3) + (t15 + t16) * t47 + (-t14 + t11) * t46 + (-t34 + t188) * t73 + (-0.2e1 * pkin(2) * t75 + (Ifges(4,1) * t166 + Ifges(4,5) * t131) * t135 + (Ifges(4,6) * t131 + (0.2e1 * Ifges(4,4) * t135 + (Ifges(4,2) + Ifges(5,3)) * t139) * t129 + t178) * t139) * t129 + m(4) * (pkin(2) ^ 2 * t129 ^ 2 + t76 ^ 2 + t77 ^ 2) + m(5) * (t32 ^ 2 + t33 ^ 2 + t58 ^ 2) + m(6) * (t26 ^ 2 + t3 ^ 2 + t4 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2); mrSges(5,3) * t167 - t40 * mrSges(4,2) + t176 * t8 + (-mrSges(4,1) + t92) * t38 + t177 * t10 + (t134 * mrSges(5,3) + t78 + t79) * t22 + m(6) * (pkin(10) * t169 + t10 * t56 - t55 * t8) + m(7) * (t10 * t50 + t22 * t67 + t51 * t8) + m(5) * (-pkin(3) * t38 + (t167 + t169) * pkin(10)); (t35 / 0.2e1 - t32 * mrSges(5,3) + t157 * t137 + t158 * t133 + t179 * pkin(10)) * t134 - t160 * t165 / 0.2e1 + t155 * t47 + t156 * t46 + m(5) * (-pkin(3) * t58 + (-t32 * t134 + t33 * t138) * pkin(10)) + t74 * t100 / 0.2e1 + t4 * t83 + t3 * t84 + t2 * t85 + t1 * t86 + t58 * t92 + t76 * mrSges(4,1) - t77 * mrSges(4,2) + t5 * t78 + t26 * t79 + t55 * t30 + t56 * t29 + t67 * t17 + t50 * t28 + t51 * t31 - pkin(3) * t36 + (t34 / 0.2e1 - t12 / 0.2e1 - t13 / 0.2e1 + pkin(10) * t48 + t33 * mrSges(5,3)) * t138 + m(6) * (pkin(10) * t134 * t26 + t3 * t55 + t4 * t56) + m(7) * (t1 * t50 + t2 * t51 + t5 * t67) + (-t97 / 0.2e1 + t62 / 0.2e1 + t63 / 0.2e1) * t73 + t150; -0.2e1 * pkin(3) * t92 + 0.2e1 * t50 * t86 + 0.2e1 * t51 * t85 + 0.2e1 * t55 * t84 + 0.2e1 * t56 * t83 + 0.2e1 * t67 * t78 + Ifges(4,3) + (t126 + t128) * mrSges(5,3) * t190 + (-t63 - t62 + t97) * t138 + m(7) * (t50 ^ 2 + t51 ^ 2 + t67 ^ 2) + m(6) * (t55 ^ 2 + t56 ^ 2 + t122) + m(5) * (pkin(3) ^ 2 + t128 * t142 + t122) + (t79 * t190 + t100 + (t65 + t66) * t137 + (t61 - t64) * t133) * t134; -t24 * mrSges(5,2) + (t90 + t175) * t22 + m(6) * (-pkin(4) * t22 + t149) + m(7) * (t22 * t88 + t149) + t189 * t191; -Ifges(5,3) * t165 + t32 * mrSges(5,1) - t33 * mrSges(5,2) - pkin(4) * t18 + t88 * t17 + t26 * t91 + t5 * t90 + t153 * t73 + t152 * t47 + t154 * t46 + (t1 * mrSges(7,2) + t4 * mrSges(6,3) + pkin(11) * t181 - t158) * t137 + (t2 * mrSges(7,2) - t3 * mrSges(6,3) + pkin(11) * t180 + t157) * t133 + m(7) * (t5 * t88 + (t1 * t137 + t2 * t133) * pkin(11)) + m(6) * (-pkin(4) * t26 + (-t3 * t133 + t4 * t137) * pkin(11)) - t178; -pkin(4) * t79 + t88 * t78 + (m(7) * t88 + t90) * t67 + (-pkin(10) * mrSges(5,2) - t153) * t138 + (t50 * mrSges(7,2) + t56 * mrSges(6,3) - t156) * t137 + (t51 * mrSges(7,2) - t55 * mrSges(6,3) + t155) * t133 + (t177 * t137 + t176 * t133 + m(6) * (-t133 * t55 + t137 * t56) + m(7) * (t133 * t51 + t137 * t50)) * pkin(11) + (t152 * t137 + t154 * t133 + (-m(6) * pkin(4) + t175) * pkin(10)) * t134 + t160; -0.2e1 * pkin(4) * t91 + 0.2e1 * t88 * t90 + Ifges(5,3) + (-t93 + t96) * t137 + (t99 + t98) * t133 + m(7) * (t88 ^ 2 + t159) + m(6) * (pkin(4) ^ 2 + t159) + 0.2e1 * t189 * pkin(11) * t192; (-mrSges(6,1) + t187) * t8 + (m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3)) * t10; t1 * mrSges(7,3) + qJ(6) * t28 + m(7) * (-pkin(5) * t2 + qJ(6) * t1) + t3 * mrSges(6,1) - t2 * mrSges(7,1) - t4 * mrSges(6,2) - pkin(5) * t31 + t188; qJ(6) * t86 + m(7) * (-pkin(5) * t51 + qJ(6) * t50) + t50 * mrSges(7,3) - t51 * mrSges(7,1) - t56 * mrSges(6,2) + t55 * mrSges(6,1) - pkin(5) * t85 + t182 * t138 + t147 + t161; t143 * mrSges(7,2) + (m(7) * t143 - t144 - t145) * pkin(11) + t95 + t94; 0.2e1 * pkin(5) * mrSges(7,1) + 0.2e1 * qJ(6) * mrSges(7,3) + m(7) * (pkin(5) ^ 2 + qJ(6) ^ 2) - t182; m(7) * t8; m(7) * t2 + t31; m(7) * t51 + t85; (m(7) * pkin(11) + mrSges(7,2)) * t133; t187; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t6(1) t6(2) t6(4) t6(7) t6(11) t6(16); t6(2) t6(3) t6(5) t6(8) t6(12) t6(17); t6(4) t6(5) t6(6) t6(9) t6(13) t6(18); t6(7) t6(8) t6(9) t6(10) t6(14) t6(19); t6(11) t6(12) t6(13) t6(14) t6(15) t6(20); t6(16) t6(17) t6(18) t6(19) t6(20) t6(21);];
Mq  = res;
