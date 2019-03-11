% Calculate joint inertia matrix for
% S6RRPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 12:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRP5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP5_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP5_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP5_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP5_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP5_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:57:57
% EndTime: 2019-03-09 11:58:01
% DurationCPUTime: 1.77s
% Computational Cost: add. (2816->383), mult. (6733->520), div. (0->0), fcn. (7234->10), ass. (0->141)
t186 = Ifges(3,3) + Ifges(4,3);
t185 = -2 * mrSges(7,3);
t132 = sin(pkin(11));
t116 = pkin(2) * t132 + pkin(9);
t184 = 0.2e1 * t116;
t136 = sin(qJ(5));
t139 = cos(qJ(5));
t135 = cos(pkin(6));
t137 = sin(qJ(4));
t140 = cos(qJ(4));
t133 = sin(pkin(6));
t134 = cos(pkin(11));
t138 = sin(qJ(2));
t141 = cos(qJ(2));
t68 = (t132 * t141 + t134 * t138) * t133;
t56 = t135 * t137 + t140 * t68;
t160 = t133 * t141;
t161 = t133 * t138;
t67 = t132 * t161 - t134 * t160;
t41 = -t136 * t56 + t139 * t67;
t42 = t136 * t67 + t139 * t56;
t55 = -t135 * t140 + t137 * t68;
t6 = Ifges(7,5) * t42 + Ifges(7,6) * t41 + Ifges(7,3) * t55;
t7 = Ifges(6,5) * t42 + Ifges(6,6) * t41 + Ifges(6,3) * t55;
t183 = t6 + t7;
t182 = Ifges(4,5) * t68 - Ifges(4,6) * t67;
t77 = (pkin(5) * t136 + t116) * t137;
t158 = t137 * t139;
t159 = t136 * t137;
t81 = mrSges(7,1) * t159 + mrSges(7,2) * t158;
t181 = m(7) * t77 + t81;
t84 = (-pkin(2) * t141 - pkin(1)) * t133;
t180 = 0.2e1 * t84;
t179 = m(6) * pkin(4);
t178 = m(7) * pkin(5);
t176 = pkin(1) * t135;
t175 = pkin(4) * t140;
t174 = pkin(10) * t137;
t109 = t141 * t176;
t79 = -pkin(8) * t161 + t109;
t173 = t79 * mrSges(3,1);
t80 = pkin(8) * t160 + t138 * t176;
t172 = t80 * mrSges(3,2);
t92 = -mrSges(6,1) * t139 + mrSges(6,2) * t136;
t171 = mrSges(5,1) - t92;
t170 = -Ifges(6,3) - Ifges(7,3);
t169 = -qJ(6) - pkin(10);
t57 = pkin(2) * t135 + t109 + (-pkin(8) - qJ(3)) * t161;
t62 = qJ(3) * t160 + t80;
t39 = t132 * t57 + t134 * t62;
t32 = pkin(9) * t135 + t39;
t43 = pkin(3) * t67 - pkin(9) * t68 + t84;
t19 = t137 * t43 + t140 * t32;
t14 = pkin(10) * t67 + t19;
t38 = -t132 * t62 + t134 * t57;
t31 = -pkin(3) * t135 - t38;
t17 = pkin(4) * t55 - pkin(10) * t56 + t31;
t4 = t136 * t17 + t139 * t14;
t21 = -mrSges(6,1) * t41 + mrSges(6,2) * t42;
t45 = mrSges(5,1) * t67 - mrSges(5,3) * t56;
t168 = -t45 + t21;
t162 = t116 * t140;
t117 = -pkin(2) * t134 - pkin(3);
t83 = t117 - t174 - t175;
t54 = t136 * t83 + t139 * t162;
t167 = mrSges(6,2) * t139;
t166 = Ifges(6,4) * t136;
t165 = Ifges(6,4) * t139;
t164 = Ifges(7,4) * t136;
t163 = Ifges(7,4) * t139;
t95 = Ifges(7,5) * t136 + Ifges(7,6) * t139;
t96 = Ifges(6,5) * t136 + Ifges(6,6) * t139;
t157 = Ifges(5,5) * t137 + Ifges(5,6) * t140;
t156 = t136 ^ 2 + t139 ^ 2;
t129 = t137 ^ 2;
t131 = t140 ^ 2;
t155 = t129 + t131;
t8 = Ifges(7,4) * t42 + Ifges(7,2) * t41 + Ifges(7,6) * t55;
t9 = Ifges(6,4) * t42 + Ifges(6,2) * t41 + Ifges(6,6) * t55;
t154 = -t8 / 0.2e1 - t9 / 0.2e1;
t153 = Ifges(5,5) * t56 - Ifges(5,6) * t55 + Ifges(5,3) * t67;
t10 = Ifges(7,1) * t42 + Ifges(7,4) * t41 + Ifges(7,5) * t55;
t11 = Ifges(6,1) * t42 + Ifges(6,4) * t41 + Ifges(6,5) * t55;
t152 = t10 / 0.2e1 + t11 / 0.2e1;
t71 = -Ifges(7,6) * t140 + (-Ifges(7,2) * t136 + t163) * t137;
t72 = -Ifges(6,6) * t140 + (-Ifges(6,2) * t136 + t165) * t137;
t151 = t71 / 0.2e1 + t72 / 0.2e1;
t73 = -Ifges(7,5) * t140 + (Ifges(7,1) * t139 - t164) * t137;
t74 = -Ifges(6,5) * t140 + (Ifges(6,1) * t139 - t166) * t137;
t150 = t73 / 0.2e1 + t74 / 0.2e1;
t149 = t95 / 0.2e1 + t96 / 0.2e1;
t97 = Ifges(7,2) * t139 + t164;
t98 = Ifges(6,2) * t139 + t166;
t148 = t97 / 0.2e1 + t98 / 0.2e1;
t100 = Ifges(7,1) * t136 + t163;
t101 = Ifges(6,1) * t136 + t165;
t147 = t100 / 0.2e1 + t101 / 0.2e1;
t20 = -t41 * mrSges(7,1) + t42 * mrSges(7,2);
t3 = -t136 * t14 + t139 * t17;
t18 = -t137 * t32 + t140 * t43;
t146 = t156 * mrSges(6,3);
t91 = -t139 * mrSges(7,1) + t136 * mrSges(7,2);
t145 = mrSges(6,1) * t136 + t167;
t13 = -pkin(4) * t67 - t18;
t144 = Ifges(3,5) * t161 + Ifges(3,6) * t160 + t135 * t186 + t182;
t118 = -pkin(5) * t139 - pkin(4);
t114 = t116 ^ 2;
t113 = Ifges(6,5) * t158;
t112 = Ifges(7,5) * t158;
t103 = t129 * t114;
t102 = Ifges(5,1) * t137 + Ifges(5,4) * t140;
t99 = Ifges(5,4) * t137 + Ifges(5,2) * t140;
t94 = t169 * t139;
t93 = -t140 * mrSges(5,1) + t137 * mrSges(5,2);
t90 = t169 * t136;
t88 = -mrSges(6,1) * t140 - mrSges(6,3) * t158;
t87 = -mrSges(7,1) * t140 - mrSges(7,3) * t158;
t86 = mrSges(6,2) * t140 - mrSges(6,3) * t159;
t85 = mrSges(7,2) * t140 - mrSges(7,3) * t159;
t82 = t145 * t137;
t76 = t139 * t83;
t70 = -Ifges(6,6) * t159 - Ifges(6,3) * t140 + t113;
t69 = -Ifges(7,6) * t159 - Ifges(7,3) * t140 + t112;
t63 = t68 * mrSges(4,2);
t59 = mrSges(4,1) * t135 - mrSges(4,3) * t68;
t58 = -mrSges(4,2) * t135 - mrSges(4,3) * t67;
t53 = -t136 * t162 + t76;
t47 = -qJ(6) * t159 + t54;
t46 = -qJ(6) * t158 + t76 + (-t116 * t136 - pkin(5)) * t140;
t44 = -mrSges(5,2) * t67 - mrSges(5,3) * t55;
t30 = mrSges(5,1) * t55 + mrSges(5,2) * t56;
t27 = Ifges(5,1) * t56 - Ifges(5,4) * t55 + Ifges(5,5) * t67;
t26 = Ifges(5,4) * t56 - Ifges(5,2) * t55 + Ifges(5,6) * t67;
t25 = mrSges(6,1) * t55 - mrSges(6,3) * t42;
t24 = mrSges(7,1) * t55 - mrSges(7,3) * t42;
t23 = -mrSges(6,2) * t55 + mrSges(6,3) * t41;
t22 = -mrSges(7,2) * t55 + mrSges(7,3) * t41;
t5 = -pkin(5) * t41 + t13;
t2 = qJ(6) * t41 + t4;
t1 = pkin(5) * t55 - qJ(6) * t42 + t3;
t12 = [t56 * t27 + 0.2e1 * t39 * t58 + 0.2e1 * t38 * t59 + 0.2e1 * t19 * t44 + 0.2e1 * t18 * t45 + 0.2e1 * t2 * t22 + 0.2e1 * t4 * t23 + 0.2e1 * t1 * t24 + 0.2e1 * t3 * t25 + 0.2e1 * t31 * t30 + 0.2e1 * t5 * t20 + 0.2e1 * t13 * t21 + ((Ifges(3,5) * t138 + Ifges(3,6) * t141) * t135 + 0.2e1 * (-t138 * t79 + t141 * t80) * mrSges(3,3) + (t141 * (Ifges(3,4) * t138 + Ifges(3,2) * t141) - 0.2e1 * pkin(1) * (-mrSges(3,1) * t141 + mrSges(3,2) * t138) + t138 * (Ifges(3,1) * t138 + Ifges(3,4) * t141) + m(3) * pkin(1) ^ 2) * t133) * t133 + (mrSges(4,1) * t180 - 0.2e1 * Ifges(4,4) * t68 + Ifges(4,2) * t67 + t153) * t67 + m(5) * (t18 ^ 2 + t19 ^ 2 + t31 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) + m(6) * (t13 ^ 2 + t3 ^ 2 + t4 ^ 2) + m(4) * (t38 ^ 2 + t39 ^ 2 + t84 ^ 2) + (t10 + t11) * t42 + (t8 + t9) * t41 + (-t26 + t183) * t55 + (t144 - 0.2e1 * t172 + 0.2e1 * t173 + t182) * t135 + m(3) * (t79 ^ 2 + t80 ^ 2) + Ifges(2,3) + t63 * t180 + Ifges(4,1) * t68 ^ 2; t144 + t56 * t102 / 0.2e1 + t117 * t30 + t2 * t85 + t4 * t86 + t1 * t87 + t3 * t88 + t31 * t93 + t77 * t20 + t5 * t81 + t13 * t82 + t53 * t25 + t54 * t23 + t46 * t24 + t47 * t22 + t38 * mrSges(4,1) - t39 * mrSges(4,2) - t172 + t173 + (t26 / 0.2e1 - t6 / 0.2e1 - t7 / 0.2e1 + t116 * t44 + t19 * mrSges(5,3)) * t140 + (t27 / 0.2e1 - t18 * mrSges(5,3) + t152 * t139 + t154 * t136 + t168 * t116) * t137 + t67 * t157 / 0.2e1 + t150 * t42 + t151 * t41 + m(6) * (t116 * t13 * t137 + t3 * t53 + t4 * t54) + m(7) * (t1 * t46 + t2 * t47 + t5 * t77) + (t132 * t58 + t134 * t59 + m(4) * (t132 * t39 + t134 * t38)) * pkin(2) + (-t99 / 0.2e1 + t69 / 0.2e1 + t70 / 0.2e1) * t55 + m(5) * (t117 * t31 + (-t137 * t18 + t140 * t19) * t116); 0.2e1 * t117 * t93 + 0.2e1 * t46 * t87 + 0.2e1 * t47 * t85 + 0.2e1 * t53 * t88 + 0.2e1 * t54 * t86 + 0.2e1 * t77 * t81 + (t99 - t69 - t70) * t140 + m(5) * (t114 * t131 + t117 ^ 2 + t103) + m(6) * (t53 ^ 2 + t54 ^ 2 + t103) + m(7) * (t46 ^ 2 + t47 ^ 2 + t77 ^ 2) + m(4) * (t132 ^ 2 + t134 ^ 2) * pkin(2) ^ 2 + (t82 * t184 + t102 + (t73 + t74) * t139 + (-t71 - t72) * t136) * t137 + 0.2e1 * (mrSges(4,1) * t134 - mrSges(4,2) * t132) * pkin(2) + t155 * mrSges(5,3) * t184 + t186; t67 * mrSges(4,1) + t63 + (-t20 - t168) * t140 + (t44 + (t22 + t23) * t139 + (-t24 - t25) * t136) * t137 + m(7) * (-t140 * t5 + (-t1 * t136 + t139 * t2) * t137) + m(6) * (-t13 * t140 + (-t136 * t3 + t139 * t4) * t137) + m(5) * (t137 * t19 + t140 * t18) + m(4) * t84; (-t82 - t181) * t140 + ((t85 + t86) * t139 + (-t87 - t88) * t136 + m(7) * (-t136 * t46 + t139 * t47) + m(6) * (-t136 * t53 + t139 * t54 - t162)) * t137; m(4) + m(5) * t155 + 0.2e1 * (m(7) / 0.2e1 + m(6) / 0.2e1) * (t156 * t129 + t131); t18 * mrSges(5,1) - t19 * mrSges(5,2) - pkin(4) * t21 + t118 * t20 - t94 * t22 + t90 * t24 + t5 * t91 + t149 * t55 + t147 * t42 + t148 * t41 + m(7) * (t1 * t90 + t118 * t5 - t2 * t94) + (t2 * mrSges(7,3) + t4 * mrSges(6,3) + (m(6) * t4 + t23) * pkin(10) - t154) * t139 + (-t1 * mrSges(7,3) - t3 * mrSges(6,3) + (-m(6) * t3 - t25) * pkin(10) + t152) * t136 + t153 + (t92 - t179) * t13; t118 * t81 + t90 * t87 + t77 * t91 - t94 * t85 - pkin(4) * t82 + m(7) * (t118 * t77 + t46 * t90 - t47 * t94) - t149 * t140 + (-t140 * mrSges(5,2) + (-t171 - t179) * t137) * t116 + (t54 * mrSges(6,3) + t47 * mrSges(7,3) + t147 * t137 + (m(6) * t54 + t86) * pkin(10) + t151) * t139 + (-t53 * mrSges(6,3) - t46 * mrSges(7,3) - t148 * t137 + (-m(6) * t53 - t88) * pkin(10) + t150) * t136 + t157; (-t91 + t171) * t140 + (t156 * mrSges(7,3) - mrSges(5,2) + t146) * t137 + m(7) * (-t118 * t140 + (-t136 * t90 - t139 * t94) * t137) + m(6) * (t156 * t174 + t175); -0.2e1 * pkin(4) * t92 + 0.2e1 * t118 * t91 + Ifges(5,3) + 0.2e1 * pkin(10) * t146 + m(7) * (t118 ^ 2 + t90 ^ 2 + t94 ^ 2) + m(6) * (t156 * pkin(10) ^ 2 + pkin(4) ^ 2) + (t185 * t94 + t97 + t98) * t139 + (t185 * t90 + t100 + t101) * t136; mrSges(6,1) * t3 + mrSges(7,1) * t1 - mrSges(6,2) * t4 - mrSges(7,2) * t2 + (m(7) * t1 + t24) * pkin(5) + t183; mrSges(6,1) * t53 + mrSges(7,1) * t46 - mrSges(6,2) * t54 - mrSges(7,2) * t47 + t112 + t113 + t170 * t140 + (-Ifges(6,6) - Ifges(7,6)) * t159 + (m(7) * t46 + t87) * pkin(5); (-t167 + (-mrSges(6,1) - t178) * t136) * t137 - t81; mrSges(7,1) * t90 + mrSges(7,2) * t94 - t145 * pkin(10) + (m(7) * t90 - t136 * mrSges(7,3)) * pkin(5) + t96 + t95; (0.2e1 * mrSges(7,1) + t178) * pkin(5) - t170; m(7) * t5 + t20; t181; -m(7) * t140; m(7) * t118 + t91; 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t12(1) t12(2) t12(4) t12(7) t12(11) t12(16); t12(2) t12(3) t12(5) t12(8) t12(12) t12(17); t12(4) t12(5) t12(6) t12(9) t12(13) t12(18); t12(7) t12(8) t12(9) t12(10) t12(14) t12(19); t12(11) t12(12) t12(13) t12(14) t12(15) t12(20); t12(16) t12(17) t12(18) t12(19) t12(20) t12(21);];
Mq  = res;
