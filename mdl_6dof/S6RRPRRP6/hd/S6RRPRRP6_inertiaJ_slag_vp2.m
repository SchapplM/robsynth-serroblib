% Calculate joint inertia matrix for
% S6RRPRRP6
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 17:14
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRP6_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP6_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP6_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP6_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP6_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP6_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:14:20
% EndTime: 2018-11-23 17:14:22
% DurationCPUTime: 1.71s
% Computational Cost: add. (2845->370), mult. (6794->509), div. (0->0), fcn. (7258->10), ass. (0->144)
t194 = Ifges(3,3) + Ifges(4,3);
t136 = sin(qJ(5));
t139 = cos(qJ(5));
t164 = t136 ^ 2 + t139 ^ 2;
t132 = sin(pkin(11));
t116 = pkin(2) * t132 + pkin(9);
t193 = 0.2e1 * t116;
t135 = cos(pkin(6));
t137 = sin(qJ(4));
t140 = cos(qJ(4));
t133 = sin(pkin(6));
t134 = cos(pkin(11));
t138 = sin(qJ(2));
t141 = cos(qJ(2));
t69 = (t132 * t141 + t134 * t138) * t133;
t55 = t135 * t137 + t140 * t69;
t171 = t133 * t141;
t172 = t133 * t138;
t68 = t132 * t172 - t134 * t171;
t40 = t136 * t55 - t139 * t68;
t41 = t136 * t68 + t139 * t55;
t54 = -t135 * t140 + t137 * t69;
t7 = Ifges(6,5) * t41 - Ifges(6,6) * t40 + Ifges(6,3) * t54;
t8 = Ifges(7,4) * t41 + Ifges(7,2) * t54 + Ifges(7,6) * t40;
t192 = t7 + t8;
t191 = Ifges(4,5) * t69 - Ifges(4,6) * t68;
t190 = (mrSges(7,2) + mrSges(6,3)) * t164;
t83 = (-pkin(2) * t141 - pkin(1)) * t133;
t189 = 0.2e1 * t83;
t188 = pkin(1) * t135;
t187 = pkin(4) * t140;
t186 = pkin(10) * t137;
t107 = t141 * t188;
t78 = -pkin(8) * t172 + t107;
t185 = t78 * mrSges(3,1);
t79 = pkin(8) * t171 + t138 * t188;
t184 = t79 * mrSges(3,2);
t91 = -mrSges(6,1) * t139 + mrSges(6,2) * t136;
t183 = mrSges(5,1) - t91;
t182 = -Ifges(6,3) - Ifges(7,2);
t56 = pkin(2) * t135 + t107 + (-pkin(8) - qJ(3)) * t172;
t62 = qJ(3) * t171 + t79;
t38 = t132 * t56 + t134 * t62;
t31 = pkin(9) * t135 + t38;
t42 = pkin(3) * t68 - pkin(9) * t69 + t83;
t18 = t137 * t42 + t140 * t31;
t14 = pkin(10) * t68 + t18;
t37 = -t132 * t62 + t134 * t56;
t30 = -pkin(3) * t135 - t37;
t16 = pkin(4) * t54 - pkin(10) * t55 + t30;
t4 = t136 * t16 + t139 * t14;
t21 = -mrSges(6,2) * t54 - mrSges(6,3) * t40;
t22 = -mrSges(7,2) * t40 + mrSges(7,3) * t54;
t181 = t21 + t22;
t23 = mrSges(6,1) * t54 - mrSges(6,3) * t41;
t24 = -t54 * mrSges(7,1) + t41 * mrSges(7,2);
t180 = -t23 + t24;
t20 = mrSges(6,1) * t40 + mrSges(6,2) * t41;
t44 = mrSges(5,1) * t68 - mrSges(5,3) * t55;
t179 = -t44 + t20;
t173 = t116 * t140;
t117 = -pkin(2) * t134 - pkin(3);
t82 = t117 - t186 - t187;
t53 = t136 * t82 + t139 * t173;
t178 = Ifges(6,4) * t136;
t177 = Ifges(6,4) * t139;
t176 = Ifges(7,5) * t136;
t175 = Ifges(7,5) * t139;
t174 = t139 * t82;
t170 = t136 * t137;
t169 = t137 * t139;
t86 = t140 * mrSges(7,1) + mrSges(7,2) * t169;
t168 = Ifges(7,4) * t169 + Ifges(7,6) * t170;
t167 = t164 * t186;
t94 = Ifges(6,5) * t136 + Ifges(6,6) * t139;
t166 = Ifges(5,5) * t137 + Ifges(5,6) * t140;
t165 = t164 * pkin(10) ^ 2;
t129 = t137 ^ 2;
t131 = t140 ^ 2;
t163 = t129 + t131;
t6 = Ifges(7,5) * t41 + Ifges(7,6) * t54 + Ifges(7,3) * t40;
t9 = Ifges(6,4) * t41 - Ifges(6,2) * t40 + Ifges(6,6) * t54;
t162 = t6 / 0.2e1 - t9 / 0.2e1;
t161 = Ifges(5,5) * t55 - Ifges(5,6) * t54 + Ifges(5,3) * t68;
t10 = Ifges(7,1) * t41 + Ifges(7,4) * t54 + Ifges(7,5) * t40;
t11 = Ifges(6,1) * t41 - Ifges(6,4) * t40 + Ifges(6,5) * t54;
t160 = t10 / 0.2e1 + t11 / 0.2e1;
t70 = -Ifges(7,6) * t140 + (Ifges(7,3) * t136 + t175) * t137;
t73 = -Ifges(6,6) * t140 + (-Ifges(6,2) * t136 + t177) * t137;
t159 = t70 / 0.2e1 - t73 / 0.2e1;
t74 = -Ifges(7,4) * t140 + (Ifges(7,1) * t139 + t176) * t137;
t75 = -Ifges(6,5) * t140 + (Ifges(6,1) * t139 - t178) * t137;
t158 = t74 / 0.2e1 + t75 / 0.2e1;
t93 = -Ifges(7,3) * t139 + t176;
t96 = Ifges(6,2) * t139 + t178;
t157 = t93 / 0.2e1 - t96 / 0.2e1;
t95 = Ifges(7,4) * t136 - Ifges(7,6) * t139;
t156 = t94 / 0.2e1 + t95 / 0.2e1;
t98 = Ifges(7,1) * t136 - t175;
t99 = Ifges(6,1) * t136 + t177;
t155 = t98 / 0.2e1 + t99 / 0.2e1;
t17 = -t137 * t31 + t140 * t42;
t153 = Ifges(6,5) * t169 - Ifges(6,6) * t170;
t1 = qJ(6) * t54 + t4;
t3 = -t136 * t14 + t139 * t16;
t2 = -pkin(5) * t54 - t3;
t152 = t1 * t139 + t136 * t2;
t151 = -t136 * t3 + t139 * t4;
t150 = t136 * mrSges(6,1) + t139 * mrSges(6,2);
t149 = t136 * mrSges(7,1) - t139 * mrSges(7,3);
t148 = -pkin(5) * t136 + qJ(6) * t139;
t52 = -t136 * t173 + t174;
t147 = -t136 * t52 + t139 * t53;
t13 = -pkin(4) * t68 - t17;
t146 = Ifges(3,5) * t172 + Ifges(3,6) * t171 + t194 * t135 + t191;
t45 = -qJ(6) * t140 + t53;
t46 = -t174 + (t116 * t136 + pkin(5)) * t140;
t84 = mrSges(6,2) * t140 - mrSges(6,3) * t170;
t85 = -mrSges(6,1) * t140 - mrSges(6,3) * t169;
t87 = -mrSges(7,2) * t170 - mrSges(7,3) * t140;
t145 = m(7) * (t136 * t46 + t139 * t45) + (t84 + t87) * t139 + (-t85 + t86) * t136;
t144 = m(7) * t148 - t149 - t150;
t114 = t116 ^ 2;
t101 = t129 * t114;
t100 = Ifges(5,1) * t137 + Ifges(5,4) * t140;
t97 = Ifges(5,4) * t137 + Ifges(5,2) * t140;
t92 = -t140 * mrSges(5,1) + t137 * mrSges(5,2);
t90 = -mrSges(7,1) * t139 - mrSges(7,3) * t136;
t89 = -pkin(5) * t139 - qJ(6) * t136 - pkin(4);
t81 = t150 * t137;
t80 = t149 * t137;
t72 = -Ifges(7,2) * t140 + t168;
t71 = -Ifges(6,3) * t140 + t153;
t64 = t69 * mrSges(4,2);
t61 = (t116 - t148) * t137;
t58 = mrSges(4,1) * t135 - mrSges(4,3) * t69;
t57 = -mrSges(4,2) * t135 - mrSges(4,3) * t68;
t43 = -mrSges(5,2) * t68 - mrSges(5,3) * t54;
t29 = mrSges(5,1) * t54 + mrSges(5,2) * t55;
t26 = Ifges(5,1) * t55 - Ifges(5,4) * t54 + Ifges(5,5) * t68;
t25 = Ifges(5,4) * t55 - Ifges(5,2) * t54 + Ifges(5,6) * t68;
t19 = mrSges(7,1) * t40 - mrSges(7,3) * t41;
t5 = pkin(5) * t40 - qJ(6) * t41 + t13;
t12 = [(-t25 + t192) * t54 + (t146 - 0.2e1 * t184 + 0.2e1 * t185 + t191) * t135 + m(3) * (t78 ^ 2 + t79 ^ 2) + 0.2e1 * t38 * t57 + 0.2e1 * t37 * t58 + t55 * t26 + 0.2e1 * t18 * t43 + 0.2e1 * t17 * t44 + 0.2e1 * t4 * t21 + 0.2e1 * t1 * t22 + 0.2e1 * t3 * t23 + 0.2e1 * t2 * t24 + 0.2e1 * t30 * t29 + 0.2e1 * t5 * t19 + 0.2e1 * t13 * t20 + m(4) * (t37 ^ 2 + t38 ^ 2 + t83 ^ 2) + m(5) * (t17 ^ 2 + t18 ^ 2 + t30 ^ 2) + m(6) * (t13 ^ 2 + t3 ^ 2 + t4 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) + (t10 + t11) * t41 + (t6 - t9) * t40 + ((Ifges(3,5) * t138 + Ifges(3,6) * t141) * t135 + 0.2e1 * (-t138 * t78 + t141 * t79) * mrSges(3,3) + (m(3) * pkin(1) ^ 2 - 0.2e1 * pkin(1) * (-mrSges(3,1) * t141 + mrSges(3,2) * t138) + t138 * (Ifges(3,1) * t138 + Ifges(3,4) * t141) + t141 * (Ifges(3,4) * t138 + Ifges(3,2) * t141)) * t133) * t133 + Ifges(4,1) * t69 ^ 2 + (mrSges(4,1) * t189 - 0.2e1 * Ifges(4,4) * t69 + Ifges(4,2) * t68 + t161) * t68 + Ifges(2,3) + t64 * t189; t55 * t100 / 0.2e1 + t117 * t29 + t4 * t84 + t3 * t85 + t2 * t86 + t1 * t87 + t30 * t92 + t5 * t80 + t13 * t81 + t61 * t19 + t52 * t23 + t53 * t21 + t45 * t22 + t46 * t24 + t37 * mrSges(4,1) - t38 * mrSges(4,2) + (t116 * t43 + t18 * mrSges(5,3) + t25 / 0.2e1 - t7 / 0.2e1 - t8 / 0.2e1) * t140 + m(6) * (t116 * t13 * t137 + t3 * t52 + t4 * t53) + m(7) * (t1 * t45 + t2 * t46 + t5 * t61) + (-t17 * mrSges(5,3) + t26 / 0.2e1 + t160 * t139 + t162 * t136 + t179 * t116) * t137 + t68 * t166 / 0.2e1 + m(5) * (t117 * t30 + (-t137 * t17 + t140 * t18) * t116) + (-t97 / 0.2e1 + t71 / 0.2e1 + t72 / 0.2e1) * t54 + t158 * t41 + t159 * t40 + t146 + (t132 * t57 + t134 * t58 + m(4) * (t132 * t38 + t134 * t37)) * pkin(2) - t184 + t185; 0.2e1 * t117 * t92 + 0.2e1 * t45 * t87 + 0.2e1 * t46 * t86 + 0.2e1 * t52 * t85 + 0.2e1 * t53 * t84 + 0.2e1 * t61 * t80 + (t97 - t71 - t72) * t140 + m(5) * (t114 * t131 + t117 ^ 2 + t101) + m(6) * (t52 ^ 2 + t53 ^ 2 + t101) + m(7) * (t45 ^ 2 + t46 ^ 2 + t61 ^ 2) + m(4) * (t132 ^ 2 + t134 ^ 2) * pkin(2) ^ 2 + (t81 * t193 + t100 + (t74 + t75) * t139 + (t70 - t73) * t136) * t137 + 0.2e1 * (mrSges(4,1) * t134 - mrSges(4,2) * t132) * pkin(2) + t163 * mrSges(5,3) * t193 + t194; t68 * mrSges(4,1) + t64 + (-t19 - t179) * t140 + (t180 * t136 + t181 * t139 + t43) * t137 + m(6) * (-t13 * t140 + t151 * t137) + m(7) * (t152 * t137 - t140 * t5) + m(5) * (t137 * t18 + t140 * t17) + m(4) * t83; (-m(7) * t61 - t80 - t81) * t140 + (m(6) * (t147 - t173) + t145) * t137; m(4) + m(5) * t163 + 0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * (t164 * t129 + t131); t17 * mrSges(5,1) - t18 * mrSges(5,2) - pkin(4) * t20 + t13 * t91 + t89 * t19 + t5 * t90 + t156 * t54 + t155 * t41 + t157 * t40 + (t1 * mrSges(7,2) + t4 * mrSges(6,3) + t181 * pkin(10) - t162) * t139 + (t2 * mrSges(7,2) - t3 * mrSges(6,3) + t180 * pkin(10) + t160) * t136 + m(7) * (t152 * pkin(10) + t5 * t89) + m(6) * (-pkin(4) * t13 + t151 * pkin(10)) + t161; -pkin(4) * t81 + t89 * t80 + (m(7) * t89 + t90) * t61 + (-t116 * mrSges(5,2) - t156) * t140 + (t45 * mrSges(7,2) + t53 * mrSges(6,3) - t159) * t139 + (t46 * mrSges(7,2) - t52 * mrSges(6,3) + t158) * t136 + (m(6) * t147 + t145) * pkin(10) + (t155 * t139 + t157 * t136 + (-m(6) * pkin(4) - t183) * t116) * t137 + t166; (-t90 + t183) * t140 + m(6) * (t167 + t187) + m(7) * (-t140 * t89 + t167) + (-mrSges(5,2) + t190) * t137; -0.2e1 * pkin(4) * t91 + 0.2e1 * t89 * t90 + Ifges(5,3) + (t96 - t93) * t139 + (t98 + t99) * t136 + m(7) * (t89 ^ 2 + t165) + m(6) * (pkin(4) ^ 2 + t165) + 0.2e1 * pkin(10) * t190; t3 * mrSges(6,1) - pkin(5) * t24 + qJ(6) * t22 + m(7) * (-pkin(5) * t2 + qJ(6) * t1) + t1 * mrSges(7,3) - t2 * mrSges(7,1) - t4 * mrSges(6,2) + t192; -t53 * mrSges(6,2) + t52 * mrSges(6,1) - t46 * mrSges(7,1) - pkin(5) * t86 + t45 * mrSges(7,3) + qJ(6) * t87 + m(7) * (-pkin(5) * t46 + qJ(6) * t45) + t182 * t140 + t153 + t168; t144 * t137; t148 * mrSges(7,2) + t144 * pkin(10) + t94 + t95; 0.2e1 * pkin(5) * mrSges(7,1) + 0.2e1 * qJ(6) * mrSges(7,3) + m(7) * (pkin(5) ^ 2 + qJ(6) ^ 2) - t182; m(7) * t2 + t24; m(7) * t46 + t86; m(7) * t170; (m(7) * pkin(10) + mrSges(7,2)) * t136; -m(7) * pkin(5) - mrSges(7,1); m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t12(1) t12(2) t12(4) t12(7) t12(11) t12(16); t12(2) t12(3) t12(5) t12(8) t12(12) t12(17); t12(4) t12(5) t12(6) t12(9) t12(13) t12(18); t12(7) t12(8) t12(9) t12(10) t12(14) t12(19); t12(11) t12(12) t12(13) t12(14) t12(15) t12(20); t12(16) t12(17) t12(18) t12(19) t12(20) t12(21);];
Mq  = res;
