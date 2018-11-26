% Calculate joint inertia matrix for
% S6RRRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
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
% Datum: 2018-11-23 17:37
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRRPPR6_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR6_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR6_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR6_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR6_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR6_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:36:41
% EndTime: 2018-11-23 17:36:43
% DurationCPUTime: 1.75s
% Computational Cost: add. (2774->355), mult. (6105->496), div. (0->0), fcn. (6716->10), ass. (0->133)
t181 = Ifges(6,4) - Ifges(5,5);
t180 = Ifges(6,5) - Ifges(5,6);
t131 = sin(pkin(6));
t139 = cos(qJ(2));
t156 = t131 * t139;
t135 = sin(qJ(3));
t138 = cos(qJ(3));
t130 = sin(pkin(11));
t132 = cos(pkin(11));
t97 = t130 * t135 - t132 * t138;
t98 = t130 * t138 + t132 * t135;
t179 = Ifges(4,5) * t135 + Ifges(4,6) * t138 + t180 * t97 - t181 * t98;
t117 = pkin(3) * t130 + qJ(5);
t178 = t117 ^ 2;
t134 = sin(qJ(6));
t137 = cos(qJ(6));
t133 = cos(pkin(6));
t136 = sin(qJ(2));
t157 = t131 * t136;
t86 = t133 * t138 - t135 * t157;
t87 = t133 * t135 + t138 * t157;
t55 = t130 * t87 - t132 * t86;
t32 = t134 * t156 + t137 * t55;
t177 = t32 / 0.2e1;
t160 = Ifges(7,4) * t137;
t36 = Ifges(7,5) * t98 + (Ifges(7,1) * t134 + t160) * t97;
t176 = t36 / 0.2e1;
t175 = pkin(4) + pkin(10);
t104 = Ifges(7,5) * t137 - Ifges(7,6) * t134;
t174 = t104 / 0.2e1;
t173 = t134 / 0.2e1;
t172 = t137 / 0.2e1;
t170 = pkin(1) * t139;
t112 = pkin(8) * t157;
t88 = t133 * t170 - t112;
t169 = t88 * mrSges(3,1);
t89 = t133 * t136 * pkin(1) + pkin(8) * t156;
t168 = t89 * mrSges(3,2);
t153 = t134 ^ 2 + t137 ^ 2;
t100 = m(7) * t153;
t167 = m(6) + t100;
t165 = mrSges(6,2) - mrSges(5,1);
t164 = -qJ(4) - pkin(9);
t79 = pkin(9) * t133 + t89;
t80 = (-pkin(2) * t139 - pkin(9) * t136 - pkin(1)) * t131;
t42 = -t135 * t79 + t138 * t80;
t26 = -pkin(3) * t156 - qJ(4) * t87 + t42;
t43 = t135 * t80 + t138 * t79;
t29 = qJ(4) * t86 + t43;
t12 = t130 * t26 + t132 * t29;
t161 = Ifges(7,4) * t134;
t159 = t134 * t97;
t158 = t137 * t97;
t154 = t135 ^ 2 + t138 ^ 2;
t152 = Ifges(6,1) + Ifges(4,3) + Ifges(5,3);
t33 = t134 * t55 - t137 * t156;
t56 = t130 * t86 + t132 * t87;
t7 = Ifges(7,5) * t33 + Ifges(7,6) * t32 + Ifges(7,3) * t56;
t34 = Ifges(7,5) * t159 + Ifges(7,6) * t158 + Ifges(7,3) * t98;
t103 = t164 * t138;
t149 = t164 * t135;
t69 = -t103 * t130 - t132 * t149;
t71 = -t132 * t103 + t130 * t149;
t151 = t69 ^ 2 + t71 ^ 2;
t150 = Ifges(3,5) * t157 + Ifges(3,6) * t156 + Ifges(3,3) * t133;
t121 = -pkin(3) * t138 - pkin(2);
t120 = -pkin(3) * t132 - pkin(4);
t11 = -t130 * t29 + t132 * t26;
t148 = t153 * mrSges(7,3);
t10 = pkin(4) * t156 - t11;
t3 = pkin(5) * t56 + pkin(10) * t156 + t10;
t78 = t112 + (-pkin(2) - t170) * t133;
t57 = -pkin(3) * t86 + t78;
t142 = -qJ(5) * t56 + t57;
t5 = t175 * t55 + t142;
t1 = -t134 * t5 + t137 * t3;
t2 = t134 * t3 + t137 * t5;
t147 = t1 * t137 + t134 * t2;
t146 = mrSges(7,1) * t137 - t134 * mrSges(7,2);
t143 = -qJ(5) * t98 + t121;
t37 = t175 * t97 + t143;
t46 = pkin(5) * t98 + t69;
t15 = -t134 * t37 + t137 * t46;
t16 = t134 * t46 + t137 * t37;
t145 = t134 * t16 + t137 * t15;
t144 = -Ifges(4,5) * t87 - Ifges(4,6) * t86 - t180 * t55 + t181 * t56;
t6 = qJ(5) * t156 - t12;
t39 = t56 * mrSges(6,1) - mrSges(6,2) * t156;
t116 = -pkin(10) + t120;
t108 = Ifges(4,1) * t135 + Ifges(4,4) * t138;
t107 = Ifges(7,1) * t137 - t161;
t106 = Ifges(4,4) * t135 + Ifges(4,2) * t138;
t105 = -Ifges(7,2) * t134 + t160;
t102 = mrSges(7,1) * t134 + mrSges(7,2) * t137;
t101 = -mrSges(4,1) * t138 + mrSges(4,2) * t135;
t91 = t98 * mrSges(6,3);
t90 = t98 * mrSges(5,2);
t74 = -mrSges(4,1) * t156 - mrSges(4,3) * t87;
t73 = mrSges(4,2) * t156 + mrSges(4,3) * t86;
t68 = Ifges(5,1) * t98 - Ifges(5,4) * t97;
t67 = Ifges(5,4) * t98 - Ifges(5,2) * t97;
t66 = -Ifges(6,2) * t98 + Ifges(6,6) * t97;
t65 = -Ifges(6,6) * t98 + Ifges(6,3) * t97;
t64 = -t97 * mrSges(6,2) - t91;
t63 = t97 * mrSges(5,1) + t90;
t62 = -mrSges(7,2) * t98 + mrSges(7,3) * t158;
t61 = mrSges(7,1) * t98 - mrSges(7,3) * t159;
t60 = pkin(4) * t97 + t143;
t59 = -mrSges(4,1) * t86 + mrSges(4,2) * t87;
t58 = t146 * t97;
t49 = t56 * mrSges(6,3);
t48 = t56 * mrSges(5,2);
t47 = -t97 * pkin(5) + t71;
t45 = Ifges(4,1) * t87 + Ifges(4,4) * t86 - Ifges(4,5) * t156;
t44 = Ifges(4,4) * t87 + Ifges(4,2) * t86 - Ifges(4,6) * t156;
t41 = -mrSges(5,1) * t156 - mrSges(5,3) * t56;
t40 = mrSges(5,2) * t156 - mrSges(5,3) * t55;
t38 = mrSges(6,1) * t55 + mrSges(6,3) * t156;
t35 = Ifges(7,6) * t98 + (Ifges(7,2) * t137 + t161) * t97;
t24 = -t55 * mrSges(6,2) - t49;
t23 = t55 * mrSges(5,1) + t48;
t22 = Ifges(5,1) * t56 - Ifges(5,4) * t55 - Ifges(5,5) * t156;
t21 = Ifges(5,4) * t56 - Ifges(5,2) * t55 - Ifges(5,6) * t156;
t20 = -Ifges(6,4) * t156 - Ifges(6,2) * t56 + Ifges(6,6) * t55;
t19 = -Ifges(6,5) * t156 - Ifges(6,6) * t56 + Ifges(6,3) * t55;
t18 = mrSges(7,1) * t56 - mrSges(7,3) * t33;
t17 = -mrSges(7,2) * t56 + mrSges(7,3) * t32;
t14 = -mrSges(7,1) * t32 + mrSges(7,2) * t33;
t13 = pkin(4) * t55 + t142;
t9 = Ifges(7,1) * t33 + Ifges(7,4) * t32 + Ifges(7,5) * t56;
t8 = Ifges(7,4) * t33 + Ifges(7,2) * t32 + Ifges(7,6) * t56;
t4 = -pkin(5) * t55 - t6;
t25 = [t86 * t44 + t87 * t45 + m(3) * (pkin(1) ^ 2 * t131 ^ 2 + t88 ^ 2 + t89 ^ 2) + m(4) * (t42 ^ 2 + t43 ^ 2 + t78 ^ 2) + m(5) * (t11 ^ 2 + t12 ^ 2 + t57 ^ 2) + m(6) * (t10 ^ 2 + t13 ^ 2 + t6 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t4 ^ 2) + 0.2e1 * t43 * t73 + 0.2e1 * t42 * t74 + 0.2e1 * t78 * t59 + 0.2e1 * t57 * t23 + t32 * t8 + t33 * t9 + 0.2e1 * t6 * t38 + 0.2e1 * t10 * t39 + 0.2e1 * t12 * t40 + 0.2e1 * t11 * t41 + 0.2e1 * t2 * t17 + 0.2e1 * t1 * t18 + 0.2e1 * t13 * t24 + 0.2e1 * t4 * t14 + (t22 + t7 - t20) * t56 + (t19 - t21) * t55 + (t150 - 0.2e1 * t168 + 0.2e1 * t169) * t133 + ((-0.2e1 * t88 * mrSges(3,3) + Ifges(3,5) * t133 + (-0.2e1 * mrSges(3,2) * pkin(1) + Ifges(3,1) * t136) * t131) * t136 + (0.2e1 * t89 * mrSges(3,3) + Ifges(3,6) * t133 + (0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * Ifges(3,4) * t136 + (Ifges(3,2) + t152) * t139) * t131 + t144) * t139) * t131 + Ifges(2,3); -t179 * t156 / 0.2e1 + (t40 - t38) * t71 + (t39 - t41) * t69 + t86 * t106 / 0.2e1 + t87 * t108 / 0.2e1 + t121 * t23 + t78 * t101 + (pkin(9) * t73 + t43 * mrSges(4,3) + t44 / 0.2e1) * t138 + (-pkin(9) * t74 - t42 * mrSges(4,3) + t45 / 0.2e1) * t135 + t150 - t4 * t58 - pkin(2) * t59 + t60 * t24 + t1 * t61 + t2 * t62 + t57 * t63 + t13 * t64 + t47 * t14 + t16 * t17 + t15 * t18 + m(5) * (-t11 * t69 + t12 * t71 + t121 * t57) + m(6) * (t10 * t69 + t13 * t60 - t6 * t71) + m(7) * (t1 * t15 + t16 * t2 + t4 * t47) + (t22 / 0.2e1 + t7 / 0.2e1 - t20 / 0.2e1 - t11 * mrSges(5,3) + t10 * mrSges(6,1)) * t98 + (-t66 / 0.2e1 + t68 / 0.2e1 + t34 / 0.2e1) * t56 + (t65 / 0.2e1 - t67 / 0.2e1) * t55 + t33 * t176 + t35 * t177 - t168 + t169 + (t8 * t172 + t9 * t173 - t21 / 0.2e1 + t19 / 0.2e1 + t6 * mrSges(6,1) - t12 * mrSges(5,3)) * t97 + m(4) * (-pkin(2) * t78 + (-t135 * t42 + t138 * t43) * pkin(9)); -0.2e1 * pkin(2) * t101 + t138 * t106 + t135 * t108 + 0.2e1 * t121 * t63 + 0.2e1 * t15 * t61 + 0.2e1 * t16 * t62 - 0.2e1 * t47 * t58 + 0.2e1 * t60 * t64 + Ifges(3,3) + (t34 - t66 + t68) * t98 + m(4) * (t154 * pkin(9) ^ 2 + pkin(2) ^ 2) + m(5) * (t121 ^ 2 + t151) + m(6) * (t60 ^ 2 + t151) + m(7) * (t15 ^ 2 + t16 ^ 2 + t47 ^ 2) + (t134 * t36 + t137 * t35 + t65 - t67) * t97 + 0.2e1 * t154 * pkin(9) * mrSges(4,3) + 0.2e1 * (t69 * t98 - t71 * t97) * (mrSges(6,1) + mrSges(5,3)); t56 * t174 + t105 * t177 + t33 * t107 / 0.2e1 + t120 * t39 + t4 * t102 + t42 * mrSges(4,1) - t43 * mrSges(4,2) - t12 * mrSges(5,2) - t6 * mrSges(6,3) + t10 * mrSges(6,2) + t11 * mrSges(5,1) - t144 + (t116 * t18 - t1 * mrSges(7,3) + t9 / 0.2e1) * t137 + (t116 * t17 - t2 * mrSges(7,3) - t8 / 0.2e1) * t134 + (t14 - t38) * t117 + m(6) * (t10 * t120 - t117 * t6) + m(7) * (t147 * t116 + t117 * t4) + (m(5) * (t11 * t132 + t12 * t130) + t132 * t41 + t130 * t40) * pkin(3) - t152 * t156; t47 * t102 - t117 * t58 + (t120 * mrSges(6,1) + t174) * t98 + (mrSges(6,3) - mrSges(5,2)) * t71 + t165 * t69 + (-mrSges(4,1) * t135 - mrSges(4,2) * t138) * pkin(9) + (-t15 * mrSges(7,3) + t116 * t61 + t176) * t137 + (t116 * t62 - t16 * mrSges(7,3) - t35 / 0.2e1) * t134 + m(7) * (t145 * t116 + t117 * t47) + m(6) * (t117 * t71 + t120 * t69) + (-t117 * mrSges(6,1) + t105 * t172 + t107 * t173) * t97 + (m(5) * (t130 * t71 - t132 * t69) + (-t130 * t97 - t132 * t98) * mrSges(5,3)) * pkin(3) + t179; 0.2e1 * t120 * mrSges(6,2) - t134 * t105 + t137 * t107 + m(7) * (t153 * t116 ^ 2 + t178) + m(6) * (t120 ^ 2 + t178) + m(5) * (t130 ^ 2 + t132 ^ 2) * pkin(3) ^ 2 + t152 + 0.2e1 * (t102 + mrSges(6,3)) * t117 + 0.2e1 * (mrSges(5,1) * t132 - mrSges(5,2) * t130) * pkin(3) - 0.2e1 * t116 * t148; -t134 * t18 + t137 * t17 + t48 - t49 - t165 * t55 + m(7) * (-t1 * t134 + t137 * t2) + m(6) * t13 + m(5) * t57; -t134 * t61 + t137 * t62 + t90 - t91 - t165 * t97 + m(7) * (-t134 * t15 + t137 * t16) + m(6) * t60 + m(5) * t121; 0; m(5) + t167; m(6) * t10 + m(7) * t147 + t134 * t17 + t137 * t18 + t39; m(6) * t69 + m(7) * t145 + t98 * mrSges(6,1) + t134 * t62 + t137 * t61; m(6) * t120 + t116 * t100 + mrSges(6,2) - t148; 0; t167; mrSges(7,1) * t1 - mrSges(7,2) * t2 + t7; mrSges(7,1) * t15 - mrSges(7,2) * t16 + t34; t146 * t116 + t104; -t102; t146; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t25(1) t25(2) t25(4) t25(7) t25(11) t25(16); t25(2) t25(3) t25(5) t25(8) t25(12) t25(17); t25(4) t25(5) t25(6) t25(9) t25(13) t25(18); t25(7) t25(8) t25(9) t25(10) t25(14) t25(19); t25(11) t25(12) t25(13) t25(14) t25(15) t25(20); t25(16) t25(17) t25(18) t25(19) t25(20) t25(21);];
Mq  = res;
