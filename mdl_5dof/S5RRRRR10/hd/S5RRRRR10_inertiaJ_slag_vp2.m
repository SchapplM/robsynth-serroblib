% Calculate joint inertia matrix for
% S5RRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR10_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR10_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR10_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR10_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR10_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR10_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:32:01
% EndTime: 2019-12-31 22:32:05
% DurationCPUTime: 1.32s
% Computational Cost: add. (2086->285), mult. (4684->425), div. (0->0), fcn. (5043->10), ass. (0->117)
t122 = sin(qJ(5));
t126 = cos(qJ(5));
t128 = cos(qJ(3));
t110 = -pkin(3) * t128 - pkin(2);
t123 = sin(qJ(4));
t124 = sin(qJ(3));
t127 = cos(qJ(4));
t90 = t123 * t124 - t127 * t128;
t91 = t123 * t128 + t124 * t127;
t56 = pkin(4) * t90 - pkin(10) * t91 + t110;
t170 = -pkin(9) - pkin(8);
t100 = t170 * t128;
t142 = t170 * t124;
t69 = -t127 * t100 + t123 * t142;
t26 = -t122 * t69 + t126 * t56;
t27 = t122 * t56 + t126 * t69;
t137 = -t122 * t26 + t126 * t27;
t121 = cos(pkin(5));
t120 = sin(pkin(5));
t125 = sin(qJ(2));
t148 = t120 * t125;
t81 = t121 * t128 - t124 * t148;
t82 = t121 * t124 + t128 * t148;
t50 = t123 * t82 - t127 * t81;
t51 = t123 * t81 + t127 * t82;
t103 = pkin(7) * t148;
t129 = cos(qJ(2));
t166 = pkin(1) * t129;
t74 = t103 + (-pkin(2) - t166) * t121;
t52 = -pkin(3) * t81 + t74;
t15 = pkin(4) * t50 - pkin(10) * t51 + t52;
t147 = t120 * t129;
t84 = t121 * t125 * pkin(1) + pkin(7) * t147;
t75 = pkin(8) * t121 + t84;
t76 = (-pkin(2) * t129 - pkin(8) * t125 - pkin(1)) * t120;
t43 = -t124 * t75 + t128 * t76;
t24 = -pkin(3) * t147 - pkin(9) * t82 + t43;
t44 = t124 * t76 + t128 * t75;
t28 = pkin(9) * t81 + t44;
t14 = t123 * t24 + t127 * t28;
t8 = -pkin(10) * t147 + t14;
t2 = -t122 * t8 + t126 * t15;
t3 = t122 * t15 + t126 * t8;
t139 = -t122 * t2 + t126 * t3;
t176 = -Ifges(4,5) * t82 - Ifges(4,6) * t81;
t67 = -t100 * t123 - t127 * t142;
t175 = t67 ^ 2;
t174 = 0.2e1 * t67;
t35 = -t122 * t51 - t126 * t147;
t173 = t35 / 0.2e1;
t95 = Ifges(6,5) * t122 + Ifges(6,6) * t126;
t172 = t95 / 0.2e1;
t155 = Ifges(6,4) * t126;
t98 = Ifges(6,1) * t122 + t155;
t171 = t98 / 0.2e1;
t169 = t122 / 0.2e1;
t168 = t126 / 0.2e1;
t165 = pkin(10) * t122;
t164 = pkin(10) * t126;
t83 = t121 * t166 - t103;
t161 = t83 * mrSges(3,1);
t160 = t84 * mrSges(3,2);
t159 = Ifges(4,3) + Ifges(5,3);
t158 = -Ifges(5,5) * t51 + Ifges(5,6) * t50;
t157 = Ifges(5,5) * t91 - Ifges(5,6) * t90;
t156 = Ifges(6,4) * t122;
t153 = t122 * t91;
t151 = t126 * t91;
t108 = pkin(3) * t123 + pkin(10);
t150 = t108 * t122;
t149 = t108 * t126;
t146 = Ifges(4,5) * t124 + Ifges(4,6) * t128;
t145 = t122 ^ 2 + t126 ^ 2;
t144 = t124 ^ 2 + t128 ^ 2;
t96 = Ifges(6,2) * t126 + t156;
t143 = t122 * t98 + t126 * t96 + Ifges(5,3);
t36 = -t122 * t147 + t126 * t51;
t9 = Ifges(6,5) * t36 + Ifges(6,6) * t35 + Ifges(6,3) * t50;
t141 = Ifges(3,5) * t148 + Ifges(3,6) * t147 + Ifges(3,3) * t121;
t140 = t145 * t108;
t138 = mrSges(6,1) * t122 + mrSges(6,2) * t126;
t13 = -t123 * t28 + t127 * t24;
t136 = 0.2e1 * t145 * mrSges(6,3);
t38 = Ifges(6,5) * t151 - Ifges(6,6) * t153 + Ifges(6,3) * t90;
t135 = (mrSges(5,1) * t127 - mrSges(5,2) * t123) * pkin(3);
t10 = Ifges(6,4) * t36 + Ifges(6,2) * t35 + Ifges(6,6) * t50;
t11 = Ifges(6,1) * t36 + Ifges(6,4) * t35 + Ifges(6,5) * t50;
t7 = pkin(4) * t147 - t13;
t93 = -mrSges(6,1) * t126 + mrSges(6,2) * t122;
t134 = t13 * mrSges(5,1) - t14 * mrSges(5,2) + t139 * mrSges(6,3) + t10 * t168 + t11 * t169 + t36 * t171 + t50 * t172 + t96 * t173 + t7 * t93 - t158;
t39 = Ifges(6,6) * t90 + (-Ifges(6,2) * t122 + t155) * t91;
t40 = Ifges(6,5) * t90 + (Ifges(6,1) * t126 - t156) * t91;
t133 = -t69 * mrSges(5,2) - t96 * t153 / 0.2e1 + t39 * t168 + t40 * t169 + t157 + t151 * t171 + t90 * t172 + (-mrSges(5,1) + t93) * t67 + t137 * mrSges(6,3);
t109 = -pkin(3) * t127 - pkin(4);
t99 = Ifges(4,1) * t124 + Ifges(4,4) * t128;
t97 = Ifges(4,4) * t124 + Ifges(4,2) * t128;
t94 = -mrSges(4,1) * t128 + mrSges(4,2) * t124;
t66 = -mrSges(4,1) * t147 - mrSges(4,3) * t82;
t65 = mrSges(4,2) * t147 + mrSges(4,3) * t81;
t61 = Ifges(5,1) * t91 - Ifges(5,4) * t90;
t60 = Ifges(5,4) * t91 - Ifges(5,2) * t90;
t59 = mrSges(5,1) * t90 + mrSges(5,2) * t91;
t58 = mrSges(6,1) * t90 - mrSges(6,3) * t151;
t57 = -mrSges(6,2) * t90 - mrSges(6,3) * t153;
t55 = t138 * t91;
t53 = -mrSges(4,1) * t81 + mrSges(4,2) * t82;
t46 = Ifges(4,1) * t82 + Ifges(4,4) * t81 - Ifges(4,5) * t147;
t45 = Ifges(4,4) * t82 + Ifges(4,2) * t81 - Ifges(4,6) * t147;
t42 = -mrSges(5,1) * t147 - mrSges(5,3) * t51;
t41 = mrSges(5,2) * t147 - mrSges(5,3) * t50;
t21 = mrSges(5,1) * t50 + mrSges(5,2) * t51;
t20 = Ifges(5,1) * t51 - Ifges(5,4) * t50 - Ifges(5,5) * t147;
t19 = Ifges(5,4) * t51 - Ifges(5,2) * t50 - Ifges(5,6) * t147;
t18 = mrSges(6,1) * t50 - mrSges(6,3) * t36;
t17 = -mrSges(6,2) * t50 + mrSges(6,3) * t35;
t16 = -mrSges(6,1) * t35 + mrSges(6,2) * t36;
t1 = [t35 * t10 + t36 * t11 + 0.2e1 * t13 * t42 + 0.2e1 * t14 * t41 + 0.2e1 * t7 * t16 + 0.2e1 * t3 * t17 + 0.2e1 * t2 * t18 + t51 * t20 + 0.2e1 * t52 * t21 + 0.2e1 * t43 * t66 + 0.2e1 * t44 * t65 + t81 * t45 + t82 * t46 + 0.2e1 * t74 * t53 + Ifges(2,3) + (t9 - t19) * t50 + (t141 - 0.2e1 * t160 + 0.2e1 * t161) * t121 + ((-0.2e1 * t83 * mrSges(3,3) + Ifges(3,5) * t121 + (-0.2e1 * pkin(1) * mrSges(3,2) + Ifges(3,1) * t125) * t120) * t125 + (0.2e1 * t84 * mrSges(3,3) + Ifges(3,6) * t121 + (0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * Ifges(3,4) * t125 + (Ifges(3,2) + t159) * t129) * t120 + t158 + t176) * t129) * t120 + m(3) * (pkin(1) ^ 2 * t120 ^ 2 + t83 ^ 2 + t84 ^ 2) + m(4) * (t43 ^ 2 + t44 ^ 2 + t74 ^ 2) + m(5) * (t13 ^ 2 + t14 ^ 2 + t52 ^ 2) + m(6) * (t2 ^ 2 + t3 ^ 2 + t7 ^ 2); (t16 - t42) * t67 + (t44 * mrSges(4,3) + pkin(8) * t65 + t45 / 0.2e1) * t128 + (t46 / 0.2e1 - pkin(8) * t66 - t43 * mrSges(4,3)) * t124 + m(5) * (t110 * t52 - t13 * t67 + t14 * t69) + m(6) * (t2 * t26 + t27 * t3 + t67 * t7) + t141 + (t9 / 0.2e1 - t19 / 0.2e1 - t14 * mrSges(5,3)) * t90 + t74 * t94 + t81 * t97 / 0.2e1 + t82 * t99 / 0.2e1 + t110 * t21 + t7 * t55 + t3 * t57 + t2 * t58 + t52 * t59 + t51 * t61 / 0.2e1 + t69 * t41 - pkin(2) * t53 - (t146 + t157) * t147 / 0.2e1 + t36 * t40 / 0.2e1 + t26 * t18 + t27 * t17 + (-t60 / 0.2e1 + t38 / 0.2e1) * t50 + t39 * t173 + m(4) * (-pkin(2) * t74 + (-t43 * t124 + t44 * t128) * pkin(8)) + (t11 * t168 + t20 / 0.2e1 - t122 * t10 / 0.2e1 - t13 * mrSges(5,3)) * t91 - t160 + t161; -0.2e1 * pkin(2) * t94 + 0.2e1 * t110 * t59 + t124 * t99 + t128 * t97 + 0.2e1 * t26 * t58 + 0.2e1 * t27 * t57 + t55 * t174 + Ifges(3,3) + 0.2e1 * t144 * pkin(8) * mrSges(4,3) + (-0.2e1 * mrSges(5,3) * t69 + t38 - t60) * t90 + m(6) * (t26 ^ 2 + t27 ^ 2 + t175) + m(5) * (t110 ^ 2 + t69 ^ 2 + t175) + m(4) * (pkin(8) ^ 2 * t144 + pkin(2) ^ 2) + (mrSges(5,3) * t174 - t122 * t39 + t126 * t40 + t61) * t91; -t159 * t147 + t17 * t149 + t109 * t16 - t18 * t150 + t43 * mrSges(4,1) - t44 * mrSges(4,2) + t134 + m(6) * (t108 * t139 + t109 * t7) + (m(5) * (t123 * t14 + t127 * t13) + t127 * t42 + t123 * t41) * pkin(3) - t176; (-mrSges(4,1) * t124 - mrSges(4,2) * t128) * pkin(8) + t57 * t149 + t109 * t55 - t58 * t150 + t133 + m(6) * (t108 * t137 + t109 * t67) + (m(5) * (t123 * t69 - t127 * t67) + (-t123 * t90 - t127 * t91) * mrSges(5,3)) * pkin(3) + t146; 0.2e1 * t109 * t93 + Ifges(4,3) + 0.2e1 * t135 + t108 * t136 + m(6) * (t108 ^ 2 * t145 + t109 ^ 2) + m(5) * (t123 ^ 2 + t127 ^ 2) * pkin(3) ^ 2 + t143; -Ifges(5,3) * t147 + t17 * t164 - t18 * t165 - pkin(4) * t16 + t134 + m(6) * (-pkin(4) * t7 + pkin(10) * t139); t57 * t164 - pkin(4) * t55 - t58 * t165 + t133 + m(6) * (-pkin(4) * t67 + pkin(10) * t137); m(6) * (-pkin(4) * t109 + pkin(10) * t140) + (-pkin(4) + t109) * t93 + t135 + (pkin(10) * t145 + t140) * mrSges(6,3) + t143; -0.2e1 * pkin(4) * t93 + m(6) * (pkin(10) ^ 2 * t145 + pkin(4) ^ 2) + pkin(10) * t136 + t143; mrSges(6,1) * t2 - mrSges(6,2) * t3 + t9; mrSges(6,1) * t26 - mrSges(6,2) * t27 + t38; -t138 * t108 + t95; -t138 * pkin(10) + t95; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
