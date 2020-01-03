% Calculate joint inertia matrix for
% S5RRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-31 21:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR10_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR10_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR10_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR10_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR10_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR10_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:26:41
% EndTime: 2019-12-31 21:26:44
% DurationCPUTime: 1.03s
% Computational Cost: add. (1808->271), mult. (4070->409), div. (0->0), fcn. (4410->10), ass. (0->108)
t105 = sin(pkin(5));
t113 = cos(qJ(2));
t124 = t105 * t113;
t109 = sin(qJ(3));
t112 = cos(qJ(3));
t104 = sin(pkin(10));
t106 = cos(pkin(10));
t74 = t104 * t109 - t106 * t112;
t75 = t104 * t112 + t106 * t109;
t148 = Ifges(4,5) * t109 + Ifges(5,5) * t75 + Ifges(4,6) * t112 - Ifges(5,6) * t74;
t134 = -qJ(4) - pkin(8);
t120 = t134 * t109;
t79 = t134 * t112;
t51 = -t104 * t79 - t106 * t120;
t147 = t51 ^ 2;
t146 = 0.2e1 * t51;
t108 = sin(qJ(5));
t111 = cos(qJ(5));
t107 = cos(pkin(5));
t110 = sin(qJ(2));
t125 = t105 * t110;
t66 = t107 * t112 - t109 * t125;
t67 = t107 * t109 + t112 * t125;
t41 = t104 * t66 + t106 * t67;
t25 = -t108 * t41 - t111 * t124;
t145 = t25 / 0.2e1;
t26 = -t108 * t124 + t111 * t41;
t144 = t26 / 0.2e1;
t80 = Ifges(6,5) * t108 + Ifges(6,6) * t111;
t143 = t80 / 0.2e1;
t142 = -t108 / 0.2e1;
t141 = t108 / 0.2e1;
t140 = t111 / 0.2e1;
t138 = pkin(1) * t113;
t87 = pkin(7) * t125;
t68 = t107 * t138 - t87;
t137 = t68 * mrSges(3,1);
t69 = t107 * t110 * pkin(1) + pkin(7) * t124;
t136 = t69 * mrSges(3,2);
t135 = Ifges(4,3) + Ifges(5,3);
t60 = pkin(8) * t107 + t69;
t61 = (-pkin(2) * t113 - pkin(8) * t110 - pkin(1)) * t105;
t32 = -t109 * t60 + t112 * t61;
t18 = -pkin(3) * t124 - qJ(4) * t67 + t32;
t33 = t109 * t61 + t112 * t60;
t22 = qJ(4) * t66 + t33;
t9 = t104 * t18 + t106 * t22;
t131 = Ifges(6,4) * t108;
t130 = Ifges(6,4) * t111;
t129 = t108 * t75;
t92 = pkin(3) * t104 + pkin(9);
t128 = t108 * t92;
t127 = t111 * t75;
t126 = t111 * t92;
t123 = t108 ^ 2 + t111 ^ 2;
t122 = t109 ^ 2 + t112 ^ 2;
t40 = t104 * t67 - t106 * t66;
t3 = Ifges(6,5) * t26 + Ifges(6,6) * t25 + Ifges(6,3) * t40;
t121 = Ifges(3,5) * t125 + Ifges(3,6) * t124 + Ifges(3,3) * t107;
t94 = -pkin(3) * t112 - pkin(2);
t16 = t40 * mrSges(5,1) + t41 * mrSges(5,2);
t48 = t74 * mrSges(5,1) + t75 * mrSges(5,2);
t119 = -Ifges(4,5) * t67 - Ifges(5,5) * t41 - Ifges(4,6) * t66 + Ifges(5,6) * t40;
t59 = t87 + (-pkin(2) - t138) * t107;
t42 = -pkin(3) * t66 + t59;
t10 = pkin(4) * t40 - pkin(9) * t41 + t42;
t7 = -pkin(9) * t124 + t9;
t1 = t10 * t111 - t108 * t7;
t2 = t10 * t108 + t111 * t7;
t118 = -t1 * t108 + t111 * t2;
t117 = mrSges(6,1) * t108 + mrSges(6,2) * t111;
t8 = -t104 * t22 + t106 * t18;
t45 = pkin(4) * t74 - pkin(9) * t75 + t94;
t53 = t104 * t120 - t106 * t79;
t19 = -t108 * t53 + t111 * t45;
t20 = t108 * t45 + t111 * t53;
t116 = -t108 * t19 + t111 * t20;
t27 = Ifges(6,5) * t127 - Ifges(6,6) * t129 + Ifges(6,3) * t74;
t93 = -pkin(3) * t106 - pkin(4);
t84 = Ifges(4,1) * t109 + Ifges(4,4) * t112;
t83 = Ifges(6,1) * t108 + t130;
t82 = Ifges(4,4) * t109 + Ifges(4,2) * t112;
t81 = Ifges(6,2) * t111 + t131;
t78 = -mrSges(4,1) * t112 + mrSges(4,2) * t109;
t77 = -mrSges(6,1) * t111 + mrSges(6,2) * t108;
t55 = -mrSges(4,1) * t124 - mrSges(4,3) * t67;
t54 = mrSges(4,2) * t124 + mrSges(4,3) * t66;
t50 = Ifges(5,1) * t75 - Ifges(5,4) * t74;
t49 = Ifges(5,4) * t75 - Ifges(5,2) * t74;
t47 = mrSges(6,1) * t74 - mrSges(6,3) * t127;
t46 = -mrSges(6,2) * t74 - mrSges(6,3) * t129;
t44 = -mrSges(4,1) * t66 + mrSges(4,2) * t67;
t43 = t117 * t75;
t35 = Ifges(4,1) * t67 + Ifges(4,4) * t66 - Ifges(4,5) * t124;
t34 = Ifges(4,4) * t67 + Ifges(4,2) * t66 - Ifges(4,6) * t124;
t31 = -mrSges(5,1) * t124 - mrSges(5,3) * t41;
t30 = mrSges(5,2) * t124 - mrSges(5,3) * t40;
t29 = Ifges(6,5) * t74 + (Ifges(6,1) * t111 - t131) * t75;
t28 = Ifges(6,6) * t74 + (-Ifges(6,2) * t108 + t130) * t75;
t15 = Ifges(5,1) * t41 - Ifges(5,4) * t40 - Ifges(5,5) * t124;
t14 = Ifges(5,4) * t41 - Ifges(5,2) * t40 - Ifges(5,6) * t124;
t13 = mrSges(6,1) * t40 - mrSges(6,3) * t26;
t12 = -mrSges(6,2) * t40 + mrSges(6,3) * t25;
t11 = -mrSges(6,1) * t25 + mrSges(6,2) * t26;
t6 = pkin(4) * t124 - t8;
t5 = Ifges(6,1) * t26 + Ifges(6,4) * t25 + Ifges(6,5) * t40;
t4 = Ifges(6,4) * t26 + Ifges(6,2) * t25 + Ifges(6,6) * t40;
t17 = [0.2e1 * t1 * t13 + 0.2e1 * t6 * t11 + 0.2e1 * t2 * t12 + t41 * t15 + 0.2e1 * t42 * t16 + t25 * t4 + t26 * t5 + 0.2e1 * t9 * t30 + 0.2e1 * t8 * t31 + 0.2e1 * t32 * t55 + 0.2e1 * t33 * t54 + t66 * t34 + t67 * t35 + 0.2e1 * t59 * t44 + Ifges(2,3) + (-t14 + t3) * t40 + (t121 - 0.2e1 * t136 + 0.2e1 * t137) * t107 + ((-0.2e1 * t68 * mrSges(3,3) + Ifges(3,5) * t107 + (-0.2e1 * pkin(1) * mrSges(3,2) + Ifges(3,1) * t110) * t105) * t110 + (0.2e1 * t69 * mrSges(3,3) + Ifges(3,6) * t107 + (0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * Ifges(3,4) * t110 + (Ifges(3,2) + t135) * t113) * t105 + t119) * t113) * t105 + m(3) * (pkin(1) ^ 2 * t105 ^ 2 + t68 ^ 2 + t69 ^ 2) + m(4) * (t32 ^ 2 + t33 ^ 2 + t59 ^ 2) + m(5) * (t42 ^ 2 + t8 ^ 2 + t9 ^ 2) + m(6) * (t1 ^ 2 + t2 ^ 2 + t6 ^ 2); m(4) * (-pkin(2) * t59 + (-t32 * t109 + t33 * t112) * pkin(8)) - t148 * t124 / 0.2e1 + (t11 - t31) * t51 - t136 + t137 + (pkin(8) * t54 + t33 * mrSges(4,3) + t34 / 0.2e1) * t112 + (t35 / 0.2e1 - pkin(8) * t55 - t32 * mrSges(4,3)) * t109 + t94 * t16 + t59 * t78 + t66 * t82 / 0.2e1 + t67 * t84 / 0.2e1 + t53 * t30 + t6 * t43 - pkin(2) * t44 + t2 * t46 + t1 * t47 + t42 * t48 + t41 * t50 / 0.2e1 + t19 * t13 + t20 * t12 + t121 + t28 * t145 + t29 * t144 + (t5 * t140 + t15 / 0.2e1 + t4 * t142 - t8 * mrSges(5,3)) * t75 + (-t49 / 0.2e1 + t27 / 0.2e1) * t40 + (t3 / 0.2e1 - t14 / 0.2e1 - t9 * mrSges(5,3)) * t74 + m(5) * (t42 * t94 - t51 * t8 + t53 * t9) + m(6) * (t1 * t19 + t2 * t20 + t51 * t6); -0.2e1 * pkin(2) * t78 + t109 * t84 + t112 * t82 + 0.2e1 * t19 * t47 + 0.2e1 * t20 * t46 + t43 * t146 + 0.2e1 * t94 * t48 + Ifges(3,3) + 0.2e1 * t122 * pkin(8) * mrSges(4,3) + (-0.2e1 * mrSges(5,3) * t53 + t27 - t49) * t74 + m(6) * (t19 ^ 2 + t20 ^ 2 + t147) + m(5) * (t53 ^ 2 + t94 ^ 2 + t147) + m(4) * (t122 * pkin(8) ^ 2 + pkin(2) ^ 2) + (mrSges(5,3) * t146 - t108 * t28 + t111 * t29 + t50) * t75; (m(5) * (t104 * t9 + t106 * t8) + t106 * t31 + t104 * t30) * pkin(3) - t135 * t124 + t12 * t126 + t4 * t140 + t5 * t141 + t93 * t11 + t6 * t77 + t40 * t143 + t81 * t145 + t83 * t144 + t32 * mrSges(4,1) - t33 * mrSges(4,2) + t8 * mrSges(5,1) - t9 * mrSges(5,2) - t13 * t128 - t119 + m(6) * (t118 * t92 + t6 * t93) + t118 * mrSges(6,3); t46 * t126 - t47 * t128 + t51 * t77 + t74 * t143 + t29 * t141 + t28 * t140 + t93 * t43 + m(6) * (t116 * t92 + t51 * t93) - t53 * mrSges(5,2) - t51 * mrSges(5,1) + (t83 * t140 + t81 * t142) * t75 + (-mrSges(4,1) * t109 - mrSges(4,2) * t112) * pkin(8) + t116 * mrSges(6,3) + (m(5) * (t104 * t53 - t106 * t51) + (-t104 * t74 - t106 * t75) * mrSges(5,3)) * pkin(3) + t148; t108 * t83 + t111 * t81 + 0.2e1 * t93 * t77 + m(6) * (t123 * t92 ^ 2 + t93 ^ 2) + m(5) * (t104 ^ 2 + t106 ^ 2) * pkin(3) ^ 2 + t135 + 0.2e1 * (mrSges(5,1) * t106 - mrSges(5,2) * t104) * pkin(3) + 0.2e1 * t123 * t92 * mrSges(6,3); t108 * t12 + t111 * t13 + m(6) * (t1 * t111 + t108 * t2) + m(5) * t42 + t16; t108 * t46 + t111 * t47 + m(6) * (t108 * t20 + t111 * t19) + m(5) * t94 + t48; 0; m(6) * t123 + m(5); mrSges(6,1) * t1 - mrSges(6,2) * t2 + t3; mrSges(6,1) * t19 - mrSges(6,2) * t20 + t27; -t117 * t92 + t80; -t77; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t17(1), t17(2), t17(4), t17(7), t17(11); t17(2), t17(3), t17(5), t17(8), t17(12); t17(4), t17(5), t17(6), t17(9), t17(13); t17(7), t17(8), t17(9), t17(10), t17(14); t17(11), t17(12), t17(13), t17(14), t17(15);];
Mq = res;
