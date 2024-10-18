% Calculate joint inertia matrix for
% S5RRRRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d1,d2,d3,d4,d5]';
% m [6x1]
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
% Datum: 2024-09-27 22:28
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR15_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR15_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR15_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR15_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR15_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR15_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 22:23:59
% EndTime: 2024-09-27 22:23:59
% DurationCPUTime: 0.34s
% Computational Cost: add. (2896->243), mult. (7362->362), div. (0->0), fcn. (8077->12), ass. (0->114)
t111 = sin(qJ(4));
t115 = cos(qJ(4));
t107 = sin(pkin(5));
t112 = sin(qJ(3));
t113 = sin(qJ(2));
t116 = cos(qJ(3));
t117 = cos(qJ(2));
t72 = (-t112 * t113 + t116 * t117) * t107;
t73 = (t112 * t117 + t113 * t116) * t107;
t47 = -t111 * t73 + t115 * t72;
t48 = t111 * t72 + t115 * t73;
t152 = Ifges(5,5) * t48 + Ifges(5,6) * t47;
t108 = cos(pkin(6));
t114 = cos(qJ(5));
t131 = t108 * t114;
t106 = sin(pkin(6));
t110 = sin(qJ(5));
t136 = t106 * t110;
t79 = pkin(4) * t131 - pkin(11) * t136;
t83 = mrSges(6,1) * t108 - mrSges(6,3) * t136;
t54 = t79 * t83;
t132 = t108 * t110;
t135 = t106 * t114;
t81 = pkin(4) * t132 + pkin(11) * t135;
t84 = -mrSges(6,2) * t108 + mrSges(6,3) * t135;
t55 = t81 * t84;
t151 = t54 + t55;
t104 = t106 * pkin(11);
t86 = pkin(3) * t111 + t104;
t146 = pkin(3) * t115;
t97 = pkin(4) + t146;
t61 = -t110 * t86 + t97 * t131;
t49 = t61 * t83;
t62 = t114 * t86 + t97 * t132;
t50 = t62 * t84;
t99 = mrSges(5,1) * t146;
t150 = t49 + t50 + t99;
t109 = cos(pkin(5));
t133 = t107 * t117;
t134 = t107 * t113;
t149 = Ifges(3,5) * t134 + Ifges(3,6) * t133 + Ifges(3,3) * t109;
t68 = Ifges(4,5) * t73;
t67 = Ifges(4,6) * t72;
t148 = pkin(1) * t109;
t147 = pkin(2) * t112;
t145 = pkin(4) * t106;
t98 = pkin(2) * t116 + pkin(3);
t77 = t111 * t98 + t115 * t147;
t144 = t77 * mrSges(5,2);
t95 = t117 * t148;
t80 = -pkin(8) * t134 + t95;
t143 = t80 * mrSges(3,1);
t82 = pkin(8) * t133 + t113 * t148;
t142 = t82 * mrSges(3,2);
t56 = pkin(2) * t109 + t95 + (-pkin(8) - pkin(9)) * t134;
t65 = pkin(9) * t133 + t82;
t38 = -t112 * t65 + t116 * t56;
t30 = pkin(3) * t109 - pkin(10) * t73 + t38;
t39 = t112 * t56 + t116 * t65;
t33 = pkin(10) * t72 + t39;
t20 = t111 * t30 + t115 * t33;
t105 = t106 ^ 2;
t141 = t105 * t97;
t76 = -t111 * t147 + t115 * t98;
t75 = pkin(4) + t76;
t140 = t106 * t75;
t78 = (-mrSges(6,1) * t114 + mrSges(6,2) * t110) * t106;
t139 = t106 * t78;
t138 = t106 * t97;
t137 = t111 * mrSges(5,2);
t130 = -0.2e1 * t139;
t129 = pkin(3) * t137;
t124 = t106 * t109 + t108 * t47;
t26 = -t110 * t48 + t124 * t114;
t27 = t124 * t110 + t114 * t48;
t40 = -t106 * t47 + t108 * t109;
t10 = Ifges(6,5) * t27 + Ifges(6,6) * t26 + Ifges(6,3) * t40;
t128 = Ifges(5,3) * t109 + t152;
t69 = Ifges(6,5) * t136 + Ifges(6,6) * t135 + Ifges(6,3) * t108;
t19 = -t111 * t33 + t115 * t30;
t70 = Ifges(6,6) * t108 + (Ifges(6,4) * t110 + Ifges(6,2) * t114) * t106;
t71 = Ifges(6,5) * t108 + (Ifges(6,1) * t110 + Ifges(6,4) * t114) * t106;
t127 = t108 * t69 + t70 * t135 + t71 * t136 + Ifges(5,3);
t126 = Ifges(4,3) + t127;
t85 = (-pkin(2) * t117 - pkin(1)) * t107;
t15 = -pkin(11) * t108 * t48 + pkin(4) * t109 + t19;
t51 = -pkin(3) * t72 + t85;
t23 = -pkin(4) * t47 - t48 * t104 + t51;
t125 = t106 * t23 + t108 * t15;
t123 = (mrSges(4,1) * t116 - mrSges(4,2) * t112) * pkin(2);
t66 = t104 + t77;
t43 = -t110 * t66 + t75 * t131;
t34 = t43 * t83;
t44 = t114 * t66 + t75 * t132;
t35 = t44 * t84;
t74 = t76 * mrSges(5,1);
t122 = t127 + t34 + t35 + t74 - t144;
t11 = Ifges(6,4) * t27 + Ifges(6,2) * t26 + Ifges(6,6) * t40;
t12 = Ifges(6,1) * t27 + Ifges(6,4) * t26 + Ifges(6,5) * t40;
t13 = t124 * pkin(11) + t20;
t3 = -t110 * t13 + t125 * t114;
t4 = t125 * t110 + t114 * t13;
t6 = -t106 * t15 + t108 * t23;
t121 = -t20 * mrSges(5,2) + t12 * t136 / 0.2e1 + t3 * t83 + t4 * t84 + t6 * t78 + t128 + t19 * mrSges(5,1) + t26 * t70 / 0.2e1 + t27 * t71 / 0.2e1 + t40 * t69 / 0.2e1 + t11 * t135 / 0.2e1 + t108 * t10 / 0.2e1;
t102 = Ifges(4,3) * t109;
t120 = t38 * mrSges(4,1) - t39 * mrSges(4,2) + t102 + t121 + t67 + t68;
t60 = mrSges(4,1) * t109 - mrSges(4,3) * t73;
t59 = -mrSges(4,2) * t109 + mrSges(4,3) * t72;
t42 = mrSges(5,1) * t109 - mrSges(5,3) * t48;
t41 = -mrSges(5,2) * t109 + mrSges(5,3) * t47;
t17 = mrSges(6,1) * t40 - mrSges(6,3) * t27;
t16 = -mrSges(6,2) * t40 + mrSges(6,3) * t26;
t14 = -mrSges(6,1) * t26 + mrSges(6,2) * t27;
t1 = [Ifges(2,3) + t73 * (Ifges(4,1) * t73 + Ifges(4,4) * t72) + t72 * (Ifges(4,4) * t73 + Ifges(4,2) * t72) + t48 * (Ifges(5,1) * t48 + Ifges(5,4) * t47) + t47 * (Ifges(5,4) * t48 + Ifges(5,2) * t47) + 0.2e1 * t85 * (-mrSges(4,1) * t72 + mrSges(4,2) * t73) + 0.2e1 * t51 * (-mrSges(5,1) * t47 + mrSges(5,2) * t48) + 0.2e1 * t39 * t59 + 0.2e1 * t38 * t60 + t40 * t10 + 0.2e1 * t20 * t41 + 0.2e1 * t19 * t42 + t26 * t11 + t27 * t12 + 0.2e1 * t6 * t14 + 0.2e1 * t4 * t16 + 0.2e1 * t3 * t17 + (t102 + t128 - 0.2e1 * t142 + 0.2e1 * t143 + 0.2e1 * t67 + 0.2e1 * t68 + t149 + t152) * t109 + m(3) * (t80 ^ 2 + t82 ^ 2) + m(4) * (t38 ^ 2 + t39 ^ 2 + t85 ^ 2) + m(5) * (t19 ^ 2 + t20 ^ 2 + t51 ^ 2) + m(6) * (t3 ^ 2 + t4 ^ 2 + t6 ^ 2) + ((t113 * Ifges(3,5) + t117 * Ifges(3,6)) * t109 + 0.2e1 * (-t80 * t113 + t82 * t117) * mrSges(3,3) + (t113 * (Ifges(3,1) * t113 + Ifges(3,4) * t117) - 0.2e1 * pkin(1) * (-mrSges(3,1) * t117 + mrSges(3,2) * t113) + t117 * (Ifges(3,4) * t113 + Ifges(3,2) * t117) + m(3) * pkin(1) ^ 2) * t107) * t107; t143 - t142 + t76 * t42 + t77 * t41 + t43 * t17 + t44 * t16 + t120 - t14 * t140 + (t112 * t59 + m(4) * (t112 * t39 + t116 * t38) + t116 * t60) * pkin(2) + m(6) * (-t6 * t140 + t3 * t43 + t4 * t44) + m(5) * (t19 * t76 + t20 * t77) + t149; t75 * t130 - 0.2e1 * t144 + Ifges(3,3) + 0.2e1 * t34 + 0.2e1 * t35 + 0.2e1 * t74 + 0.2e1 * t123 + m(6) * (t105 * t75 ^ 2 + t43 ^ 2 + t44 ^ 2) + m(5) * (t76 ^ 2 + t77 ^ 2) + m(4) * (t112 ^ 2 + t116 ^ 2) * pkin(2) ^ 2 + t126; m(6) * (-t6 * t138 + t3 * t61 + t4 * t62) + t61 * t17 + t62 * t16 + t120 - t14 * t138 + (m(5) * (t111 * t20 + t115 * t19) + t115 * t42 + t111 * t41) * pkin(3); Ifges(4,3) + m(6) * (t75 * t141 + t43 * t61 + t44 * t62) + t123 + (-t75 - t97) * t139 + t122 + (m(5) * (t111 * t77 + t115 * t76) - t137) * pkin(3) + t150; -0.2e1 * t129 + t97 * t130 + 0.2e1 * t49 + 0.2e1 * t50 + 0.2e1 * t99 + m(6) * (t105 * t97 ^ 2 + t61 ^ 2 + t62 ^ 2) + m(5) * (t111 ^ 2 + t115 ^ 2) * pkin(3) ^ 2 + t126; t79 * t17 + t81 * t16 + t121 - t14 * t145 + m(6) * (-t6 * t145 + t3 * t79 + t4 * t81); m(6) * (pkin(4) * t105 * t75 + t43 * t79 + t44 * t81) + (-pkin(4) - t75) * t139 + t122 + t151; -t129 + m(6) * (pkin(4) * t141 + t61 * t79 + t62 * t81) + (-pkin(4) - t97) * t139 + t127 + t150 + t151; pkin(4) * t130 + 0.2e1 * t55 + 0.2e1 * t54 + m(6) * (pkin(4) ^ 2 * t105 + t79 ^ 2 + t81 ^ 2) + t127; mrSges(6,1) * t3 - mrSges(6,2) * t4 + t10; mrSges(6,1) * t43 - mrSges(6,2) * t44 + t69; mrSges(6,1) * t61 - mrSges(6,2) * t62 + t69; mrSges(6,1) * t79 - mrSges(6,2) * t81 + t69; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
