% Calculate joint inertia matrix for
% S6PPPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
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
% Datum: 2019-03-08 18:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PPPRRR1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_inertiaJ_slag_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPPRRR1_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPPRRR1_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPPRRR1_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:39:24
% EndTime: 2019-03-08 18:39:26
% DurationCPUTime: 0.63s
% Computational Cost: add. (1197->205), mult. (3482->320), div. (0->0), fcn. (4201->16), ass. (0->94)
t113 = 2 * pkin(10);
t112 = m(7) * pkin(11) + mrSges(7,3);
t68 = sin(qJ(6));
t71 = cos(qJ(6));
t40 = -t71 * mrSges(7,1) + t68 * mrSges(7,2);
t111 = -m(7) * pkin(5) - mrSges(6,1) + t40;
t58 = sin(pkin(14));
t59 = sin(pkin(13));
t62 = sin(pkin(6));
t63 = cos(pkin(14));
t64 = cos(pkin(13));
t66 = cos(pkin(7));
t89 = t64 * t66;
t61 = sin(pkin(7));
t67 = cos(pkin(6));
t90 = t61 * t67;
t17 = t63 * t90 + (-t58 * t59 + t63 * t89) * t62;
t31 = -t62 * t64 * t61 + t67 * t66;
t60 = sin(pkin(8));
t65 = cos(pkin(8));
t12 = -t17 * t60 + t31 * t65;
t69 = sin(qJ(5));
t72 = cos(qJ(5));
t18 = t62 * t59 * t63 + (t62 * t89 + t90) * t58;
t70 = sin(qJ(4));
t73 = cos(qJ(4));
t8 = t18 * t73 + (t17 * t65 + t31 * t60) * t70;
t3 = -t12 * t72 + t8 * t69;
t110 = t3 ^ 2;
t88 = t65 * t73;
t92 = t60 * t73;
t6 = -t17 * t88 + t18 * t70 - t31 * t92;
t109 = t6 ^ 2;
t93 = t60 * t70;
t21 = t66 * t93 + (t63 * t65 * t70 + t58 * t73) * t61;
t91 = t61 * t63;
t30 = -t60 * t91 + t65 * t66;
t13 = t69 * t21 - t72 * t30;
t108 = t13 ^ 2;
t19 = t70 * t61 * t58 - t66 * t92 - t88 * t91;
t107 = t19 ^ 2;
t32 = -t72 * t65 + t69 * t93;
t106 = t32 ^ 2;
t105 = t71 / 0.2e1;
t104 = pkin(10) * t72;
t103 = t13 * t3;
t102 = t19 * t6;
t101 = t3 * t69;
t100 = t32 * t3;
t99 = Ifges(7,4) * t68;
t98 = Ifges(7,4) * t71;
t97 = Ifges(7,6) * t72;
t96 = t13 * t69;
t95 = t32 * t13;
t94 = t32 * t69;
t87 = t68 * t69;
t86 = t69 * t71;
t41 = -t72 * mrSges(6,1) + t69 * mrSges(6,2);
t85 = mrSges(5,1) - t41;
t84 = t68 ^ 2 + t71 ^ 2;
t5 = t12 * t69 + t8 * t72;
t82 = t5 * t72 + t101;
t80 = mrSges(7,1) * t68 + mrSges(7,2) * t71;
t15 = t72 * t21 + t69 * t30;
t79 = t15 * t72 + t96;
t39 = -t72 * pkin(5) - t69 * pkin(11) - pkin(4);
t25 = -t68 * t104 + t71 * t39;
t26 = t71 * t104 + t68 * t39;
t77 = -t25 * t68 + t26 * t71;
t34 = t69 * t65 + t72 * t93;
t76 = t34 * t72 + t94;
t75 = pkin(10) ^ 2;
t57 = t72 ^ 2;
t55 = t69 ^ 2;
t51 = t60 ^ 2;
t50 = t55 * t75;
t49 = Ifges(7,5) * t68;
t48 = Ifges(7,6) * t71;
t46 = t51 * t73 ^ 2;
t45 = Ifges(7,5) * t86;
t43 = Ifges(7,1) * t68 + t98;
t42 = Ifges(7,2) * t71 + t99;
t38 = -t72 * mrSges(7,1) - mrSges(7,3) * t86;
t37 = t72 * mrSges(7,2) - mrSges(7,3) * t87;
t35 = t80 * t69;
t29 = -Ifges(7,5) * t72 + (Ifges(7,1) * t71 - t99) * t69;
t28 = -t97 + (-Ifges(7,2) * t68 + t98) * t69;
t23 = t71 * t34 - t68 * t92;
t22 = -t68 * t34 - t71 * t92;
t10 = t71 * t15 + t68 * t19;
t9 = -t68 * t15 + t71 * t19;
t2 = t5 * t71 + t6 * t68;
t1 = -t5 * t68 + t6 * t71;
t4 = [m(2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t110) + m(6) * (t5 ^ 2 + t109 + t110) + m(5) * (t12 ^ 2 + t8 ^ 2 + t109) + m(4) * (t17 ^ 2 + t18 ^ 2 + t31 ^ 2) + m(3) * (t67 ^ 2 + (t59 ^ 2 + t64 ^ 2) * t62 ^ 2); m(3) * t67 + m(7) * (t9 * t1 + t10 * t2 + t103) + m(6) * (t15 * t5 + t102 + t103) + m(5) * (t30 * t12 + t21 * t8 + t102) + m(4) * (t66 * t31 + (t17 * t63 + t18 * t58) * t61); m(3) + m(7) * (t10 ^ 2 + t9 ^ 2 + t108) + m(6) * (t15 ^ 2 + t107 + t108) + m(5) * (t21 ^ 2 + t30 ^ 2 + t107) + m(4) * (t66 ^ 2 + (t58 ^ 2 + t63 ^ 2) * t61 ^ 2); m(7) * (t22 * t1 + t23 * t2 + t100) + m(6) * (t34 * t5 - t6 * t92 + t100) + m(5) * (t65 * t12 + (-t6 * t73 + t70 * t8) * t60) + m(4) * t31; m(4) * t66 + m(7) * (t23 * t10 + t22 * t9 + t95) + m(6) * (t34 * t15 - t19 * t92 + t95) + m(5) * (t65 * t30 + (-t19 * t73 + t21 * t70) * t60); m(4) + m(7) * (t22 ^ 2 + t23 ^ 2 + t106) + m(6) * (t34 ^ 2 + t106 + t46) + m(5) * (t51 * t70 ^ 2 + t65 ^ 2 + t46); -t8 * mrSges(5,2) + t1 * t38 + t2 * t37 + t3 * t35 - t85 * t6 + t82 * mrSges(6,3) + m(7) * (pkin(10) * t101 + t25 * t1 + t26 * t2) + m(6) * (-pkin(4) * t6 + t82 * pkin(10)); -t21 * mrSges(5,2) + t10 * t37 + t13 * t35 + t9 * t38 - t85 * t19 + t79 * mrSges(6,3) + m(7) * (pkin(10) * t96 + t26 * t10 + t25 * t9) + m(6) * (-pkin(4) * t19 + t79 * pkin(10)); t22 * t38 + t23 * t37 + t32 * t35 + t76 * mrSges(6,3) + (-t70 * mrSges(5,2) + t85 * t73) * t60 + m(7) * (pkin(10) * t94 + t25 * t22 + t26 * t23) + m(6) * (pkin(4) * t92 + t76 * pkin(10)); -0.2e1 * pkin(4) * t41 + 0.2e1 * t25 * t38 + 0.2e1 * t26 * t37 + Ifges(5,3) + (t55 + t57) * mrSges(6,3) * t113 + m(7) * (t25 ^ 2 + t26 ^ 2 + t50) + m(6) * (pkin(4) ^ 2 + t57 * t75 + t50) + (-t45 + (Ifges(7,3) + Ifges(6,2)) * t72) * t72 + (Ifges(6,1) * t69 + 0.2e1 * Ifges(6,4) * t72 + t35 * t113 + t71 * t29 + (-t28 + t97) * t68) * t69; -t5 * mrSges(6,2) + t112 * (-t1 * t68 + t2 * t71) + t111 * t3; -t15 * mrSges(6,2) + t112 * (t10 * t71 - t9 * t68) + t111 * t13; -t34 * mrSges(6,2) + t112 * (-t22 * t68 + t23 * t71) + t111 * t32; -pkin(5) * t35 + t68 * t29 / 0.2e1 + t28 * t105 + (-t49 / 0.2e1 - t48 / 0.2e1 + Ifges(6,6) - pkin(10) * mrSges(6,2)) * t72 + t77 * mrSges(7,3) + (m(7) * t77 + t71 * t37 - t68 * t38) * pkin(11) + (Ifges(6,5) + t43 * t105 - t68 * t42 / 0.2e1 + t111 * pkin(10)) * t69; Ifges(6,3) + m(7) * (t84 * pkin(11) ^ 2 + pkin(5) ^ 2) - 0.2e1 * pkin(5) * t40 + t68 * t43 + t71 * t42 + 0.2e1 * t84 * pkin(11) * mrSges(7,3); t1 * mrSges(7,1) - t2 * mrSges(7,2); t9 * mrSges(7,1) - t10 * mrSges(7,2); t22 * mrSges(7,1) - t23 * mrSges(7,2); t25 * mrSges(7,1) - t26 * mrSges(7,2) - Ifges(7,6) * t87 - Ifges(7,3) * t72 + t45; -t80 * pkin(11) + t48 + t49; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t4(1) t4(2) t4(4) t4(7) t4(11) t4(16); t4(2) t4(3) t4(5) t4(8) t4(12) t4(17); t4(4) t4(5) t4(6) t4(9) t4(13) t4(18); t4(7) t4(8) t4(9) t4(10) t4(14) t4(19); t4(11) t4(12) t4(13) t4(14) t4(15) t4(20); t4(16) t4(17) t4(18) t4(19) t4(20) t4(21);];
Mq  = res;
