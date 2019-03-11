% Calculate joint inertia matrix for
% S6PRPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-03-08 19:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRP1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP1_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP1_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP1_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:55:18
% EndTime: 2019-03-08 19:55:20
% DurationCPUTime: 0.82s
% Computational Cost: add. (664->214), mult. (1536->287), div. (0->0), fcn. (1439->10), ass. (0->91)
t117 = Ifges(6,5) + Ifges(7,5);
t92 = Ifges(6,6) + Ifges(7,6);
t116 = -2 * mrSges(7,3);
t62 = sin(pkin(11));
t48 = t62 * pkin(2) + pkin(8);
t115 = 0.2e1 * t48;
t114 = m(6) + m(7);
t113 = m(6) * pkin(4);
t69 = cos(qJ(5));
t50 = -t69 * pkin(5) - pkin(4);
t112 = m(7) * t50;
t66 = sin(qJ(5));
t111 = t117 * t66 + t92 * t69;
t63 = sin(pkin(6));
t64 = cos(pkin(11));
t68 = sin(qJ(2));
t71 = cos(qJ(2));
t14 = (t62 * t71 + t64 * t68) * t63;
t65 = cos(pkin(6));
t67 = sin(qJ(4));
t70 = cos(qJ(4));
t9 = t14 * t67 - t65 * t70;
t110 = t9 ^ 2;
t12 = (t62 * t68 - t64 * t71) * t63;
t109 = t12 ^ 2;
t108 = m(7) * pkin(5);
t21 = (pkin(5) * t66 + t48) * t67;
t106 = m(7) * t21;
t105 = t67 * pkin(9);
t104 = t70 * pkin(4);
t103 = t70 * t9;
t102 = t9 * t67;
t101 = mrSges(6,2) * t69;
t100 = Ifges(6,4) * t66;
t99 = Ifges(6,4) * t69;
t98 = Ifges(7,4) * t66;
t97 = Ifges(7,4) * t69;
t11 = t14 * t70 + t65 * t67;
t96 = t11 * t70;
t95 = t48 * t70;
t94 = t66 * t67;
t93 = t67 * t69;
t91 = Ifges(6,3) + Ifges(7,3);
t90 = -qJ(6) - pkin(9);
t23 = mrSges(7,1) * t94 + mrSges(7,2) * t93;
t74 = mrSges(6,1) * t66 + t101;
t24 = t74 * t67;
t89 = -t23 - t24;
t26 = t70 * mrSges(7,2) - mrSges(7,3) * t94;
t27 = t70 * mrSges(6,2) - mrSges(6,3) * t94;
t88 = t26 + t27;
t28 = -t70 * mrSges(7,1) - mrSges(7,3) * t93;
t29 = -t70 * mrSges(6,1) - mrSges(6,3) * t93;
t87 = t28 + t29;
t49 = -t64 * pkin(2) - pkin(3);
t25 = -t104 + t49 - t105;
t8 = t66 * t25 + t69 * t95;
t33 = -t69 * mrSges(6,1) + t66 * mrSges(6,2);
t86 = t33 - mrSges(5,1);
t85 = t117 * t93;
t82 = t66 ^ 2 + t69 ^ 2;
t59 = t67 ^ 2;
t61 = t70 ^ 2;
t81 = t59 + t61;
t80 = qJ(6) * t67;
t32 = -t69 * mrSges(7,1) + t66 * mrSges(7,2);
t78 = -t32 - t86;
t77 = -mrSges(6,1) - t108;
t76 = t82 * mrSges(6,3);
t3 = -t11 * t66 + t12 * t69;
t4 = t11 * t69 + t12 * t66;
t75 = -t3 * t66 + t4 * t69;
t57 = t65 ^ 2;
t46 = t48 ^ 2;
t40 = t59 * t46;
t39 = Ifges(6,1) * t66 + t99;
t38 = Ifges(7,1) * t66 + t97;
t37 = Ifges(6,2) * t69 + t100;
t36 = Ifges(7,2) * t69 + t98;
t35 = t90 * t69;
t34 = -t70 * mrSges(5,1) + t67 * mrSges(5,2);
t31 = t90 * t66;
t20 = t69 * t25;
t18 = -Ifges(6,5) * t70 + (Ifges(6,1) * t69 - t100) * t67;
t17 = -Ifges(7,5) * t70 + (Ifges(7,1) * t69 - t98) * t67;
t16 = -Ifges(6,6) * t70 + (-Ifges(6,2) * t66 + t99) * t67;
t15 = -Ifges(7,6) * t70 + (-Ifges(7,2) * t66 + t97) * t67;
t7 = -t66 * t95 + t20;
t6 = -t66 * t80 + t8;
t5 = -t69 * t80 + t20 + (-t48 * t66 - pkin(5)) * t70;
t1 = [m(2) + m(5) * (t11 ^ 2 + t109 + t110) + m(4) * (t14 ^ 2 + t109 + t57) + m(3) * (t57 + (t68 ^ 2 + t71 ^ 2) * t63 ^ 2) + (t3 ^ 2 + t4 ^ 2 + t110) * t114; mrSges(5,3) * t96 - t14 * mrSges(4,2) + (t71 * mrSges(3,1) - t68 * mrSges(3,2)) * t63 + t88 * t4 + t87 * t3 + (-mrSges(4,1) + t34) * t12 + (t67 * mrSges(5,3) - t89) * t9 + m(6) * (t48 * t102 + t7 * t3 + t8 * t4) + m(7) * (t21 * t9 + t5 * t3 + t6 * t4) + m(5) * (t49 * t12 + (t96 + t102) * t48) + m(4) * (-t12 * t64 + t14 * t62) * pkin(2); 0.2e1 * t21 * t23 + 0.2e1 * t6 * t26 + 0.2e1 * t8 * t27 + 0.2e1 * t5 * t28 + 0.2e1 * t7 * t29 + 0.2e1 * t49 * t34 + Ifges(3,3) + Ifges(4,3) + ((Ifges(5,2) + t91) * t70 - t85) * t70 + m(5) * (t61 * t46 + t49 ^ 2 + t40) + m(6) * (t7 ^ 2 + t8 ^ 2 + t40) + m(7) * (t21 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(4) * (t62 ^ 2 + t64 ^ 2) * pkin(2) ^ 2 + (Ifges(5,1) * t67 + 0.2e1 * Ifges(5,4) * t70 + t24 * t115 + (t17 + t18) * t69 + (t92 * t70 - t15 - t16) * t66) * t67 + 0.2e1 * (t64 * mrSges(4,1) - t62 * mrSges(4,2)) * pkin(2) + t81 * mrSges(5,3) * t115; m(4) * t65 + m(5) * (t67 * t11 - t103) + (t75 * t67 - t103) * t114; (t89 - t106) * t70 + (t88 * t69 - t87 * t66 + m(7) * (-t5 * t66 + t6 * t69) + m(6) * (-t66 * t7 + t69 * t8 - t95)) * t67; m(4) + m(5) * t81 + (t82 * t59 + t61) * t114; -t11 * mrSges(5,2) + m(7) * (t31 * t3 - t35 * t4) + (t112 - t78 - t113) * t9 + (m(6) * pkin(9) + mrSges(6,3) + mrSges(7,3)) * t75; t31 * t28 + t21 * t32 - t35 * t26 + t50 * t23 - pkin(4) * t24 + m(7) * (t50 * t21 + t31 * t5 - t35 * t6) + (t15 / 0.2e1 + t16 / 0.2e1 + t8 * mrSges(6,3) + t6 * mrSges(7,3) + (m(6) * t8 + t27) * pkin(9)) * t69 + (t17 / 0.2e1 + t18 / 0.2e1 - t5 * mrSges(7,3) - t7 * mrSges(6,3) + (-m(6) * t7 - t29) * pkin(9)) * t66 + (Ifges(5,5) + (t86 - t113) * t48 + (t38 / 0.2e1 + t39 / 0.2e1) * t69 + (-t36 / 0.2e1 - t37 / 0.2e1) * t66) * t67 + (Ifges(5,6) - t48 * mrSges(5,2) - t111 / 0.2e1) * t70; t78 * t70 + (t82 * mrSges(7,3) - mrSges(5,2) + t76) * t67 + m(6) * (t82 * t105 + t104) + m(7) * (-t50 * t70 + (-t31 * t66 - t35 * t69) * t67); -0.2e1 * pkin(4) * t33 + 0.2e1 * t50 * t32 + Ifges(5,3) + 0.2e1 * pkin(9) * t76 + m(6) * (t82 * pkin(9) ^ 2 + pkin(4) ^ 2) + m(7) * (t31 ^ 2 + t35 ^ 2 + t50 ^ 2) + (t35 * t116 + t36 + t37) * t69 + (t31 * t116 + t38 + t39) * t66; (-mrSges(6,2) - mrSges(7,2)) * t4 + (mrSges(7,1) - t77) * t3; t7 * mrSges(6,1) + t5 * mrSges(7,1) - t8 * mrSges(6,2) - t6 * mrSges(7,2) - t91 * t70 - t92 * t94 + (m(7) * t5 + t28) * pkin(5) + t85; (t77 * t66 - t101) * t67 - t23; t31 * mrSges(7,1) + t35 * mrSges(7,2) - t74 * pkin(9) + (m(7) * t31 - t66 * mrSges(7,3)) * pkin(5) + t111; (0.2e1 * mrSges(7,1) + t108) * pkin(5) + t91; m(7) * t9; t23 + t106; -m(7) * t70; t32 + t112; 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
