% Calculate joint inertia matrix for
% S5RRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2022-01-20 12:09
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR6_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR6_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR6_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR6_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 12:07:39
% EndTime: 2022-01-20 12:07:39
% DurationCPUTime: 0.53s
% Computational Cost: add. (1090->169), mult. (2000->229), div. (0->0), fcn. (2012->8), ass. (0->79)
t71 = cos(qJ(3));
t114 = t71 ^ 2;
t66 = sin(qJ(4));
t67 = sin(qJ(3));
t70 = cos(qJ(4));
t46 = t66 * t71 + t70 * t67;
t103 = t46 * pkin(9);
t68 = sin(qJ(2));
t55 = t68 * pkin(1) + pkin(7);
t41 = (-pkin(8) - t55) * t67;
t62 = t71 * pkin(8);
t42 = t71 * t55 + t62;
t19 = t70 * t41 - t66 * t42;
t11 = t19 - t103;
t20 = t66 * t41 + t70 * t42;
t45 = -t66 * t67 + t70 * t71;
t40 = t45 * pkin(9);
t12 = t40 + t20;
t65 = sin(qJ(5));
t69 = cos(qJ(5));
t2 = t69 * t11 - t65 * t12;
t3 = t65 * t11 + t69 * t12;
t113 = t2 * mrSges(6,1) - t3 * mrSges(6,2);
t52 = (-pkin(8) - pkin(7)) * t67;
t53 = t71 * pkin(7) + t62;
t27 = t70 * t52 - t66 * t53;
t28 = t66 * t52 + t70 * t53;
t112 = t27 * mrSges(5,1) - t28 * mrSges(5,2);
t111 = Ifges(5,5) * t46 + Ifges(5,6) * t45;
t15 = t27 - t103;
t16 = t40 + t28;
t7 = t69 * t15 - t65 * t16;
t8 = t65 * t15 + t69 * t16;
t110 = t7 * mrSges(6,1) - t8 * mrSges(6,2);
t109 = t19 * mrSges(5,1) - t20 * mrSges(5,2);
t23 = t69 * t45 - t65 * t46;
t99 = t23 * mrSges(6,3);
t108 = t65 * pkin(4) * t99 + t111;
t24 = t65 * t45 + t69 * t46;
t9 = -t23 * mrSges(6,1) + t24 * mrSges(6,2);
t107 = 0.2e1 * t9;
t25 = -t45 * mrSges(5,1) + t46 * mrSges(5,2);
t106 = 0.2e1 * t25;
t105 = pkin(3) * t66;
t72 = cos(qJ(2));
t102 = t72 * pkin(1);
t98 = t24 * mrSges(6,3);
t56 = t70 * pkin(3) + pkin(4);
t37 = t69 * t105 + t65 * t56;
t96 = t37 * mrSges(6,2);
t95 = t65 * mrSges(6,2);
t94 = Ifges(5,3) + Ifges(6,3);
t93 = Ifges(6,5) * t24 + Ifges(6,6) * t23;
t92 = t67 ^ 2 + t114;
t91 = 2 * mrSges(5,3);
t90 = 0.2e1 * mrSges(6,3);
t89 = pkin(4) * t95;
t88 = t69 * t98;
t87 = t70 * t46 * mrSges(5,3);
t58 = -t71 * pkin(3) - pkin(2);
t86 = t92 * t55;
t36 = -t65 * t105 + t69 * t56;
t31 = t36 * mrSges(6,1);
t85 = Ifges(6,3) + t31 - t96;
t84 = -t67 * mrSges(4,1) - t71 * mrSges(4,2);
t83 = t93 + t113;
t82 = t93 + t110;
t81 = 0.2e1 * t92 * mrSges(4,3);
t30 = -t45 * pkin(4) + t58;
t80 = Ifges(5,1) * t46 ^ 2 + Ifges(6,1) * t24 ^ 2 + Ifges(4,2) * t114 + Ifges(3,3) + (Ifges(4,1) * t67 + 0.2e1 * Ifges(4,4) * t71) * t67 + (0.2e1 * Ifges(5,4) * t46 + Ifges(5,2) * t45) * t45 + (0.2e1 * Ifges(6,4) * t24 + Ifges(6,2) * t23) * t23;
t79 = (t72 * mrSges(3,1) - t68 * mrSges(3,2)) * pkin(1);
t78 = (t70 * mrSges(5,1) - t66 * mrSges(5,2)) * pkin(3);
t77 = t45 * mrSges(5,3) * t105 + Ifges(4,5) * t67 + Ifges(4,6) * t71 - t36 * t98 + t37 * t99 + t111 + t93;
t59 = t69 * pkin(4) * mrSges(6,1);
t57 = -pkin(2) - t102;
t51 = -t71 * mrSges(4,1) + t67 * mrSges(4,2);
t50 = t58 - t102;
t29 = t30 - t102;
t1 = [t80 + (-t19 * t46 + t20 * t45) * t91 + (-t2 * t24 + t3 * t23) * t90 + m(3) * (t68 ^ 2 + t72 ^ 2) * pkin(1) ^ 2 + m(4) * (t92 * t55 ^ 2 + t57 ^ 2) + t55 * t81 + m(6) * (t2 ^ 2 + t29 ^ 2 + t3 ^ 2) + m(5) * (t19 ^ 2 + t20 ^ 2 + t50 ^ 2) + 0.2e1 * t57 * t51 + t50 * t106 + t29 * t107 + Ifges(2,3) + 0.2e1 * t79; t80 + t79 + m(4) * (-pkin(2) * t57 + pkin(7) * t86) + (t92 * pkin(7) + t86) * mrSges(4,3) + ((-t19 - t27) * t46 + (t20 + t28) * t45) * mrSges(5,3) + ((-t2 - t7) * t24 + (t3 + t8) * t23) * mrSges(6,3) + m(6) * (t7 * t2 + t30 * t29 + t8 * t3) + m(5) * (t27 * t19 + t28 * t20 + t58 * t50) + (t29 + t30) * t9 + (-pkin(2) + t57) * t51 + (t58 + t50) * t25; -0.2e1 * pkin(2) * t51 + t58 * t106 + t30 * t107 + (t8 * t23 - t7 * t24) * t90 + (-t27 * t46 + t28 * t45) * t91 + pkin(7) * t81 + m(6) * (t30 ^ 2 + t7 ^ 2 + t8 ^ 2) + m(5) * (t27 ^ 2 + t28 ^ 2 + t58 ^ 2) + m(4) * (t92 * pkin(7) ^ 2 + pkin(2) ^ 2) + t80; t77 + (-t87 + m(5) * (t19 * t70 + t20 * t66)) * pkin(3) + m(6) * (t36 * t2 + t37 * t3) + t84 * t55 + t109 + t113; t77 + (-t87 + m(5) * (t27 * t70 + t28 * t66)) * pkin(3) + m(6) * (t36 * t7 + t37 * t8) + t84 * pkin(7) + t110 + t112; -0.2e1 * t96 + Ifges(4,3) + 0.2e1 * t31 + 0.2e1 * t78 + m(6) * (t36 ^ 2 + t37 ^ 2) + m(5) * (t66 ^ 2 + t70 ^ 2) * pkin(3) ^ 2 + t94; (-t88 + m(6) * (t2 * t69 + t3 * t65)) * pkin(4) + t83 + t108 + t109; (-t88 + m(6) * (t65 * t8 + t69 * t7)) * pkin(4) + t82 + t108 + t112; Ifges(5,3) + t59 + t78 + (-t95 + m(6) * (t36 * t69 + t37 * t65)) * pkin(4) + t85; -0.2e1 * t89 + 0.2e1 * t59 + m(6) * (t65 ^ 2 + t69 ^ 2) * pkin(4) ^ 2 + t94; t83; t82; t85; Ifges(6,3) + t59 - t89; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
