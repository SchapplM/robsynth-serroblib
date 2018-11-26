% Calculate joint inertia matrix for
% S6RPPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2]';
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
% Datum: 2018-11-23 15:48
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRR3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR3_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR3_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR3_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR3_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR3_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:48:37
% EndTime: 2018-11-23 15:48:37
% DurationCPUTime: 0.64s
% Computational Cost: add. (869->218), mult. (1661->310), div. (0->0), fcn. (1479->8), ass. (0->91)
t65 = cos(pkin(10));
t54 = -t65 * pkin(1) - pkin(2);
t50 = -pkin(7) + t54;
t108 = -0.2e1 * t50;
t71 = cos(qJ(4));
t67 = sin(qJ(5));
t70 = cos(qJ(5));
t79 = t67 * mrSges(6,1) + t70 * mrSges(6,2);
t33 = t79 * t71;
t66 = sin(qJ(6));
t69 = cos(qJ(6));
t38 = t66 * t70 + t69 * t67;
t28 = t38 * t71;
t37 = -t66 * t67 + t69 * t70;
t30 = t37 * t71;
t7 = t28 * mrSges(7,1) + t30 * mrSges(7,2);
t107 = -t33 - t7;
t68 = sin(qJ(4));
t106 = t71 * t68;
t92 = t70 * t71;
t105 = Ifges(6,5) * t92 + Ifges(6,3) * t68;
t64 = sin(pkin(10));
t52 = t64 * pkin(1) + qJ(3);
t104 = t52 ^ 2;
t103 = 0.2e1 * t52;
t102 = m(7) * pkin(5);
t101 = t70 / 0.2e1;
t100 = -pkin(9) - pkin(8);
t99 = t68 * pkin(4);
t61 = t68 ^ 2;
t63 = t71 ^ 2;
t88 = t63 + t61;
t98 = m(5) * t88 + m(4);
t97 = Ifges(6,4) * t67;
t96 = Ifges(6,4) * t70;
t95 = Ifges(6,6) * t68;
t94 = t50 * t68;
t93 = t67 * t71;
t91 = t71 * mrSges(6,3);
t32 = -t71 * pkin(8) + t52 + t99;
t10 = t67 * t32 + t70 * t94;
t44 = -t70 * mrSges(6,1) + t67 * mrSges(6,2);
t90 = t44 - mrSges(5,1);
t89 = t67 ^ 2 + t70 ^ 2;
t11 = -t37 * mrSges(7,1) + t38 * mrSges(7,2);
t87 = -t11 - t90;
t86 = Ifges(7,5) * t30 - Ifges(7,6) * t28 + Ifges(7,3) * t68;
t85 = pkin(8) * t89;
t84 = t89 * mrSges(6,3);
t83 = t88 * mrSges(5,3);
t27 = t38 * t68;
t29 = t37 * t68;
t82 = -t27 * mrSges(7,1) - t29 * mrSges(7,2);
t26 = t70 * t32;
t9 = -t67 * t94 + t26;
t81 = t10 * t70 - t67 * t9;
t80 = -mrSges(5,2) + t84;
t40 = -t68 * mrSges(6,2) - t67 * t91;
t41 = t68 * mrSges(6,1) - t70 * t91;
t78 = t70 * t40 - t67 * t41;
t47 = t100 * t67;
t48 = t100 * t70;
t15 = t69 * t47 + t66 * t48;
t16 = t66 * t47 - t69 * t48;
t34 = Ifges(7,6) * t37;
t35 = Ifges(7,5) * t38;
t77 = t15 * mrSges(7,1) - t16 * mrSges(7,2) + t34 + t35;
t6 = -pkin(9) * t92 + t26 + (-t50 * t67 + pkin(5)) * t68;
t8 = -pkin(9) * t93 + t10;
t2 = t69 * t6 - t66 * t8;
t3 = t66 * t6 + t69 * t8;
t76 = t2 * mrSges(7,1) - t3 * mrSges(7,2) + t86;
t75 = (t69 * mrSges(7,1) - t66 * mrSges(7,2)) * pkin(5);
t59 = Ifges(6,5) * t67;
t58 = Ifges(6,6) * t70;
t55 = -t70 * pkin(5) - pkin(4);
t49 = t50 ^ 2;
t46 = Ifges(6,1) * t67 + t96;
t45 = Ifges(6,2) * t70 + t97;
t43 = t63 * t50;
t42 = t63 * t49;
t31 = (pkin(5) * t67 - t50) * t71;
t24 = Ifges(6,5) * t68 + (Ifges(6,1) * t70 - t97) * t71;
t23 = t95 + (-Ifges(6,2) * t67 + t96) * t71;
t18 = t68 * mrSges(7,1) - t30 * mrSges(7,3);
t17 = -t68 * mrSges(7,2) - t28 * mrSges(7,3);
t13 = Ifges(7,1) * t38 + Ifges(7,4) * t37;
t12 = Ifges(7,4) * t38 + Ifges(7,2) * t37;
t5 = Ifges(7,1) * t30 - Ifges(7,4) * t28 + Ifges(7,5) * t68;
t4 = Ifges(7,4) * t30 - Ifges(7,2) * t28 + Ifges(7,6) * t68;
t1 = [0.2e1 * t54 * mrSges(4,2) + mrSges(4,3) * t103 + 0.2e1 * t10 * t40 + 0.2e1 * t3 * t17 + 0.2e1 * t2 * t18 - t28 * t4 + t30 * t5 + 0.2e1 * t31 * t7 + 0.2e1 * t9 * t41 + Ifges(4,1) + Ifges(2,3) + Ifges(3,3) + (mrSges(5,1) * t103 + Ifges(5,2) * t68 + t105 + t86) * t68 + (mrSges(5,2) * t103 + Ifges(5,1) * t71 - 0.2e1 * Ifges(5,4) * t68 + t70 * t24 + t33 * t108 + (-t23 - t95) * t67) * t71 + m(5) * (t61 * t49 + t104 + t42) + m(4) * (t54 ^ 2 + t104) + m(6) * (t10 ^ 2 + t9 ^ 2 + t42) + m(7) * (t2 ^ 2 + t3 ^ 2 + t31 ^ 2) + m(3) * (t64 ^ 2 + t65 ^ 2) * pkin(1) ^ 2 + 0.2e1 * (t65 * mrSges(3,1) - t64 * mrSges(3,2)) * pkin(1) + t83 * t108; t30 * t17 - t28 * t18 + m(7) * (-t2 * t28 + t3 * t30) + (m(6) * (t81 - t94) + t78) * t71 + (m(7) * t31 - t107) * t68; m(3) + m(7) * (t28 ^ 2 + t30 ^ 2 + t61) + m(6) * (t89 * t63 + t61) + t98; t29 * t17 - t27 * t18 + mrSges(4,2) + t107 * t71 + t78 * t68 - t83 + m(7) * (-t27 * t2 + t29 * t3 - t71 * t31) + m(6) * (t81 * t68 + t43) + m(5) * (t61 * t50 + t43) + m(4) * t54; m(6) * (-0.1e1 + t89) * t106 + m(7) * (t27 * t28 + t29 * t30 - t106); m(6) * (t89 * t61 + t63) + m(7) * (t27 ^ 2 + t29 ^ 2 + t63) + t98; t23 * t101 + t67 * t24 / 0.2e1 + t55 * t7 + t30 * t13 / 0.2e1 + t31 * t11 - pkin(4) * t33 + t37 * t4 / 0.2e1 + t38 * t5 / 0.2e1 - t28 * t12 / 0.2e1 + t16 * t17 + t15 * t18 + m(7) * (t15 * t2 + t16 * t3 + t55 * t31) + (t59 / 0.2e1 + t58 / 0.2e1 + t35 / 0.2e1 + t34 / 0.2e1 - Ifges(5,6) - t50 * mrSges(5,2)) * t68 + (-t2 * t38 + t3 * t37) * mrSges(7,3) + t81 * mrSges(6,3) + (m(6) * t81 + t78) * pkin(8) + (t46 * t101 - t67 * t45 / 0.2e1 + Ifges(5,5) + (m(6) * pkin(4) - t90) * t50) * t71; (t28 * t38 + t30 * t37) * mrSges(7,3) + t80 * t71 - t87 * t68 + m(7) * (-t15 * t28 + t16 * t30 + t55 * t68) + m(6) * (t71 * t85 - t99); (t27 * t38 + t29 * t37) * mrSges(7,3) + t87 * t71 + t80 * t68 + m(6) * (pkin(4) * t71 + t68 * t85) + m(7) * (-t15 * t27 + t16 * t29 - t55 * t71); -0.2e1 * pkin(4) * t44 + 0.2e1 * t55 * t11 + t37 * t12 + t38 * t13 + t70 * t45 + t67 * t46 + Ifges(5,3) + m(7) * (t15 ^ 2 + t16 ^ 2 + t55 ^ 2) + m(6) * (t89 * pkin(8) ^ 2 + pkin(4) ^ 2) + 0.2e1 * (-t15 * t38 + t16 * t37) * mrSges(7,3) + 0.2e1 * pkin(8) * t84; -Ifges(6,6) * t93 + t9 * mrSges(6,1) - t10 * mrSges(6,2) + (m(7) * (t2 * t69 + t3 * t66) + t66 * t17 + t69 * t18) * pkin(5) + t76 + t105; (-t28 * t69 + t30 * t66) * t102 + t107; -t79 * t68 + (-t27 * t69 + t29 * t66) * t102 + t82; t58 + t59 - t79 * pkin(8) + (m(7) * (t15 * t69 + t16 * t66) + (t66 * t37 - t69 * t38) * mrSges(7,3)) * pkin(5) + t77; Ifges(6,3) + Ifges(7,3) + m(7) * (t66 ^ 2 + t69 ^ 2) * pkin(5) ^ 2 + 0.2e1 * t75; t76; -t7; t82; t77; Ifges(7,3) + t75; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
