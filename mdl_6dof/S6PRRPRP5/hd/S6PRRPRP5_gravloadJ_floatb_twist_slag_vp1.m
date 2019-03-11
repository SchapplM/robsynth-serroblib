% Calculate Gravitation load on the joints for
% S6PRRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRPRP5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP5_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRP5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:45:43
% EndTime: 2019-03-08 21:45:45
% DurationCPUTime: 0.79s
% Computational Cost: add. (476->147), mult. (1181->215), div. (0->0), fcn. (1410->10), ass. (0->73)
t102 = pkin(4) + pkin(8);
t56 = sin(qJ(2));
t59 = cos(qJ(2));
t81 = cos(pkin(10));
t82 = cos(pkin(6));
t66 = t82 * t81;
t80 = sin(pkin(10));
t36 = t56 * t66 + t80 * t59;
t55 = sin(qJ(3));
t58 = cos(qJ(3));
t53 = sin(pkin(6));
t71 = t53 * t81;
t15 = t36 * t58 - t55 * t71;
t65 = t82 * t80;
t38 = -t56 * t65 + t81 * t59;
t70 = t53 * t80;
t17 = t38 * t58 + t55 * t70;
t91 = t53 * t56;
t40 = t82 * t55 + t58 * t91;
t101 = g(1) * t17 + g(2) * t15 + g(3) * t40;
t97 = g(3) * t53;
t96 = rSges(5,1) + pkin(8);
t95 = rSges(7,1) + pkin(5);
t94 = rSges(4,3) + pkin(8);
t35 = t80 * t56 - t59 * t66;
t93 = t35 * t58;
t37 = t81 * t56 + t59 * t65;
t92 = t37 * t58;
t90 = t53 * t59;
t54 = sin(qJ(5));
t89 = t54 * t55;
t88 = t54 * t59;
t57 = cos(qJ(5));
t87 = t55 * t57;
t86 = pkin(2) * t90 + pkin(8) * t91;
t85 = qJ(4) * t55;
t84 = rSges(5,3) + qJ(4);
t83 = rSges(7,3) + qJ(6);
t79 = -m(5) - m(6) - m(7);
t78 = t57 * t90;
t77 = t58 * t90;
t30 = t35 * pkin(2);
t76 = -pkin(3) * t93 - t35 * t85 - t30;
t32 = t37 * pkin(2);
t75 = -pkin(3) * t92 - t37 * t85 - t32;
t14 = t36 * t55 + t58 * t71;
t12 = t14 * pkin(3);
t74 = -pkin(9) * t14 - t12;
t16 = t38 * t55 - t58 * t70;
t13 = t16 * pkin(3);
t73 = -pkin(9) * t16 - t13;
t39 = t55 * t91 - t82 * t58;
t34 = t39 * pkin(3);
t72 = -pkin(9) * t39 - t34;
t69 = pkin(3) * t77 + t85 * t90 + t86;
t68 = rSges(4,1) * t58 - rSges(4,2) * t55;
t67 = rSges(5,2) * t58 - rSges(5,3) * t55;
t64 = pkin(4) * t91 + pkin(9) * t77 + t69;
t62 = -pkin(9) * t93 + t102 * t36 + t76;
t61 = -pkin(9) * t92 + t102 * t38 + t75;
t23 = (t55 * t88 + t56 * t57) * t53;
t22 = t54 * t91 - t55 * t78;
t19 = -t39 * t54 + t78;
t18 = t39 * t57 + t53 * t88;
t9 = -t37 * t89 + t38 * t57;
t8 = t37 * t87 + t38 * t54;
t7 = -t35 * t89 + t36 * t57;
t6 = t35 * t87 + t36 * t54;
t5 = t16 * t54 + t37 * t57;
t4 = -t16 * t57 + t37 * t54;
t3 = t14 * t54 + t35 * t57;
t2 = -t14 * t57 + t35 * t54;
t1 = [(-m(2) - m(3) - m(4) + t79) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t37 - rSges(3,2) * t38) + g(2) * (-rSges(3,1) * t35 - rSges(3,2) * t36) + (rSges(3,1) * t59 - rSges(3,2) * t56) * t97) - m(4) * (g(1) * (-t68 * t37 + t94 * t38 - t32) + g(2) * (-t68 * t35 + t94 * t36 - t30) + g(3) * t86 + (rSges(4,3) * t56 + t68 * t59) * t97) - m(5) * (g(1) * (t67 * t37 + t96 * t38 + t75) + g(2) * (t67 * t35 + t96 * t36 + t76) + g(3) * t69 + (rSges(5,1) * t56 - t67 * t59) * t97) - m(6) * (g(1) * (rSges(6,1) * t9 - rSges(6,2) * t8 - rSges(6,3) * t92 + t61) + g(2) * (rSges(6,1) * t7 - rSges(6,2) * t6 - rSges(6,3) * t93 + t62) + g(3) * (t23 * rSges(6,1) - t22 * rSges(6,2) + rSges(6,3) * t77 + t64)) - m(7) * (g(1) * (-rSges(7,2) * t92 + t83 * t8 + t95 * t9 + t61) + g(2) * (-rSges(7,2) * t93 + t83 * t6 + t95 * t7 + t62) + g(3) * (rSges(7,2) * t77 + t83 * t22 + t95 * t23 + t64)) -m(4) * (g(1) * (-rSges(4,1) * t16 - rSges(4,2) * t17) + g(2) * (-rSges(4,1) * t14 - rSges(4,2) * t15) + g(3) * (-rSges(4,1) * t39 - rSges(4,2) * t40)) - m(5) * (g(1) * (t16 * rSges(5,2) + t84 * t17 - t13) + g(2) * (t14 * rSges(5,2) + t84 * t15 - t12) + g(3) * (t39 * rSges(5,2) + t84 * t40 - t34)) + (-g(1) * (-rSges(7,2) * t16 + t73) - g(2) * (-rSges(7,2) * t14 + t74) - g(3) * (-rSges(7,2) * t39 + t72) - t101 * (t95 * t54 - t83 * t57 + qJ(4))) * m(7) + (-g(1) * (-rSges(6,3) * t16 + t73) - g(2) * (-rSges(6,3) * t14 + t74) - g(3) * (-rSges(6,3) * t39 + t72) - t101 * (rSges(6,1) * t54 + rSges(6,2) * t57 + qJ(4))) * m(6), t79 * (g(1) * t16 + g(2) * t14 + g(3) * t39) -m(6) * (g(1) * (-rSges(6,1) * t4 - rSges(6,2) * t5) + g(2) * (-rSges(6,1) * t2 - rSges(6,2) * t3) + g(3) * (rSges(6,1) * t18 + rSges(6,2) * t19)) - m(7) * (g(1) * (-t95 * t4 + t83 * t5) + g(2) * (-t95 * t2 + t83 * t3) + g(3) * (t95 * t18 - t83 * t19)) -m(7) * (g(1) * t4 + g(2) * t2 - g(3) * t18)];
taug  = t1(:);
