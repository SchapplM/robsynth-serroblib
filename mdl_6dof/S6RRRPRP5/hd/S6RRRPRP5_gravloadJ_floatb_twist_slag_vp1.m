% Calculate Gravitation load on the joints for
% S6RRRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-03-09 16:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRP5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP5_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:47:46
% EndTime: 2019-03-09 16:47:48
% DurationCPUTime: 0.89s
% Computational Cost: add. (602->160), mult. (650->219), div. (0->0), fcn. (638->10), ass. (0->74)
t91 = rSges(7,1) + pkin(5);
t78 = rSges(7,3) + qJ(6);
t53 = -qJ(4) - pkin(8);
t82 = rSges(5,3) - t53;
t59 = cos(qJ(1));
t56 = sin(qJ(1));
t95 = g(2) * t56;
t103 = g(1) * t59 + t95;
t55 = sin(qJ(2));
t58 = cos(qJ(2));
t90 = rSges(4,3) + pkin(8);
t102 = t58 * pkin(2) + t90 * t55;
t52 = qJ(3) + pkin(10);
t46 = qJ(5) + t52;
t41 = sin(t46);
t42 = cos(t46);
t101 = t78 * t41 + t91 * t42;
t86 = t56 * t58;
t11 = t41 * t86 + t42 * t59;
t84 = t59 * t41;
t12 = t42 * t86 - t84;
t100 = -t11 * rSges(6,1) - t12 * rSges(6,2);
t13 = -t56 * t42 + t58 * t84;
t85 = t58 * t59;
t14 = t56 * t41 + t42 * t85;
t99 = -t13 * rSges(6,1) - t14 * rSges(6,2);
t54 = sin(qJ(3));
t98 = pkin(3) * t54;
t97 = g(1) * t56;
t45 = cos(t52);
t57 = cos(qJ(3));
t48 = t57 * pkin(3);
t34 = pkin(4) * t45 + t48;
t32 = pkin(2) + t34;
t23 = t58 * t32;
t94 = g(3) * t23;
t93 = g(3) * t55;
t89 = rSges(3,2) * t55;
t88 = t41 * t55;
t44 = sin(t52);
t33 = pkin(4) * t44 + t98;
t31 = t59 * t33;
t51 = -pkin(9) + t53;
t83 = rSges(7,2) - t51;
t81 = rSges(6,3) - t51;
t80 = -t58 * t31 + t56 * t34;
t79 = t59 * pkin(1) + t56 * pkin(7);
t49 = t59 * pkin(7);
t77 = t56 * t55 * t51 + t31 + t49;
t76 = -pkin(1) - t23;
t75 = t83 * t59;
t74 = t81 * t59;
t73 = -t33 * t86 - t34 * t59;
t72 = t32 * t85 + t56 * t33 + t79;
t71 = rSges(3,1) * t58 - t89;
t69 = rSges(6,1) * t42 - rSges(6,2) * t41;
t68 = -t91 * t11 + t78 * t12;
t67 = -t91 * t13 + t78 * t14;
t66 = rSges(4,1) * t57 - rSges(4,2) * t54 + pkin(2);
t26 = -t54 * t85 + t56 * t57;
t24 = t54 * t86 + t57 * t59;
t43 = t48 + pkin(2);
t65 = rSges(5,1) * t45 - rSges(5,2) * t44 + t43;
t63 = (-rSges(6,1) * t41 - rSges(6,2) * t42) * t55;
t62 = t78 * t42 * t55 - t91 * t88;
t61 = t43 * t58 + t82 * t55;
t28 = t55 * t33;
t27 = t56 * t54 + t57 * t85;
t25 = t54 * t59 - t57 * t86;
t18 = t56 * t44 + t45 * t85;
t17 = -t44 * t85 + t56 * t45;
t16 = t44 * t59 - t45 * t86;
t15 = t44 * t86 + t45 * t59;
t1 = [-m(2) * (g(1) * (-t56 * rSges(2,1) - rSges(2,2) * t59) + g(2) * (rSges(2,1) * t59 - t56 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) * t59 + t49) + g(2) * (rSges(3,1) * t85 - t59 * t89 + t79) + (g(1) * (-pkin(1) - t71) + g(2) * rSges(3,3)) * t56) - m(4) * (g(1) * (rSges(4,1) * t25 + rSges(4,2) * t24 + t49) + (-pkin(1) - t102) * t97 + (t27 * rSges(4,1) + t26 * rSges(4,2) + t102 * t59 + t79) * g(2)) - m(5) * (g(1) * (t16 * rSges(5,1) + t15 * rSges(5,2) + t49) + g(2) * (t18 * rSges(5,1) + t17 * rSges(5,2) + t79) + (g(1) * t98 + g(2) * t61) * t59 + (g(1) * (-pkin(1) - t61) + g(2) * t98) * t56) - m(6) * (g(1) * (-rSges(6,1) * t12 + rSges(6,2) * t11 + t77) + g(2) * (t14 * rSges(6,1) - t13 * rSges(6,2) + t55 * t74 + t72) + (-rSges(6,3) * t55 + t76) * t97) - m(7) * (g(1) * (-t78 * t11 - t91 * t12 + t77) + g(2) * (t78 * t13 + t91 * t14 + t55 * t75 + t72) + (-rSges(7,2) * t55 + t76) * t97) -m(3) * (g(3) * t71 + t103 * (-rSges(3,1) * t55 - rSges(3,2) * t58)) - m(4) * ((g(3) * t66 + t103 * t90) * t58 + (g(3) * t90 - t103 * t66) * t55) - m(5) * ((g(3) * t65 + t103 * t82) * t58 + (g(3) * t82 - t103 * t65) * t55) - m(6) * (t94 + (g(1) * t74 + g(3) * t69 + t81 * t95) * t58 + (g(3) * t81 + t103 * (-t32 - t69)) * t55) - m(7) * (t94 + (g(1) * t75 + g(3) * t101 + t83 * t95) * t58 + (g(3) * t83 + t103 * (-t101 - t32)) * t55) -m(4) * (g(1) * (rSges(4,1) * t26 - rSges(4,2) * t27) + g(2) * (-rSges(4,1) * t24 + rSges(4,2) * t25) + (-rSges(4,1) * t54 - rSges(4,2) * t57) * t93) - m(5) * (g(1) * (t17 * rSges(5,1) - t18 * rSges(5,2)) + g(2) * (-t15 * rSges(5,1) + t16 * rSges(5,2)) + (g(1) * t26 - g(2) * t24) * pkin(3) + (-rSges(5,1) * t44 - rSges(5,2) * t45 - t98) * t93) - m(6) * (g(1) * (t80 + t99) + g(2) * (t73 + t100) + g(3) * (-t28 + t63)) - m(7) * (g(1) * (t67 + t80) + g(2) * (t68 + t73) + g(3) * (-t28 + t62)) (-m(5) - m(6) - m(7)) * (-g(3) * t58 + t103 * t55) -m(6) * (g(1) * t99 + g(2) * t100 + g(3) * t63) - m(7) * (g(1) * t67 + g(2) * t68 + g(3) * t62) -m(7) * (g(1) * t13 + g(2) * t11 + g(3) * t88)];
taug  = t1(:);
