% Calculate Gravitation load on the joints for
% S6RRRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 19:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRR12_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR12_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR12_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR12_gravloadJ_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR12_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR12_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:38:12
% EndTime: 2019-03-09 19:38:16
% DurationCPUTime: 1.31s
% Computational Cost: add. (840->198), mult. (1575->283), div. (0->0), fcn. (1881->14), ass. (0->75)
t52 = -pkin(10) - qJ(4);
t85 = pkin(11) - t52 + rSges(7,3);
t84 = -t52 + rSges(6,3);
t83 = qJ(4) + rSges(5,3);
t54 = sin(qJ(2));
t82 = cos(pkin(6));
t91 = cos(qJ(1));
t71 = t82 * t91;
t89 = sin(qJ(1));
t90 = cos(qJ(2));
t25 = t54 * t71 + t89 * t90;
t53 = sin(qJ(3));
t55 = cos(qJ(3));
t50 = sin(pkin(6));
t79 = t50 * t91;
t12 = t25 * t53 + t55 * t79;
t70 = t82 * t89;
t27 = -t54 * t70 + t90 * t91;
t77 = t50 * t89;
t16 = t27 * t53 - t55 * t77;
t88 = t50 * t54;
t22 = t53 * t88 - t55 * t82;
t101 = g(1) * t16 + g(2) * t12 + g(3) * t22;
t13 = t25 * t55 - t53 * t79;
t17 = t27 * t55 + t53 * t77;
t23 = t53 * t82 + t55 * t88;
t100 = g(1) * t17 + g(2) * t13 + g(3) * t23;
t49 = sin(pkin(12));
t99 = pkin(4) * t49;
t94 = g(3) * t50;
t93 = rSges(4,3) + pkin(9);
t48 = pkin(12) + qJ(5);
t42 = sin(t48);
t29 = pkin(5) * t42 + t99;
t92 = pkin(9) + t29;
t51 = cos(pkin(12));
t41 = t51 * pkin(4) + pkin(3);
t87 = t91 * pkin(1) + pkin(8) * t77;
t78 = t50 * t90;
t86 = pkin(2) * t78 + pkin(9) * t88;
t81 = t27 * pkin(2) + t87;
t80 = pkin(9) + t99;
t76 = t53 * t90;
t75 = t55 * t90;
t74 = t50 * t76;
t73 = t50 * t75;
t72 = -pkin(1) * t89 + pkin(8) * t79;
t69 = -rSges(4,1) * t55 + rSges(4,2) * t53;
t24 = t54 * t89 - t71 * t90;
t43 = cos(t48);
t68 = -t13 * t42 + t24 * t43;
t26 = t54 * t91 + t70 * t90;
t8 = -t17 * t42 + t26 * t43;
t67 = -t25 * pkin(2) + t72;
t66 = -rSges(5,1) * t51 + rSges(5,2) * t49 - pkin(3);
t65 = rSges(5,1) * t49 + rSges(5,2) * t51 + pkin(9);
t64 = rSges(6,1) * t43 - rSges(6,2) * t42 + t41;
t28 = pkin(5) * t43 + t41;
t44 = qJ(6) + t48;
t39 = sin(t44);
t40 = cos(t44);
t63 = rSges(7,1) * t40 - rSges(7,2) * t39 + t28;
t62 = rSges(7,1) * t39 + rSges(7,2) * t40 + t92;
t61 = -t23 * t42 - t43 * t78;
t60 = rSges(6,1) * t42 + rSges(6,2) * t43 + t80;
t59 = -t53 * t83 + t55 * t66;
t58 = -t53 * t84 - t55 * t64;
t57 = -t53 * t85 - t55 * t63;
t6 = -t17 * t39 + t26 * t40;
t7 = t17 * t40 + t26 * t39;
t56 = m(7) * (g(1) * (t6 * rSges(7,1) - t7 * rSges(7,2)) + g(2) * ((-t13 * t39 + t24 * t40) * rSges(7,1) + (-t13 * t40 - t24 * t39) * rSges(7,2)) + g(3) * ((-t23 * t39 - t40 * t78) * rSges(7,1) + (-t23 * t40 + t39 * t78) * rSges(7,2)));
t20 = t26 * pkin(2);
t18 = t24 * pkin(2);
t9 = t17 * t43 + t26 * t42;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t89 - rSges(2,2) * t91) + g(2) * (rSges(2,1) * t91 - rSges(2,2) * t89)) - m(3) * (g(1) * (-t25 * rSges(3,1) + t24 * rSges(3,2) + rSges(3,3) * t79 + t72) + g(2) * (t27 * rSges(3,1) - t26 * rSges(3,2) + rSges(3,3) * t77 + t87)) - m(4) * (g(1) * (-rSges(4,1) * t13 + rSges(4,2) * t12 - t24 * t93 + t67) + g(2) * (rSges(4,1) * t17 - rSges(4,2) * t16 + t26 * t93 + t81)) - m(5) * (g(1) * (-t13 * pkin(3) - t24 * pkin(9) + (-t13 * t51 - t24 * t49) * rSges(5,1) + (t13 * t49 - t24 * t51) * rSges(5,2) - t83 * t12 + t67) + g(2) * (t26 * pkin(9) + t17 * pkin(3) + (t17 * t51 + t26 * t49) * rSges(5,1) + (-t17 * t49 + t26 * t51) * rSges(5,2) + t83 * t16 + t81)) - m(6) * (g(1) * (-t12 * t84 - t13 * t64 - t24 * t60 + t67) + g(2) * (rSges(6,1) * t9 + rSges(6,2) * t8 + t16 * t84 + t17 * t41 + t26 * t80 + t81)) - m(7) * (g(1) * (-t12 * t85 - t13 * t63 - t24 * t62 + t67) + g(2) * (rSges(7,1) * t7 + rSges(7,2) * t6 + t16 * t85 + t17 * t28 + t26 * t92 + t81)) -m(3) * (g(1) * (-rSges(3,1) * t26 - rSges(3,2) * t27) + g(2) * (-rSges(3,1) * t24 - rSges(3,2) * t25) + (rSges(3,1) * t90 - rSges(3,2) * t54) * t94) - m(4) * (g(1) * (t26 * t69 + t27 * t93 - t20) + g(2) * (t24 * t69 + t25 * t93 - t18) + g(3) * t86 + (rSges(4,1) * t75 - rSges(4,2) * t76 + rSges(4,3) * t54) * t94) - m(5) * (g(1) * (t26 * t59 + t27 * t65 - t20) + g(2) * (t24 * t59 + t25 * t65 - t18) + g(3) * (pkin(3) * t73 + t86 + t83 * t74 + ((t49 * t54 + t51 * t75) * rSges(5,1) + (-t49 * t75 + t51 * t54) * rSges(5,2)) * t50)) - m(6) * (g(1) * (t26 * t58 + t27 * t60 - t20) + g(2) * (t24 * t58 + t25 * t60 - t18) + g(3) * (t41 * t73 + t88 * t99 + t86 + t84 * t74 + ((t42 * t54 + t43 * t75) * rSges(6,1) + (-t42 * t75 + t43 * t54) * rSges(6,2)) * t50)) - m(7) * (g(1) * (t26 * t57 + t27 * t62 - t20) + g(2) * (t24 * t57 + t25 * t62 - t18) + g(3) * (t28 * t73 + t29 * t88 + t86 + t85 * t74 + ((t39 * t54 + t40 * t75) * rSges(7,1) + (-t39 * t75 + t40 * t54) * rSges(7,2)) * t50)) -m(4) * (g(1) * (-rSges(4,1) * t16 - rSges(4,2) * t17) + g(2) * (-rSges(4,1) * t12 - rSges(4,2) * t13) + g(3) * (-rSges(4,1) * t22 - rSges(4,2) * t23)) - m(5) * (t100 * t83 + t101 * t66) - m(6) * (t100 * t84 - t101 * t64) - m(7) * (t100 * t85 - t101 * t63) (-m(5) - m(6) - m(7)) * t101, -m(6) * (g(1) * (rSges(6,1) * t8 - rSges(6,2) * t9) + g(2) * (t68 * rSges(6,1) + (-t13 * t43 - t24 * t42) * rSges(6,2)) + g(3) * (t61 * rSges(6,1) + (-t23 * t43 + t42 * t78) * rSges(6,2))) - t56 - m(7) * (g(1) * t8 + g(2) * t68 + g(3) * t61) * pkin(5), -t56];
taug  = t1(:);
