% Calculate Gravitation load on the joints for
% S5RRRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha4,d1,d2,d3,d4,d5]';
% m [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 17:33
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRR13_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR13_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR13_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR13_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR13_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR13_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 17:30:24
% EndTime: 2024-09-27 17:30:25
% DurationCPUTime: 0.31s
% Computational Cost: add. (585->91), mult. (362->125), div. (0->0), fcn. (326->16), ass. (0->59)
t85 = rSges(5,3) + pkin(9);
t50 = qJ(1) + qJ(2);
t44 = sin(t50);
t84 = pkin(2) * t44;
t51 = sin(pkin(5));
t83 = g(3) * t51;
t54 = sin(qJ(1));
t82 = t54 * pkin(1);
t47 = qJ(3) + t50;
t39 = sin(t47);
t81 = t39 * t51;
t40 = cos(t47);
t80 = t40 * t51;
t52 = cos(pkin(5));
t53 = sin(qJ(4));
t79 = t52 * t53;
t55 = cos(qJ(4));
t78 = t52 * t55;
t49 = qJ(4) + qJ(5);
t42 = pkin(5) + t49;
t77 = sin(t42) / 0.2e1;
t46 = cos(t50);
t76 = t46 * rSges(3,1) - t44 * rSges(3,2);
t75 = t40 * rSges(4,1) - t39 * rSges(4,2);
t74 = pkin(5) - t49;
t38 = pkin(2) * t46;
t73 = t38 + t75;
t17 = -t39 * t78 - t40 * t53;
t18 = -t39 * t79 + t40 * t55;
t72 = t18 * rSges(5,1) + t17 * rSges(5,2) + t40 * pkin(3) + t85 * t81;
t71 = sin(t74);
t70 = -t44 * rSges(3,1) - t46 * rSges(3,2);
t69 = -t39 * rSges(4,1) - t40 * rSges(4,2);
t21 = t77 - t71 / 0.2e1;
t45 = cos(t49);
t68 = -t40 * t21 - t39 * t45;
t67 = t39 * t21 - t40 * t45;
t66 = t38 + t72;
t65 = -t39 * t53 + t40 * t78;
t33 = cos(t42) / 0.2e1;
t37 = cos(t74);
t22 = t37 / 0.2e1 + t33;
t43 = sin(t49);
t10 = -t39 * t22 - t40 * t43;
t23 = pkin(4) * t79 - t51 * (pkin(9) + pkin(10));
t41 = t55 * pkin(4) + pkin(3);
t64 = -t67 * rSges(6,1) + t10 * rSges(6,2) + rSges(6,3) * t81 - t39 * t23 + t40 * t41;
t16 = -t39 * t55 - t40 * t79;
t63 = t16 * rSges(5,1) - t65 * rSges(5,2) - t39 * pkin(3) + t85 * t80;
t62 = t38 + t64;
t61 = t69 - t84;
t9 = -t40 * t22 + t39 * t43;
t60 = t68 * rSges(6,1) + t9 * rSges(6,2) + rSges(6,3) * t80 - t40 * t23 - t39 * t41;
t59 = t63 - t84;
t58 = t60 - t84;
t57 = g(1) * (t10 * rSges(6,1) + t67 * rSges(6,2)) + g(2) * (-t9 * rSges(6,1) + t68 * rSges(6,2)) + g(3) * ((t77 + t71 / 0.2e1) * rSges(6,1) + (t33 - t37 / 0.2e1) * rSges(6,2));
t56 = cos(qJ(1));
t48 = t56 * pkin(1);
t1 = [-m(2) * (g(1) * (-t54 * rSges(2,1) - t56 * rSges(2,2)) + g(2) * (t56 * rSges(2,1) - t54 * rSges(2,2))) - m(3) * (g(1) * (t70 - t82) + g(2) * (t48 + t76)) - m(4) * (g(1) * (t61 - t82) + g(2) * (t48 + t73)) - m(5) * (g(1) * (t59 - t82) + g(2) * (t48 + t66)) - m(6) * (g(1) * (t58 - t82) + g(2) * (t48 + t62)), -m(3) * (g(1) * t70 + g(2) * t76) - m(4) * (g(1) * t61 + g(2) * t73) - m(5) * (g(1) * t59 + g(2) * t66) - m(6) * (g(1) * t58 + g(2) * t62), -m(4) * (g(1) * t69 + g(2) * t75) - m(5) * (g(1) * t63 + g(2) * t72) - m(6) * (g(1) * t60 + g(2) * t64), -m(5) * (g(1) * (t17 * rSges(5,1) - t18 * rSges(5,2)) + g(2) * (rSges(5,1) * t65 + t16 * rSges(5,2)) + (rSges(5,1) * t55 - rSges(5,2) * t53) * t83) - m(6) * ((g(1) * t17 + g(2) * t65 + t55 * t83) * pkin(4) + t57), -m(6) * t57];
taug = t1(:);
