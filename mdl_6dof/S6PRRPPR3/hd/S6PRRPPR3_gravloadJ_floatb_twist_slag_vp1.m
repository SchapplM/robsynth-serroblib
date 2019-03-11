% Calculate Gravitation load on the joints for
% S6PRRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
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
% Datum: 2019-03-08 21:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRPPR3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR3_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPPR3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:08:15
% EndTime: 2019-03-08 21:08:17
% DurationCPUTime: 0.67s
% Computational Cost: add. (393->134), mult. (958->189), div. (0->0), fcn. (1119->10), ass. (0->62)
t89 = rSges(7,3) + pkin(9);
t40 = sin(qJ(2));
t43 = cos(qJ(2));
t37 = cos(pkin(10));
t64 = cos(pkin(6));
t56 = t37 * t64;
t63 = sin(pkin(10));
t21 = -t63 * t40 + t43 * t56;
t50 = t64 * t63;
t23 = -t37 * t40 - t43 * t50;
t88 = g(1) * t23 + g(2) * t21;
t22 = t40 * t56 + t63 * t43;
t24 = t37 * t43 - t40 * t50;
t87 = g(1) * t24 + g(2) * t22;
t86 = -m(6) - m(7);
t39 = sin(qJ(3));
t36 = sin(pkin(6));
t42 = cos(qJ(3));
t72 = t36 * t42;
t7 = t22 * t39 + t37 * t72;
t4 = t7 * pkin(3);
t85 = -t7 * pkin(4) - t4;
t57 = t36 * t63;
t9 = t24 * t39 - t42 * t57;
t6 = t9 * pkin(3);
t84 = -t9 * pkin(4) - t6;
t79 = g(3) * t36;
t78 = rSges(5,2) + pkin(8);
t77 = rSges(4,3) + pkin(8);
t76 = rSges(6,1) * t39;
t75 = t21 * t42;
t74 = t23 * t42;
t73 = t36 * t40;
t71 = t36 * t43;
t25 = t39 * t73 - t64 * t42;
t20 = t25 * pkin(3);
t70 = -t25 * pkin(4) - t20;
t69 = pkin(2) * t71 + pkin(8) * t73;
t68 = qJ(4) * t39;
t67 = rSges(6,1) + qJ(4);
t66 = rSges(5,3) + qJ(4);
t65 = -rSges(6,3) - qJ(5);
t62 = -m(5) + t86;
t61 = t42 * t71;
t17 = t21 * pkin(2);
t59 = pkin(3) * t75 + t21 * t68 + t17;
t18 = t23 * pkin(2);
t58 = pkin(3) * t74 + t23 * t68 + t18;
t55 = pkin(4) * t75 + t59;
t54 = pkin(4) * t74 + t58;
t53 = pkin(3) * t61 + t68 * t71 + t69;
t52 = rSges(4,1) * t42 - rSges(4,2) * t39;
t51 = rSges(5,1) * t42 + rSges(5,3) * t39;
t38 = sin(qJ(6));
t41 = cos(qJ(6));
t49 = rSges(7,1) * t41 - rSges(7,2) * t38 + pkin(5);
t48 = g(3) * (pkin(4) * t61 + t53);
t47 = -t38 * rSges(7,1) - t41 * rSges(7,2) - qJ(5);
t26 = t64 * t39 + t40 * t72;
t10 = t24 * t42 + t39 * t57;
t8 = -t37 * t36 * t39 + t22 * t42;
t1 = [(-m(2) - m(3) - m(4) + t62) * g(3), -m(3) * (g(1) * (rSges(3,1) * t23 - rSges(3,2) * t24) + g(2) * (rSges(3,1) * t21 - rSges(3,2) * t22) + (rSges(3,1) * t43 - rSges(3,2) * t40) * t79) - m(4) * (g(1) * (t52 * t23 + t77 * t24 + t18) + g(2) * (t52 * t21 + t77 * t22 + t17) + g(3) * t69 + (rSges(4,3) * t40 + t52 * t43) * t79) - m(5) * (g(1) * (t51 * t23 + t78 * t24 + t58) + g(2) * (t51 * t21 + t78 * t22 + t59) + g(3) * t53 + (rSges(5,2) * t40 + t51 * t43) * t79) - m(6) * (g(1) * (-rSges(6,2) * t74 + t23 * t76 + t54) + g(2) * (-rSges(6,2) * t75 + t21 * t76 + t55) + t48 + ((-rSges(6,2) * t42 + t76) * t43 + t65 * t40) * t79 + t87 * (pkin(8) + t65)) - m(7) * (t47 * t40 * t79 + g(1) * t54 + g(2) * t55 + t48 + t87 * (pkin(8) + t47) + (t43 * t79 + t88) * (t49 * t39 + t42 * t89)) -m(4) * (g(1) * (-rSges(4,1) * t9 - rSges(4,2) * t10) + g(2) * (-rSges(4,1) * t7 - rSges(4,2) * t8) + g(3) * (-rSges(4,1) * t25 - rSges(4,2) * t26)) - m(5) * (g(1) * (-rSges(5,1) * t9 + t66 * t10 - t6) + g(2) * (-rSges(5,1) * t7 + t66 * t8 - t4) + g(3) * (-rSges(5,1) * t25 + t66 * t26 - t20)) - m(6) * (g(1) * (rSges(6,2) * t9 + t67 * t10 + t84) + g(2) * (rSges(6,2) * t7 + t67 * t8 + t85) + g(3) * (rSges(6,2) * t25 + t67 * t26 + t70)) + (-g(1) * (-t89 * t9 + t84) - g(2) * (-t7 * t89 + t85) - g(3) * (-t89 * t25 + t70) - (g(1) * t10 + g(2) * t8 + g(3) * t26) * (qJ(4) + t49)) * m(7), t62 * (g(1) * t9 + g(2) * t7 + g(3) * t25) t86 * (g(3) * t71 + t88) -m(7) * (g(1) * ((t23 * t41 - t38 * t9) * rSges(7,1) + (-t23 * t38 - t41 * t9) * rSges(7,2)) + g(2) * ((t21 * t41 - t38 * t7) * rSges(7,1) + (-t21 * t38 - t41 * t7) * rSges(7,2)) + g(3) * ((-t25 * t38 + t41 * t71) * rSges(7,1) + (-t25 * t41 - t38 * t71) * rSges(7,2)))];
taug  = t1(:);
