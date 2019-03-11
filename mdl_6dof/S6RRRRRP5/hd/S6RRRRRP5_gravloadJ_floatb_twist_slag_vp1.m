% Calculate Gravitation load on the joints for
% S6RRRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 01:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRP5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP5_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:18:11
% EndTime: 2019-03-10 01:18:12
% DurationCPUTime: 0.84s
% Computational Cost: add. (650->176), mult. (680->235), div. (0->0), fcn. (657->10), ass. (0->70)
t54 = -pkin(9) - pkin(8);
t71 = rSges(5,3) - t54;
t46 = -pkin(10) + t54;
t70 = rSges(6,3) - t46;
t69 = rSges(7,3) + qJ(6) - t46;
t50 = sin(qJ(1));
t53 = cos(qJ(1));
t88 = g(1) * t53 + g(2) * t50;
t49 = sin(qJ(2));
t52 = cos(qJ(2));
t75 = rSges(4,3) + pkin(8);
t87 = t52 * pkin(2) + t75 * t49;
t47 = qJ(3) + qJ(4);
t41 = qJ(5) + t47;
t35 = sin(t41);
t36 = cos(t41);
t73 = t50 * t52;
t10 = t35 * t53 - t36 * t73;
t9 = t35 * t73 + t36 * t53;
t86 = -t9 * rSges(6,1) + t10 * rSges(6,2);
t85 = -t9 * rSges(7,1) + t10 * rSges(7,2);
t72 = t52 * t53;
t11 = -t35 * t72 + t36 * t50;
t12 = t35 * t50 + t36 * t72;
t84 = t11 * rSges(7,1) - t12 * rSges(7,2);
t83 = t11 * rSges(6,1) - t12 * rSges(6,2);
t48 = sin(qJ(3));
t82 = pkin(3) * t48;
t38 = sin(t47);
t81 = pkin(4) * t38;
t80 = pkin(5) * t35;
t77 = g(3) * t49;
t74 = rSges(3,2) * t49;
t39 = cos(t47);
t34 = pkin(4) * t39;
t51 = cos(qJ(3));
t43 = t51 * pkin(3);
t32 = t34 + t43;
t68 = t53 * pkin(1) + t50 * pkin(7);
t33 = pkin(5) * t36;
t23 = t33 + t32;
t28 = -t80 - t81;
t67 = rSges(3,1) * t52 - t74;
t65 = -rSges(5,1) * t38 - rSges(5,2) * t39;
t64 = -rSges(6,1) * t35 - rSges(6,2) * t36;
t63 = -rSges(7,1) * t35 - rSges(7,2) * t36;
t62 = rSges(4,1) * t51 - rSges(4,2) * t48 + pkin(2);
t20 = -t38 * t72 + t39 * t50;
t18 = t38 * t73 + t39 * t53;
t26 = -t48 * t72 + t50 * t51;
t24 = t48 * t73 + t51 * t53;
t37 = t43 + pkin(2);
t61 = rSges(5,1) * t39 - rSges(5,2) * t38 + t37;
t30 = pkin(2) + t32;
t60 = rSges(6,1) * t36 - rSges(6,2) * t35 + t30;
t17 = pkin(2) + t23;
t59 = rSges(7,1) * t36 - rSges(7,2) * t35 + t17;
t58 = t52 * t37 + t71 * t49;
t57 = t30 * t52 + t70 * t49;
t56 = t17 * t52 + t69 * t49;
t19 = t38 * t53 - t39 * t73;
t21 = t38 * t50 + t39 * t72;
t55 = g(1) * (t20 * rSges(5,1) - t21 * rSges(5,2)) + g(2) * (-t18 * rSges(5,1) + t19 * rSges(5,2));
t44 = t53 * pkin(7);
t31 = t81 + t82;
t29 = t33 + t34;
t27 = t48 * t50 + t51 * t72;
t25 = t48 * t53 - t51 * t73;
t22 = -t28 + t82;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t50 - rSges(2,2) * t53) + g(2) * (rSges(2,1) * t53 - rSges(2,2) * t50)) - m(3) * (g(1) * (rSges(3,3) * t53 + t44) + g(2) * (rSges(3,1) * t72 - t53 * t74 + t68) + (g(1) * (-pkin(1) - t67) + g(2) * rSges(3,3)) * t50) - m(4) * ((rSges(4,1) * t27 + rSges(4,2) * t26 + t87 * t53 + t68) * g(2) + (rSges(4,1) * t25 + rSges(4,2) * t24 + t44 + (-pkin(1) - t87) * t50) * g(1)) - m(5) * (g(1) * (t19 * rSges(5,1) + t18 * rSges(5,2) + t44) + g(2) * (t21 * rSges(5,1) + t20 * rSges(5,2) + t68) + (g(1) * t82 + g(2) * t58) * t53 + (g(1) * (-pkin(1) - t58) + g(2) * t82) * t50) - m(6) * (g(1) * (rSges(6,1) * t10 + rSges(6,2) * t9 + t44) + g(2) * (rSges(6,1) * t12 + rSges(6,2) * t11 + t68) + (g(1) * t31 + g(2) * t57) * t53 + (g(1) * (-pkin(1) - t57) + g(2) * t31) * t50) - m(7) * (g(1) * (rSges(7,1) * t10 + rSges(7,2) * t9 + t44) + g(2) * (rSges(7,1) * t12 + rSges(7,2) * t11 + t68) + (g(1) * t22 + g(2) * t56) * t53 + (g(1) * (-pkin(1) - t56) + g(2) * t22) * t50) -m(3) * (g(3) * t67 + t88 * (-rSges(3,1) * t49 - rSges(3,2) * t52)) - m(4) * ((g(3) * t62 + t88 * t75) * t52 + (g(3) * t75 - t88 * t62) * t49) - m(5) * ((g(3) * t61 + t88 * t71) * t52 + (g(3) * t71 - t88 * t61) * t49) - m(6) * ((g(3) * t60 + t88 * t70) * t52 + (g(3) * t70 - t88 * t60) * t49) - m(7) * ((g(3) * t59 + t88 * t69) * t52 + (g(3) * t69 - t88 * t59) * t49) -m(4) * (g(1) * (rSges(4,1) * t26 - rSges(4,2) * t27) + g(2) * (-rSges(4,1) * t24 + rSges(4,2) * t25) + (-rSges(4,1) * t48 - rSges(4,2) * t51) * t77) - m(5) * ((g(1) * t26 - g(2) * t24) * pkin(3) + (t65 - t82) * t77 + t55) - m(6) * (g(1) * (-t31 * t72 + t32 * t50 + t83) + g(2) * (-t31 * t73 - t32 * t53 + t86) + (-t31 + t64) * t77) - m(7) * (g(1) * (-t22 * t72 + t23 * t50 + t84) + g(2) * (-t22 * t73 - t23 * t53 + t85) + (-t22 + t63) * t77) -m(5) * t55 - m(6) * (g(1) * (t20 * pkin(4) + t83) + g(2) * (-t18 * pkin(4) + t86)) - m(7) * (g(1) * (t28 * t72 + t29 * t50 + t84) + g(2) * (t28 * t73 - t29 * t53 + t85)) + (-m(5) * t65 - m(6) * (t64 - t81) - m(7) * (t28 + t63)) * t77, -m(6) * (g(1) * t83 + g(2) * t86) - m(7) * (g(1) * (t11 * pkin(5) + t84) + g(2) * (-t9 * pkin(5) + t85)) + (-m(6) * t64 - m(7) * (t63 - t80)) * t77, -m(7) * (-g(3) * t52 + t88 * t49)];
taug  = t1(:);
