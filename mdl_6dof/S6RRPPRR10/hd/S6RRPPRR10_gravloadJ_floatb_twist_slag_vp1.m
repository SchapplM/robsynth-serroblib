% Calculate Gravitation load on the joints for
% S6RRPPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
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
% Datum: 2019-03-09 09:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPPRR10_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR10_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR10_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR10_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR10_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR10_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:34:11
% EndTime: 2019-03-09 09:34:13
% DurationCPUTime: 0.77s
% Computational Cost: add. (357->142), mult. (495->194), div. (0->0), fcn. (462->10), ass. (0->67)
t41 = cos(qJ(2));
t79 = g(3) * t41;
t38 = -pkin(8) - qJ(4);
t67 = rSges(6,3) - t38;
t66 = rSges(7,3) + pkin(9) - t38;
t42 = cos(qJ(1));
t40 = sin(qJ(1));
t80 = g(2) * t40;
t51 = g(1) * t42 + t80;
t62 = rSges(5,3) + qJ(4);
t36 = sin(pkin(10));
t82 = pkin(4) * t36;
t31 = t41 * pkin(2);
t37 = cos(pkin(10));
t24 = t37 * pkin(4) + pkin(3);
t77 = rSges(4,2) * t41;
t76 = rSges(5,2) * t37;
t35 = pkin(10) + qJ(5);
t25 = sin(t35);
t16 = pkin(5) * t25 + t82;
t39 = sin(qJ(2));
t75 = t16 * t39;
t74 = t36 * t39;
t73 = t39 * t42;
t27 = qJ(6) + t35;
t22 = sin(t27);
t72 = t40 * t22;
t23 = cos(t27);
t71 = t40 * t23;
t70 = t40 * t25;
t26 = cos(t35);
t69 = t40 * t26;
t68 = t41 * t42;
t28 = t39 * qJ(3);
t65 = t28 + t31;
t64 = t42 * pkin(1) + t40 * pkin(7);
t63 = qJ(3) * t41;
t61 = -m(5) - m(6) - m(7);
t60 = pkin(4) * t74;
t59 = -pkin(2) - t67;
t58 = -pkin(2) - t66;
t57 = -pkin(2) - t62;
t56 = -pkin(1) - t28;
t55 = pkin(2) * t68 + t42 * t28 + t64;
t54 = g(1) * t59;
t53 = g(1) * t58;
t52 = g(1) * t57;
t50 = rSges(3,1) * t41 - rSges(3,2) * t39;
t48 = rSges(5,1) * t36 + t76;
t47 = t37 * rSges(5,1) - t36 * rSges(5,2) + pkin(3);
t9 = t26 * t73 - t70;
t11 = t25 * t42 + t39 * t69;
t46 = rSges(7,1) * t22 + rSges(7,2) * t23 + t16;
t45 = rSges(6,1) * t25 + rSges(6,2) * t26 + t82;
t18 = t40 * t63;
t20 = t42 * t63;
t44 = g(1) * t20 + g(2) * t18 + g(3) * t65;
t5 = t23 * t73 - t72;
t6 = t22 * t73 + t71;
t7 = t22 * t42 + t39 * t71;
t8 = t23 * t42 - t39 * t72;
t43 = g(1) * (t5 * rSges(7,1) - t6 * rSges(7,2)) + g(2) * (t7 * rSges(7,1) + t8 * rSges(7,2)) + (-rSges(7,1) * t23 + rSges(7,2) * t22) * t79;
t32 = t42 * pkin(7);
t15 = pkin(5) * t26 + t24;
t12 = t26 * t42 - t39 * t70;
t10 = t25 * t73 + t69;
t1 = [-m(2) * (g(1) * (-t40 * rSges(2,1) - rSges(2,2) * t42) + g(2) * (rSges(2,1) * t42 - t40 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) * t42 + t32) + g(2) * (rSges(3,1) * t68 - rSges(3,2) * t73 + t64) + (g(1) * (-pkin(1) - t50) + g(2) * rSges(3,3)) * t40) - m(4) * (g(1) * (rSges(4,1) * t42 + t32) + g(2) * (-rSges(4,2) * t68 + rSges(4,3) * t73 + t55) + (g(1) * (-rSges(4,3) * t39 - t31 + t56 + t77) + g(2) * rSges(4,1)) * t40) - m(5) * (g(1) * t32 + g(2) * t55 + (g(1) * t47 + g(2) * (rSges(5,1) * t74 + t39 * t76 + t62 * t41)) * t42 + (g(2) * t47 + t41 * t52 + (-pkin(1) + (-qJ(3) - t48) * t39) * g(1)) * t40) - m(6) * (g(1) * (t12 * rSges(6,1) - t11 * rSges(6,2) + t32) + g(2) * (t10 * rSges(6,1) + t9 * rSges(6,2) + t55) + (g(1) * t24 + g(2) * (t67 * t41 + t60)) * t42 + (g(1) * (t56 - t60) + g(2) * t24 + t41 * t54) * t40) - m(7) * (g(1) * (t8 * rSges(7,1) - t7 * rSges(7,2) + t32) + g(2) * (t6 * rSges(7,1) + t5 * rSges(7,2) + t55) + (g(1) * t15 + g(2) * (t66 * t41 + t75)) * t42 + (g(1) * (t56 - t75) + g(2) * t15 + t41 * t53) * t40) -m(3) * (g(3) * t50 + t51 * (-rSges(3,1) * t39 - rSges(3,2) * t41)) - m(4) * (g(1) * (rSges(4,3) * t68 + t20) + g(2) * (rSges(4,3) * t40 * t41 + t18) + g(3) * (t65 - t77) + (g(3) * rSges(4,3) + t51 * (rSges(4,2) - pkin(2))) * t39) - m(5) * ((g(3) * t62 + t51 * t48) * t41 + (g(3) * t48 + t42 * t52 + t57 * t80) * t39 + t44) - m(6) * ((g(3) * t67 + t51 * t45) * t41 + (g(3) * t45 + t42 * t54 + t59 * t80) * t39 + t44) - m(7) * ((g(3) * t66 + t51 * t46) * t41 + (g(3) * t46 + t42 * t53 + t58 * t80) * t39 + t44) (-m(4) + t61) * (t39 * t51 - t79) t61 * (g(3) * t39 + t41 * t51) -m(6) * (g(1) * (rSges(6,1) * t9 - rSges(6,2) * t10) + g(2) * (rSges(6,1) * t11 + rSges(6,2) * t12) + (-rSges(6,1) * t26 + rSges(6,2) * t25) * t79) - m(7) * ((g(1) * t9 + g(2) * t11 - t26 * t79) * pkin(5) + t43) -m(7) * t43];
taug  = t1(:);
