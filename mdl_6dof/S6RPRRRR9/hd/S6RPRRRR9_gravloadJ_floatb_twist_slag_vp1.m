% Calculate Gravitation load on the joints for
% S6RPRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
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
% Datum: 2019-03-09 07:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRR9_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR9_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR9_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR9_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR9_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR9_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:23:14
% EndTime: 2019-03-09 07:23:15
% DurationCPUTime: 0.77s
% Computational Cost: add. (395->139), mult. (492->186), div. (0->0), fcn. (465->10), ass. (0->63)
t40 = sin(qJ(3));
t43 = cos(qJ(3));
t66 = rSges(5,3) + pkin(8);
t80 = t66 * t43;
t83 = t40 * pkin(3) - t80;
t41 = sin(qJ(1));
t44 = cos(qJ(1));
t79 = g(1) * t41 - g(2) * t44;
t45 = -pkin(9) - pkin(8);
t58 = rSges(7,3) + pkin(10) - t45;
t82 = t58 * t43;
t59 = rSges(6,3) - t45;
t81 = t59 * t43;
t75 = -pkin(1) - pkin(7);
t38 = qJ(4) + qJ(5);
t33 = qJ(6) + t38;
t26 = sin(t33);
t27 = cos(t33);
t63 = t40 * t41;
t5 = -t26 * t63 + t27 * t44;
t6 = t26 * t44 + t27 * t63;
t74 = t5 * rSges(7,1) - t6 * rSges(7,2);
t62 = t40 * t44;
t7 = t26 * t62 + t27 * t41;
t8 = -t26 * t41 + t27 * t62;
t73 = t7 * rSges(7,1) + t8 * rSges(7,2);
t39 = sin(qJ(4));
t72 = pkin(4) * t39;
t29 = sin(t38);
t71 = pkin(5) * t29;
t68 = g(3) * t43;
t30 = cos(t38);
t13 = -t29 * t63 + t30 * t44;
t14 = t29 * t44 + t30 * t63;
t65 = t13 * rSges(6,1) - t14 * rSges(6,2);
t64 = rSges(4,2) * t43;
t42 = cos(qJ(4));
t61 = t41 * t42;
t60 = t42 * t44;
t15 = t29 * t62 + t30 * t41;
t16 = -t29 * t41 + t30 * t62;
t57 = t15 * rSges(6,1) + t16 * rSges(6,2);
t34 = t42 * pkin(4);
t23 = pkin(5) * t30 + t34;
t56 = t44 * pkin(1) + t41 * qJ(2);
t55 = t44 * pkin(7) + t56;
t53 = rSges(4,1) * t40 + t64;
t52 = -rSges(6,1) * t29 - rSges(6,2) * t30;
t51 = -rSges(7,1) * t26 - rSges(7,2) * t27;
t50 = rSges(5,1) * t42 - rSges(5,2) * t39 + pkin(3);
t19 = t39 * t62 + t61;
t17 = -t39 * t63 + t60;
t28 = t34 + pkin(3);
t49 = rSges(6,1) * t30 - rSges(6,2) * t29 + t28;
t21 = pkin(3) + t23;
t48 = rSges(7,1) * t27 - rSges(7,2) * t26 + t21;
t47 = t40 * t28 - t81;
t46 = t40 * t21 - t82;
t32 = t44 * qJ(2);
t22 = t71 + t72;
t20 = -t39 * t41 + t40 * t60;
t18 = t39 * t44 + t40 * t61;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t41 - rSges(2,2) * t44) + g(2) * (rSges(2,1) * t44 - rSges(2,2) * t41)) - m(3) * (g(1) * (t44 * rSges(3,3) + t32 + (rSges(3,2) - pkin(1)) * t41) + g(2) * (-rSges(3,2) * t44 + rSges(3,3) * t41 + t56)) - m(4) * (g(1) * (rSges(4,1) * t62 + t44 * t64 + t32) + g(2) * (t44 * rSges(4,3) + t55) + (g(1) * (-rSges(4,3) + t75) + g(2) * t53) * t41) - m(5) * ((rSges(5,1) * t18 + rSges(5,2) * t17 + t83 * t41 + t55) * g(2) + (rSges(5,1) * t20 - rSges(5,2) * t19 + t75 * t41 + t83 * t44 + t32) * g(1)) - m(6) * (g(1) * (t16 * rSges(6,1) - t15 * rSges(6,2) + t32) + g(2) * (t14 * rSges(6,1) + t13 * rSges(6,2) + t55) + (g(1) * t47 + g(2) * t72) * t44 + (g(1) * (-t72 + t75) + g(2) * t47) * t41) - m(7) * (g(1) * (rSges(7,1) * t8 - rSges(7,2) * t7 + t32) + g(2) * (rSges(7,1) * t6 + rSges(7,2) * t5 + t55) + (g(1) * t46 + g(2) * t22) * t44 + (g(1) * (-t22 + t75) + g(2) * t46) * t41) (-m(3) - m(4) - m(5) - m(6) - m(7)) * t79, -m(4) * (-g(3) * t53 + t79 * (rSges(4,1) * t43 - rSges(4,2) * t40)) - m(5) * (g(3) * (-t50 * t40 + t80) + t79 * (t66 * t40 + t50 * t43)) - m(6) * (g(3) * (-t49 * t40 + t81) + t79 * (t59 * t40 + t49 * t43)) - m(7) * (g(3) * (-t48 * t40 + t82) + t79 * (t58 * t40 + t48 * t43)) -m(5) * (g(1) * (rSges(5,1) * t17 - rSges(5,2) * t18) + g(2) * (rSges(5,1) * t19 + rSges(5,2) * t20)) - m(6) * (g(1) * (t17 * pkin(4) + t65) + g(2) * (t19 * pkin(4) + t57)) - m(7) * (g(1) * (-t22 * t63 + t23 * t44 + t74) + g(2) * (t22 * t62 + t23 * t41 + t73)) + (-m(5) * (-rSges(5,1) * t39 - rSges(5,2) * t42) - m(6) * (t52 - t72) - m(7) * (-t22 + t51)) * t68, -m(6) * (g(1) * t65 + g(2) * t57) - m(7) * (g(1) * (t13 * pkin(5) + t74) + g(2) * (t15 * pkin(5) + t73)) + (-m(6) * t52 - m(7) * (t51 - t71)) * t68, -m(7) * (g(1) * t74 + g(2) * t73 + t51 * t68)];
taug  = t1(:);
