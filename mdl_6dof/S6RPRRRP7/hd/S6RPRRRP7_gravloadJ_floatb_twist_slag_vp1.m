% Calculate Gravitation load on the joints for
% S6RPRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 06:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRP7_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP7_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP7_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP7_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP7_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP7_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:18:44
% EndTime: 2019-03-09 06:18:46
% DurationCPUTime: 0.76s
% Computational Cost: add. (526->134), mult. (533->178), div. (0->0), fcn. (519->10), ass. (0->66)
t83 = rSges(7,1) + pkin(5);
t37 = pkin(10) + qJ(3);
t34 = cos(t37);
t43 = sin(qJ(1));
t44 = cos(qJ(4));
t72 = t43 * t44;
t42 = sin(qJ(4));
t45 = cos(qJ(1));
t74 = t42 * t45;
t17 = -t34 * t74 + t72;
t64 = rSges(7,3) + qJ(6);
t82 = rSges(5,3) + pkin(8);
t93 = g(1) * t45 + g(2) * t43;
t38 = qJ(4) + qJ(5);
t35 = sin(t38);
t36 = cos(t38);
t92 = t64 * t35 + t83 * t36;
t76 = t36 * t45;
t77 = t35 * t43;
t11 = t34 * t77 + t76;
t70 = t45 * t35;
t73 = t43 * t36;
t12 = t34 * t73 - t70;
t91 = -t11 * rSges(6,1) - t12 * rSges(6,2);
t13 = t34 * t70 - t73;
t14 = t34 * t76 + t77;
t90 = -t13 * rSges(6,1) - t14 * rSges(6,2);
t89 = pkin(4) * t42;
t41 = -pkin(7) - qJ(2);
t87 = g(2) * t41;
t32 = pkin(4) * t44 + pkin(3);
t21 = t34 * t32;
t85 = g(3) * t21;
t33 = sin(t37);
t84 = g(3) * t33;
t81 = t33 * t35;
t79 = t33 * t45;
t46 = -pkin(9) - pkin(8);
t78 = t33 * t46;
t75 = t42 * t43;
t71 = t44 * t45;
t69 = rSges(7,2) - t46;
t68 = rSges(4,3) - t41;
t67 = rSges(6,3) - t46;
t66 = t64 * t33 * t36;
t65 = rSges(3,3) + qJ(2);
t40 = cos(pkin(10));
t31 = pkin(2) * t40 + pkin(1);
t62 = -t31 - t21;
t61 = pkin(4) * t74 - t45 * t41 + t43 * t78;
t60 = rSges(4,1) * t34 - rSges(4,2) * t33;
t58 = rSges(6,1) * t36 - rSges(6,2) * t35;
t57 = -rSges(6,1) * t35 - rSges(6,2) * t36;
t56 = -t83 * t11 + t64 * t12;
t55 = t17 * pkin(4);
t54 = -t83 * t13 + t64 * t14;
t53 = rSges(3,1) * t40 - rSges(3,2) * sin(pkin(10)) + pkin(1);
t52 = rSges(5,1) * t44 - rSges(5,2) * t42 + pkin(3);
t15 = t34 * t75 + t71;
t24 = t45 * t31;
t50 = pkin(4) * t75 + t24 + (t21 - t78) * t45;
t49 = pkin(3) * t34 + t33 * t82;
t48 = t15 * pkin(4);
t18 = t34 * t71 + t75;
t16 = -t34 * t72 + t74;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t43 - rSges(2,2) * t45) + g(2) * (rSges(2,1) * t45 - rSges(2,2) * t43)) - m(3) * ((g(1) * t65 + g(2) * t53) * t45 + (-g(1) * t53 + g(2) * t65) * t43) - m(4) * (g(2) * t24 + (g(1) * t68 + g(2) * t60) * t45 + (g(1) * (-t31 - t60) + g(2) * t68) * t43) - m(5) * (g(1) * (t16 * rSges(5,1) + t15 * rSges(5,2)) + g(2) * (t18 * rSges(5,1) + t17 * rSges(5,2) + t24) + (-g(1) * t41 + g(2) * t49) * t45 + (g(1) * (-t31 - t49) - t87) * t43) - m(6) * (g(1) * (-t12 * rSges(6,1) + t11 * rSges(6,2) + t61) + g(2) * (t14 * rSges(6,1) - t13 * rSges(6,2) + rSges(6,3) * t79 + t50) + (g(1) * (-t33 * rSges(6,3) + t62) - t87) * t43) - m(7) * (g(1) * (-t64 * t11 - t12 * t83 + t61) + g(2) * (rSges(7,2) * t79 + t64 * t13 + t14 * t83 + t50) + (g(1) * (-t33 * rSges(7,2) + t62) - t87) * t43) (-m(3) - m(4) - m(5) - m(6) - m(7)) * (g(1) * t43 - g(2) * t45) -m(4) * (g(3) * t60 + t93 * (-rSges(4,1) * t33 - rSges(4,2) * t34)) - m(5) * ((g(3) * t52 + t82 * t93) * t34 + (g(3) * t82 - t52 * t93) * t33) - m(6) * (t85 + (g(3) * t58 + t67 * t93) * t34 + (g(3) * t67 + t93 * (-t32 - t58)) * t33) - m(7) * (t85 + (g(3) * t92 + t69 * t93) * t34 + (g(3) * t69 + t93 * (-t32 - t92)) * t33) -m(5) * (g(1) * (rSges(5,1) * t17 - rSges(5,2) * t18) + g(2) * (-rSges(5,1) * t15 + rSges(5,2) * t16)) - m(6) * (g(1) * (t55 + t90) + g(2) * (-t48 + t91)) - m(7) * (g(1) * (t54 + t55) + g(2) * (-t48 + t56) + g(3) * t66) + (-m(5) * (-rSges(5,1) * t42 - rSges(5,2) * t44) - m(6) * (t57 - t89) - m(7) * (-t35 * t83 - t89)) * t84, -m(6) * (g(1) * t90 + g(2) * t91 + t57 * t84) - m(7) * (g(1) * t54 + g(2) * t56 + g(3) * (-t83 * t81 + t66)) -m(7) * (g(1) * t13 + g(2) * t11 + g(3) * t81)];
taug  = t1(:);
