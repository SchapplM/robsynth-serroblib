% Calculate Gravitation load on the joints for
% S5PRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
% m_mdh [6x1]
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
% Datum: 2019-12-05 17:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRP8_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP8_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP8_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP8_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP8_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRP8_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:58:14
% EndTime: 2019-12-05 16:58:16
% DurationCPUTime: 0.49s
% Computational Cost: add. (349->118), mult. (879->185), div. (0->0), fcn. (1053->10), ass. (0->64)
t81 = rSges(6,2) + pkin(8);
t51 = cos(qJ(3));
t80 = pkin(3) * t51;
t46 = sin(pkin(5));
t79 = g(3) * t46;
t78 = rSges(6,1) + pkin(4);
t77 = rSges(4,3) + pkin(7);
t76 = rSges(5,3) + pkin(8);
t45 = sin(pkin(9));
t49 = sin(qJ(2));
t52 = cos(qJ(2));
t64 = cos(pkin(9));
t65 = cos(pkin(5));
t54 = t65 * t64;
t30 = t45 * t49 - t52 * t54;
t48 = sin(qJ(3));
t75 = t30 * t48;
t61 = t45 * t65;
t32 = t49 * t64 + t52 * t61;
t74 = t32 * t48;
t73 = t46 * t49;
t72 = t46 * t51;
t71 = t46 * t52;
t47 = sin(qJ(4));
t70 = t47 * t51;
t50 = cos(qJ(4));
t69 = t50 * t51;
t68 = t51 * t52;
t67 = pkin(2) * t71 + pkin(7) * t73;
t66 = rSges(6,3) + qJ(5);
t63 = t48 * t71;
t62 = t47 * t71;
t60 = t46 * t64;
t59 = t46 * pkin(3) * t68 + pkin(8) * t63 + t67;
t58 = rSges(4,1) * t51 - rSges(4,2) * t48;
t57 = rSges(5,1) * t50 - rSges(5,2) * t47;
t27 = t30 * pkin(2);
t31 = t45 * t52 + t49 * t54;
t56 = t31 * pkin(7) - pkin(8) * t75 - t30 * t80 - t27;
t28 = t32 * pkin(2);
t33 = -t49 * t61 + t52 * t64;
t55 = t33 * pkin(7) - pkin(8) * t74 - t32 * t80 - t28;
t35 = t48 * t65 + t49 * t72;
t34 = -t48 * t73 + t51 * t65;
t29 = t34 * pkin(3);
t18 = (t47 * t49 + t50 * t68) * t46;
t17 = -t50 * t73 + t51 * t62;
t16 = t35 * t50 - t62;
t15 = t35 * t47 + t50 * t71;
t14 = t45 * t46 * t48 + t33 * t51;
t13 = -t33 * t48 + t45 * t72;
t12 = t31 * t51 - t48 * t60;
t11 = -t31 * t48 - t51 * t60;
t10 = t13 * pkin(3);
t9 = t11 * pkin(3);
t8 = -t32 * t69 + t33 * t47;
t7 = -t32 * t70 - t33 * t50;
t6 = -t30 * t69 + t31 * t47;
t5 = -t30 * t70 - t31 * t50;
t4 = t14 * t50 + t32 * t47;
t3 = t14 * t47 - t32 * t50;
t2 = t12 * t50 + t30 * t47;
t1 = t12 * t47 - t30 * t50;
t19 = [(-m(2) - m(3) - m(4) - m(5) - m(6)) * g(3), -m(3) * (g(1) * (-t32 * rSges(3,1) - t33 * rSges(3,2)) + g(2) * (-t30 * rSges(3,1) - t31 * rSges(3,2)) + (rSges(3,1) * t52 - rSges(3,2) * t49) * t79) - m(4) * (g(1) * (-t32 * t58 + t33 * t77 - t28) + g(2) * (-t30 * t58 + t31 * t77 - t27) + g(3) * t67 + (rSges(4,3) * t49 + t52 * t58) * t79) - m(5) * (g(1) * (t8 * rSges(5,1) - t7 * rSges(5,2) - rSges(5,3) * t74 + t55) + g(2) * (t6 * rSges(5,1) - t5 * rSges(5,2) - rSges(5,3) * t75 + t56) + g(3) * (t18 * rSges(5,1) - t17 * rSges(5,2) + rSges(5,3) * t63 + t59)) - m(6) * (g(1) * (-rSges(6,2) * t74 + t66 * t7 + t78 * t8 + t55) + g(2) * (-rSges(6,2) * t75 + t5 * t66 + t6 * t78 + t56) + g(3) * (rSges(6,2) * t63 + t17 * t66 + t18 * t78 + t59)), -m(4) * (g(1) * (t13 * rSges(4,1) - t14 * rSges(4,2)) + g(2) * (t11 * rSges(4,1) - t12 * rSges(4,2)) + g(3) * (t34 * rSges(4,1) - t35 * rSges(4,2))) - m(5) * (g(1) * (t13 * t57 + t14 * t76 + t10) + g(2) * (t11 * t57 + t12 * t76 + t9) + g(3) * (t34 * t57 + t35 * t76 + t29)) + (-g(1) * (t81 * t14 + t10) - g(2) * (t81 * t12 + t9) - g(3) * (t81 * t35 + t29) - (g(1) * t13 + g(2) * t11 + g(3) * t34) * (t47 * t66 + t50 * t78)) * m(6), -m(5) * (g(1) * (-t3 * rSges(5,1) - t4 * rSges(5,2)) + g(2) * (-t1 * rSges(5,1) - t2 * rSges(5,2)) + g(3) * (-t15 * rSges(5,1) - t16 * rSges(5,2))) - m(6) * (g(1) * (-t3 * t78 + t4 * t66) + g(2) * (-t1 * t78 + t2 * t66) + g(3) * (-t15 * t78 + t16 * t66)), -m(6) * (g(1) * t3 + g(2) * t1 + g(3) * t15)];
taug = t19(:);
