% Calculate potential energy for
% S6PRPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRRP6_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP6_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRPRRP6_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP6_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP6_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP6_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRP6_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:18:07
% EndTime: 2019-03-08 20:18:07
% DurationCPUTime: 0.43s
% Computational Cost: add. (296->120), mult. (582->146), div. (0->0), fcn. (686->10), ass. (0->55)
t36 = sin(pkin(10));
t73 = g(1) * t36;
t38 = cos(pkin(10));
t72 = g(2) * t38;
t71 = rSges(7,1) + pkin(5);
t70 = rSges(7,2) + pkin(9);
t69 = rSges(6,3) + pkin(9);
t37 = sin(pkin(6));
t68 = t36 * t37;
t41 = sin(qJ(4));
t67 = t37 * t41;
t42 = sin(qJ(2));
t66 = t37 * t42;
t44 = cos(qJ(4));
t65 = t37 * t44;
t45 = cos(qJ(2));
t64 = t37 * t45;
t39 = cos(pkin(6));
t63 = t39 * t42;
t62 = t39 * t45;
t61 = qJ(3) * t45;
t60 = rSges(7,3) + qJ(6);
t59 = t36 * pkin(1) + r_base(2);
t58 = qJ(1) + r_base(3);
t57 = (-pkin(3) - pkin(7)) * t38;
t56 = t38 * pkin(1) + pkin(7) * t68 + r_base(1);
t55 = t39 * pkin(7) + t58;
t54 = pkin(2) * t66 + t55;
t19 = t36 * t42 - t38 * t62;
t20 = t36 * t45 + t38 * t63;
t53 = t20 * pkin(2) + qJ(3) * t19 + t59;
t52 = t39 * pkin(3) + pkin(8) * t66 + t54;
t21 = t36 * t62 + t38 * t42;
t22 = -t36 * t63 + t38 * t45;
t51 = t22 * pkin(2) + qJ(3) * t21 + t56;
t50 = pkin(3) * t68 + t51;
t49 = pkin(8) * t20 + t53;
t24 = t39 * t44 - t41 * t64;
t48 = t24 * pkin(4) - t37 * t61 + t52;
t8 = t21 * t41 + t36 * t65;
t47 = t8 * pkin(4) + pkin(8) * t22 + t50;
t10 = t19 * t41 - t38 * t65;
t46 = t10 * pkin(4) + t37 * t57 + t49;
t43 = cos(qJ(5));
t40 = sin(qJ(5));
t23 = t39 * t41 + t44 * t64;
t12 = t24 * t43 + t40 * t66;
t11 = t24 * t40 - t43 * t66;
t9 = t19 * t44 + t38 * t67;
t7 = -t21 * t44 + t36 * t67;
t4 = t10 * t43 + t20 * t40;
t3 = t10 * t40 - t20 * t43;
t2 = t22 * t40 + t43 * t8;
t1 = -t22 * t43 + t40 * t8;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t38 - rSges(2,2) * t36 + r_base(1)) + g(2) * (rSges(2,1) * t36 + rSges(2,2) * t38 + r_base(2)) + g(3) * (rSges(2,3) + t58)) - m(3) * (g(1) * (rSges(3,1) * t22 - rSges(3,2) * t21 + t56) + g(2) * (rSges(3,1) * t20 - rSges(3,2) * t19 + t59) + g(3) * (t39 * rSges(3,3) + t55) + (rSges(3,3) * t73 + g(3) * (rSges(3,1) * t42 + rSges(3,2) * t45) + (-rSges(3,3) - pkin(7)) * t72) * t37) - m(4) * (g(1) * (-rSges(4,2) * t22 + rSges(4,3) * t21 + t51) + g(2) * (-rSges(4,2) * t20 + rSges(4,3) * t19 + t53) + g(3) * (t39 * rSges(4,1) + t54) + (rSges(4,1) * t73 + g(3) * (-rSges(4,2) * t42 - rSges(4,3) * t45 - t61) + (-rSges(4,1) - pkin(7)) * t72) * t37) - m(5) * (g(1) * (rSges(5,1) * t8 - rSges(5,2) * t7 + (rSges(5,3) + pkin(8)) * t22 + t50) + g(2) * (rSges(5,1) * t10 + rSges(5,2) * t9 + rSges(5,3) * t20 + t49) + g(3) * (t24 * rSges(5,1) - t23 * rSges(5,2) + t52) + (g(3) * (rSges(5,3) * t42 - t61) + g(2) * t57) * t37) - m(6) * (g(1) * (rSges(6,1) * t2 - rSges(6,2) * t1 + t69 * t7 + t47) + g(2) * (rSges(6,1) * t4 - rSges(6,2) * t3 - t69 * t9 + t46) + g(3) * (t12 * rSges(6,1) - t11 * rSges(6,2) + t69 * t23 + t48)) - m(7) * (g(1) * (t60 * t1 + t71 * t2 + t70 * t7 + t47) + g(2) * (t60 * t3 + t71 * t4 - t70 * t9 + t46) + g(3) * (t60 * t11 + t71 * t12 + t70 * t23 + t48));
U  = t5;
