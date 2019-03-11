% Calculate potential energy for
% S6RRPRRP14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
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
% Datum: 2019-03-09 13:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRP14_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP14_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPRRP14_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP14_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP14_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP14_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP14_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:03:39
% EndTime: 2019-03-09 13:03:39
% DurationCPUTime: 0.45s
% Computational Cost: add. (296->120), mult. (582->143), div. (0->0), fcn. (686->10), ass. (0->56)
t41 = sin(qJ(1));
t74 = g(1) * t41;
t45 = cos(qJ(1));
t73 = g(2) * t45;
t72 = rSges(7,1) + pkin(5);
t71 = rSges(7,2) + pkin(10);
t70 = rSges(6,3) + pkin(10);
t36 = sin(pkin(6));
t40 = sin(qJ(2));
t69 = t36 * t40;
t68 = t36 * t41;
t44 = cos(qJ(2));
t67 = t36 * t44;
t66 = t36 * t45;
t65 = t40 * t41;
t64 = t40 * t45;
t63 = t41 * t44;
t62 = t44 * t45;
t61 = qJ(3) * t44;
t60 = rSges(7,3) + qJ(6);
t59 = pkin(7) + r_base(3);
t58 = t41 * pkin(1) + r_base(2);
t57 = (-pkin(3) - pkin(8)) * t45;
t37 = cos(pkin(6));
t56 = t37 * pkin(8) + t59;
t55 = t45 * pkin(1) + pkin(8) * t68 + r_base(1);
t54 = pkin(2) * t69 + t56;
t53 = t37 * pkin(3) + pkin(9) * t69 + t54;
t21 = -t37 * t62 + t65;
t22 = t37 * t64 + t63;
t52 = t22 * pkin(2) + t21 * qJ(3) + t58;
t23 = t37 * t63 + t64;
t24 = -t37 * t65 + t62;
t51 = t24 * pkin(2) + qJ(3) * t23 + t55;
t50 = pkin(3) * t68 + t51;
t49 = t22 * pkin(9) + t52;
t39 = sin(qJ(4));
t43 = cos(qJ(4));
t20 = t37 * t43 - t39 * t67;
t48 = t20 * pkin(4) - t36 * t61 + t53;
t10 = t23 * t39 + t43 * t68;
t47 = t10 * pkin(4) + pkin(9) * t24 + t50;
t12 = t21 * t39 - t43 * t66;
t46 = t12 * pkin(4) + t36 * t57 + t49;
t42 = cos(qJ(5));
t38 = sin(qJ(5));
t19 = t37 * t39 + t43 * t67;
t11 = t21 * t43 + t39 * t66;
t9 = -t23 * t43 + t39 * t68;
t8 = t20 * t42 + t38 * t69;
t7 = t20 * t38 - t42 * t69;
t4 = t12 * t42 + t22 * t38;
t3 = t12 * t38 - t22 * t42;
t2 = t10 * t42 + t24 * t38;
t1 = t10 * t38 - t24 * t42;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t45 - t41 * rSges(2,2) + r_base(1)) + g(2) * (t41 * rSges(2,1) + rSges(2,2) * t45 + r_base(2)) + g(3) * (rSges(2,3) + t59)) - m(3) * (g(1) * (rSges(3,1) * t24 - rSges(3,2) * t23 + t55) + g(2) * (t22 * rSges(3,1) - t21 * rSges(3,2) + t58) + g(3) * (rSges(3,3) * t37 + t56) + (rSges(3,3) * t74 + g(3) * (rSges(3,1) * t40 + rSges(3,2) * t44) + (-rSges(3,3) - pkin(8)) * t73) * t36) - m(4) * (g(1) * (-rSges(4,2) * t24 + rSges(4,3) * t23 + t51) + g(2) * (-t22 * rSges(4,2) + t21 * rSges(4,3) + t52) + g(3) * (rSges(4,1) * t37 + t54) + (rSges(4,1) * t74 + g(3) * (-rSges(4,2) * t40 - rSges(4,3) * t44 - t61) + (-rSges(4,1) - pkin(8)) * t73) * t36) - m(5) * (g(1) * (rSges(5,1) * t10 - rSges(5,2) * t9 + (rSges(5,3) + pkin(9)) * t24 + t50) + g(2) * (t12 * rSges(5,1) + t11 * rSges(5,2) + t22 * rSges(5,3) + t49) + g(3) * (rSges(5,1) * t20 - rSges(5,2) * t19 + t53) + (g(3) * (rSges(5,3) * t40 - t61) + g(2) * t57) * t36) - m(6) * (g(1) * (rSges(6,1) * t2 - rSges(6,2) * t1 + t70 * t9 + t47) + g(2) * (t4 * rSges(6,1) - t3 * rSges(6,2) - t70 * t11 + t46) + g(3) * (rSges(6,1) * t8 - rSges(6,2) * t7 + t70 * t19 + t48)) - m(7) * (g(1) * (t60 * t1 + t72 * t2 + t71 * t9 + t47) + g(2) * (-t71 * t11 + t60 * t3 + t72 * t4 + t46) + g(3) * (t71 * t19 + t60 * t7 + t72 * t8 + t48));
U  = t5;
