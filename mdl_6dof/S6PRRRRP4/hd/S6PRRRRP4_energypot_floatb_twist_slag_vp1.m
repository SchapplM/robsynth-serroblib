% Calculate potential energy for
% S6PRRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-03-09 00:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRRRP4_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP4_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRRRRP4_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP4_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP4_energypot_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP4_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRP4_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:13:01
% EndTime: 2019-03-09 00:13:01
% DurationCPUTime: 0.38s
% Computational Cost: add. (367->128), mult. (660->157), div. (0->0), fcn. (793->12), ass. (0->57)
t42 = sin(pkin(6));
t76 = pkin(7) * t42;
t75 = rSges(7,1) + pkin(5);
t74 = rSges(4,3) + pkin(8);
t73 = pkin(9) + rSges(5,3);
t41 = sin(pkin(11));
t43 = cos(pkin(11));
t47 = sin(qJ(2));
t44 = cos(pkin(6));
t50 = cos(qJ(2));
t65 = t44 * t50;
t23 = t41 * t47 - t43 * t65;
t45 = sin(qJ(4));
t72 = t23 * t45;
t25 = t41 * t65 + t43 * t47;
t71 = t25 * t45;
t46 = sin(qJ(3));
t70 = t42 * t46;
t69 = t42 * t47;
t49 = cos(qJ(3));
t68 = t42 * t49;
t67 = t42 * t50;
t66 = t44 * t47;
t64 = rSges(7,3) + qJ(6);
t63 = t41 * pkin(1) + r_base(2);
t62 = qJ(1) + r_base(3);
t61 = t43 * pkin(1) + t41 * t76 + r_base(1);
t60 = t44 * pkin(7) + t62;
t26 = -t41 * t66 + t43 * t50;
t59 = t26 * pkin(2) + t61;
t58 = pkin(2) * t69 + t60;
t24 = t41 * t50 + t43 * t66;
t57 = t24 * pkin(2) - t43 * t76 + t63;
t56 = t25 * pkin(8) + t59;
t55 = t23 * pkin(8) + t57;
t13 = t26 * t46 - t41 * t68;
t14 = t26 * t49 + t41 * t70;
t48 = cos(qJ(4));
t34 = t48 * pkin(4) + pkin(3);
t51 = -pkin(10) - pkin(9);
t54 = pkin(4) * t71 - t13 * t51 + t14 * t34 + t56;
t11 = t24 * t46 + t43 * t68;
t12 = t24 * t49 - t43 * t70;
t53 = pkin(4) * t72 - t11 * t51 + t12 * t34 + t55;
t27 = -t44 * t49 + t46 * t69;
t28 = t44 * t46 + t47 * t68;
t52 = t28 * t34 + (-pkin(4) * t45 - pkin(8)) * t67 - t27 * t51 + t58;
t40 = qJ(4) + qJ(5);
t36 = cos(t40);
t35 = sin(t40);
t8 = t28 * t36 - t35 * t67;
t7 = t28 * t35 + t36 * t67;
t4 = t14 * t36 + t25 * t35;
t3 = t14 * t35 - t25 * t36;
t2 = t12 * t36 + t23 * t35;
t1 = t12 * t35 - t23 * t36;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t43 * rSges(2,1) - t41 * rSges(2,2) + r_base(1)) + g(2) * (t41 * rSges(2,1) + t43 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t62)) - m(3) * (g(1) * (t26 * rSges(3,1) - t25 * rSges(3,2) + t61) + g(2) * (t24 * rSges(3,1) - t23 * rSges(3,2) + t63) + g(3) * (t44 * rSges(3,3) + t60) + (g(1) * rSges(3,3) * t41 + g(3) * (rSges(3,1) * t47 + rSges(3,2) * t50) + g(2) * (-rSges(3,3) - pkin(7)) * t43) * t42) - m(4) * (g(1) * (t14 * rSges(4,1) - t13 * rSges(4,2) + t25 * t74 + t59) + g(2) * (t12 * rSges(4,1) - t11 * rSges(4,2) + t23 * t74 + t57) + g(3) * (t28 * rSges(4,1) - t27 * rSges(4,2) - t67 * t74 + t58)) - m(5) * (g(1) * (t14 * pkin(3) + (t14 * t48 + t71) * rSges(5,1) + (-t14 * t45 + t25 * t48) * rSges(5,2) + t73 * t13 + t56) + g(2) * (t12 * pkin(3) + (t12 * t48 + t72) * rSges(5,1) + (-t12 * t45 + t23 * t48) * rSges(5,2) + t73 * t11 + t55) + g(3) * (t28 * pkin(3) - pkin(8) * t67 + (t28 * t48 - t45 * t67) * rSges(5,1) + (-t28 * t45 - t48 * t67) * rSges(5,2) + t73 * t27 + t58)) - m(6) * (g(1) * (t4 * rSges(6,1) - t3 * rSges(6,2) + t13 * rSges(6,3) + t54) + g(2) * (t2 * rSges(6,1) - t1 * rSges(6,2) + t11 * rSges(6,3) + t53) + g(3) * (t8 * rSges(6,1) - t7 * rSges(6,2) + t27 * rSges(6,3) + t52)) - m(7) * (g(1) * (t13 * rSges(7,2) + t3 * t64 + t4 * t75 + t54) + g(2) * (t11 * rSges(7,2) + t1 * t64 + t2 * t75 + t53) + g(3) * (t27 * rSges(7,2) + t64 * t7 + t75 * t8 + t52));
U  = t5;
