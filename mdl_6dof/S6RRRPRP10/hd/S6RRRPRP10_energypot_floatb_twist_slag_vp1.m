% Calculate potential energy for
% S6RRRPRP10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-03-09 17:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPRP10_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP10_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRPRP10_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP10_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP10_energypot_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP10_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP10_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:29:22
% EndTime: 2019-03-09 17:29:22
% DurationCPUTime: 0.38s
% Computational Cost: add. (367->128), mult. (660->155), div. (0->0), fcn. (793->12), ass. (0->59)
t78 = rSges(7,1) + pkin(5);
t77 = rSges(4,3) + pkin(9);
t44 = cos(pkin(6));
t50 = cos(qJ(2));
t51 = cos(qJ(1));
t66 = t51 * t50;
t47 = sin(qJ(2));
t48 = sin(qJ(1));
t69 = t48 * t47;
t25 = -t44 * t66 + t69;
t41 = sin(pkin(11));
t76 = t25 * t41;
t67 = t51 * t47;
t68 = t48 * t50;
t27 = t44 * t68 + t67;
t75 = t27 * t41;
t42 = sin(pkin(6));
t74 = t42 * t47;
t73 = t42 * t48;
t49 = cos(qJ(3));
t72 = t42 * t49;
t71 = t42 * t50;
t70 = t42 * t51;
t65 = rSges(7,3) + qJ(6);
t64 = qJ(4) + rSges(5,3);
t63 = pkin(7) + r_base(3);
t62 = t48 * pkin(1) + r_base(2);
t61 = t44 * pkin(8) + t63;
t60 = t51 * pkin(1) + pkin(8) * t73 + r_base(1);
t59 = pkin(2) * t74 + t61;
t28 = -t44 * t69 + t66;
t58 = t28 * pkin(2) + t60;
t26 = t44 * t67 + t68;
t57 = t26 * pkin(2) - pkin(8) * t70 + t62;
t56 = t27 * pkin(9) + t58;
t55 = t25 * pkin(9) + t57;
t46 = sin(qJ(3));
t13 = t28 * t46 - t48 * t72;
t14 = t28 * t49 + t46 * t73;
t43 = cos(pkin(11));
t34 = t43 * pkin(4) + pkin(3);
t45 = -pkin(10) - qJ(4);
t54 = pkin(4) * t75 - t13 * t45 + t14 * t34 + t56;
t11 = t26 * t46 + t49 * t70;
t12 = t26 * t49 - t46 * t70;
t53 = pkin(4) * t76 - t11 * t45 + t12 * t34 + t55;
t23 = -t44 * t49 + t46 * t74;
t24 = t44 * t46 + t47 * t72;
t52 = t24 * t34 + (-pkin(4) * t41 - pkin(9)) * t71 - t23 * t45 + t59;
t40 = pkin(11) + qJ(5);
t36 = cos(t40);
t35 = sin(t40);
t8 = t24 * t36 - t35 * t71;
t7 = t24 * t35 + t36 * t71;
t4 = t14 * t36 + t27 * t35;
t3 = t14 * t35 - t27 * t36;
t2 = t12 * t36 + t25 * t35;
t1 = t12 * t35 - t25 * t36;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t51 * rSges(2,1) - t48 * rSges(2,2) + r_base(1)) + g(2) * (t48 * rSges(2,1) + t51 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t63)) - m(3) * (g(1) * (t28 * rSges(3,1) - t27 * rSges(3,2) + t60) + g(2) * (t26 * rSges(3,1) - t25 * rSges(3,2) + t62) + g(3) * (t44 * rSges(3,3) + t61) + (g(1) * rSges(3,3) * t48 + g(3) * (rSges(3,1) * t47 + rSges(3,2) * t50) + g(2) * (-rSges(3,3) - pkin(8)) * t51) * t42) - m(4) * (g(1) * (t14 * rSges(4,1) - t13 * rSges(4,2) + t77 * t27 + t58) + g(2) * (t12 * rSges(4,1) - t11 * rSges(4,2) + t77 * t25 + t57) + g(3) * (t24 * rSges(4,1) - t23 * rSges(4,2) - t77 * t71 + t59)) - m(5) * (g(1) * (t14 * pkin(3) + (t14 * t43 + t75) * rSges(5,1) + (-t14 * t41 + t27 * t43) * rSges(5,2) + t64 * t13 + t56) + g(2) * (t12 * pkin(3) + (t12 * t43 + t76) * rSges(5,1) + (-t12 * t41 + t25 * t43) * rSges(5,2) + t64 * t11 + t55) + g(3) * (t24 * pkin(3) - pkin(9) * t71 + (t24 * t43 - t41 * t71) * rSges(5,1) + (-t24 * t41 - t43 * t71) * rSges(5,2) + t64 * t23 + t59)) - m(6) * (g(1) * (t4 * rSges(6,1) - t3 * rSges(6,2) + t13 * rSges(6,3) + t54) + g(2) * (t2 * rSges(6,1) - t1 * rSges(6,2) + t11 * rSges(6,3) + t53) + g(3) * (t8 * rSges(6,1) - t7 * rSges(6,2) + t23 * rSges(6,3) + t52)) - m(7) * (g(1) * (t13 * rSges(7,2) + t65 * t3 + t78 * t4 + t54) + g(2) * (t11 * rSges(7,2) + t65 * t1 + t78 * t2 + t53) + g(3) * (t23 * rSges(7,2) + t65 * t7 + t78 * t8 + t52));
U  = t5;
