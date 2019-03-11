% Calculate potential energy for
% S6RRRPRP7
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
% Datum: 2019-03-09 17:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPRP7_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP7_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRPRP7_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP7_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP7_energypot_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP7_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP7_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:03:44
% EndTime: 2019-03-09 17:03:44
% DurationCPUTime: 0.38s
% Computational Cost: add. (391->126), mult. (591->151), div. (0->0), fcn. (697->12), ass. (0->57)
t45 = sin(qJ(3));
t76 = pkin(3) * t45;
t75 = rSges(7,1) + pkin(5);
t74 = rSges(7,2) + pkin(10);
t73 = rSges(6,3) + pkin(10);
t72 = pkin(9) + rSges(4,3);
t41 = sin(pkin(6));
t46 = sin(qJ(2));
t71 = t41 * t46;
t47 = sin(qJ(1));
t70 = t41 * t47;
t50 = cos(qJ(2));
t69 = t41 * t50;
t51 = cos(qJ(1));
t68 = t41 * t51;
t67 = t47 * t46;
t66 = t47 * t50;
t65 = t51 * t46;
t64 = t51 * t50;
t63 = rSges(7,3) + qJ(6);
t62 = pkin(7) + r_base(3);
t61 = t47 * pkin(1) + r_base(2);
t60 = t45 * t70;
t42 = cos(pkin(6));
t59 = t42 * pkin(8) + t62;
t58 = t51 * pkin(1) + pkin(8) * t70 + r_base(1);
t49 = cos(qJ(3));
t34 = t49 * pkin(3) + pkin(2);
t43 = -qJ(4) - pkin(9);
t57 = t34 * t71 + t42 * t76 + t43 * t69 + t59;
t24 = t42 * t66 + t65;
t25 = -t42 * t67 + t64;
t56 = pkin(3) * t60 - t24 * t43 + t25 * t34 + t58;
t40 = qJ(3) + pkin(11);
t35 = sin(t40);
t36 = cos(t40);
t17 = t42 * t35 + t36 * t71;
t55 = t17 * pkin(4) + t57;
t12 = t25 * t36 + t35 * t70;
t54 = t12 * pkin(4) + t56;
t22 = -t42 * t64 + t67;
t23 = t42 * t65 + t66;
t53 = t23 * t34 + (-pkin(8) - t76) * t68 - t22 * t43 + t61;
t10 = t23 * t36 - t35 * t68;
t52 = t10 * pkin(4) + t53;
t48 = cos(qJ(5));
t44 = sin(qJ(5));
t16 = t35 * t71 - t42 * t36;
t11 = t25 * t35 - t36 * t70;
t9 = t23 * t35 + t36 * t68;
t8 = t17 * t48 - t44 * t69;
t7 = t17 * t44 + t48 * t69;
t4 = t12 * t48 + t24 * t44;
t3 = t12 * t44 - t24 * t48;
t2 = t10 * t48 + t22 * t44;
t1 = t10 * t44 - t22 * t48;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t51 * rSges(2,1) - t47 * rSges(2,2) + r_base(1)) + g(2) * (t47 * rSges(2,1) + t51 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t62)) - m(3) * (g(1) * (t25 * rSges(3,1) - t24 * rSges(3,2) + t58) + g(2) * (t23 * rSges(3,1) - t22 * rSges(3,2) + t61) + g(3) * (t42 * rSges(3,3) + t59) + (g(1) * rSges(3,3) * t47 + g(3) * (rSges(3,1) * t46 + rSges(3,2) * t50) + g(2) * (-rSges(3,3) - pkin(8)) * t51) * t41) - m(4) * (g(1) * (t25 * pkin(2) + (t25 * t49 + t60) * rSges(4,1) + (-t25 * t45 + t49 * t70) * rSges(4,2) + t72 * t24 + t58) + g(2) * (t23 * pkin(2) - pkin(8) * t68 + (t23 * t49 - t45 * t68) * rSges(4,1) + (-t23 * t45 - t49 * t68) * rSges(4,2) + t72 * t22 + t61) + g(3) * ((t45 * rSges(4,1) + t49 * rSges(4,2)) * t42 + (-t72 * t50 + (t49 * rSges(4,1) - t45 * rSges(4,2) + pkin(2)) * t46) * t41 + t59)) - m(5) * (g(1) * (t12 * rSges(5,1) - t11 * rSges(5,2) + t24 * rSges(5,3) + t56) + g(2) * (t10 * rSges(5,1) - t9 * rSges(5,2) + t22 * rSges(5,3) + t53) + g(3) * (t17 * rSges(5,1) - t16 * rSges(5,2) - rSges(5,3) * t69 + t57)) - m(6) * (g(1) * (t4 * rSges(6,1) - t3 * rSges(6,2) + t11 * t73 + t54) + g(2) * (t2 * rSges(6,1) - t1 * rSges(6,2) + t73 * t9 + t52) + g(3) * (t8 * rSges(6,1) - t7 * rSges(6,2) + t16 * t73 + t55)) - m(7) * (g(1) * (t11 * t74 + t3 * t63 + t4 * t75 + t54) + g(2) * (t1 * t63 + t2 * t75 + t74 * t9 + t52) + g(3) * (t16 * t74 + t63 * t7 + t75 * t8 + t55));
U  = t5;
