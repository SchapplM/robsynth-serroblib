% Calculate potential energy for
% S6PRRPRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-03-08 21:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRPRP2_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP2_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRRPRP2_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP2_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP2_energypot_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP2_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRP2_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:28:17
% EndTime: 2019-03-08 21:28:17
% DurationCPUTime: 0.37s
% Computational Cost: add. (391->126), mult. (591->155), div. (0->0), fcn. (697->12), ass. (0->57)
t47 = sin(qJ(3));
t76 = pkin(3) * t47;
t75 = rSges(7,1) + pkin(5);
t74 = rSges(7,2) + pkin(9);
t73 = rSges(6,3) + pkin(9);
t72 = pkin(8) + rSges(4,3);
t41 = sin(pkin(10));
t42 = sin(pkin(6));
t71 = t41 * t42;
t43 = cos(pkin(10));
t70 = t42 * t43;
t69 = t42 * t47;
t48 = sin(qJ(2));
t68 = t42 * t48;
t50 = cos(qJ(3));
t67 = t42 * t50;
t51 = cos(qJ(2));
t66 = t42 * t51;
t44 = cos(pkin(6));
t65 = t44 * t48;
t64 = t44 * t51;
t63 = rSges(7,3) + qJ(6);
t62 = t41 * pkin(1) + r_base(2);
t61 = t41 * t69;
t60 = qJ(1) + r_base(3);
t59 = t43 * pkin(1) + pkin(7) * t71 + r_base(1);
t58 = t44 * pkin(7) + t60;
t24 = t41 * t64 + t43 * t48;
t25 = -t41 * t65 + t43 * t51;
t34 = t50 * pkin(3) + pkin(2);
t45 = -qJ(4) - pkin(8);
t57 = pkin(3) * t61 - t24 * t45 + t25 * t34 + t59;
t56 = t34 * t68 + t44 * t76 + t45 * t66 + t58;
t40 = qJ(3) + pkin(11);
t35 = sin(t40);
t36 = cos(t40);
t10 = t25 * t36 + t35 * t71;
t55 = t10 * pkin(4) + t57;
t17 = t44 * t35 + t36 * t68;
t54 = t17 * pkin(4) + t56;
t22 = t41 * t48 - t43 * t64;
t23 = t41 * t51 + t43 * t65;
t53 = t23 * t34 + (-pkin(7) - t76) * t70 - t22 * t45 + t62;
t8 = t23 * t36 - t35 * t70;
t52 = t8 * pkin(4) + t53;
t49 = cos(qJ(5));
t46 = sin(qJ(5));
t16 = t35 * t68 - t44 * t36;
t12 = t17 * t49 - t46 * t66;
t11 = t17 * t46 + t49 * t66;
t9 = t25 * t35 - t36 * t71;
t7 = t23 * t35 + t36 * t70;
t4 = t10 * t49 + t24 * t46;
t3 = t10 * t46 - t24 * t49;
t2 = t22 * t46 + t8 * t49;
t1 = -t22 * t49 + t8 * t46;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t43 * rSges(2,1) - t41 * rSges(2,2) + r_base(1)) + g(2) * (t41 * rSges(2,1) + t43 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t60)) - m(3) * (g(1) * (t25 * rSges(3,1) - t24 * rSges(3,2) + t59) + g(2) * (t23 * rSges(3,1) - t22 * rSges(3,2) + t62) + g(3) * (t44 * rSges(3,3) + t58) + (g(1) * rSges(3,3) * t41 + g(3) * (rSges(3,1) * t48 + rSges(3,2) * t51) + g(2) * (-rSges(3,3) - pkin(7)) * t43) * t42) - m(4) * (g(1) * (t25 * pkin(2) + (t25 * t50 + t61) * rSges(4,1) + (-t25 * t47 + t41 * t67) * rSges(4,2) + t72 * t24 + t59) + g(2) * (t23 * pkin(2) - pkin(7) * t70 + (t23 * t50 - t43 * t69) * rSges(4,1) + (-t23 * t47 - t43 * t67) * rSges(4,2) + t72 * t22 + t62) + g(3) * ((t47 * rSges(4,1) + t50 * rSges(4,2)) * t44 + (-t72 * t51 + (t50 * rSges(4,1) - t47 * rSges(4,2) + pkin(2)) * t48) * t42 + t58)) - m(5) * (g(1) * (t10 * rSges(5,1) - t9 * rSges(5,2) + t24 * rSges(5,3) + t57) + g(2) * (t8 * rSges(5,1) - t7 * rSges(5,2) + t22 * rSges(5,3) + t53) + g(3) * (t17 * rSges(5,1) - t16 * rSges(5,2) - rSges(5,3) * t66 + t56)) - m(6) * (g(1) * (t4 * rSges(6,1) - t3 * rSges(6,2) + t73 * t9 + t55) + g(2) * (t2 * rSges(6,1) - t1 * rSges(6,2) + t73 * t7 + t52) + g(3) * (t12 * rSges(6,1) - t11 * rSges(6,2) + t73 * t16 + t54)) - m(7) * (g(1) * (t63 * t3 + t75 * t4 + t74 * t9 + t55) + g(2) * (t63 * t1 + t75 * t2 + t74 * t7 + t52) + g(3) * (t63 * t11 + t75 * t12 + t74 * t16 + t54));
U  = t5;
