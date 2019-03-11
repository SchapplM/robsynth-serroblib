% Calculate potential energy for
% S6PRRPRP3
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
% Datum: 2019-03-08 21:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRPRP3_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP3_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRRPRP3_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP3_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP3_energypot_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP3_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRP3_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:34:06
% EndTime: 2019-03-08 21:34:06
% DurationCPUTime: 0.36s
% Computational Cost: add. (367->128), mult. (660->157), div. (0->0), fcn. (793->12), ass. (0->57)
t43 = sin(pkin(6));
t76 = pkin(7) * t43;
t75 = rSges(7,1) + pkin(5);
t74 = rSges(4,3) + pkin(8);
t42 = sin(pkin(10));
t45 = cos(pkin(10));
t49 = sin(qJ(2));
t46 = cos(pkin(6));
t51 = cos(qJ(2));
t66 = t46 * t51;
t23 = t42 * t49 - t45 * t66;
t41 = sin(pkin(11));
t73 = t23 * t41;
t25 = t42 * t66 + t45 * t49;
t72 = t25 * t41;
t48 = sin(qJ(3));
t71 = t43 * t48;
t70 = t43 * t49;
t50 = cos(qJ(3));
t69 = t43 * t50;
t68 = t43 * t51;
t67 = t46 * t49;
t65 = rSges(7,3) + qJ(6);
t64 = qJ(4) + rSges(5,3);
t63 = t42 * pkin(1) + r_base(2);
t62 = qJ(1) + r_base(3);
t61 = t45 * pkin(1) + t42 * t76 + r_base(1);
t60 = t46 * pkin(7) + t62;
t26 = -t42 * t67 + t45 * t51;
t59 = t26 * pkin(2) + t61;
t58 = pkin(2) * t70 + t60;
t24 = t42 * t51 + t45 * t67;
t57 = t24 * pkin(2) - t45 * t76 + t63;
t56 = t25 * pkin(8) + t59;
t55 = t23 * pkin(8) + t57;
t13 = t26 * t48 - t42 * t69;
t14 = t26 * t50 + t42 * t71;
t44 = cos(pkin(11));
t34 = t44 * pkin(4) + pkin(3);
t47 = -pkin(9) - qJ(4);
t54 = pkin(4) * t72 - t13 * t47 + t14 * t34 + t56;
t11 = t24 * t48 + t45 * t69;
t12 = t24 * t50 - t45 * t71;
t53 = pkin(4) * t73 - t11 * t47 + t12 * t34 + t55;
t27 = -t46 * t50 + t48 * t70;
t28 = t46 * t48 + t49 * t69;
t52 = t28 * t34 + (-pkin(4) * t41 - pkin(8)) * t68 - t27 * t47 + t58;
t40 = pkin(11) + qJ(5);
t36 = cos(t40);
t35 = sin(t40);
t8 = t28 * t36 - t35 * t68;
t7 = t28 * t35 + t36 * t68;
t4 = t14 * t36 + t25 * t35;
t3 = t14 * t35 - t25 * t36;
t2 = t12 * t36 + t23 * t35;
t1 = t12 * t35 - t23 * t36;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t45 * rSges(2,1) - t42 * rSges(2,2) + r_base(1)) + g(2) * (t42 * rSges(2,1) + t45 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t62)) - m(3) * (g(1) * (t26 * rSges(3,1) - t25 * rSges(3,2) + t61) + g(2) * (t24 * rSges(3,1) - t23 * rSges(3,2) + t63) + g(3) * (t46 * rSges(3,3) + t60) + (g(1) * rSges(3,3) * t42 + g(3) * (rSges(3,1) * t49 + rSges(3,2) * t51) + g(2) * (-rSges(3,3) - pkin(7)) * t45) * t43) - m(4) * (g(1) * (t14 * rSges(4,1) - t13 * rSges(4,2) + t74 * t25 + t59) + g(2) * (t12 * rSges(4,1) - t11 * rSges(4,2) + t74 * t23 + t57) + g(3) * (t28 * rSges(4,1) - t27 * rSges(4,2) - t74 * t68 + t58)) - m(5) * (g(1) * (t14 * pkin(3) + (t14 * t44 + t72) * rSges(5,1) + (-t14 * t41 + t25 * t44) * rSges(5,2) + t64 * t13 + t56) + g(2) * (t12 * pkin(3) + (t12 * t44 + t73) * rSges(5,1) + (-t12 * t41 + t23 * t44) * rSges(5,2) + t64 * t11 + t55) + g(3) * (t28 * pkin(3) - pkin(8) * t68 + (t28 * t44 - t41 * t68) * rSges(5,1) + (-t28 * t41 - t44 * t68) * rSges(5,2) + t64 * t27 + t58)) - m(6) * (g(1) * (t4 * rSges(6,1) - t3 * rSges(6,2) + t13 * rSges(6,3) + t54) + g(2) * (t2 * rSges(6,1) - t1 * rSges(6,2) + t11 * rSges(6,3) + t53) + g(3) * (t8 * rSges(6,1) - t7 * rSges(6,2) + t27 * rSges(6,3) + t52)) - m(7) * (g(1) * (t13 * rSges(7,2) + t65 * t3 + t75 * t4 + t54) + g(2) * (t11 * rSges(7,2) + t65 * t1 + t75 * t2 + t53) + g(3) * (t27 * rSges(7,2) + t65 * t7 + t75 * t8 + t52));
U  = t5;
