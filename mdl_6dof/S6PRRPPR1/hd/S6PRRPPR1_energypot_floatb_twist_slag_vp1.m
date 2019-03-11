% Calculate potential energy for
% S6PRRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
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
% Datum: 2019-03-08 21:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRPPR1_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRRPPR1_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR1_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_energypot_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR1_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPPR1_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:58:24
% EndTime: 2019-03-08 20:58:25
% DurationCPUTime: 0.47s
% Computational Cost: add. (376->138), mult. (545->174), div. (0->0), fcn. (633->14), ass. (0->53)
t42 = sin(qJ(3));
t68 = pkin(3) * t42;
t67 = pkin(8) + rSges(4,3);
t35 = sin(pkin(10));
t38 = cos(pkin(10));
t43 = sin(qJ(2));
t39 = cos(pkin(6));
t45 = cos(qJ(2));
t57 = t39 * t45;
t11 = t35 * t43 - t38 * t57;
t34 = sin(pkin(12));
t66 = t11 * t34;
t13 = t35 * t57 + t38 * t43;
t65 = t13 * t34;
t36 = sin(pkin(6));
t64 = t35 * t36;
t63 = t36 * t38;
t62 = t36 * t42;
t61 = t36 * t43;
t44 = cos(qJ(3));
t60 = t36 * t44;
t59 = t36 * t45;
t58 = t39 * t43;
t56 = pkin(9) + qJ(5) + rSges(7,3);
t55 = qJ(5) + rSges(6,3);
t54 = t35 * pkin(1) + r_base(2);
t53 = t35 * t62;
t52 = t34 * t59;
t51 = qJ(1) + r_base(3);
t50 = t38 * pkin(1) + pkin(7) * t64 + r_base(1);
t49 = t39 * pkin(7) + t51;
t14 = -t35 * t58 + t38 * t45;
t24 = t44 * pkin(3) + pkin(2);
t40 = -qJ(4) - pkin(8);
t48 = pkin(3) * t53 - t13 * t40 + t14 * t24 + t50;
t47 = t24 * t61 + t39 * t68 + t40 * t59 + t49;
t12 = t35 * t45 + t38 * t58;
t46 = t12 * t24 + (-pkin(7) - t68) * t63 - t11 * t40 + t54;
t37 = cos(pkin(12));
t33 = qJ(3) + pkin(11);
t32 = pkin(12) + qJ(6);
t28 = cos(t33);
t27 = cos(t32);
t26 = sin(t33);
t25 = sin(t32);
t23 = t37 * pkin(5) + pkin(4);
t8 = t39 * t26 + t28 * t61;
t7 = t26 * t61 - t39 * t28;
t4 = t14 * t28 + t26 * t64;
t3 = t14 * t26 - t28 * t64;
t2 = t12 * t28 - t26 * t63;
t1 = t12 * t26 + t28 * t63;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t38 * rSges(2,1) - t35 * rSges(2,2) + r_base(1)) + g(2) * (t35 * rSges(2,1) + t38 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t51)) - m(3) * (g(1) * (t14 * rSges(3,1) - t13 * rSges(3,2) + t50) + g(2) * (t12 * rSges(3,1) - t11 * rSges(3,2) + t54) + g(3) * (t39 * rSges(3,3) + t49) + (g(1) * rSges(3,3) * t35 + g(3) * (rSges(3,1) * t43 + rSges(3,2) * t45) + g(2) * (-rSges(3,3) - pkin(7)) * t38) * t36) - m(4) * (g(1) * (t14 * pkin(2) + (t14 * t44 + t53) * rSges(4,1) + (-t14 * t42 + t35 * t60) * rSges(4,2) + t67 * t13 + t50) + g(2) * (t12 * pkin(2) - pkin(7) * t63 + (t12 * t44 - t38 * t62) * rSges(4,1) + (-t12 * t42 - t38 * t60) * rSges(4,2) + t67 * t11 + t54) + g(3) * ((t42 * rSges(4,1) + t44 * rSges(4,2)) * t39 + (-t67 * t45 + (t44 * rSges(4,1) - t42 * rSges(4,2) + pkin(2)) * t43) * t36 + t49)) - m(5) * (g(1) * (t4 * rSges(5,1) - t3 * rSges(5,2) + t13 * rSges(5,3) + t48) + g(2) * (t2 * rSges(5,1) - t1 * rSges(5,2) + t11 * rSges(5,3) + t46) + g(3) * (t8 * rSges(5,1) - t7 * rSges(5,2) - rSges(5,3) * t59 + t47)) - m(6) * (g(1) * (t4 * pkin(4) + (t4 * t37 + t65) * rSges(6,1) + (t13 * t37 - t4 * t34) * rSges(6,2) + t55 * t3 + t48) + g(2) * (t2 * pkin(4) + (t2 * t37 + t66) * rSges(6,1) + (t11 * t37 - t2 * t34) * rSges(6,2) + t55 * t1 + t46) + g(3) * (t8 * pkin(4) + (t8 * t37 - t52) * rSges(6,1) + (-t8 * t34 - t37 * t59) * rSges(6,2) + t55 * t7 + t47)) - m(7) * (g(1) * (t4 * t23 + pkin(5) * t65 + (t13 * t25 + t4 * t27) * rSges(7,1) + (t13 * t27 - t4 * t25) * rSges(7,2) + t56 * t3 + t48) + g(2) * (t2 * t23 + pkin(5) * t66 + (t11 * t25 + t2 * t27) * rSges(7,1) + (t11 * t27 - t2 * t25) * rSges(7,2) + t56 * t1 + t46) + g(3) * (t8 * t23 - pkin(5) * t52 + (-t25 * t59 + t8 * t27) * rSges(7,1) + (-t8 * t25 - t27 * t59) * rSges(7,2) + t56 * t7 + t47));
U  = t5;
