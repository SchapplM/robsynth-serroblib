% Calculate potential energy for
% S6RRRRRP9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 02:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRRP9_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP9_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRRRP9_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP9_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP9_energypot_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP9_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP9_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:00:10
% EndTime: 2019-03-10 02:00:11
% DurationCPUTime: 0.44s
% Computational Cost: add. (344->133), mult. (604->162), div. (0->0), fcn. (717->12), ass. (0->53)
t36 = sin(qJ(4));
t52 = t36 * pkin(4) + pkin(9);
t43 = -pkin(11) - pkin(10);
t65 = rSges(4,3) + pkin(9);
t34 = qJ(4) + qJ(5);
t27 = sin(t34);
t64 = pkin(5) * t27 + t52;
t63 = pkin(10) + rSges(5,3);
t40 = cos(qJ(4));
t26 = t40 * pkin(4) + pkin(3);
t62 = cos(qJ(3));
t35 = sin(pkin(6));
t38 = sin(qJ(2));
t61 = t35 * t38;
t39 = sin(qJ(1));
t60 = t35 * t39;
t41 = cos(qJ(2));
t59 = t35 * t41;
t42 = cos(qJ(1));
t58 = t35 * t42;
t57 = rSges(6,3) - t43;
t56 = rSges(7,3) + qJ(6) - t43;
t55 = cos(pkin(6));
t54 = pkin(7) + r_base(3);
t53 = t39 * pkin(1) + r_base(2);
t51 = t35 * t62;
t50 = t55 * pkin(8) + t54;
t49 = t39 * t55;
t48 = t42 * t55;
t47 = t42 * pkin(1) + pkin(8) * t60 + r_base(1);
t46 = pkin(2) * t61 + t50;
t18 = -t38 * t49 + t42 * t41;
t45 = t18 * pkin(2) + t47;
t16 = t38 * t48 + t39 * t41;
t44 = t16 * pkin(2) - pkin(8) * t58 + t53;
t37 = sin(qJ(3));
t28 = cos(t34);
t19 = pkin(5) * t28 + t26;
t17 = t42 * t38 + t41 * t49;
t15 = t39 * t38 - t41 * t48;
t14 = t55 * t37 + t38 * t51;
t13 = t37 * t61 - t55 * t62;
t10 = t18 * t62 + t37 * t60;
t9 = t18 * t37 - t39 * t51;
t8 = t16 * t62 - t37 * t58;
t7 = t16 * t37 + t42 * t51;
t6 = t14 * t28 - t27 * t59;
t5 = -t14 * t27 - t28 * t59;
t4 = t10 * t28 + t17 * t27;
t3 = -t10 * t27 + t17 * t28;
t2 = t15 * t27 + t8 * t28;
t1 = t15 * t28 - t8 * t27;
t11 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t42 * rSges(2,1) - t39 * rSges(2,2) + r_base(1)) + g(2) * (t39 * rSges(2,1) + t42 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t54)) - m(3) * (g(1) * (t18 * rSges(3,1) - t17 * rSges(3,2) + t47) + g(2) * (t16 * rSges(3,1) - t15 * rSges(3,2) + t53) + g(3) * (t55 * rSges(3,3) + t50) + (g(1) * rSges(3,3) * t39 + g(3) * (rSges(3,1) * t38 + rSges(3,2) * t41) + g(2) * (-rSges(3,3) - pkin(8)) * t42) * t35) - m(4) * (g(1) * (t10 * rSges(4,1) - t9 * rSges(4,2) + t65 * t17 + t45) + g(2) * (t8 * rSges(4,1) - t7 * rSges(4,2) + t65 * t15 + t44) + g(3) * (t14 * rSges(4,1) - t13 * rSges(4,2) - t65 * t59 + t46)) - m(5) * (g(1) * (t10 * pkin(3) + t17 * pkin(9) + (t10 * t40 + t17 * t36) * rSges(5,1) + (-t10 * t36 + t17 * t40) * rSges(5,2) + t63 * t9 + t45) + g(2) * (t8 * pkin(3) + t15 * pkin(9) + (t15 * t36 + t8 * t40) * rSges(5,1) + (t15 * t40 - t8 * t36) * rSges(5,2) + t63 * t7 + t44) + g(3) * (t14 * pkin(3) - pkin(9) * t59 + (t14 * t40 - t36 * t59) * rSges(5,1) + (-t14 * t36 - t40 * t59) * rSges(5,2) + t63 * t13 + t46)) - m(6) * (g(1) * (t4 * rSges(6,1) + t3 * rSges(6,2) + t10 * t26 + t52 * t17 + t57 * t9 + t45) + g(2) * (t2 * rSges(6,1) + t1 * rSges(6,2) + t52 * t15 + t8 * t26 + t57 * t7 + t44) + g(3) * (t6 * rSges(6,1) + t5 * rSges(6,2) + t57 * t13 + t14 * t26 - t52 * t59 + t46)) - m(7) * (g(1) * (t4 * rSges(7,1) + t3 * rSges(7,2) + t10 * t19 + t64 * t17 + t56 * t9 + t45) + g(2) * (t2 * rSges(7,1) + t1 * rSges(7,2) + t64 * t15 + t8 * t19 + t56 * t7 + t44) + g(3) * (t6 * rSges(7,1) + t5 * rSges(7,2) + t56 * t13 + t14 * t19 - t64 * t59 + t46));
U  = t11;
