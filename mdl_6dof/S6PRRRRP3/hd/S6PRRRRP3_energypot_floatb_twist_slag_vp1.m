% Calculate potential energy for
% S6PRRRRP3
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
% Datum: 2019-03-09 00:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRRRP3_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP3_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRRRRP3_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP3_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP3_energypot_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP3_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRP3_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:06:15
% EndTime: 2019-03-09 00:06:15
% DurationCPUTime: 0.45s
% Computational Cost: add. (344->133), mult. (604->162), div. (0->0), fcn. (717->12), ass. (0->53)
t38 = sin(qJ(4));
t52 = t38 * pkin(4) + pkin(8);
t43 = -pkin(10) - pkin(9);
t36 = sin(pkin(6));
t66 = pkin(7) * t36;
t64 = rSges(4,3) + pkin(8);
t34 = qJ(4) + qJ(5);
t27 = sin(t34);
t63 = pkin(5) * t27 + t52;
t62 = pkin(9) + rSges(5,3);
t41 = cos(qJ(4));
t26 = t41 * pkin(4) + pkin(3);
t61 = cos(qJ(3));
t39 = sin(qJ(3));
t60 = t36 * t39;
t40 = sin(qJ(2));
t59 = t36 * t40;
t42 = cos(qJ(2));
t58 = t36 * t42;
t57 = rSges(6,3) - t43;
t56 = rSges(7,3) + qJ(6) - t43;
t55 = cos(pkin(6));
t35 = sin(pkin(11));
t54 = t35 * pkin(1) + r_base(2);
t53 = qJ(1) + r_base(3);
t51 = t36 * t61;
t50 = t40 * t55;
t49 = t42 * t55;
t37 = cos(pkin(11));
t48 = t37 * pkin(1) + t35 * t66 + r_base(1);
t47 = t55 * pkin(7) + t53;
t16 = -t35 * t50 + t37 * t42;
t46 = t16 * pkin(2) + t48;
t45 = pkin(2) * t59 + t47;
t14 = t35 * t42 + t37 * t50;
t44 = t14 * pkin(2) - t37 * t66 + t54;
t28 = cos(t34);
t19 = pkin(5) * t28 + t26;
t18 = t55 * t39 + t40 * t51;
t17 = t39 * t59 - t55 * t61;
t15 = t35 * t49 + t37 * t40;
t13 = t35 * t40 - t37 * t49;
t10 = t16 * t61 + t35 * t60;
t9 = t16 * t39 - t35 * t51;
t8 = t14 * t61 - t37 * t60;
t7 = t14 * t39 + t37 * t51;
t6 = t18 * t28 - t27 * t58;
t5 = -t18 * t27 - t28 * t58;
t4 = t10 * t28 + t15 * t27;
t3 = -t10 * t27 + t15 * t28;
t2 = t13 * t27 + t8 * t28;
t1 = t13 * t28 - t8 * t27;
t11 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t37 * rSges(2,1) - t35 * rSges(2,2) + r_base(1)) + g(2) * (t35 * rSges(2,1) + t37 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t53)) - m(3) * (g(1) * (t16 * rSges(3,1) - t15 * rSges(3,2) + t48) + g(2) * (t14 * rSges(3,1) - t13 * rSges(3,2) + t54) + g(3) * (t55 * rSges(3,3) + t47) + (g(1) * rSges(3,3) * t35 + g(3) * (rSges(3,1) * t40 + rSges(3,2) * t42) + g(2) * (-rSges(3,3) - pkin(7)) * t37) * t36) - m(4) * (g(1) * (t10 * rSges(4,1) - t9 * rSges(4,2) + t64 * t15 + t46) + g(2) * (t8 * rSges(4,1) - t7 * rSges(4,2) + t64 * t13 + t44) + g(3) * (t18 * rSges(4,1) - t17 * rSges(4,2) - t64 * t58 + t45)) - m(5) * (g(1) * (t10 * pkin(3) + t15 * pkin(8) + (t10 * t41 + t15 * t38) * rSges(5,1) + (-t10 * t38 + t15 * t41) * rSges(5,2) + t62 * t9 + t46) + g(2) * (t8 * pkin(3) + t13 * pkin(8) + (t13 * t38 + t8 * t41) * rSges(5,1) + (t13 * t41 - t8 * t38) * rSges(5,2) + t62 * t7 + t44) + g(3) * (t18 * pkin(3) - pkin(8) * t58 + (t18 * t41 - t38 * t58) * rSges(5,1) + (-t18 * t38 - t41 * t58) * rSges(5,2) + t62 * t17 + t45)) - m(6) * (g(1) * (t4 * rSges(6,1) + t3 * rSges(6,2) + t10 * t26 + t52 * t15 + t57 * t9 + t46) + g(2) * (t2 * rSges(6,1) + t1 * rSges(6,2) + t52 * t13 + t8 * t26 + t57 * t7 + t44) + g(3) * (t6 * rSges(6,1) + t5 * rSges(6,2) + t57 * t17 + t18 * t26 - t52 * t58 + t45)) - m(7) * (g(1) * (t4 * rSges(7,1) + t3 * rSges(7,2) + t10 * t19 + t63 * t15 + t56 * t9 + t46) + g(2) * (t2 * rSges(7,1) + t1 * rSges(7,2) + t63 * t13 + t8 * t19 + t56 * t7 + t44) + g(3) * (t6 * rSges(7,1) + t5 * rSges(7,2) + t56 * t17 + t18 * t19 - t63 * t58 + t45));
U  = t11;
