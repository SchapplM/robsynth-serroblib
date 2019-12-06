% Calculate potential energy for
% S5PRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRRP7_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP7_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PRRRP7_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP7_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP7_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP7_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRP7_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:54:01
% EndTime: 2019-12-05 16:54:01
% DurationCPUTime: 0.32s
% Computational Cost: add. (225->104), mult. (435->129), div. (0->0), fcn. (511->10), ass. (0->46)
t29 = sin(pkin(5));
t56 = pkin(6) * t29;
t55 = rSges(4,3) + pkin(7);
t54 = rSges(5,3) + pkin(8);
t34 = sin(qJ(3));
t53 = t29 * t34;
t35 = sin(qJ(2));
t52 = t29 * t35;
t37 = cos(qJ(3));
t51 = t29 * t37;
t38 = cos(qJ(2));
t50 = t29 * t38;
t31 = cos(pkin(5));
t49 = t31 * t35;
t48 = t31 * t38;
t47 = rSges(6,3) + qJ(5) + pkin(8);
t28 = sin(pkin(9));
t46 = t28 * pkin(1) + r_base(2);
t45 = qJ(1) + r_base(3);
t33 = sin(qJ(4));
t44 = pkin(4) * t33 + pkin(7);
t30 = cos(pkin(9));
t43 = t30 * pkin(1) + t28 * t56 + r_base(1);
t42 = t31 * pkin(6) + t45;
t16 = -t28 * t49 + t30 * t38;
t41 = t16 * pkin(2) + t43;
t40 = pkin(2) * t52 + t42;
t14 = t28 * t38 + t30 * t49;
t39 = t14 * pkin(2) - t30 * t56 + t46;
t36 = cos(qJ(4));
t24 = t36 * pkin(4) + pkin(3);
t18 = t31 * t34 + t35 * t51;
t17 = -t31 * t37 + t34 * t52;
t15 = t28 * t48 + t30 * t35;
t13 = t28 * t35 - t30 * t48;
t10 = t18 * t36 - t33 * t50;
t9 = -t18 * t33 - t36 * t50;
t8 = t16 * t37 + t28 * t53;
t7 = t16 * t34 - t28 * t51;
t6 = t14 * t37 - t30 * t53;
t5 = t14 * t34 + t30 * t51;
t4 = t15 * t33 + t8 * t36;
t3 = t15 * t36 - t8 * t33;
t2 = t13 * t33 + t6 * t36;
t1 = t13 * t36 - t6 * t33;
t11 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t30 * rSges(2,1) - t28 * rSges(2,2) + r_base(1)) + g(2) * (t28 * rSges(2,1) + t30 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t45)) - m(3) * (g(1) * (t16 * rSges(3,1) - t15 * rSges(3,2) + t43) + g(2) * (t14 * rSges(3,1) - t13 * rSges(3,2) + t46) + g(3) * (t31 * rSges(3,3) + t42) + (g(1) * rSges(3,3) * t28 + g(3) * (rSges(3,1) * t35 + rSges(3,2) * t38) + g(2) * (-rSges(3,3) - pkin(6)) * t30) * t29) - m(4) * (g(1) * (t8 * rSges(4,1) - t7 * rSges(4,2) + t55 * t15 + t41) + g(2) * (t6 * rSges(4,1) - t5 * rSges(4,2) + t55 * t13 + t39) + g(3) * (t18 * rSges(4,1) - t17 * rSges(4,2) - t55 * t50 + t40)) - m(5) * (g(1) * (t4 * rSges(5,1) + t3 * rSges(5,2) + t8 * pkin(3) + t15 * pkin(7) + t54 * t7 + t41) + g(2) * (t2 * rSges(5,1) + t1 * rSges(5,2) + t6 * pkin(3) + t13 * pkin(7) + t54 * t5 + t39) + g(3) * (t10 * rSges(5,1) + t9 * rSges(5,2) + t18 * pkin(3) - pkin(7) * t50 + t54 * t17 + t40)) - m(6) * (g(1) * (t4 * rSges(6,1) + t3 * rSges(6,2) + t44 * t15 + t8 * t24 + t47 * t7 + t41) + g(2) * (t2 * rSges(6,1) + t1 * rSges(6,2) + t44 * t13 + t6 * t24 + t47 * t5 + t39) + g(3) * (t10 * rSges(6,1) + t9 * rSges(6,2) + t47 * t17 + t18 * t24 - t44 * t50 + t40));
U = t11;
