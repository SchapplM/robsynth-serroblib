% Calculate potential energy for
% S6RRPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
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
% Datum: 2019-03-09 12:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRP7_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP7_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPRRP7_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP7_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP7_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP7_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP7_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:16:39
% EndTime: 2019-03-09 12:16:39
% DurationCPUTime: 0.37s
% Computational Cost: add. (207->106), mult. (351->120), div. (0->0), fcn. (377->8), ass. (0->40)
t56 = rSges(7,2) + pkin(9);
t57 = rSges(7,1) + pkin(5);
t55 = rSges(6,3) + pkin(9);
t31 = sin(qJ(2));
t32 = sin(qJ(1));
t54 = t31 * t32;
t34 = cos(qJ(4));
t53 = t31 * t34;
t35 = cos(qJ(2));
t52 = t32 * t35;
t36 = cos(qJ(1));
t51 = t35 * t36;
t50 = qJ(3) * t31;
t49 = rSges(7,3) + qJ(6);
t48 = pkin(6) + r_base(3);
t47 = t32 * pkin(1) + r_base(2);
t46 = t31 * pkin(2) + t48;
t45 = t36 * pkin(1) + t32 * pkin(7) + r_base(1);
t44 = pkin(2) * t52 + t32 * t50 + t47;
t30 = sin(qJ(4));
t12 = t30 * t31 + t34 * t35;
t43 = pkin(2) * t51 + t36 * t50 + t45;
t42 = pkin(3) * t52 + t36 * pkin(8) + t44;
t41 = pkin(3) * t51 + t43;
t40 = t31 * pkin(3) - t35 * qJ(3) + t46;
t13 = -t30 * t35 + t53;
t39 = t13 * pkin(4) + t40;
t8 = t12 * t32;
t38 = t8 * pkin(4) - t36 * pkin(7) + t42;
t10 = t12 * t36;
t37 = t10 * pkin(4) - t32 * pkin(8) + t41;
t33 = cos(qJ(5));
t29 = sin(qJ(5));
t9 = t30 * t51 - t36 * t53;
t7 = t30 * t52 - t32 * t53;
t4 = t10 * t33 - t29 * t32;
t3 = t10 * t29 + t32 * t33;
t2 = t29 * t36 + t8 * t33;
t1 = t29 * t8 - t36 * t33;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t36 - t32 * rSges(2,2) + r_base(1)) + g(2) * (t32 * rSges(2,1) + rSges(2,2) * t36 + r_base(2)) + g(3) * (rSges(2,3) + t48)) - m(3) * (g(1) * (t32 * rSges(3,3) + t45) + g(2) * (rSges(3,1) * t52 - rSges(3,2) * t54 + t47) + g(3) * (rSges(3,1) * t31 + rSges(3,2) * t35 + t48) + (g(1) * (rSges(3,1) * t35 - rSges(3,2) * t31) + g(2) * (-rSges(3,3) - pkin(7))) * t36) - m(4) * (g(1) * (t32 * rSges(4,2) + t43) + g(2) * (rSges(4,1) * t52 + rSges(4,3) * t54 + t44) + g(3) * (t31 * rSges(4,1) + (-rSges(4,3) - qJ(3)) * t35 + t46) + (g(1) * (rSges(4,1) * t35 + rSges(4,3) * t31) + g(2) * (-rSges(4,2) - pkin(7))) * t36) - m(5) * (g(1) * (t10 * rSges(5,1) - t9 * rSges(5,2) + (-rSges(5,3) - pkin(8)) * t32 + t41) + g(2) * (t8 * rSges(5,1) - t7 * rSges(5,2) + (rSges(5,3) - pkin(7)) * t36 + t42) + g(3) * (rSges(5,1) * t13 - rSges(5,2) * t12 + t40)) - m(6) * (g(1) * (t4 * rSges(6,1) - t3 * rSges(6,2) + t55 * t9 + t37) + g(2) * (t2 * rSges(6,1) - t1 * rSges(6,2) + t55 * t7 + t38) + g(3) * ((rSges(6,1) * t33 - rSges(6,2) * t29) * t13 + t55 * t12 + t39)) - m(7) * (g(1) * (t49 * t3 + t57 * t4 + t56 * t9 + t37) + g(2) * (t49 * t1 + t57 * t2 + t56 * t7 + t38) + (t39 + (t49 * t29 + t57 * t33) * t13 + t56 * t12) * g(3));
U  = t5;
