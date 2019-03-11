% Calculate potential energy for
% S6RRPRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
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
% Datum: 2019-03-09 11:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRPR11_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR11_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPRPR11_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR11_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR11_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR11_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR11_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:11:35
% EndTime: 2019-03-09 11:11:35
% DurationCPUTime: 0.58s
% Computational Cost: add. (207->118), mult. (237->135), div. (0->0), fcn. (221->10), ass. (0->40)
t54 = rSges(5,3) + pkin(8);
t20 = -qJ(5) - pkin(8);
t53 = rSges(6,3) - t20;
t52 = rSges(7,3) + pkin(9) - t20;
t23 = sin(qJ(1));
t26 = cos(qJ(1));
t51 = g(1) * t26 + g(2) * t23;
t21 = sin(qJ(4));
t50 = pkin(4) * t21;
t24 = cos(qJ(4));
t9 = t24 * pkin(4) + pkin(3);
t22 = sin(qJ(2));
t46 = t22 * t23;
t45 = t22 * t26;
t19 = qJ(4) + pkin(10);
t10 = sin(t19);
t44 = t23 * t10;
t11 = cos(t19);
t43 = t23 * t11;
t42 = t23 * t21;
t41 = t23 * t24;
t25 = cos(qJ(2));
t40 = t23 * t25;
t39 = t24 * t26;
t36 = qJ(3) * t22;
t35 = pkin(6) + r_base(3);
t34 = t23 * pkin(1) + r_base(2);
t33 = t21 * t45;
t32 = t22 * t42;
t31 = t22 * pkin(2) + t35;
t30 = t26 * pkin(1) + t23 * pkin(7) + r_base(1);
t29 = pkin(2) * t40 + t23 * t36 + t34;
t28 = t30 + (pkin(2) * t25 + t36) * t26;
t27 = -t26 * pkin(7) + t29;
t12 = qJ(6) + t19;
t8 = cos(t12);
t7 = sin(t12);
t2 = pkin(5) * t10 + t50;
t1 = pkin(5) * t11 + t9;
t3 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t26 - t23 * rSges(2,2) + r_base(1)) + g(2) * (t23 * rSges(2,1) + rSges(2,2) * t26 + r_base(2)) + g(3) * (rSges(2,3) + t35)) - m(3) * (g(1) * (t23 * rSges(3,3) + t30) + g(2) * (rSges(3,1) * t40 - rSges(3,2) * t46 + t34) + g(3) * (rSges(3,1) * t22 + rSges(3,2) * t25 + t35) + (g(1) * (rSges(3,1) * t25 - rSges(3,2) * t22) + g(2) * (-rSges(3,3) - pkin(7))) * t26) - m(4) * (g(1) * (t23 * rSges(4,1) + t28) + g(2) * (-rSges(4,2) * t40 + rSges(4,3) * t46 + t29) + g(3) * (-rSges(4,2) * t22 + (-rSges(4,3) - qJ(3)) * t25 + t31) + (g(1) * (-rSges(4,2) * t25 + rSges(4,3) * t22) + g(2) * (-rSges(4,1) - pkin(7))) * t26) - m(5) * (g(1) * (t23 * pkin(3) + (t33 + t41) * rSges(5,1) + (t22 * t39 - t42) * rSges(5,2) + t28) + g(2) * (-t26 * pkin(3) + (t32 - t39) * rSges(5,1) + (t21 * t26 + t22 * t41) * rSges(5,2) + t27) + g(3) * (t54 * t22 + t31) + (g(3) * (-rSges(5,1) * t21 - rSges(5,2) * t24 - qJ(3)) + t51 * t54) * t25) - m(6) * (g(1) * (t23 * t9 + pkin(4) * t33 + (t10 * t45 + t43) * rSges(6,1) + (t11 * t45 - t44) * rSges(6,2) + t28) + g(2) * (-t26 * t9 + pkin(4) * t32 + (-t11 * t26 + t22 * t44) * rSges(6,1) + (t10 * t26 + t22 * t43) * rSges(6,2) + t27) + g(3) * (t53 * t22 + t31) + (g(3) * (-rSges(6,1) * t10 - rSges(6,2) * t11 - qJ(3) - t50) + t51 * t53) * t25) - m(7) * (g(1) * (t23 * t1 + t2 * t45 + (t23 * t8 + t7 * t45) * rSges(7,1) + (-t23 * t7 + t8 * t45) * rSges(7,2) + t28) + g(2) * (-t26 * t1 + t2 * t46 + (-t26 * t8 + t7 * t46) * rSges(7,1) + (t26 * t7 + t8 * t46) * rSges(7,2) + t27) + g(3) * (t52 * t22 + t31) + (g(3) * (-rSges(7,1) * t7 - rSges(7,2) * t8 - qJ(3) - t2) + t51 * t52) * t25);
U  = t3;
