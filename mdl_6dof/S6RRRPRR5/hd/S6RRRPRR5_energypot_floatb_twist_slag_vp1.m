% Calculate potential energy for
% S6RRRPRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 18:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPRR5_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR5_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRPRR5_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR5_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR5_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR5_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR5_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:20:29
% EndTime: 2019-03-09 18:20:30
% DurationCPUTime: 0.47s
% Computational Cost: add. (227->108), mult. (204->122), div. (0->0), fcn. (184->10), ass. (0->42)
t55 = rSges(6,3) + pkin(9);
t54 = rSges(7,3) + pkin(10) + pkin(9);
t20 = sin(qJ(1));
t23 = cos(qJ(1));
t53 = g(1) * t23 + g(2) * t20;
t50 = rSges(3,3) + pkin(7);
t17 = qJ(2) + qJ(3);
t12 = sin(t17);
t48 = t12 * t23;
t16 = qJ(5) + qJ(6);
t11 = sin(t16);
t47 = t20 * t11;
t13 = cos(t16);
t46 = t20 * t13;
t18 = sin(qJ(5));
t45 = t20 * t18;
t21 = cos(qJ(5));
t44 = t20 * t21;
t43 = t23 * t11;
t42 = t23 * t13;
t14 = cos(t17);
t41 = t23 * t14;
t40 = t23 * t18;
t39 = t23 * t21;
t37 = qJ(4) * t12;
t36 = pkin(6) + r_base(3);
t22 = cos(qJ(2));
t9 = t22 * pkin(2) + pkin(1);
t35 = t23 * t9 + r_base(1);
t34 = t12 * t45;
t33 = t12 * t40;
t19 = sin(qJ(2));
t32 = t19 * pkin(2) + t36;
t25 = -pkin(8) - pkin(7);
t31 = t20 * t9 + t23 * t25 + r_base(2);
t30 = pkin(3) * t41 + t23 * t37 + t35;
t29 = t12 * pkin(3) + t32;
t28 = t31 + (pkin(3) * t14 + t37) * t20;
t27 = rSges(3,1) * t22 - rSges(3,2) * t19 + pkin(1);
t26 = -t20 * t25 + t30;
t8 = t21 * pkin(5) + pkin(4);
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t23 * rSges(2,1) - t20 * rSges(2,2) + r_base(1)) + g(2) * (t20 * rSges(2,1) + t23 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t36)) - m(3) * (g(1) * r_base(1) + g(2) * r_base(2) + g(3) * (t19 * rSges(3,1) + t22 * rSges(3,2) + t36) + (g(1) * t27 - g(2) * t50) * t23 + (g(1) * t50 + g(2) * t27) * t20) - m(4) * (g(1) * (rSges(4,1) * t41 - rSges(4,2) * t48 + t35) + g(2) * (-t23 * rSges(4,3) + t31) + g(3) * (t12 * rSges(4,1) + t14 * rSges(4,2) + t32) + (g(1) * (rSges(4,3) - t25) + g(2) * (rSges(4,1) * t14 - rSges(4,2) * t12)) * t20) - m(5) * (g(1) * (-rSges(5,2) * t41 + rSges(5,3) * t48 + t30) + g(2) * (-t23 * rSges(5,1) + t28) + g(3) * (-t12 * rSges(5,2) + (-rSges(5,3) - qJ(4)) * t14 + t29) + (g(1) * (rSges(5,1) - t25) + g(2) * (-rSges(5,2) * t14 + rSges(5,3) * t12)) * t20) - m(6) * (g(1) * (t20 * pkin(4) + (t33 + t44) * rSges(6,1) + (t12 * t39 - t45) * rSges(6,2) + t26) + g(2) * (-t23 * pkin(4) + (t34 - t39) * rSges(6,1) + (t12 * t44 + t40) * rSges(6,2) + t28) + g(3) * (t55 * t12 + t29) + (g(3) * (-rSges(6,1) * t18 - rSges(6,2) * t21 - qJ(4)) + t53 * t55) * t14) - m(7) * (g(1) * (t20 * t8 + pkin(5) * t33 + (t12 * t43 + t46) * rSges(7,1) + (t12 * t42 - t47) * rSges(7,2) + t26) + g(2) * (-t23 * t8 + pkin(5) * t34 + (t12 * t47 - t42) * rSges(7,1) + (t12 * t46 + t43) * rSges(7,2) + t28) + g(3) * (t54 * t12 + t29) + (g(3) * (-rSges(7,1) * t11 - rSges(7,2) * t13 - pkin(5) * t18 - qJ(4)) + t53 * t54) * t14);
U  = t1;
