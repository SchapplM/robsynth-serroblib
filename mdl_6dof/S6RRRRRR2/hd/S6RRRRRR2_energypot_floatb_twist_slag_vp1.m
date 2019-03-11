% Calculate potential energy for
% S6RRRRRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 03:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRRR2_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR2_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRRRR2_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR2_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR2_energypot_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR2_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRR2_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:33:05
% EndTime: 2019-03-10 03:33:06
% DurationCPUTime: 0.51s
% Computational Cost: add. (256->104), mult. (186->115), div. (0->0), fcn. (166->12), ass. (0->46)
t56 = rSges(6,3) + pkin(10);
t55 = rSges(7,3) + pkin(11) + pkin(10);
t22 = sin(qJ(1));
t25 = cos(qJ(1));
t54 = g(1) * t25 + g(2) * t22;
t27 = -pkin(8) - pkin(7);
t19 = qJ(2) + qJ(3);
t14 = qJ(4) + t19;
t6 = sin(t14);
t53 = rSges(5,2) * t6;
t7 = cos(t14);
t50 = t22 * t7;
t49 = t25 * t7;
t48 = rSges(3,3) + pkin(7);
t24 = cos(qJ(2));
t9 = t24 * pkin(2) + pkin(1);
t18 = qJ(5) + qJ(6);
t10 = sin(t18);
t46 = t22 * t10;
t12 = cos(t18);
t45 = t22 * t12;
t20 = sin(qJ(5));
t44 = t22 * t20;
t23 = cos(qJ(5));
t43 = t22 * t23;
t42 = t25 * t10;
t41 = t25 * t12;
t40 = t25 * t20;
t39 = t25 * t23;
t38 = rSges(4,3) - t27;
t36 = pkin(6) + r_base(3);
t13 = cos(t19);
t3 = pkin(3) * t13 + t9;
t35 = t25 * t3 + r_base(1);
t17 = -pkin(9) + t27;
t34 = t25 * t17 + t22 * t3 + r_base(2);
t21 = sin(qJ(2));
t33 = t21 * pkin(2) + t36;
t11 = sin(t19);
t32 = pkin(3) * t11 + t33;
t31 = -t22 * t17 + t35;
t30 = rSges(3,1) * t24 - rSges(3,2) * t21 + pkin(1);
t29 = rSges(4,1) * t13 - rSges(4,2) * t11 + t9;
t28 = g(1) * r_base(1) + g(2) * r_base(2);
t8 = t23 * pkin(5) + pkin(4);
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t25 * rSges(2,1) - t22 * rSges(2,2) + r_base(1)) + g(2) * (t22 * rSges(2,1) + t25 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t36)) - m(3) * (g(3) * (t21 * rSges(3,1) + t24 * rSges(3,2) + t36) + (g(1) * t30 - g(2) * t48) * t25 + (g(1) * t48 + g(2) * t30) * t22 + t28) - m(4) * (g(3) * (t11 * rSges(4,1) + t13 * rSges(4,2) + t33) + (g(1) * t29 - g(2) * t38) * t25 + (g(1) * t38 + g(2) * t29) * t22 + t28) - m(5) * (g(1) * (rSges(5,1) * t49 - t25 * t53 + t35) + g(2) * (-t25 * rSges(5,3) + t34) + g(3) * (t6 * rSges(5,1) + t7 * rSges(5,2) + t32) + (g(1) * (rSges(5,3) - t17) + g(2) * (rSges(5,1) * t7 - t53)) * t22) - m(6) * (g(1) * (pkin(4) * t49 + (t7 * t39 + t44) * rSges(6,1) + (-t7 * t40 + t43) * rSges(6,2) + t31) + g(2) * (pkin(4) * t50 + (t7 * t43 - t40) * rSges(6,1) + (-t7 * t44 - t39) * rSges(6,2) + t34) + g(3) * (-t56 * t7 + t32) + (g(3) * (rSges(6,1) * t23 - rSges(6,2) * t20 + pkin(4)) + t54 * t56) * t6) - m(7) * (g(1) * (t8 * t49 + pkin(5) * t44 + (t7 * t41 + t46) * rSges(7,1) + (-t7 * t42 + t45) * rSges(7,2) + t31) + g(2) * (t8 * t50 - pkin(5) * t40 + (t7 * t45 - t42) * rSges(7,1) + (-t7 * t46 - t41) * rSges(7,2) + t34) + g(3) * (-t55 * t7 + t32) + (g(3) * (rSges(7,1) * t12 - rSges(7,2) * t10 + t8) + t54 * t55) * t6);
U  = t1;
