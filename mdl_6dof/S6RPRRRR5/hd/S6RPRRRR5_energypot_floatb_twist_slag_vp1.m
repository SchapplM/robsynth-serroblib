% Calculate potential energy for
% S6RPRRRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 07:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRR5_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR5_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRRRR5_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR5_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR5_energypot_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR5_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR5_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:08:06
% EndTime: 2019-03-09 07:08:07
% DurationCPUTime: 0.48s
% Computational Cost: add. (256->104), mult. (186->115), div. (0->0), fcn. (166->12), ass. (0->46)
t56 = rSges(6,3) + pkin(9);
t55 = rSges(7,3) + pkin(10) + pkin(9);
t24 = sin(qJ(1));
t26 = cos(qJ(1));
t54 = g(1) * t26 + g(2) * t24;
t18 = pkin(11) + qJ(3);
t12 = qJ(4) + t18;
t6 = sin(t12);
t53 = rSges(5,2) * t6;
t7 = cos(t12);
t50 = t24 * t7;
t49 = t26 * t7;
t21 = cos(pkin(11));
t8 = t21 * pkin(2) + pkin(1);
t19 = qJ(5) + qJ(6);
t13 = sin(t19);
t47 = t24 * t13;
t14 = cos(t19);
t46 = t24 * t14;
t23 = sin(qJ(5));
t45 = t24 * t23;
t25 = cos(qJ(5));
t44 = t24 * t25;
t43 = t26 * t13;
t42 = t26 * t14;
t41 = t26 * t23;
t40 = t26 * t25;
t22 = -pkin(7) - qJ(2);
t39 = rSges(4,3) - t22;
t37 = rSges(3,3) + qJ(2);
t36 = pkin(6) + r_base(3);
t11 = cos(t18);
t3 = pkin(3) * t11 + t8;
t35 = t26 * t3 + r_base(1);
t17 = -pkin(8) + t22;
t34 = t26 * t17 + t24 * t3 + r_base(2);
t20 = sin(pkin(11));
t33 = t20 * pkin(2) + t36;
t10 = sin(t18);
t32 = pkin(3) * t10 + t33;
t31 = -t24 * t17 + t35;
t30 = rSges(3,1) * t21 - rSges(3,2) * t20 + pkin(1);
t29 = rSges(4,1) * t11 - rSges(4,2) * t10 + t8;
t28 = g(1) * r_base(1) + g(2) * r_base(2);
t9 = t25 * pkin(5) + pkin(4);
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t26 * rSges(2,1) - t24 * rSges(2,2) + r_base(1)) + g(2) * (t24 * rSges(2,1) + t26 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t36)) - m(3) * (g(3) * (t20 * rSges(3,1) + t21 * rSges(3,2) + t36) + (g(1) * t30 - g(2) * t37) * t26 + (g(1) * t37 + g(2) * t30) * t24 + t28) - m(4) * (g(3) * (t10 * rSges(4,1) + t11 * rSges(4,2) + t33) + (g(1) * t29 - g(2) * t39) * t26 + (g(1) * t39 + g(2) * t29) * t24 + t28) - m(5) * (g(1) * (rSges(5,1) * t49 - t26 * t53 + t35) + g(2) * (-t26 * rSges(5,3) + t34) + g(3) * (t6 * rSges(5,1) + t7 * rSges(5,2) + t32) + (g(1) * (rSges(5,3) - t17) + g(2) * (rSges(5,1) * t7 - t53)) * t24) - m(6) * (g(1) * (pkin(4) * t49 + (t7 * t40 + t45) * rSges(6,1) + (-t7 * t41 + t44) * rSges(6,2) + t31) + g(2) * (pkin(4) * t50 + (t7 * t44 - t41) * rSges(6,1) + (-t7 * t45 - t40) * rSges(6,2) + t34) + g(3) * (-t56 * t7 + t32) + (g(3) * (rSges(6,1) * t25 - rSges(6,2) * t23 + pkin(4)) + t54 * t56) * t6) - m(7) * (g(1) * (t9 * t49 + pkin(5) * t45 + (t7 * t42 + t47) * rSges(7,1) + (-t7 * t43 + t46) * rSges(7,2) + t31) + g(2) * (t9 * t50 - pkin(5) * t41 + (t7 * t46 - t43) * rSges(7,1) + (-t7 * t47 - t42) * rSges(7,2) + t34) + g(3) * (-t55 * t7 + t32) + (g(3) * (rSges(7,1) * t14 - rSges(7,2) * t13 + t9) + t54 * t55) * t6);
U  = t1;
