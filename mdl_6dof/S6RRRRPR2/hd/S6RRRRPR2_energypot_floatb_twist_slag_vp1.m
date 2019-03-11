% Calculate potential energy for
% S6RRRRPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 22:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRPR2_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR2_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRRPR2_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR2_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR2_energypot_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR2_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR2_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:57:18
% EndTime: 2019-03-09 21:57:19
% DurationCPUTime: 0.52s
% Computational Cost: add. (256->104), mult. (186->115), div. (0->0), fcn. (166->12), ass. (0->46)
t56 = rSges(7,3) + pkin(10) + qJ(5);
t24 = sin(qJ(1));
t26 = cos(qJ(1));
t55 = g(1) * t26 + g(2) * t24;
t54 = rSges(6,3) + qJ(5);
t27 = -pkin(8) - pkin(7);
t19 = qJ(2) + qJ(3);
t14 = qJ(4) + t19;
t7 = sin(t14);
t53 = rSges(5,2) * t7;
t8 = cos(t14);
t50 = t24 * t8;
t49 = t26 * t8;
t48 = rSges(3,3) + pkin(7);
t25 = cos(qJ(2));
t9 = t25 * pkin(2) + pkin(1);
t17 = pkin(11) + qJ(6);
t10 = sin(t17);
t47 = t24 * t10;
t11 = cos(t17);
t46 = t24 * t11;
t20 = sin(pkin(11));
t45 = t24 * t20;
t21 = cos(pkin(11));
t44 = t24 * t21;
t43 = t26 * t10;
t42 = t26 * t11;
t41 = t26 * t20;
t40 = t26 * t21;
t39 = rSges(4,3) - t27;
t36 = pkin(6) + r_base(3);
t13 = cos(t19);
t3 = pkin(3) * t13 + t9;
t35 = t26 * t3 + r_base(1);
t18 = -pkin(9) + t27;
t34 = t26 * t18 + t24 * t3 + r_base(2);
t23 = sin(qJ(2));
t33 = t23 * pkin(2) + t36;
t12 = sin(t19);
t32 = pkin(3) * t12 + t33;
t31 = -t24 * t18 + t35;
t30 = rSges(3,1) * t25 - rSges(3,2) * t23 + pkin(1);
t29 = rSges(4,1) * t13 - rSges(4,2) * t12 + t9;
t28 = g(1) * r_base(1) + g(2) * r_base(2);
t5 = t21 * pkin(5) + pkin(4);
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t26 * rSges(2,1) - t24 * rSges(2,2) + r_base(1)) + g(2) * (t24 * rSges(2,1) + t26 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t36)) - m(3) * (g(3) * (t23 * rSges(3,1) + t25 * rSges(3,2) + t36) + (g(1) * t30 - g(2) * t48) * t26 + (g(1) * t48 + g(2) * t30) * t24 + t28) - m(4) * (g(3) * (t12 * rSges(4,1) + t13 * rSges(4,2) + t33) + (g(1) * t29 - g(2) * t39) * t26 + (g(1) * t39 + g(2) * t29) * t24 + t28) - m(5) * (g(1) * (rSges(5,1) * t49 - t26 * t53 + t35) + g(2) * (-t26 * rSges(5,3) + t34) + g(3) * (t7 * rSges(5,1) + t8 * rSges(5,2) + t32) + (g(1) * (rSges(5,3) - t18) + g(2) * (rSges(5,1) * t8 - t53)) * t24) - m(6) * (g(1) * (pkin(4) * t49 + (t8 * t40 + t45) * rSges(6,1) + (-t8 * t41 + t44) * rSges(6,2) + t31) + g(2) * (pkin(4) * t50 + (t8 * t44 - t41) * rSges(6,1) + (-t8 * t45 - t40) * rSges(6,2) + t34) + g(3) * (-t54 * t8 + t32) + (g(3) * (rSges(6,1) * t21 - rSges(6,2) * t20 + pkin(4)) + t55 * t54) * t7) - m(7) * (g(1) * (t5 * t49 + pkin(5) * t45 + (t8 * t42 + t47) * rSges(7,1) + (-t8 * t43 + t46) * rSges(7,2) + t31) + g(2) * (t5 * t50 - pkin(5) * t41 + (t8 * t46 - t43) * rSges(7,1) + (-t8 * t47 - t42) * rSges(7,2) + t34) + g(3) * (-t56 * t8 + t32) + (g(3) * (rSges(7,1) * t11 - rSges(7,2) * t10 + t5) + t55 * t56) * t7);
U  = t1;
