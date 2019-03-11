% Calculate potential energy for
% S6RPRRRR4
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
% Datum: 2019-03-09 07:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRR4_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR4_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRRRR4_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR4_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR4_energypot_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR4_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR4_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:04:42
% EndTime: 2019-03-09 07:04:43
% DurationCPUTime: 0.38s
% Computational Cost: add. (252->96), mult. (162->102), div. (0->0), fcn. (138->12), ass. (0->42)
t50 = rSges(7,3) + pkin(10);
t22 = pkin(11) + qJ(3);
t17 = qJ(4) + t22;
t14 = qJ(5) + t17;
t7 = sin(t14);
t49 = rSges(6,2) * t7;
t29 = cos(qJ(1));
t8 = cos(t14);
t48 = t29 * t8;
t24 = cos(pkin(11));
t13 = t24 * pkin(2) + pkin(1);
t26 = sin(qJ(6));
t27 = sin(qJ(1));
t46 = t27 * t26;
t28 = cos(qJ(6));
t45 = t27 * t28;
t44 = t29 * t26;
t43 = t29 * t28;
t25 = -pkin(7) - qJ(2);
t42 = rSges(4,3) - t25;
t21 = -pkin(8) + t25;
t41 = rSges(5,3) - t21;
t40 = rSges(3,3) + qJ(2);
t39 = pkin(6) + r_base(3);
t12 = cos(t17);
t16 = cos(t22);
t4 = pkin(3) * t16 + t13;
t3 = pkin(4) * t12 + t4;
t38 = t29 * t3 + r_base(1);
t18 = -pkin(9) + t21;
t37 = t29 * t18 + t27 * t3 + r_base(2);
t23 = sin(pkin(11));
t36 = t23 * pkin(2) + t39;
t15 = sin(t22);
t35 = pkin(3) * t15 + t36;
t11 = sin(t17);
t34 = pkin(4) * t11 + t35;
t33 = rSges(3,1) * t24 - rSges(3,2) * t23 + pkin(1);
t32 = rSges(5,1) * t12 - rSges(5,2) * t11 + t4;
t31 = rSges(4,1) * t16 - rSges(4,2) * t15 + t13;
t30 = g(1) * r_base(1) + g(2) * r_base(2);
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t29 * rSges(2,1) - t27 * rSges(2,2) + r_base(1)) + g(2) * (t27 * rSges(2,1) + t29 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t39)) - m(3) * (g(3) * (t23 * rSges(3,1) + t24 * rSges(3,2) + t39) + (g(1) * t33 - g(2) * t40) * t29 + (g(1) * t40 + g(2) * t33) * t27 + t30) - m(4) * (g(3) * (t15 * rSges(4,1) + t16 * rSges(4,2) + t36) + (g(1) * t31 - g(2) * t42) * t29 + (g(1) * t42 + g(2) * t31) * t27 + t30) - m(5) * (g(3) * (t11 * rSges(5,1) + t12 * rSges(5,2) + t35) + (g(1) * t32 - g(2) * t41) * t29 + (g(1) * t41 + g(2) * t32) * t27 + t30) - m(6) * (g(1) * (rSges(6,1) * t48 - t29 * t49 + t38) + g(2) * (-t29 * rSges(6,3) + t37) + g(3) * (t7 * rSges(6,1) + t8 * rSges(6,2) + t34) + (g(1) * (rSges(6,3) - t18) + g(2) * (rSges(6,1) * t8 - t49)) * t27) - m(7) * (g(1) * (pkin(5) * t48 - t27 * t18 + (t43 * t8 + t46) * rSges(7,1) + (-t44 * t8 + t45) * rSges(7,2) + t38) + g(2) * (t27 * t8 * pkin(5) + (t45 * t8 - t44) * rSges(7,1) + (-t46 * t8 - t43) * rSges(7,2) + t37) + g(3) * (-t50 * t8 + t34) + (g(3) * (rSges(7,1) * t28 - rSges(7,2) * t26 + pkin(5)) + (g(1) * t29 + g(2) * t27) * t50) * t7);
U  = t1;
