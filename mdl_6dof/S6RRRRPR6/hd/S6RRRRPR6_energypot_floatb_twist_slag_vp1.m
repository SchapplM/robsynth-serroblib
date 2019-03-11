% Calculate potential energy for
% S6RRRRPR6
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
% Datum: 2019-03-09 22:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRPR6_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR6_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRRPR6_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR6_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR6_energypot_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR6_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR6_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:18:32
% EndTime: 2019-03-09 22:18:33
% DurationCPUTime: 0.62s
% Computational Cost: add. (266->123), mult. (240->142), div. (0->0), fcn. (228->12), ass. (0->37)
t50 = rSges(4,3) + pkin(8);
t30 = -pkin(9) - pkin(8);
t49 = rSges(5,3) - t30;
t22 = -qJ(5) + t30;
t48 = rSges(6,3) - t22;
t47 = rSges(7,3) + pkin(10) - t22;
t26 = sin(qJ(1));
t29 = cos(qJ(1));
t46 = g(1) * t29 + g(2) * t26;
t27 = cos(qJ(3));
t12 = t27 * pkin(3) + pkin(2);
t23 = qJ(3) + qJ(4);
t14 = sin(t23);
t24 = sin(qJ(3));
t4 = t24 * pkin(3) + pkin(4) * t14;
t25 = sin(qJ(2));
t42 = rSges(3,2) * t25;
t41 = t24 * t26;
t40 = t24 * t29;
t28 = cos(qJ(2));
t39 = t26 * t28;
t38 = t28 * t29;
t34 = pkin(6) + r_base(3);
t33 = t26 * pkin(1) + r_base(2);
t15 = cos(t23);
t3 = pkin(4) * t15 + t12;
t13 = pkin(11) + t23;
t32 = t29 * pkin(1) + t26 * pkin(7) + r_base(1);
t31 = -t29 * pkin(7) + t33;
t11 = qJ(6) + t13;
t8 = cos(t13);
t7 = sin(t13);
t6 = cos(t11);
t5 = sin(t11);
t2 = pkin(5) * t7 + t4;
t1 = pkin(5) * t8 + t3;
t9 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t29 - rSges(2,2) * t26 + r_base(1)) + g(2) * (rSges(2,1) * t26 + rSges(2,2) * t29 + r_base(2)) + g(3) * (rSges(2,3) + t34)) - m(3) * (g(1) * (rSges(3,3) * t26 + t32) + g(2) * (rSges(3,1) * t39 - t26 * t42 + t33) + g(3) * (rSges(3,1) * t25 + rSges(3,2) * t28 + t34) + (g(1) * (rSges(3,1) * t28 - t42) + g(2) * (-rSges(3,3) - pkin(7))) * t29) - m(4) * (g(1) * (pkin(2) * t38 + (t27 * t38 + t41) * rSges(4,1) + (-t24 * t38 + t26 * t27) * rSges(4,2) + t32) + g(2) * (pkin(2) * t39 + (t27 * t39 - t40) * rSges(4,1) + (-t24 * t39 - t27 * t29) * rSges(4,2) + t31) + g(3) * (-t50 * t28 + t34) + (g(3) * (rSges(4,1) * t27 - rSges(4,2) * t24 + pkin(2)) + t46 * t50) * t25) - m(5) * (g(1) * (t12 * t38 + pkin(3) * t41 + (t14 * t26 + t15 * t38) * rSges(5,1) + (-t14 * t38 + t15 * t26) * rSges(5,2) + t32) + g(2) * (t12 * t39 - pkin(3) * t40 + (-t14 * t29 + t15 * t39) * rSges(5,1) + (-t14 * t39 - t15 * t29) * rSges(5,2) + t31) + g(3) * (-t49 * t28 + t34) + (g(3) * (rSges(5,1) * t15 - rSges(5,2) * t14 + t12) + t46 * t49) * t25) - m(6) * (g(1) * (t3 * t38 + t26 * t4 + (t26 * t7 + t8 * t38) * rSges(6,1) + (t26 * t8 - t7 * t38) * rSges(6,2) + t32) + g(2) * (t3 * t39 - t29 * t4 + (-t29 * t7 + t8 * t39) * rSges(6,1) + (-t29 * t8 - t7 * t39) * rSges(6,2) + t31) + g(3) * (-t48 * t28 + t34) + (g(3) * (rSges(6,1) * t8 - rSges(6,2) * t7 + t3) + t46 * t48) * t25) - m(7) * (g(1) * (t1 * t38 + t26 * t2 + (t26 * t5 + t6 * t38) * rSges(7,1) + (t26 * t6 - t5 * t38) * rSges(7,2) + t32) + g(2) * (t1 * t39 - t29 * t2 + (-t29 * t5 + t6 * t39) * rSges(7,1) + (-t29 * t6 - t5 * t39) * rSges(7,2) + t31) + g(3) * (-t47 * t28 + t34) + (g(3) * (rSges(7,1) * t6 - rSges(7,2) * t5 + t1) + t46 * t47) * t25);
U  = t9;
