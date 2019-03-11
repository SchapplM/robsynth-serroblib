% Calculate potential energy for
% S6RRRPRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 18:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPRR6_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR6_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRPRR6_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR6_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR6_energypot_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR6_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR6_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:24:59
% EndTime: 2019-03-09 18:25:00
% DurationCPUTime: 0.62s
% Computational Cost: add. (266->123), mult. (240->142), div. (0->0), fcn. (228->12), ass. (0->37)
t50 = rSges(4,3) + pkin(8);
t24 = -qJ(4) - pkin(8);
t49 = rSges(5,3) - t24;
t22 = -pkin(9) + t24;
t48 = rSges(6,3) - t22;
t47 = rSges(7,3) + pkin(10) - t22;
t27 = sin(qJ(1));
t30 = cos(qJ(1));
t46 = g(1) * t30 + g(2) * t27;
t28 = cos(qJ(3));
t12 = t28 * pkin(3) + pkin(2);
t23 = qJ(3) + pkin(11);
t13 = sin(t23);
t25 = sin(qJ(3));
t4 = t25 * pkin(3) + pkin(4) * t13;
t26 = sin(qJ(2));
t42 = rSges(3,2) * t26;
t41 = t25 * t30;
t40 = t27 * t25;
t29 = cos(qJ(2));
t39 = t27 * t29;
t38 = t29 * t30;
t34 = pkin(6) + r_base(3);
t14 = cos(t23);
t3 = pkin(4) * t14 + t12;
t33 = t27 * pkin(1) + r_base(2);
t15 = qJ(5) + t23;
t32 = t30 * pkin(1) + t27 * pkin(7) + r_base(1);
t31 = -t30 * pkin(7) + t33;
t11 = qJ(6) + t15;
t10 = cos(t15);
t9 = sin(t15);
t6 = cos(t11);
t5 = sin(t11);
t2 = pkin(5) * t9 + t4;
t1 = pkin(5) * t10 + t3;
t7 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t30 - t27 * rSges(2,2) + r_base(1)) + g(2) * (t27 * rSges(2,1) + rSges(2,2) * t30 + r_base(2)) + g(3) * (rSges(2,3) + t34)) - m(3) * (g(1) * (t27 * rSges(3,3) + t32) + g(2) * (rSges(3,1) * t39 - t27 * t42 + t33) + g(3) * (rSges(3,1) * t26 + rSges(3,2) * t29 + t34) + (g(1) * (rSges(3,1) * t29 - t42) + g(2) * (-rSges(3,3) - pkin(7))) * t30) - m(4) * (g(1) * (pkin(2) * t38 + (t28 * t38 + t40) * rSges(4,1) + (-t25 * t38 + t27 * t28) * rSges(4,2) + t32) + g(2) * (pkin(2) * t39 + (t28 * t39 - t41) * rSges(4,1) + (-t25 * t39 - t28 * t30) * rSges(4,2) + t31) + g(3) * (-t50 * t29 + t34) + (g(3) * (rSges(4,1) * t28 - rSges(4,2) * t25 + pkin(2)) + t46 * t50) * t26) - m(5) * (g(1) * (t12 * t38 + pkin(3) * t40 + (t27 * t13 + t14 * t38) * rSges(5,1) + (-t13 * t38 + t27 * t14) * rSges(5,2) + t32) + g(2) * (t12 * t39 - pkin(3) * t41 + (-t13 * t30 + t14 * t39) * rSges(5,1) + (-t13 * t39 - t14 * t30) * rSges(5,2) + t31) + g(3) * (-t49 * t29 + t34) + (g(3) * (rSges(5,1) * t14 - rSges(5,2) * t13 + t12) + t46 * t49) * t26) - m(6) * (g(1) * (t3 * t38 + t27 * t4 + (t10 * t38 + t27 * t9) * rSges(6,1) + (t27 * t10 - t9 * t38) * rSges(6,2) + t32) + g(2) * (t3 * t39 - t30 * t4 + (t10 * t39 - t30 * t9) * rSges(6,1) + (-t10 * t30 - t9 * t39) * rSges(6,2) + t31) + g(3) * (-t48 * t29 + t34) + (g(3) * (rSges(6,1) * t10 - rSges(6,2) * t9 + t3) + t46 * t48) * t26) - m(7) * (g(1) * (t1 * t38 + t27 * t2 + (t27 * t5 + t6 * t38) * rSges(7,1) + (t27 * t6 - t5 * t38) * rSges(7,2) + t32) + g(2) * (t1 * t39 - t30 * t2 + (-t30 * t5 + t6 * t39) * rSges(7,1) + (-t30 * t6 - t5 * t39) * rSges(7,2) + t31) + g(3) * (-t47 * t29 + t34) + (g(3) * (rSges(7,1) * t6 - rSges(7,2) * t5 + t1) + t46 * t47) * t26);
U  = t7;
