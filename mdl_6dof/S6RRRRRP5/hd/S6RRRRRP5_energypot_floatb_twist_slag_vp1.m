% Calculate potential energy for
% S6RRRRRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 01:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRRP5_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP5_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRRRP5_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP5_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP5_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP5_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP5_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:18:02
% EndTime: 2019-03-10 01:18:02
% DurationCPUTime: 0.55s
% Computational Cost: add. (256->118), mult. (240->134), div. (0->0), fcn. (228->10), ass. (0->38)
t51 = rSges(4,3) + pkin(8);
t31 = -pkin(9) - pkin(8);
t50 = rSges(5,3) - t31;
t23 = -pkin(10) + t31;
t49 = rSges(6,3) - t23;
t48 = rSges(7,3) + qJ(6) - t23;
t27 = sin(qJ(1));
t30 = cos(qJ(1));
t47 = g(1) * t30 + g(2) * t27;
t28 = cos(qJ(3));
t13 = t28 * pkin(3) + pkin(2);
t24 = qJ(3) + qJ(4);
t14 = sin(t24);
t25 = sin(qJ(3));
t8 = t25 * pkin(3) + pkin(4) * t14;
t26 = sin(qJ(2));
t43 = rSges(3,2) * t26;
t42 = t25 * t27;
t41 = t25 * t30;
t29 = cos(qJ(2));
t40 = t27 * t29;
t39 = t29 * t30;
t35 = pkin(6) + r_base(3);
t34 = t27 * pkin(1) + r_base(2);
t15 = cos(t24);
t7 = pkin(4) * t15 + t13;
t33 = t30 * pkin(1) + t27 * pkin(7) + r_base(1);
t32 = -t30 * pkin(7) + t34;
t17 = qJ(5) + t24;
t12 = cos(t17);
t11 = sin(t17);
t6 = pkin(5) * t11 + t8;
t5 = pkin(5) * t12 + t7;
t4 = t11 * t27 + t12 * t39;
t3 = -t11 * t39 + t12 * t27;
t2 = -t11 * t30 + t12 * t40;
t1 = -t11 * t40 - t12 * t30;
t9 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t30 - rSges(2,2) * t27 + r_base(1)) + g(2) * (rSges(2,1) * t27 + rSges(2,2) * t30 + r_base(2)) + g(3) * (rSges(2,3) + t35)) - m(3) * (g(1) * (t27 * rSges(3,3) + t33) + g(2) * (rSges(3,1) * t40 - t27 * t43 + t34) + g(3) * (rSges(3,1) * t26 + rSges(3,2) * t29 + t35) + (g(1) * (rSges(3,1) * t29 - t43) + g(2) * (-rSges(3,3) - pkin(7))) * t30) - m(4) * (g(1) * (pkin(2) * t39 + (t28 * t39 + t42) * rSges(4,1) + (-t25 * t39 + t27 * t28) * rSges(4,2) + t33) + g(2) * (pkin(2) * t40 + (t28 * t40 - t41) * rSges(4,1) + (-t25 * t40 - t28 * t30) * rSges(4,2) + t32) + g(3) * (-t51 * t29 + t35) + (g(3) * (rSges(4,1) * t28 - rSges(4,2) * t25 + pkin(2)) + t47 * t51) * t26) - m(5) * (g(1) * (t13 * t39 + pkin(3) * t42 + (t14 * t27 + t15 * t39) * rSges(5,1) + (-t14 * t39 + t15 * t27) * rSges(5,2) + t33) + g(2) * (t13 * t40 - pkin(3) * t41 + (-t14 * t30 + t15 * t40) * rSges(5,1) + (-t14 * t40 - t15 * t30) * rSges(5,2) + t32) + g(3) * (-t50 * t29 + t35) + (g(3) * (rSges(5,1) * t15 - rSges(5,2) * t14 + t13) + t47 * t50) * t26) - m(6) * (g(1) * (t4 * rSges(6,1) + t3 * rSges(6,2) + t27 * t8 + t7 * t39 + t33) + g(2) * (t2 * rSges(6,1) + t1 * rSges(6,2) - t30 * t8 + t7 * t40 + t32) + g(3) * (-t49 * t29 + t35) + (g(3) * (rSges(6,1) * t12 - rSges(6,2) * t11 + t7) + t47 * t49) * t26) - m(7) * (g(1) * (t4 * rSges(7,1) + t3 * rSges(7,2) + t27 * t6 + t5 * t39 + t33) + g(2) * (t2 * rSges(7,1) + t1 * rSges(7,2) - t30 * t6 + t5 * t40 + t32) + g(3) * (-t48 * t29 + t35) + (g(3) * (rSges(7,1) * t12 - rSges(7,2) * t11 + t5) + t47 * t48) * t26);
U  = t9;
