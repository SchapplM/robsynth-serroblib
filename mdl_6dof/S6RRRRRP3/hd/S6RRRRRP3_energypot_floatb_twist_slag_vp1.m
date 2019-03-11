% Calculate potential energy for
% S6RRRRRP3
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
% Datum: 2019-03-10 01:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRRP3_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP3_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRRRP3_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP3_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP3_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP3_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP3_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:05:45
% EndTime: 2019-03-10 01:05:45
% DurationCPUTime: 0.48s
% Computational Cost: add. (247->109), mult. (212->122), div. (0->0), fcn. (196->10), ass. (0->41)
t51 = rSges(5,3) + pkin(9);
t27 = -pkin(10) - pkin(9);
t50 = rSges(6,3) - t27;
t49 = rSges(7,3) + qJ(6) - t27;
t23 = sin(qJ(1));
t26 = cos(qJ(1));
t48 = g(1) * t26 + g(2) * t23;
t45 = rSges(3,3) + pkin(7);
t24 = cos(qJ(4));
t9 = t24 * pkin(4) + pkin(3);
t20 = qJ(2) + qJ(3);
t13 = sin(t20);
t43 = rSges(4,2) * t13;
t15 = cos(t20);
t42 = t23 * t15;
t21 = sin(qJ(4));
t41 = t23 * t21;
t40 = t23 * t24;
t39 = t26 * t15;
t38 = t26 * t21;
t37 = t26 * t24;
t34 = pkin(6) + r_base(3);
t25 = cos(qJ(2));
t10 = t25 * pkin(2) + pkin(1);
t33 = t26 * t10 + r_base(1);
t22 = sin(qJ(2));
t32 = t22 * pkin(2) + t34;
t28 = -pkin(8) - pkin(7);
t31 = t23 * t10 + t26 * t28 + r_base(2);
t30 = -t23 * t28 + t33;
t29 = rSges(3,1) * t25 - rSges(3,2) * t22 + pkin(1);
t19 = qJ(4) + qJ(5);
t14 = cos(t19);
t12 = sin(t19);
t6 = t21 * pkin(4) + pkin(5) * t12;
t5 = pkin(5) * t14 + t9;
t4 = t23 * t12 + t14 * t39;
t3 = -t12 * t39 + t23 * t14;
t2 = -t26 * t12 + t14 * t42;
t1 = -t12 * t42 - t26 * t14;
t7 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t26 * rSges(2,1) - t23 * rSges(2,2) + r_base(1)) + g(2) * (t23 * rSges(2,1) + t26 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t34)) - m(3) * (g(1) * r_base(1) + g(2) * r_base(2) + g(3) * (t22 * rSges(3,1) + t25 * rSges(3,2) + t34) + (g(1) * t29 - g(2) * t45) * t26 + (g(1) * t45 + g(2) * t29) * t23) - m(4) * (g(1) * (rSges(4,1) * t39 - t26 * t43 + t33) + g(2) * (-t26 * rSges(4,3) + t31) + g(3) * (t13 * rSges(4,1) + t15 * rSges(4,2) + t32) + (g(1) * (rSges(4,3) - t28) + g(2) * (rSges(4,1) * t15 - t43)) * t23) - m(5) * (g(1) * (pkin(3) * t39 + (t15 * t37 + t41) * rSges(5,1) + (-t15 * t38 + t40) * rSges(5,2) + t30) + g(2) * (pkin(3) * t42 + (t15 * t40 - t38) * rSges(5,1) + (-t15 * t41 - t37) * rSges(5,2) + t31) + g(3) * (-t51 * t15 + t32) + (g(3) * (rSges(5,1) * t24 - rSges(5,2) * t21 + pkin(3)) + t48 * t51) * t13) - m(6) * (g(1) * (t4 * rSges(6,1) + t3 * rSges(6,2) + pkin(4) * t41 + t9 * t39 + t30) + g(2) * (t2 * rSges(6,1) + t1 * rSges(6,2) - pkin(4) * t38 + t9 * t42 + t31) + g(3) * (-t50 * t15 + t32) + (g(3) * (rSges(6,1) * t14 - rSges(6,2) * t12 + t9) + t48 * t50) * t13) - m(7) * (g(1) * (t4 * rSges(7,1) + t3 * rSges(7,2) + t23 * t6 + t5 * t39 + t30) + g(2) * (t2 * rSges(7,1) + t1 * rSges(7,2) - t26 * t6 + t5 * t42 + t31) + g(3) * (-t49 * t15 + t32) + (g(3) * (rSges(7,1) * t14 - rSges(7,2) * t12 + t5) + t48 * t49) * t13);
U  = t7;
