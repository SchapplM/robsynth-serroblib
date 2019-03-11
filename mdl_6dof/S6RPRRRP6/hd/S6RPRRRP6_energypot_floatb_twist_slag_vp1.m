% Calculate potential energy for
% S6RPRRRP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 06:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRP6_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP6_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRRRP6_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP6_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP6_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP6_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP6_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:14:08
% EndTime: 2019-03-09 06:14:09
% DurationCPUTime: 0.48s
% Computational Cost: add. (247->109), mult. (212->122), div. (0->0), fcn. (196->10), ass. (0->45)
t55 = rSges(5,3) + pkin(8);
t28 = -pkin(9) - pkin(8);
t54 = rSges(6,3) - t28;
t53 = rSges(7,3) + qJ(6) - t28;
t25 = sin(qJ(1));
t27 = cos(qJ(1));
t52 = g(1) * t27 + g(2) * t25;
t26 = cos(qJ(4));
t11 = t26 * pkin(4) + pkin(3);
t19 = pkin(10) + qJ(3);
t12 = sin(t19);
t48 = rSges(4,2) * t12;
t13 = cos(t19);
t47 = t25 * t13;
t20 = qJ(4) + qJ(5);
t14 = sin(t20);
t46 = t25 * t14;
t15 = cos(t20);
t45 = t25 * t15;
t24 = sin(qJ(4));
t44 = t25 * t24;
t43 = t25 * t26;
t42 = t27 * t13;
t41 = t27 * t14;
t40 = t27 * t15;
t39 = t27 * t24;
t38 = t27 * t26;
t35 = rSges(3,3) + qJ(2);
t34 = pkin(6) + r_base(3);
t22 = cos(pkin(10));
t9 = t22 * pkin(2) + pkin(1);
t33 = t27 * t9 + r_base(1);
t21 = sin(pkin(10));
t32 = t21 * pkin(2) + t34;
t23 = -pkin(7) - qJ(2);
t31 = t27 * t23 + t25 * t9 + r_base(2);
t30 = -t25 * t23 + t33;
t29 = rSges(3,1) * t22 - rSges(3,2) * t21 + pkin(1);
t6 = t24 * pkin(4) + pkin(5) * t14;
t5 = pkin(5) * t15 + t11;
t4 = t13 * t40 + t46;
t3 = -t13 * t41 + t45;
t2 = t13 * t45 - t41;
t1 = -t13 * t46 - t40;
t7 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t27 * rSges(2,1) - t25 * rSges(2,2) + r_base(1)) + g(2) * (t25 * rSges(2,1) + t27 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t34)) - m(3) * (g(1) * r_base(1) + g(2) * r_base(2) + g(3) * (t21 * rSges(3,1) + t22 * rSges(3,2) + t34) + (g(1) * t29 - g(2) * t35) * t27 + (g(1) * t35 + g(2) * t29) * t25) - m(4) * (g(1) * (rSges(4,1) * t42 - t27 * t48 + t33) + g(2) * (-t27 * rSges(4,3) + t31) + g(3) * (t12 * rSges(4,1) + t13 * rSges(4,2) + t32) + (g(1) * (rSges(4,3) - t23) + g(2) * (rSges(4,1) * t13 - t48)) * t25) - m(5) * (g(1) * (pkin(3) * t42 + (t13 * t38 + t44) * rSges(5,1) + (-t13 * t39 + t43) * rSges(5,2) + t30) + g(2) * (pkin(3) * t47 + (t13 * t43 - t39) * rSges(5,1) + (-t13 * t44 - t38) * rSges(5,2) + t31) + g(3) * (-t55 * t13 + t32) + (g(3) * (rSges(5,1) * t26 - rSges(5,2) * t24 + pkin(3)) + t52 * t55) * t12) - m(6) * (g(1) * (t4 * rSges(6,1) + t3 * rSges(6,2) + pkin(4) * t44 + t11 * t42 + t30) + g(2) * (t2 * rSges(6,1) + t1 * rSges(6,2) - pkin(4) * t39 + t11 * t47 + t31) + g(3) * (-t54 * t13 + t32) + (g(3) * (rSges(6,1) * t15 - rSges(6,2) * t14 + t11) + t52 * t54) * t12) - m(7) * (g(1) * (t4 * rSges(7,1) + t3 * rSges(7,2) + t25 * t6 + t5 * t42 + t30) + g(2) * (t2 * rSges(7,1) + t1 * rSges(7,2) - t27 * t6 + t5 * t47 + t31) + g(3) * (-t53 * t13 + t32) + (g(3) * (rSges(7,1) * t15 - rSges(7,2) * t14 + t5) + t52 * t53) * t12);
U  = t7;
