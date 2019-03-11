% Calculate potential energy for
% S6RPRRRP4
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
% Datum: 2019-03-09 06:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRP4_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP4_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRRRP4_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP4_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP4_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP4_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP4_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:06:52
% EndTime: 2019-03-09 06:06:53
% DurationCPUTime: 0.40s
% Computational Cost: add. (246->99), mult. (186->107), div. (0->0), fcn. (166->10), ass. (0->43)
t53 = rSges(6,3) + pkin(9);
t52 = rSges(7,3) + qJ(6) + pkin(9);
t26 = sin(qJ(1));
t28 = cos(qJ(1));
t51 = g(1) * t28 + g(2) * t26;
t22 = cos(pkin(10));
t12 = t22 * pkin(2) + pkin(1);
t20 = pkin(10) + qJ(3);
t16 = qJ(4) + t20;
t10 = sin(t16);
t47 = rSges(5,2) * t10;
t11 = cos(t16);
t46 = t11 * t26;
t45 = t11 * t28;
t25 = sin(qJ(5));
t44 = t26 * t25;
t27 = cos(qJ(5));
t43 = t26 * t27;
t42 = t28 * t25;
t41 = t28 * t27;
t24 = -pkin(7) - qJ(2);
t40 = rSges(4,3) - t24;
t38 = rSges(3,3) + qJ(2);
t37 = pkin(6) + r_base(3);
t15 = cos(t20);
t7 = pkin(3) * t15 + t12;
t36 = t28 * t7 + r_base(1);
t19 = -pkin(8) + t24;
t35 = t28 * t19 + t26 * t7 + r_base(2);
t21 = sin(pkin(10));
t34 = t21 * pkin(2) + t37;
t14 = sin(t20);
t33 = pkin(3) * t14 + t34;
t32 = -t26 * t19 + t36;
t31 = rSges(3,1) * t22 - rSges(3,2) * t21 + pkin(1);
t30 = rSges(4,1) * t15 - rSges(4,2) * t14 + t12;
t29 = g(1) * r_base(1) + g(2) * r_base(2);
t13 = t27 * pkin(5) + pkin(4);
t4 = t11 * t41 + t44;
t3 = -t11 * t42 + t43;
t2 = t11 * t43 - t42;
t1 = -t11 * t44 - t41;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t28 * rSges(2,1) - t26 * rSges(2,2) + r_base(1)) + g(2) * (t26 * rSges(2,1) + t28 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t37)) - m(3) * (g(3) * (t21 * rSges(3,1) + t22 * rSges(3,2) + t37) + (g(1) * t31 - g(2) * t38) * t28 + (g(1) * t38 + g(2) * t31) * t26 + t29) - m(4) * (g(3) * (t14 * rSges(4,1) + t15 * rSges(4,2) + t34) + (g(1) * t30 - g(2) * t40) * t28 + (g(1) * t40 + g(2) * t30) * t26 + t29) - m(5) * (g(1) * (rSges(5,1) * t45 - t28 * t47 + t36) + g(2) * (-t28 * rSges(5,3) + t35) + g(3) * (t10 * rSges(5,1) + t11 * rSges(5,2) + t33) + (g(1) * (rSges(5,3) - t19) + g(2) * (rSges(5,1) * t11 - t47)) * t26) - m(6) * (g(1) * (t4 * rSges(6,1) + t3 * rSges(6,2) + pkin(4) * t45 + t32) + g(2) * (t2 * rSges(6,1) + t1 * rSges(6,2) + pkin(4) * t46 + t35) + g(3) * (-t53 * t11 + t33) + (g(3) * (rSges(6,1) * t27 - rSges(6,2) * t25 + pkin(4)) + t51 * t53) * t10) - m(7) * (g(1) * (t4 * rSges(7,1) + t3 * rSges(7,2) + pkin(5) * t44 + t13 * t45 + t32) + g(2) * (t2 * rSges(7,1) + t1 * rSges(7,2) - pkin(5) * t42 + t13 * t46 + t35) + g(3) * (-t52 * t11 + t33) + (g(3) * (rSges(7,1) * t27 - rSges(7,2) * t25 + t13) + t51 * t52) * t10);
U  = t5;
