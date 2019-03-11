% Calculate potential energy for
% S6RPRRRP5
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
% Datum: 2019-03-09 06:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRP5_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP5_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRRRP5_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP5_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP5_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP5_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP5_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:10:40
% EndTime: 2019-03-09 06:10:40
% DurationCPUTime: 0.36s
% Computational Cost: add. (260->98), mult. (199->107), div. (0->0), fcn. (183->10), ass. (0->43)
t54 = rSges(7,1) + pkin(5);
t53 = rSges(7,3) + qJ(6);
t27 = cos(pkin(10));
t17 = t27 * pkin(2) + pkin(1);
t25 = pkin(10) + qJ(3);
t21 = qJ(4) + t25;
t15 = sin(t21);
t30 = sin(qJ(1));
t52 = t15 * t30;
t32 = cos(qJ(1));
t51 = t15 * t32;
t16 = cos(t21);
t50 = t16 * t32;
t29 = sin(qJ(5));
t49 = t30 * t29;
t31 = cos(qJ(5));
t48 = t30 * t31;
t47 = t32 * t29;
t46 = t32 * t31;
t28 = -pkin(7) - qJ(2);
t45 = rSges(4,3) - t28;
t44 = rSges(3,3) + qJ(2);
t43 = pkin(6) + r_base(3);
t20 = cos(t25);
t7 = pkin(3) * t20 + t17;
t42 = t32 * t7 + r_base(1);
t26 = sin(pkin(10));
t41 = t26 * pkin(2) + t43;
t24 = -pkin(8) + t28;
t40 = t32 * t24 + t30 * t7 + r_base(2);
t19 = sin(t25);
t39 = pkin(3) * t19 + t41;
t38 = t30 * t16 * pkin(4) + pkin(9) * t52 + t40;
t37 = t15 * pkin(4) + t39;
t36 = rSges(3,1) * t27 - rSges(3,2) * t26 + pkin(1);
t35 = rSges(4,1) * t20 - rSges(4,2) * t19 + t17;
t34 = g(1) * r_base(1) + g(2) * r_base(2);
t33 = pkin(4) * t50 + pkin(9) * t51 - t30 * t24 + t42;
t4 = t16 * t46 + t49;
t3 = t16 * t47 - t48;
t2 = t16 * t48 - t47;
t1 = t16 * t49 + t46;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t32 * rSges(2,1) - t30 * rSges(2,2) + r_base(1)) + g(2) * (t30 * rSges(2,1) + t32 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t43)) - m(3) * (g(3) * (t26 * rSges(3,1) + t27 * rSges(3,2) + t43) + (g(1) * t36 - g(2) * t44) * t32 + (g(1) * t44 + g(2) * t36) * t30 + t34) - m(4) * (g(3) * (t19 * rSges(4,1) + t20 * rSges(4,2) + t41) + (g(1) * t35 - g(2) * t45) * t32 + (g(1) * t45 + g(2) * t35) * t30 + t34) - m(5) * (g(1) * (rSges(5,1) * t50 - rSges(5,2) * t51 + t42) + g(2) * (-t32 * rSges(5,3) + t40) + g(3) * (t15 * rSges(5,1) + t16 * rSges(5,2) + t39) + (g(1) * (rSges(5,3) - t24) + g(2) * (rSges(5,1) * t16 - rSges(5,2) * t15)) * t30) - m(6) * (g(1) * (t4 * rSges(6,1) - t3 * rSges(6,2) + rSges(6,3) * t51 + t33) + g(2) * (t2 * rSges(6,1) - t1 * rSges(6,2) + rSges(6,3) * t52 + t38) + g(3) * ((-rSges(6,3) - pkin(9)) * t16 + (rSges(6,1) * t31 - rSges(6,2) * t29) * t15 + t37)) - m(7) * (g(1) * (t53 * t3 + t54 * t4 + t33) + g(2) * (t53 * t1 + t54 * t2 + t38) + g(3) * (t37 + (-rSges(7,2) - pkin(9)) * t16) + (g(3) * (t53 * t29 + t54 * t31) + (g(1) * t32 + g(2) * t30) * rSges(7,2)) * t15);
U  = t5;
