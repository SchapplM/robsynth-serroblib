% Calculate potential energy for
% S6RPRRRP3
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
% Datum: 2019-03-09 06:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRP3_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP3_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRRRP3_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP3_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP3_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP3_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP3_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:02:48
% EndTime: 2019-03-09 06:02:49
% DurationCPUTime: 0.46s
% Computational Cost: add. (271->101), mult. (213->115), div. (0->0), fcn. (201->10), ass. (0->42)
t57 = rSges(7,1) + pkin(5);
t56 = rSges(5,3) + pkin(8);
t21 = qJ(1) + pkin(10);
t15 = sin(t21);
t16 = cos(t21);
t55 = g(1) * t16 + g(2) * t15;
t54 = rSges(7,3) + qJ(6);
t24 = sin(qJ(3));
t50 = rSges(4,2) * t24;
t23 = sin(qJ(4));
t49 = t15 * t23;
t27 = cos(qJ(3));
t48 = t15 * t27;
t47 = t16 * t23;
t46 = t16 * t27;
t22 = qJ(4) + qJ(5);
t17 = sin(t22);
t45 = t17 * t27;
t18 = cos(t22);
t44 = t18 * t27;
t43 = t23 * t27;
t26 = cos(qJ(4));
t42 = t26 * t27;
t39 = pkin(6) + r_base(3);
t25 = sin(qJ(1));
t38 = t25 * pkin(1) + r_base(2);
t28 = cos(qJ(1));
t37 = t28 * pkin(1) + r_base(1);
t36 = qJ(2) + t39;
t35 = t15 * pkin(2) + t38;
t34 = t16 * pkin(2) + t15 * pkin(7) + t37;
t13 = t26 * pkin(4) + pkin(3);
t29 = -pkin(9) - pkin(8);
t33 = t24 * t13 + t27 * t29 + t36;
t32 = -t16 * pkin(7) + t35;
t31 = pkin(4) * t49 + t13 * t46 + t34;
t30 = -pkin(4) * t47 + t13 * t48 + t32;
t4 = t15 * t17 + t16 * t44;
t3 = -t15 * t18 + t16 * t45;
t2 = t15 * t44 - t16 * t17;
t1 = t15 * t45 + t16 * t18;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t28 * rSges(2,1) - t25 * rSges(2,2) + r_base(1)) + g(2) * (t25 * rSges(2,1) + t28 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t39)) - m(3) * (g(1) * (t16 * rSges(3,1) - t15 * rSges(3,2) + t37) + g(2) * (t15 * rSges(3,1) + t16 * rSges(3,2) + t38) + g(3) * (rSges(3,3) + t36)) - m(4) * (g(1) * (t15 * rSges(4,3) + t34) + g(2) * (rSges(4,1) * t48 - t15 * t50 + t35) + g(3) * (t24 * rSges(4,1) + t27 * rSges(4,2) + t36) + (g(1) * (rSges(4,1) * t27 - t50) + g(2) * (-rSges(4,3) - pkin(7))) * t16) - m(5) * (g(1) * (pkin(3) * t46 + (t16 * t42 + t49) * rSges(5,1) + (t15 * t26 - t16 * t43) * rSges(5,2) + t34) + g(2) * (pkin(3) * t48 + (t15 * t42 - t47) * rSges(5,1) + (-t15 * t43 - t16 * t26) * rSges(5,2) + t32) + g(3) * (-t56 * t27 + t36) + (g(3) * (rSges(5,1) * t26 - rSges(5,2) * t23 + pkin(3)) + t55 * t56) * t24) - m(6) * (g(1) * (t4 * rSges(6,1) - t3 * rSges(6,2) + t31) + g(2) * (t2 * rSges(6,1) - t1 * rSges(6,2) + t30) + g(3) * (-t27 * rSges(6,3) + t33) + (g(3) * (rSges(6,1) * t18 - rSges(6,2) * t17) + t55 * (rSges(6,3) - t29)) * t24) - m(7) * (g(1) * (t54 * t3 + t57 * t4 + t31) + g(2) * (t54 * t1 + t57 * t2 + t30) + g(3) * (-t27 * rSges(7,2) + t33) + (g(3) * (t54 * t17 + t57 * t18) + t55 * (rSges(7,2) - t29)) * t24);
U  = t5;
