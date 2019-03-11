% Calculate potential energy for
% S6RPPRRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-03-09 01:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPRRP1_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP1_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPPRRP1_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP1_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP1_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP1_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRP1_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:57:32
% EndTime: 2019-03-09 01:57:32
% DurationCPUTime: 0.39s
% Computational Cost: add. (245->96), mult. (172->103), div. (0->0), fcn. (152->10), ass. (0->40)
t50 = rSges(6,3) + pkin(8);
t49 = rSges(7,3) + qJ(6) + pkin(8);
t18 = qJ(1) + pkin(9);
t11 = sin(t18);
t13 = cos(t18);
t48 = g(1) * t13 + g(2) * t11;
t17 = pkin(10) + qJ(4);
t10 = sin(t17);
t44 = rSges(5,2) * t10;
t12 = cos(t17);
t43 = t11 * t12;
t23 = sin(qJ(5));
t42 = t11 * t23;
t25 = cos(qJ(5));
t41 = t11 * t25;
t40 = t12 * t13;
t39 = t13 * t23;
t38 = t13 * t25;
t36 = rSges(4,3) + qJ(3);
t35 = pkin(6) + r_base(3);
t24 = sin(qJ(1));
t34 = t24 * pkin(1) + r_base(2);
t26 = cos(qJ(1));
t33 = t26 * pkin(1) + r_base(1);
t20 = cos(pkin(10));
t8 = t20 * pkin(3) + pkin(2);
t32 = t13 * t8 + t33;
t31 = qJ(2) + t35;
t22 = -pkin(7) - qJ(3);
t30 = t11 * t8 + t13 * t22 + t34;
t19 = sin(pkin(10));
t29 = t19 * pkin(3) + t31;
t28 = rSges(4,1) * t20 - rSges(4,2) * t19 + pkin(2);
t27 = -t11 * t22 + t32;
t9 = t25 * pkin(5) + pkin(4);
t4 = t12 * t38 + t42;
t3 = -t12 * t39 + t41;
t2 = t12 * t41 - t39;
t1 = -t12 * t42 - t38;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t26 * rSges(2,1) - t24 * rSges(2,2) + r_base(1)) + g(2) * (t24 * rSges(2,1) + t26 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t35)) - m(3) * (g(1) * (t13 * rSges(3,1) - t11 * rSges(3,2) + t33) + g(2) * (t11 * rSges(3,1) + t13 * rSges(3,2) + t34) + g(3) * (rSges(3,3) + t31)) - m(4) * (g(1) * t33 + g(2) * t34 + g(3) * (t19 * rSges(4,1) + t20 * rSges(4,2) + t31) + (g(1) * t28 - g(2) * t36) * t13 + (g(1) * t36 + g(2) * t28) * t11) - m(5) * (g(1) * (rSges(5,1) * t40 - t13 * t44 + t32) + g(2) * (-t13 * rSges(5,3) + t30) + g(3) * (t10 * rSges(5,1) + t12 * rSges(5,2) + t29) + (g(1) * (rSges(5,3) - t22) + g(2) * (rSges(5,1) * t12 - t44)) * t11) - m(6) * (g(1) * (t4 * rSges(6,1) + t3 * rSges(6,2) + pkin(4) * t40 + t27) + g(2) * (t2 * rSges(6,1) + t1 * rSges(6,2) + pkin(4) * t43 + t30) + g(3) * (-t50 * t12 + t29) + (g(3) * (rSges(6,1) * t25 - rSges(6,2) * t23 + pkin(4)) + t48 * t50) * t10) - m(7) * (g(1) * (t4 * rSges(7,1) + t3 * rSges(7,2) + pkin(5) * t42 + t40 * t9 + t27) + g(2) * (t2 * rSges(7,1) + t1 * rSges(7,2) - pkin(5) * t39 + t43 * t9 + t30) + g(3) * (-t49 * t12 + t29) + (g(3) * (rSges(7,1) * t25 - rSges(7,2) * t23 + t9) + t48 * t49) * t10);
U  = t5;
