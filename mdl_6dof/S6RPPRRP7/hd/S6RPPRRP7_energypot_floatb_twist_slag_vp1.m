% Calculate potential energy for
% S6RPPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
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
% Datum: 2019-03-09 02:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPRRP7_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP7_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPPRRP7_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP7_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP7_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP7_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRP7_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:12:30
% EndTime: 2019-03-09 02:12:31
% DurationCPUTime: 0.42s
% Computational Cost: add. (184->99), mult. (180->105), div. (0->0), fcn. (160->8), ass. (0->38)
t43 = rSges(6,3) + pkin(8);
t42 = rSges(7,3) + qJ(6) + pkin(8);
t13 = pkin(9) + qJ(4);
t7 = sin(t13);
t41 = pkin(4) * t7;
t14 = sin(pkin(9));
t40 = pkin(3) * t14;
t18 = sin(qJ(5));
t19 = sin(qJ(1));
t39 = t19 * t18;
t20 = cos(qJ(5));
t38 = t19 * t20;
t21 = cos(qJ(1));
t37 = t21 * t18;
t36 = t21 * t20;
t17 = -pkin(7) - qJ(3);
t35 = rSges(5,3) - t17;
t34 = rSges(4,3) + qJ(3);
t33 = pkin(6) + r_base(3);
t32 = t19 * pkin(1) + r_base(2);
t31 = pkin(2) + t33;
t30 = pkin(5) * t18 - t17;
t29 = t21 * pkin(1) + t19 * qJ(2) + r_base(1);
t28 = -qJ(2) - t40;
t27 = g(2) * t32;
t15 = cos(pkin(9));
t26 = t15 * pkin(3) + t31;
t25 = t19 * t40 + t29;
t8 = cos(t13);
t24 = rSges(5,1) * t7 + rSges(5,2) * t8;
t23 = rSges(4,1) * t14 + rSges(4,2) * t15;
t6 = t20 * pkin(5) + pkin(4);
t22 = -t42 * t8 + t6 * t7;
t4 = -t36 * t7 + t39;
t3 = t37 * t7 + t38;
t2 = t38 * t7 + t37;
t1 = -t39 * t7 + t36;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t21 * rSges(2,1) - t19 * rSges(2,2) + r_base(1)) + g(2) * (t19 * rSges(2,1) + t21 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t33)) - m(3) * (g(1) * (-t21 * rSges(3,2) + t19 * rSges(3,3) + t29) + g(2) * (-t19 * rSges(3,2) + (-rSges(3,3) - qJ(2)) * t21 + t32) + g(3) * (rSges(3,1) + t33)) - m(4) * (g(1) * t29 + t27 + g(3) * (t15 * rSges(4,1) - t14 * rSges(4,2) + t31) + (g(1) * t23 + g(2) * t34) * t19 + (g(1) * t34 + g(2) * (-qJ(2) - t23)) * t21) - m(5) * (g(1) * t25 + t27 + g(3) * (t8 * rSges(5,1) - t7 * rSges(5,2) + t26) + (g(1) * t24 + g(2) * t35) * t19 + (g(1) * t35 + g(2) * (-t24 + t28)) * t21) - m(6) * (g(1) * (t2 * rSges(6,1) + t1 * rSges(6,2) + t19 * t41 + t25) + g(2) * (t4 * rSges(6,1) + t3 * rSges(6,2) - t19 * t17 + t32) + g(3) * (t43 * t7 + t26) + (g(3) * (rSges(6,1) * t20 - rSges(6,2) * t18 + pkin(4)) - g(1) * t43 * t19) * t8 + (-g(1) * t17 + g(2) * (t43 * t8 + t28 - t41)) * t21) - m(7) * (g(1) * (t2 * rSges(7,1) + t1 * rSges(7,2) + t25) + g(2) * (t4 * rSges(7,1) + t3 * rSges(7,2) + t32) + (g(1) * t22 + g(2) * t30) * t19 + (g(1) * t30 + g(2) * (-t22 + t28)) * t21 + (t26 + (rSges(7,1) * t20 - rSges(7,2) * t18 + t6) * t8 + t42 * t7) * g(3));
U  = t5;
