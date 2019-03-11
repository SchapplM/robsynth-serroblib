% Calculate potential energy for
% S6RPRRPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
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
% Datum: 2019-03-09 05:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRPR5_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR5_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRRPR5_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR5_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR5_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR5_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR5_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:12:05
% EndTime: 2019-03-09 05:12:06
% DurationCPUTime: 0.42s
% Computational Cost: add. (238->99), mult. (173->107), div. (0->0), fcn. (149->10), ass. (0->38)
t49 = rSges(7,3) + pkin(9);
t22 = cos(pkin(10));
t13 = t22 * pkin(2) + pkin(1);
t20 = pkin(10) + qJ(3);
t16 = qJ(4) + t20;
t11 = sin(t16);
t27 = cos(qJ(1));
t47 = t11 * t27;
t24 = sin(qJ(6));
t25 = sin(qJ(1));
t46 = t25 * t24;
t26 = cos(qJ(6));
t45 = t25 * t26;
t12 = cos(t16);
t44 = t27 * t12;
t43 = t27 * t24;
t42 = t27 * t26;
t23 = -pkin(7) - qJ(2);
t41 = rSges(4,3) - t23;
t40 = qJ(5) * t11;
t39 = rSges(3,3) + qJ(2);
t38 = pkin(6) + r_base(3);
t15 = cos(t20);
t3 = pkin(3) * t15 + t13;
t37 = t27 * t3 + r_base(1);
t19 = -pkin(8) + t23;
t36 = t27 * t19 + t25 * t3 + r_base(2);
t21 = sin(pkin(10));
t35 = t21 * pkin(2) + t38;
t34 = pkin(4) * t44 + t27 * t40 + t37;
t14 = sin(t20);
t33 = pkin(3) * t14 + t35;
t32 = t36 + (pkin(4) * t12 + t40) * t25;
t31 = t11 * pkin(4) + t33;
t30 = rSges(3,1) * t22 - rSges(3,2) * t21 + pkin(1);
t29 = rSges(4,1) * t15 - rSges(4,2) * t14 + t13;
t28 = g(1) * r_base(1) + g(2) * r_base(2);
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t27 * rSges(2,1) - t25 * rSges(2,2) + r_base(1)) + g(2) * (t25 * rSges(2,1) + t27 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t38)) - m(3) * (g(3) * (t21 * rSges(3,1) + t22 * rSges(3,2) + t38) + (g(1) * t30 - g(2) * t39) * t27 + (g(1) * t39 + g(2) * t30) * t25 + t28) - m(4) * (g(3) * (t14 * rSges(4,1) + t15 * rSges(4,2) + t35) + (g(1) * t29 - g(2) * t41) * t27 + (g(1) * t41 + g(2) * t29) * t25 + t28) - m(5) * (g(1) * (rSges(5,1) * t44 - rSges(5,2) * t47 + t37) + g(2) * (-t27 * rSges(5,3) + t36) + g(3) * (t11 * rSges(5,1) + t12 * rSges(5,2) + t33) + (g(1) * (rSges(5,3) - t19) + g(2) * (rSges(5,1) * t12 - rSges(5,2) * t11)) * t25) - m(6) * (g(1) * (-rSges(6,2) * t44 + rSges(6,3) * t47 + t34) + g(2) * (-t27 * rSges(6,1) + t32) + g(3) * (-t11 * rSges(6,2) + (-rSges(6,3) - qJ(5)) * t12 + t31) + (g(1) * (rSges(6,1) - t19) + g(2) * (-rSges(6,2) * t12 + rSges(6,3) * t11)) * t25) - m(7) * (g(1) * ((t11 * t43 + t45) * rSges(7,1) + (t11 * t42 - t46) * rSges(7,2) + t34 + (pkin(5) - t19) * t25) + g(2) * (-t27 * pkin(5) + (t11 * t46 - t42) * rSges(7,1) + (t11 * t45 + t43) * rSges(7,2) + t32) + g(3) * (t49 * t11 + t31) + (g(3) * (-rSges(7,1) * t24 - rSges(7,2) * t26 - qJ(5)) + (g(1) * t27 + g(2) * t25) * t49) * t12);
U  = t1;
