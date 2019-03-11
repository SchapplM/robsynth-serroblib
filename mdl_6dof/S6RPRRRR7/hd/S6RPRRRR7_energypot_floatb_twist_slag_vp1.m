% Calculate potential energy for
% S6RPRRRR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
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
% Datum: 2019-03-09 07:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRR7_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR7_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRRRR7_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR7_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR7_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR7_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR7_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:16:52
% EndTime: 2019-03-09 07:16:52
% DurationCPUTime: 0.40s
% Computational Cost: add. (193->87), mult. (156->83), div. (0->0), fcn. (132->10), ass. (0->37)
t15 = sin(qJ(6));
t18 = cos(qJ(6));
t26 = rSges(7,1) * t18 - rSges(7,2) * t15 + pkin(5);
t14 = qJ(3) + qJ(4);
t9 = qJ(5) + t14;
t4 = sin(t9);
t41 = rSges(7,3) + pkin(10);
t5 = cos(t9);
t43 = t26 * t4 - t41 * t5;
t21 = -pkin(8) - pkin(7);
t16 = sin(qJ(3));
t40 = t16 * pkin(3);
t39 = rSges(4,3) + pkin(7);
t38 = rSges(5,3) - t21;
t13 = -pkin(9) + t21;
t37 = rSges(6,3) - t13;
t6 = sin(t14);
t2 = pkin(4) * t6 + t40;
t36 = -qJ(2) - t2;
t35 = pkin(6) + r_base(3);
t17 = sin(qJ(1));
t34 = t17 * pkin(1) + r_base(2);
t33 = pkin(2) + t35;
t20 = cos(qJ(1));
t32 = t20 * pkin(1) + t17 * qJ(2) + r_base(1);
t31 = g(2) * t34;
t19 = cos(qJ(3));
t30 = t19 * pkin(3) + t33;
t29 = rSges(6,1) * t4 + rSges(6,2) * t5;
t7 = cos(t14);
t28 = pkin(4) * t7 + t30;
t27 = rSges(4,1) * t16 + rSges(4,2) * t19;
t25 = t15 * rSges(7,1) + t18 * rSges(7,2) - t13;
t24 = rSges(5,1) * t6 + rSges(5,2) * t7 + t40;
t23 = g(1) * t32 + t31;
t22 = g(1) * (t17 * t2 + t32) + t31;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t20 * rSges(2,1) - t17 * rSges(2,2) + r_base(1)) + g(2) * (t17 * rSges(2,1) + t20 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t35)) - m(3) * (g(1) * (-t20 * rSges(3,2) + t17 * rSges(3,3) + t32) + g(2) * (-t17 * rSges(3,2) + (-rSges(3,3) - qJ(2)) * t20 + t34) + g(3) * (rSges(3,1) + t35)) - m(4) * (g(3) * (t19 * rSges(4,1) - t16 * rSges(4,2) + t33) + (g(1) * t27 + g(2) * t39) * t17 + (g(1) * t39 + g(2) * (-qJ(2) - t27)) * t20 + t23) - m(5) * (g(3) * (t7 * rSges(5,1) - t6 * rSges(5,2) + t30) + (g(1) * t24 + g(2) * t38) * t17 + (g(1) * t38 + g(2) * (-qJ(2) - t24)) * t20 + t23) - m(6) * (g(3) * (t5 * rSges(6,1) - t4 * rSges(6,2) + t28) + (g(1) * t29 + g(2) * t37) * t17 + (g(1) * t37 + g(2) * (-t29 + t36)) * t20 + t22) - m(7) * ((t43 * g(1) + g(2) * t25) * t17 + (g(1) * t25 + (t36 - t43) * g(2)) * t20 + t22 + (t26 * t5 + t41 * t4 + t28) * g(3));
U  = t1;
