% Calculate potential energy for
% S6RPRRRP9
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
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
% Datum: 2019-03-09 06:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRP9_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP9_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRRRP9_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP9_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP9_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP9_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP9_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:26:03
% EndTime: 2019-03-09 06:26:03
% DurationCPUTime: 0.46s
% Computational Cost: add. (185->105), mult. (200->116), div. (0->0), fcn. (184->8), ass. (0->36)
t47 = rSges(5,3) + pkin(8);
t24 = -pkin(9) - pkin(8);
t46 = rSges(6,3) - t24;
t45 = rSges(7,3) + qJ(6) - t24;
t20 = sin(qJ(1));
t23 = cos(qJ(1));
t44 = -g(1) * t20 + g(2) * t23;
t21 = cos(qJ(4));
t7 = pkin(4) * t21 + pkin(3);
t22 = cos(qJ(3));
t40 = rSges(4,2) * t22;
t18 = sin(qJ(4));
t39 = t20 * t18;
t19 = sin(qJ(3));
t38 = t20 * t19;
t37 = t20 * t21;
t36 = t23 * t18;
t35 = t23 * t19;
t34 = t23 * t21;
t31 = pkin(6) + r_base(3);
t30 = pkin(1) * t20 + r_base(2);
t29 = pkin(2) + t31;
t28 = pkin(1) * t23 + qJ(2) * t20 + r_base(1);
t27 = pkin(7) * t20 + t30;
t26 = pkin(7) * t23 + t28;
t25 = -t23 * qJ(2) + t27;
t17 = qJ(4) + qJ(5);
t9 = cos(t17);
t8 = sin(t17);
t6 = pkin(4) * t18 + pkin(5) * t8;
t5 = pkin(5) * t9 + t7;
t4 = t20 * t8 - t35 * t9;
t3 = t20 * t9 + t35 * t8;
t2 = t23 * t8 + t38 * t9;
t1 = t23 * t9 - t38 * t8;
t10 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t23 - rSges(2,2) * t20 + r_base(1)) + g(2) * (rSges(2,1) * t20 + rSges(2,2) * t23 + r_base(2)) + g(3) * (rSges(2,3) + t31)) - m(3) * (g(1) * (-rSges(3,2) * t23 + rSges(3,3) * t20 + t28) + g(2) * (-t20 * rSges(3,2) + (-rSges(3,3) - qJ(2)) * t23 + t30) + g(3) * (rSges(3,1) + t31)) - m(4) * (g(1) * (rSges(4,1) * t38 + t20 * t40 + t26) + g(2) * (t20 * rSges(4,3) + t27) + g(3) * (rSges(4,1) * t22 - rSges(4,2) * t19 + t29) + (g(1) * rSges(4,3) + g(2) * (-rSges(4,1) * t19 - qJ(2) - t40)) * t23) - m(5) * (g(1) * (pkin(3) * t38 + (t19 * t37 + t36) * rSges(5,1) + (-t18 * t38 + t34) * rSges(5,2) + t26) + g(2) * (-pkin(3) * t35 + (-t19 * t34 + t39) * rSges(5,1) + (t18 * t35 + t37) * rSges(5,2) + t25) + g(3) * (t47 * t19 + t29) + (g(3) * (rSges(5,1) * t21 - rSges(5,2) * t18 + pkin(3)) + t44 * t47) * t22) - m(6) * (g(1) * (t2 * rSges(6,1) + t1 * rSges(6,2) + pkin(4) * t36 + t38 * t7 + t26) + g(2) * (t4 * rSges(6,1) + t3 * rSges(6,2) + pkin(4) * t39 - t35 * t7 + t25) + g(3) * (t46 * t19 + t29) + (g(3) * (rSges(6,1) * t9 - rSges(6,2) * t8 + t7) + t44 * t46) * t22) - m(7) * (g(1) * (t2 * rSges(7,1) + t1 * rSges(7,2) + t23 * t6 + t38 * t5 + t26) + g(2) * (t4 * rSges(7,1) + t3 * rSges(7,2) + t20 * t6 - t35 * t5 + t25) + g(3) * (t45 * t19 + t29) + (g(3) * (rSges(7,1) * t9 - rSges(7,2) * t8 + t5) + t44 * t45) * t22);
U  = t10;
