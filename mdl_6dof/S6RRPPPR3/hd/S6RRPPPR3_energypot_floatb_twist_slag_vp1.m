% Calculate potential energy for
% S6RRPPPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta5]';
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
% Datum: 2019-03-09 08:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPPPR3_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR3_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPPPR3_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR3_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR3_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR3_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPPR3_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:13:45
% EndTime: 2019-03-09 08:13:46
% DurationCPUTime: 0.51s
% Computational Cost: add. (173->111), mult. (226->126), div. (0->0), fcn. (206->8), ass. (0->34)
t50 = rSges(7,3) + pkin(8) + qJ(5);
t21 = sin(qJ(1));
t23 = cos(qJ(1));
t49 = g(1) * t23 + g(2) * t21;
t48 = rSges(6,3) + qJ(5);
t20 = sin(qJ(2));
t45 = t20 * t23;
t17 = sin(pkin(9));
t44 = t21 * t17;
t43 = t21 * t20;
t22 = cos(qJ(2));
t42 = t21 * t22;
t41 = t22 * t23;
t40 = t23 * t17;
t18 = cos(pkin(9));
t39 = t23 * t18;
t37 = qJ(3) * t20;
t35 = pkin(6) + r_base(3);
t34 = t21 * pkin(1) + r_base(2);
t33 = t20 * pkin(2) + t35;
t32 = t23 * pkin(1) + t21 * pkin(7) + r_base(1);
t31 = pkin(2) * t42 + t21 * t37 + t34;
t30 = t20 * pkin(3) + t33;
t29 = rSges(5,1) * t20 - rSges(5,2) * t22;
t28 = pkin(2) * t41 + t23 * t37 + t32;
t27 = pkin(3) * t42 + t23 * qJ(4) + t31;
t26 = pkin(3) * t41 + t28;
t25 = -t23 * pkin(7) + t27;
t24 = -t21 * qJ(4) + t26;
t16 = pkin(9) + qJ(6);
t9 = cos(t16);
t8 = sin(t16);
t7 = t18 * pkin(5) + pkin(4);
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t23 * rSges(2,1) - t21 * rSges(2,2) + r_base(1)) + g(2) * (t21 * rSges(2,1) + t23 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t35)) - m(3) * (g(1) * (t21 * rSges(3,3) + t32) + g(2) * (rSges(3,1) * t42 - rSges(3,2) * t43 + t34) + g(3) * (t20 * rSges(3,1) + t22 * rSges(3,2) + t35) + (g(1) * (rSges(3,1) * t22 - rSges(3,2) * t20) + g(2) * (-rSges(3,3) - pkin(7))) * t23) - m(4) * (g(1) * (t21 * rSges(4,2) + t28) + g(2) * (rSges(4,1) * t42 + rSges(4,3) * t43 + t31) + g(3) * (t20 * rSges(4,1) + (-rSges(4,3) - qJ(3)) * t22 + t33) + (g(1) * (rSges(4,1) * t22 + rSges(4,3) * t20) + g(2) * (-rSges(4,2) - pkin(7))) * t23) - m(5) * (g(1) * t26 + g(2) * t27 + g(3) * (-t20 * rSges(5,2) + (-rSges(5,1) - qJ(3)) * t22 + t30) + (g(1) * t29 + g(2) * (rSges(5,3) - pkin(7))) * t23 + (g(1) * (-rSges(5,3) - qJ(4)) + g(2) * t29) * t21) - m(6) * (g(1) * (pkin(4) * t45 + (t20 * t39 - t44) * rSges(6,1) + (-t21 * t18 - t20 * t40) * rSges(6,2) + t24) + g(2) * (pkin(4) * t43 + (t18 * t43 + t40) * rSges(6,1) + (-t17 * t43 + t39) * rSges(6,2) + t25) + g(3) * (t48 * t20 + t30) + (g(3) * (-rSges(6,1) * t18 + rSges(6,2) * t17 - pkin(4) - qJ(3)) + t49 * t48) * t22) - m(7) * (g(1) * (t7 * t45 - pkin(5) * t44 + (-t21 * t8 + t45 * t9) * rSges(7,1) + (-t21 * t9 - t45 * t8) * rSges(7,2) + t24) + g(2) * (t7 * t43 + pkin(5) * t40 + (t23 * t8 + t43 * t9) * rSges(7,1) + (t23 * t9 - t43 * t8) * rSges(7,2) + t25) + g(3) * (t50 * t20 + t30) + (g(3) * (-rSges(7,1) * t9 + rSges(7,2) * t8 - qJ(3) - t7) + t49 * t50) * t22);
U  = t1;
