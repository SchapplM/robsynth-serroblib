% Calculate potential energy for
% S6RPRRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 05:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRPR8_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR8_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRRPR8_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR8_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR8_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR8_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR8_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:22:30
% EndTime: 2019-03-09 05:22:30
% DurationCPUTime: 0.51s
% Computational Cost: add. (195->110), mult. (200->124), div. (0->0), fcn. (184->10), ass. (0->35)
t46 = rSges(5,3) + pkin(8);
t17 = -qJ(5) - pkin(8);
t45 = rSges(6,3) - t17;
t44 = rSges(7,3) + pkin(9) - t17;
t20 = sin(qJ(1));
t23 = cos(qJ(1));
t43 = -g(1) * t20 + g(2) * t23;
t21 = cos(qJ(4));
t5 = t21 * pkin(4) + pkin(3);
t22 = cos(qJ(3));
t39 = rSges(4,2) * t22;
t18 = sin(qJ(4));
t38 = t20 * t18;
t19 = sin(qJ(3));
t37 = t20 * t19;
t36 = t20 * t21;
t35 = t23 * t18;
t34 = t23 * t19;
t33 = t23 * t21;
t30 = pkin(6) + r_base(3);
t16 = qJ(4) + pkin(10);
t29 = t20 * pkin(1) + r_base(2);
t28 = pkin(2) + t30;
t27 = t23 * pkin(1) + t20 * qJ(2) + r_base(1);
t26 = t20 * pkin(7) + t29;
t25 = t23 * pkin(7) + t27;
t24 = -t23 * qJ(2) + t26;
t8 = qJ(6) + t16;
t7 = cos(t16);
t6 = sin(t16);
t4 = cos(t8);
t3 = sin(t8);
t2 = t18 * pkin(4) + pkin(5) * t6;
t1 = pkin(5) * t7 + t5;
t9 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t23 * rSges(2,1) - t20 * rSges(2,2) + r_base(1)) + g(2) * (t20 * rSges(2,1) + t23 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t30)) - m(3) * (g(1) * (-t23 * rSges(3,2) + t20 * rSges(3,3) + t27) + g(2) * (-t20 * rSges(3,2) + (-rSges(3,3) - qJ(2)) * t23 + t29) + g(3) * (rSges(3,1) + t30)) - m(4) * (g(1) * (rSges(4,1) * t37 + t20 * t39 + t25) + g(2) * (t20 * rSges(4,3) + t26) + g(3) * (t22 * rSges(4,1) - t19 * rSges(4,2) + t28) + (g(1) * rSges(4,3) + g(2) * (-rSges(4,1) * t19 - qJ(2) - t39)) * t23) - m(5) * (g(1) * (pkin(3) * t37 + (t19 * t36 + t35) * rSges(5,1) + (-t18 * t37 + t33) * rSges(5,2) + t25) + g(2) * (-pkin(3) * t34 + (-t19 * t33 + t38) * rSges(5,1) + (t18 * t34 + t36) * rSges(5,2) + t24) + g(3) * (t46 * t19 + t28) + (g(3) * (rSges(5,1) * t21 - rSges(5,2) * t18 + pkin(3)) + t43 * t46) * t22) - m(6) * (g(1) * (t5 * t37 + pkin(4) * t35 + (t23 * t6 + t37 * t7) * rSges(6,1) + (t23 * t7 - t37 * t6) * rSges(6,2) + t25) + g(2) * (-t5 * t34 + pkin(4) * t38 + (t20 * t6 - t34 * t7) * rSges(6,1) + (t20 * t7 + t34 * t6) * rSges(6,2) + t24) + g(3) * (t45 * t19 + t28) + (g(3) * (rSges(6,1) * t7 - rSges(6,2) * t6 + t5) + t43 * t45) * t22) - m(7) * (g(1) * (t1 * t37 + t23 * t2 + (t23 * t3 + t37 * t4) * rSges(7,1) + (t23 * t4 - t3 * t37) * rSges(7,2) + t25) + g(2) * (-t1 * t34 + t20 * t2 + (t20 * t3 - t34 * t4) * rSges(7,1) + (t20 * t4 + t3 * t34) * rSges(7,2) + t24) + g(3) * (t44 * t19 + t28) + (g(3) * (rSges(7,1) * t4 - rSges(7,2) * t3 + t1) + t43 * t44) * t22);
U  = t9;
