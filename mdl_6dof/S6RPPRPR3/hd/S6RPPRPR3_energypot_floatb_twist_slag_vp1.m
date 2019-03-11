% Calculate potential energy for
% S6RPPRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 01:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPRPR3_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR3_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPPRPR3_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR3_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR3_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR3_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRPR3_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:43:48
% EndTime: 2019-03-09 01:43:48
% DurationCPUTime: 0.32s
% Computational Cost: add. (210->84), mult. (141->82), div. (0->0), fcn. (117->10), ass. (0->34)
t15 = sin(qJ(6));
t18 = cos(qJ(6));
t23 = rSges(7,1) * t18 - rSges(7,2) * t15 + pkin(5);
t12 = qJ(4) + pkin(10);
t5 = sin(t12);
t41 = t23 * t5;
t40 = rSges(7,3) + pkin(8);
t16 = sin(qJ(4));
t39 = t16 * pkin(4);
t38 = rSges(5,3) + pkin(7);
t14 = -qJ(5) - pkin(7);
t36 = rSges(6,3) - t14;
t35 = pkin(6) + r_base(3);
t17 = sin(qJ(1));
t34 = t17 * pkin(1) + r_base(2);
t20 = cos(qJ(1));
t33 = t20 * pkin(1) + r_base(1);
t13 = qJ(1) + pkin(9);
t6 = sin(t13);
t32 = t6 * pkin(2) + t34;
t31 = -qJ(3) - t39;
t30 = qJ(2) + t35;
t8 = cos(t13);
t29 = t8 * pkin(2) + t6 * qJ(3) + t33;
t28 = pkin(3) + t30;
t27 = g(2) * t32;
t7 = cos(t12);
t26 = rSges(6,1) * t5 + rSges(6,2) * t7;
t19 = cos(qJ(4));
t25 = rSges(5,1) * t16 + rSges(5,2) * t19;
t24 = t19 * pkin(4) + t28;
t22 = t15 * rSges(7,1) + t18 * rSges(7,2) - t14;
t21 = g(1) * (t6 * t39 + t29) + t27;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t20 * rSges(2,1) - t17 * rSges(2,2) + r_base(1)) + g(2) * (t17 * rSges(2,1) + t20 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t35)) - m(3) * (g(1) * (t8 * rSges(3,1) - t6 * rSges(3,2) + t33) + g(2) * (t6 * rSges(3,1) + t8 * rSges(3,2) + t34) + g(3) * (rSges(3,3) + t30)) - m(4) * (g(1) * (-t8 * rSges(4,2) + t6 * rSges(4,3) + t29) + g(2) * (-t6 * rSges(4,2) + (-rSges(4,3) - qJ(3)) * t8 + t32) + g(3) * (rSges(4,1) + t30)) - m(5) * (g(1) * t29 + t27 + g(3) * (t19 * rSges(5,1) - t16 * rSges(5,2) + t28) + (g(1) * t25 + g(2) * t38) * t6 + (g(1) * t38 + g(2) * (-qJ(3) - t25)) * t8) - m(6) * (g(3) * (t7 * rSges(6,1) - t5 * rSges(6,2) + t24) + (g(1) * t26 + g(2) * t36) * t6 + (g(1) * t36 + g(2) * (-t26 + t31)) * t8 + t21) - m(7) * (g(3) * (t40 * t5 + t24) + (g(1) * t41 + g(2) * t22) * t6 + (g(1) * t22 + (t31 - t41) * g(2)) * t8 + (g(3) * t23 + (-g(1) * t6 + g(2) * t8) * t40) * t7 + t21);
U  = t1;
