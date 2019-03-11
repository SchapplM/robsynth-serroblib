% Calculate potential energy for
% S6RPPPRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2]';
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
% Datum: 2019-03-09 01:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPPRR1_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR1_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPPPRR1_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR1_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR1_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR1_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPPRR1_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:29:47
% EndTime: 2019-03-09 01:29:47
% DurationCPUTime: 0.32s
% Computational Cost: add. (180->86), mult. (127->89), div. (0->0), fcn. (103->8), ass. (0->26)
t34 = rSges(7,3) + pkin(8);
t13 = sin(qJ(5));
t33 = t13 * pkin(5);
t12 = sin(qJ(6));
t31 = t12 * t13;
t15 = cos(qJ(6));
t30 = t13 * t15;
t29 = pkin(6) + r_base(3);
t14 = sin(qJ(1));
t28 = t14 * pkin(1) + r_base(2);
t17 = cos(qJ(1));
t27 = t17 * pkin(1) + r_base(1);
t11 = qJ(1) + pkin(9);
t7 = sin(t11);
t26 = t7 * pkin(2) + t28;
t25 = qJ(2) + t29;
t24 = t7 * qJ(4) + t26;
t8 = cos(t11);
t23 = t8 * pkin(2) + t7 * qJ(3) + t27;
t22 = pkin(3) + t25;
t21 = t8 * pkin(7) + t24;
t20 = t8 * qJ(4) + t23;
t16 = cos(qJ(5));
t19 = rSges(6,1) * t13 + rSges(6,2) * t16;
t18 = pkin(4) + t22;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t17 * rSges(2,1) - t14 * rSges(2,2) + r_base(1)) + g(2) * (t14 * rSges(2,1) + t17 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t29)) - m(3) * (g(1) * (t8 * rSges(3,1) - t7 * rSges(3,2) + t27) + g(2) * (t7 * rSges(3,1) + t8 * rSges(3,2) + t28) + g(3) * (rSges(3,3) + t25)) - m(4) * (g(1) * (-t8 * rSges(4,2) + t7 * rSges(4,3) + t23) + g(2) * (-t7 * rSges(4,2) + (-rSges(4,3) - qJ(3)) * t8 + t26) + g(3) * (rSges(4,1) + t25)) - m(5) * (g(1) * (t7 * rSges(5,2) + t8 * rSges(5,3) + t20) + g(2) * (t7 * rSges(5,3) + (-rSges(5,2) - qJ(3)) * t8 + t24) + g(3) * (rSges(5,1) + t22)) - m(6) * (g(1) * t20 + g(2) * t21 + g(3) * (t16 * rSges(6,1) - t13 * rSges(6,2) + t18) + (g(1) * t19 + g(2) * (rSges(6,3) - qJ(3))) * t8 + (g(1) * (-rSges(6,3) - pkin(7)) + g(2) * t19) * t7) - m(7) * (g(1) * (t8 * t33 - t7 * pkin(7) + (-t7 * t12 + t8 * t30) * rSges(7,1) + (-t7 * t15 - t8 * t31) * rSges(7,2) + t20) + g(2) * (t7 * t33 - t8 * qJ(3) + (t8 * t12 + t7 * t30) * rSges(7,1) + (t8 * t15 - t7 * t31) * rSges(7,2) + t21) + g(3) * (t34 * t13 + t18) + (g(3) * (rSges(7,1) * t15 - rSges(7,2) * t12 + pkin(5)) - (g(1) * t8 + g(2) * t7) * t34) * t16);
U  = t1;
