% Calculate potential energy for
% S6RPPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
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
% Datum: 2019-03-09 02:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPRRP6_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP6_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPPRRP6_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP6_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP6_energypot_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP6_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRP6_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:10:07
% EndTime: 2019-03-09 02:10:07
% DurationCPUTime: 0.35s
% Computational Cost: add. (141->89), mult. (176->93), div. (0->0), fcn. (160->6), ass. (0->31)
t45 = rSges(7,1) + pkin(5);
t18 = sin(qJ(1));
t21 = cos(qJ(1));
t44 = g(1) * t21 + g(2) * t18;
t43 = rSges(7,3) + qJ(6);
t17 = sin(qJ(4));
t42 = pkin(4) * t17;
t16 = sin(qJ(5));
t37 = t18 * t16;
t19 = cos(qJ(5));
t36 = t18 * t19;
t35 = t21 * t16;
t34 = t21 * t19;
t33 = pkin(6) + r_base(3);
t32 = t18 * pkin(1) + r_base(2);
t31 = pkin(2) + t33;
t30 = t18 * qJ(3) + t32;
t29 = t21 * pkin(1) + t18 * qJ(2) + r_base(1);
t28 = pkin(3) + t31;
t27 = t21 * qJ(3) + t29;
t26 = t21 * pkin(7) + t30;
t20 = cos(qJ(4));
t25 = rSges(5,1) * t17 + rSges(5,2) * t20;
t24 = t20 * pkin(4) + t17 * pkin(8) + t28;
t23 = -t18 * pkin(7) + t21 * t42 + t27;
t22 = -t21 * qJ(2) + t18 * t42 + t26;
t4 = t17 * t34 - t37;
t3 = t17 * t35 + t36;
t2 = t17 * t36 + t35;
t1 = t17 * t37 - t34;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t21 * rSges(2,1) - t18 * rSges(2,2) + r_base(1)) + g(2) * (t18 * rSges(2,1) + t21 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t33)) - m(3) * (g(1) * (-t21 * rSges(3,2) + t18 * rSges(3,3) + t29) + g(2) * (-t18 * rSges(3,2) + (-rSges(3,3) - qJ(2)) * t21 + t32) + g(3) * (rSges(3,1) + t33)) - m(4) * (g(1) * (t18 * rSges(4,2) + t21 * rSges(4,3) + t27) + g(2) * (t18 * rSges(4,3) + (-rSges(4,2) - qJ(2)) * t21 + t30) + g(3) * (rSges(4,1) + t31)) - m(5) * (g(1) * t27 + g(2) * t26 + g(3) * (t20 * rSges(5,1) - t17 * rSges(5,2) + t28) + (g(1) * t25 + g(2) * (rSges(5,3) - qJ(2))) * t21 + (g(1) * (-rSges(5,3) - pkin(7)) + g(2) * t25) * t18) - m(6) * (g(1) * (t4 * rSges(6,1) - t3 * rSges(6,2) + t23) + g(2) * (t2 * rSges(6,1) - t1 * rSges(6,2) + t22) + g(3) * (t17 * rSges(6,3) + t24) + (g(3) * (rSges(6,1) * t19 - rSges(6,2) * t16) + t44 * (-rSges(6,3) - pkin(8))) * t20) - m(7) * (g(1) * (t43 * t3 + t45 * t4 + t23) + g(2) * (t43 * t1 + t45 * t2 + t22) + g(3) * (t17 * rSges(7,2) + t24) + (g(3) * (t43 * t16 + t45 * t19) + t44 * (-rSges(7,2) - pkin(8))) * t20);
U  = t5;
