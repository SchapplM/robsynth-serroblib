% Calculate potential energy for
% S5PRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRPRP6_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP6_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PRPRP6_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP6_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRPRP6_energypot_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP6_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRP6_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:40:16
% EndTime: 2019-12-05 15:40:16
% DurationCPUTime: 0.35s
% Computational Cost: add. (126->84), mult. (182->97), div. (0->0), fcn. (170->6), ass. (0->29)
t41 = rSges(6,1) + pkin(4);
t40 = rSges(6,3) + qJ(5);
t18 = sin(pkin(7));
t21 = sin(qJ(2));
t39 = t18 * t21;
t23 = cos(qJ(2));
t38 = t18 * t23;
t19 = cos(pkin(7));
t37 = t19 * t23;
t20 = sin(qJ(4));
t36 = t20 * t21;
t22 = cos(qJ(4));
t35 = t21 * t22;
t34 = qJ(3) * t21;
t33 = t18 * pkin(1) + r_base(2);
t32 = qJ(1) + r_base(3);
t31 = t19 * pkin(1) + t18 * pkin(5) + r_base(1);
t30 = t21 * pkin(2) + t32;
t29 = pkin(2) * t38 + t18 * t34 + t33;
t28 = t21 * pkin(6) + t30;
t27 = g(1) * t19 + g(2) * t18;
t26 = pkin(2) * t37 + t19 * t34 + t31;
t25 = t18 * pkin(3) + pkin(6) * t37 + t26;
t24 = t29 + pkin(6) * t38 + (-pkin(3) - pkin(5)) * t19;
t4 = t18 * t36 - t19 * t22;
t3 = t18 * t35 + t19 * t20;
t2 = t18 * t22 + t19 * t36;
t1 = t18 * t20 - t19 * t35;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t19 * rSges(2,1) - t18 * rSges(2,2) + r_base(1)) + g(2) * (t18 * rSges(2,1) + t19 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t32)) - m(3) * (g(1) * (t18 * rSges(3,3) + t31) + g(2) * (rSges(3,1) * t38 - rSges(3,2) * t39 + t33) + g(3) * (t21 * rSges(3,1) + t23 * rSges(3,2) + t32) + (g(1) * (rSges(3,1) * t23 - rSges(3,2) * t21) + g(2) * (-rSges(3,3) - pkin(5))) * t19) - m(4) * (g(1) * (t18 * rSges(4,1) + t26) + g(2) * (-rSges(4,2) * t38 + rSges(4,3) * t39 + t29) + g(3) * (-t21 * rSges(4,2) + (-rSges(4,3) - qJ(3)) * t23 + t30) + (g(1) * (-rSges(4,2) * t23 + rSges(4,3) * t21) + g(2) * (-rSges(4,1) - pkin(5))) * t19) - m(5) * (g(1) * (t2 * rSges(5,1) - t1 * rSges(5,2) + t25) + g(2) * (t4 * rSges(5,1) + t3 * rSges(5,2) + t24) + g(3) * (t21 * rSges(5,3) + t28) + (g(3) * (-rSges(5,1) * t20 - rSges(5,2) * t22 - qJ(3)) + t27 * rSges(5,3)) * t23) - m(6) * (g(1) * (t40 * t1 + t41 * t2 + t25) + g(2) * (-t40 * t3 + t41 * t4 + t24) + g(3) * (t21 * rSges(6,2) + t28) + (g(3) * (-t41 * t20 + t40 * t22 - qJ(3)) + t27 * rSges(6,2)) * t23);
U = t5;
