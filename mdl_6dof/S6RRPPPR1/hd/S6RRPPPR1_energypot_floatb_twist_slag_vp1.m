% Calculate potential energy for
% S6RRPPPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta4]';
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
% Datum: 2019-03-09 08:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPPPR1_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR1_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPPPR1_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR1_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR1_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR1_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPPR1_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:05:58
% EndTime: 2019-03-09 08:05:59
% DurationCPUTime: 0.56s
% Computational Cost: add. (260->111), mult. (274->129), div. (0->0), fcn. (278->10), ass. (0->40)
t53 = rSges(3,3) + pkin(7);
t52 = -rSges(7,3) - pkin(8);
t22 = qJ(2) + pkin(9);
t20 = cos(t22);
t31 = cos(qJ(1));
t51 = t20 * t31;
t19 = sin(t22);
t28 = sin(qJ(1));
t50 = t28 * t19;
t23 = sin(pkin(10));
t49 = t28 * t23;
t24 = cos(pkin(10));
t48 = t28 * t24;
t47 = t31 * t19;
t46 = t31 * t23;
t45 = t31 * t24;
t44 = qJ(4) * t19;
t43 = rSges(6,3) + qJ(5);
t42 = pkin(6) + r_base(3);
t30 = cos(qJ(2));
t18 = t30 * pkin(2) + pkin(1);
t41 = t31 * t18 + r_base(1);
t27 = sin(qJ(2));
t40 = t27 * pkin(2) + t42;
t25 = -qJ(3) - pkin(7);
t39 = t28 * t18 + t31 * t25 + r_base(2);
t38 = t19 * pkin(3) + t40;
t37 = t39 + (pkin(3) * t20 + t44) * t28;
t36 = rSges(3,1) * t30 - rSges(3,2) * t27 + pkin(1);
t35 = t38 + (pkin(4) * t24 + qJ(5) * t23) * t19;
t4 = t20 * t48 - t46;
t34 = t4 * pkin(4) + t37;
t33 = pkin(3) * t51 - t28 * t25 + t31 * t44 + t41;
t6 = t20 * t45 + t49;
t32 = t6 * pkin(4) + t33;
t29 = cos(qJ(6));
t26 = sin(qJ(6));
t5 = t20 * t46 - t48;
t3 = t20 * t49 + t45;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t31 * rSges(2,1) - t28 * rSges(2,2) + r_base(1)) + g(2) * (t28 * rSges(2,1) + t31 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t42)) - m(3) * (g(1) * r_base(1) + g(2) * r_base(2) + g(3) * (t27 * rSges(3,1) + t30 * rSges(3,2) + t42) + (g(1) * t36 - g(2) * t53) * t31 + (g(1) * t53 + g(2) * t36) * t28) - m(4) * (g(1) * (rSges(4,1) * t51 - rSges(4,2) * t47 + t41) + g(2) * (-t31 * rSges(4,3) + t39) + g(3) * (t19 * rSges(4,1) + t20 * rSges(4,2) + t40) + (g(1) * (rSges(4,3) - t25) + g(2) * (rSges(4,1) * t20 - rSges(4,2) * t19)) * t28) - m(5) * (g(1) * (t6 * rSges(5,1) - t5 * rSges(5,2) + rSges(5,3) * t47 + t33) + g(2) * (t4 * rSges(5,1) - t3 * rSges(5,2) + rSges(5,3) * t50 + t37) + g(3) * ((-rSges(5,3) - qJ(4)) * t20 + (rSges(5,1) * t24 - rSges(5,2) * t23) * t19 + t38)) - m(6) * (g(1) * (t6 * rSges(6,1) + rSges(6,2) * t47 + t43 * t5 + t32) + g(2) * (t4 * rSges(6,1) + rSges(6,2) * t50 + t43 * t3 + t34) + g(3) * ((-rSges(6,2) - qJ(4)) * t20 + (rSges(6,1) * t24 + rSges(6,3) * t23) * t19 + t35)) - m(7) * (g(1) * (t6 * pkin(5) + t5 * qJ(5) + (t5 * t26 + t6 * t29) * rSges(7,1) + (-t6 * t26 + t5 * t29) * rSges(7,2) + t32) + g(2) * (t4 * pkin(5) + t3 * qJ(5) + (t3 * t26 + t4 * t29) * rSges(7,1) + (-t4 * t26 + t3 * t29) * rSges(7,2) + t34) + (g(1) * t31 + g(2) * t28) * t19 * t52 + (t35 + (-qJ(4) - t52) * t20 + (t24 * pkin(5) + (t23 * t26 + t24 * t29) * rSges(7,1) + (t23 * t29 - t24 * t26) * rSges(7,2)) * t19) * g(3));
U  = t1;
