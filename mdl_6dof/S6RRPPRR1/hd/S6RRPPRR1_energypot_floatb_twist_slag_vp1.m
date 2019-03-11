% Calculate potential energy for
% S6RRPPRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
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
% Datum: 2019-03-09 08:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPPRR1_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR1_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPPRR1_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR1_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR1_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR1_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR1_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:45:42
% EndTime: 2019-03-09 08:45:42
% DurationCPUTime: 0.42s
% Computational Cost: add. (252->104), mult. (247->116), div. (0->0), fcn. (244->10), ass. (0->40)
t52 = pkin(9) + rSges(7,3);
t54 = rSges(3,3) + pkin(7);
t26 = -qJ(3) - pkin(7);
t53 = -pkin(8) - t26;
t25 = qJ(2) + pkin(10);
t21 = sin(t25);
t32 = cos(qJ(5));
t51 = t21 * t32;
t34 = cos(qJ(1));
t50 = t21 * t34;
t22 = cos(t25);
t30 = sin(qJ(1));
t49 = t22 * t30;
t48 = t22 * t34;
t47 = qJ(4) * t21;
t46 = pkin(6) + r_base(3);
t33 = cos(qJ(2));
t20 = t33 * pkin(2) + pkin(1);
t45 = t34 * t20 + r_base(1);
t29 = sin(qJ(2));
t44 = t29 * pkin(2) + t46;
t43 = t30 * t20 + t34 * t26 + r_base(2);
t42 = t21 * pkin(3) + t44;
t41 = pkin(3) * t48 + t34 * t47 + t45;
t28 = sin(qJ(5));
t5 = t21 * t28 + t22 * t32;
t40 = pkin(3) * t49 + t30 * t47 + t43;
t39 = pkin(4) * t48 + t41;
t38 = rSges(3,1) * t33 - rSges(3,2) * t29 + pkin(1);
t27 = sin(qJ(6));
t31 = cos(qJ(6));
t37 = t31 * rSges(7,1) - t27 * rSges(7,2) + pkin(5);
t36 = pkin(4) * t49 + t34 * pkin(8) + t40;
t35 = t21 * pkin(4) - t22 * qJ(4) + t42;
t6 = -t22 * t28 + t51;
t4 = t5 * t34;
t3 = t28 * t48 - t32 * t50;
t2 = t5 * t30;
t1 = t28 * t49 - t30 * t51;
t7 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t34 * rSges(2,1) - t30 * rSges(2,2) + r_base(1)) + g(2) * (t30 * rSges(2,1) + t34 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t46)) - m(3) * (g(1) * r_base(1) + g(2) * r_base(2) + g(3) * (t29 * rSges(3,1) + t33 * rSges(3,2) + t46) + (g(1) * t38 - g(2) * t54) * t34 + (g(1) * t54 + g(2) * t38) * t30) - m(4) * (g(1) * (rSges(4,1) * t48 - rSges(4,2) * t50 + t45) + g(2) * (-t34 * rSges(4,3) + t43) + g(3) * (t21 * rSges(4,1) + t22 * rSges(4,2) + t44) + (g(1) * (rSges(4,3) - t26) + g(2) * (rSges(4,1) * t22 - rSges(4,2) * t21)) * t30) - m(5) * (g(1) * (rSges(5,1) * t48 + rSges(5,3) * t50 + t41) + g(2) * (-t34 * rSges(5,2) + t40) + g(3) * (t21 * rSges(5,1) + (-rSges(5,3) - qJ(4)) * t22 + t42) + (g(1) * (rSges(5,2) - t26) + g(2) * (rSges(5,1) * t22 + rSges(5,3) * t21)) * t30) - m(6) * (g(2) * (t2 * rSges(6,1) - t1 * rSges(6,2) + t34 * rSges(6,3) + t36) + g(3) * (t6 * rSges(6,1) - t5 * rSges(6,2) + t35) + (t4 * rSges(6,1) - t3 * rSges(6,2) + t39 + (-rSges(6,3) + t53) * t30) * g(1)) - m(7) * (g(1) * (t37 * t4 + (-t27 * rSges(7,1) - t31 * rSges(7,2) + t53) * t30 + t52 * t3 + t39) + g(2) * (t2 * pkin(5) + (t2 * t31 + t34 * t27) * rSges(7,1) + (-t2 * t27 + t34 * t31) * rSges(7,2) + t52 * t1 + t36) + (t37 * t6 + t52 * t5 + t35) * g(3));
U  = t7;
