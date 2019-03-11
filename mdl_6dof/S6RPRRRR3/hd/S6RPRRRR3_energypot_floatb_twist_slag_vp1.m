% Calculate potential energy for
% S6RPRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 07:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRR3_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR3_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRRRR3_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR3_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR3_energypot_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR3_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR3_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:00:11
% EndTime: 2019-03-09 07:00:12
% DurationCPUTime: 0.51s
% Computational Cost: add. (266->109), mult. (198->127), div. (0->0), fcn. (182->12), ass. (0->40)
t53 = rSges(5,3) + pkin(8);
t26 = -pkin(9) - pkin(8);
t52 = rSges(6,3) - t26;
t51 = rSges(7,3) + pkin(10) - t26;
t21 = sin(qJ(3));
t24 = cos(qJ(3));
t50 = rSges(4,1) * t24 - rSges(4,2) * t21;
t17 = qJ(1) + pkin(11);
t10 = cos(t17);
t9 = sin(t17);
t49 = g(1) * t10 + g(2) * t9;
t20 = sin(qJ(4));
t46 = t9 * t20;
t45 = t9 * t24;
t23 = cos(qJ(4));
t8 = t23 * pkin(4) + pkin(3);
t41 = t10 * t20;
t40 = t10 * t24;
t19 = qJ(4) + qJ(5);
t11 = sin(t19);
t39 = t11 * t24;
t12 = cos(t19);
t38 = t12 * t24;
t37 = t20 * t24;
t36 = t23 * t24;
t33 = pkin(6) + r_base(3);
t22 = sin(qJ(1));
t32 = t22 * pkin(1) + r_base(2);
t25 = cos(qJ(1));
t31 = t25 * pkin(1) + r_base(1);
t30 = t9 * pkin(2) + t32;
t29 = qJ(2) + t33;
t28 = t10 * pkin(2) + t9 * pkin(7) + t31;
t27 = -t10 * pkin(7) + t30;
t13 = qJ(6) + t19;
t7 = cos(t13);
t6 = sin(t13);
t2 = t20 * pkin(4) + pkin(5) * t11;
t1 = pkin(5) * t12 + t8;
t3 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t25 * rSges(2,1) - t22 * rSges(2,2) + r_base(1)) + g(2) * (t22 * rSges(2,1) + t25 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t33)) - m(3) * (g(1) * (t10 * rSges(3,1) - t9 * rSges(3,2) + t31) + g(2) * (t9 * rSges(3,1) + t10 * rSges(3,2) + t32) + g(3) * (rSges(3,3) + t29)) - m(4) * (g(1) * (t9 * rSges(4,3) + t28) + g(2) * (t50 * t9 + t30) + g(3) * (t21 * rSges(4,1) + t24 * rSges(4,2) + t29) + (g(1) * t50 + g(2) * (-rSges(4,3) - pkin(7))) * t10) - m(5) * (g(1) * (pkin(3) * t40 + (t10 * t36 + t46) * rSges(5,1) + (-t10 * t37 + t9 * t23) * rSges(5,2) + t28) + g(2) * (pkin(3) * t45 + (t9 * t36 - t41) * rSges(5,1) + (-t10 * t23 - t9 * t37) * rSges(5,2) + t27) + g(3) * (-t53 * t24 + t29) + (g(3) * (rSges(5,1) * t23 - rSges(5,2) * t20 + pkin(3)) + t49 * t53) * t21) - m(6) * (g(1) * (t8 * t40 + pkin(4) * t46 + (t10 * t38 + t9 * t11) * rSges(6,1) + (-t10 * t39 + t9 * t12) * rSges(6,2) + t28) + g(2) * (t8 * t45 - pkin(4) * t41 + (-t10 * t11 + t9 * t38) * rSges(6,1) + (-t10 * t12 - t9 * t39) * rSges(6,2) + t27) + g(3) * (-t52 * t24 + t29) + (g(3) * (rSges(6,1) * t12 - rSges(6,2) * t11 + t8) + t49 * t52) * t21) - m(7) * (g(1) * (t1 * t40 + t9 * t2 + (t7 * t40 + t9 * t6) * rSges(7,1) + (-t6 * t40 + t9 * t7) * rSges(7,2) + t28) + g(2) * (t1 * t45 - t10 * t2 + (-t10 * t6 + t7 * t45) * rSges(7,1) + (-t10 * t7 - t6 * t45) * rSges(7,2) + t27) + g(3) * (-t51 * t24 + t29) + (g(3) * (rSges(7,1) * t7 - rSges(7,2) * t6 + t1) + t49 * t51) * t21);
U  = t3;
