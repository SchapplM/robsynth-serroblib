% Calculate potential energy for
% S6RPRRRR6
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
% Datum: 2019-03-09 07:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRR6_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR6_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRRRR6_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR6_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR6_energypot_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR6_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR6_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:12:09
% EndTime: 2019-03-09 07:12:09
% DurationCPUTime: 0.54s
% Computational Cost: add. (257->114), mult. (212->130), div. (0->0), fcn. (196->12), ass. (0->44)
t54 = rSges(5,3) + pkin(8);
t27 = -pkin(9) - pkin(8);
t53 = rSges(6,3) - t27;
t52 = rSges(7,3) + pkin(10) - t27;
t24 = sin(qJ(1));
t26 = cos(qJ(1));
t51 = g(1) * t26 + g(2) * t24;
t25 = cos(qJ(4));
t9 = t25 * pkin(4) + pkin(3);
t17 = pkin(11) + qJ(3);
t10 = sin(t17);
t47 = rSges(4,2) * t10;
t11 = cos(t17);
t46 = t11 * t24;
t45 = t11 * t26;
t19 = qJ(4) + qJ(5);
t12 = sin(t19);
t44 = t12 * t24;
t43 = t12 * t26;
t13 = cos(t19);
t42 = t13 * t24;
t41 = t13 * t26;
t23 = sin(qJ(4));
t40 = t23 * t24;
t39 = t23 * t26;
t38 = t24 * t25;
t37 = t25 * t26;
t34 = rSges(3,3) + qJ(2);
t33 = pkin(6) + r_base(3);
t21 = cos(pkin(11));
t5 = pkin(2) * t21 + pkin(1);
t32 = t26 * t5 + r_base(1);
t22 = -pkin(7) - qJ(2);
t31 = t26 * t22 + t24 * t5 + r_base(2);
t20 = sin(pkin(11));
t30 = t20 * pkin(2) + t33;
t29 = -t24 * t22 + t32;
t28 = rSges(3,1) * t21 - rSges(3,2) * t20 + pkin(1);
t15 = qJ(6) + t19;
t8 = cos(t15);
t7 = sin(t15);
t2 = pkin(4) * t23 + pkin(5) * t12;
t1 = pkin(5) * t13 + t9;
t3 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t26 - rSges(2,2) * t24 + r_base(1)) + g(2) * (rSges(2,1) * t24 + rSges(2,2) * t26 + r_base(2)) + g(3) * (rSges(2,3) + t33)) - m(3) * (g(1) * r_base(1) + g(2) * r_base(2) + g(3) * (rSges(3,1) * t20 + rSges(3,2) * t21 + t33) + (g(1) * t28 - g(2) * t34) * t26 + (g(1) * t34 + g(2) * t28) * t24) - m(4) * (g(1) * (rSges(4,1) * t45 - t26 * t47 + t32) + g(2) * (-t26 * rSges(4,3) + t31) + g(3) * (rSges(4,1) * t10 + rSges(4,2) * t11 + t30) + (g(1) * (rSges(4,3) - t22) + g(2) * (rSges(4,1) * t11 - t47)) * t24) - m(5) * (g(1) * (pkin(3) * t45 + (t11 * t37 + t40) * rSges(5,1) + (-t11 * t39 + t38) * rSges(5,2) + t29) + g(2) * (pkin(3) * t46 + (t11 * t38 - t39) * rSges(5,1) + (-t11 * t40 - t37) * rSges(5,2) + t31) + g(3) * (-t54 * t11 + t30) + (g(3) * (rSges(5,1) * t25 - rSges(5,2) * t23 + pkin(3)) + t51 * t54) * t10) - m(6) * (g(1) * (t9 * t45 + pkin(4) * t40 + (t11 * t41 + t44) * rSges(6,1) + (-t11 * t43 + t42) * rSges(6,2) + t29) + g(2) * (t9 * t46 - pkin(4) * t39 + (t11 * t42 - t43) * rSges(6,1) + (-t11 * t44 - t41) * rSges(6,2) + t31) + g(3) * (-t53 * t11 + t30) + (g(3) * (rSges(6,1) * t13 - rSges(6,2) * t12 + t9) + t51 * t53) * t10) - m(7) * (g(1) * (t1 * t45 + t24 * t2 + (t24 * t7 + t8 * t45) * rSges(7,1) + (t24 * t8 - t7 * t45) * rSges(7,2) + t29) + g(2) * (t1 * t46 - t26 * t2 + (-t26 * t7 + t8 * t46) * rSges(7,1) + (-t26 * t8 - t7 * t46) * rSges(7,2) + t31) + g(3) * (-t52 * t11 + t30) + (g(3) * (rSges(7,1) * t8 - rSges(7,2) * t7 + t1) + t51 * t52) * t10);
U  = t3;
