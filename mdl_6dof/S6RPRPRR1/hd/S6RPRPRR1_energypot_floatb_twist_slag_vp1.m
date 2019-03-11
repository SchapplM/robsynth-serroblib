% Calculate potential energy for
% S6RPRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-03-09 03:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPRR1_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR1_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRPRR1_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR1_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR1_energypot_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR1_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR1_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:33:47
% EndTime: 2019-03-09 03:33:48
% DurationCPUTime: 0.37s
% Computational Cost: add. (244->93), mult. (148->98), div. (0->0), fcn. (124->12), ass. (0->40)
t48 = rSges(7,3) + pkin(9);
t19 = qJ(3) + pkin(11);
t13 = qJ(5) + t19;
t6 = sin(t13);
t47 = rSges(6,2) * t6;
t20 = qJ(1) + pkin(10);
t12 = cos(t20);
t7 = cos(t13);
t46 = t12 * t7;
t45 = rSges(4,3) + pkin(7);
t26 = cos(qJ(3));
t8 = t26 * pkin(3) + pkin(2);
t10 = sin(t20);
t22 = sin(qJ(6));
t43 = t10 * t22;
t25 = cos(qJ(6));
t42 = t10 * t25;
t41 = t12 * t22;
t40 = t12 * t25;
t21 = -qJ(4) - pkin(7);
t39 = rSges(5,3) - t21;
t38 = pkin(6) + r_base(3);
t24 = sin(qJ(1));
t37 = t24 * pkin(1) + r_base(2);
t27 = cos(qJ(1));
t36 = t27 * pkin(1) + r_base(1);
t11 = cos(t19);
t3 = pkin(4) * t11 + t8;
t35 = t12 * t3 + t36;
t34 = qJ(2) + t38;
t18 = -pkin(8) + t21;
t33 = t10 * t3 + t12 * t18 + t37;
t23 = sin(qJ(3));
t32 = t23 * pkin(3) + t34;
t9 = sin(t19);
t31 = pkin(4) * t9 + t32;
t30 = rSges(5,1) * t11 - rSges(5,2) * t9 + t8;
t29 = rSges(4,1) * t26 - rSges(4,2) * t23 + pkin(2);
t28 = g(1) * t36 + g(2) * t37;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t27 * rSges(2,1) - t24 * rSges(2,2) + r_base(1)) + g(2) * (t24 * rSges(2,1) + t27 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t38)) - m(3) * (g(1) * (t12 * rSges(3,1) - t10 * rSges(3,2) + t36) + g(2) * (t10 * rSges(3,1) + t12 * rSges(3,2) + t37) + g(3) * (rSges(3,3) + t34)) - m(4) * (g(3) * (t23 * rSges(4,1) + t26 * rSges(4,2) + t34) + (g(1) * t29 - g(2) * t45) * t12 + (g(1) * t45 + g(2) * t29) * t10 + t28) - m(5) * (g(3) * (t9 * rSges(5,1) + t11 * rSges(5,2) + t32) + (g(1) * t30 - g(2) * t39) * t12 + (g(1) * t39 + g(2) * t30) * t10 + t28) - m(6) * (g(1) * (rSges(6,1) * t46 - t12 * t47 + t35) + g(2) * (-t12 * rSges(6,3) + t33) + g(3) * (t6 * rSges(6,1) + t7 * rSges(6,2) + t31) + (g(1) * (rSges(6,3) - t18) + g(2) * (rSges(6,1) * t7 - t47)) * t10) - m(7) * (g(1) * (pkin(5) * t46 - t10 * t18 + (t40 * t7 + t43) * rSges(7,1) + (-t41 * t7 + t42) * rSges(7,2) + t35) + g(2) * (t10 * t7 * pkin(5) + (t42 * t7 - t41) * rSges(7,1) + (-t43 * t7 - t40) * rSges(7,2) + t33) + g(3) * (-t48 * t7 + t31) + (g(3) * (rSges(7,1) * t25 - rSges(7,2) * t22 + pkin(5)) + (g(1) * t12 + g(2) * t10) * t48) * t6);
U  = t1;
