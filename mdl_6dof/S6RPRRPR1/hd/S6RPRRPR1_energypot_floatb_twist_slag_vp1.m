% Calculate potential energy for
% S6RPRRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 04:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRPR1_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR1_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRRPR1_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR1_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR1_energypot_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR1_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR1_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:57:27
% EndTime: 2019-03-09 04:57:27
% DurationCPUTime: 0.35s
% Computational Cost: add. (244->93), mult. (148->98), div. (0->0), fcn. (124->12), ass. (0->40)
t48 = rSges(7,3) + pkin(9);
t27 = -pkin(8) - pkin(7);
t20 = qJ(3) + qJ(4);
t11 = pkin(11) + t20;
t5 = sin(t11);
t47 = rSges(6,2) * t5;
t19 = qJ(1) + pkin(10);
t10 = cos(t19);
t6 = cos(t11);
t46 = t10 * t6;
t21 = sin(qJ(6));
t9 = sin(t19);
t45 = t9 * t21;
t24 = cos(qJ(6));
t44 = t9 * t24;
t43 = rSges(4,3) + pkin(7);
t25 = cos(qJ(3));
t8 = t25 * pkin(3) + pkin(2);
t41 = t10 * t21;
t40 = t10 * t24;
t39 = rSges(5,3) - t27;
t38 = pkin(6) + r_base(3);
t23 = sin(qJ(1));
t37 = t23 * pkin(1) + r_base(2);
t26 = cos(qJ(1));
t36 = t26 * pkin(1) + r_base(1);
t13 = cos(t20);
t3 = pkin(4) * t13 + t8;
t35 = t10 * t3 + t36;
t34 = qJ(2) + t38;
t18 = -qJ(5) + t27;
t33 = t10 * t18 + t9 * t3 + t37;
t22 = sin(qJ(3));
t32 = t22 * pkin(3) + t34;
t12 = sin(t20);
t31 = pkin(4) * t12 + t32;
t30 = rSges(4,1) * t25 - rSges(4,2) * t22 + pkin(2);
t29 = rSges(5,1) * t13 - rSges(5,2) * t12 + t8;
t28 = g(1) * t36 + g(2) * t37;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t26 * rSges(2,1) - t23 * rSges(2,2) + r_base(1)) + g(2) * (t23 * rSges(2,1) + t26 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t38)) - m(3) * (g(1) * (t10 * rSges(3,1) - t9 * rSges(3,2) + t36) + g(2) * (t9 * rSges(3,1) + t10 * rSges(3,2) + t37) + g(3) * (rSges(3,3) + t34)) - m(4) * (g(3) * (t22 * rSges(4,1) + t25 * rSges(4,2) + t34) + (g(1) * t43 + g(2) * t30) * t9 + (g(1) * t30 - g(2) * t43) * t10 + t28) - m(5) * (g(3) * (t12 * rSges(5,1) + t13 * rSges(5,2) + t32) + (g(1) * t39 + g(2) * t29) * t9 + (g(1) * t29 - g(2) * t39) * t10 + t28) - m(6) * (g(1) * (rSges(6,1) * t46 - t10 * t47 + t35) + g(2) * (-t10 * rSges(6,3) + t33) + g(3) * (t5 * rSges(6,1) + t6 * rSges(6,2) + t31) + (g(1) * (rSges(6,3) - t18) + g(2) * (rSges(6,1) * t6 - t47)) * t9) - m(7) * (g(1) * (pkin(5) * t46 - t9 * t18 + (t40 * t6 + t45) * rSges(7,1) + (-t41 * t6 + t44) * rSges(7,2) + t35) + g(2) * (t9 * t6 * pkin(5) + (t44 * t6 - t41) * rSges(7,1) + (-t45 * t6 - t40) * rSges(7,2) + t33) + g(3) * (-t48 * t6 + t31) + (g(3) * (rSges(7,1) * t24 - rSges(7,2) * t21 + pkin(5)) + (g(1) * t10 + g(2) * t9) * t48) * t5);
U  = t1;
