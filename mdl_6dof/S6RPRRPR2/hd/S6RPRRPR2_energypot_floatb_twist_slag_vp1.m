% Calculate potential energy for
% S6RPRRPR2
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
% Datum: 2019-03-09 05:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRPR2_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR2_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRRPR2_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR2_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR2_energypot_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR2_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR2_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:00:08
% EndTime: 2019-03-09 05:00:08
% DurationCPUTime: 0.53s
% Computational Cost: add. (266->110), mult. (198->126), div. (0->0), fcn. (182->12), ass. (0->38)
t49 = rSges(5,3) + pkin(8);
t20 = -qJ(5) - pkin(8);
t48 = rSges(6,3) - t20;
t47 = rSges(7,3) + pkin(9) - t20;
t19 = qJ(1) + pkin(10);
t10 = sin(t19);
t12 = cos(t19);
t46 = g(1) * t12 + g(2) * t10;
t24 = cos(qJ(4));
t8 = t24 * pkin(4) + pkin(3);
t22 = sin(qJ(3));
t42 = rSges(4,2) * t22;
t21 = sin(qJ(4));
t41 = t10 * t21;
t25 = cos(qJ(3));
t40 = t10 * t25;
t39 = t12 * t21;
t38 = t12 * t25;
t37 = t21 * t25;
t36 = t24 * t25;
t33 = pkin(6) + r_base(3);
t18 = qJ(4) + pkin(11);
t23 = sin(qJ(1));
t32 = t23 * pkin(1) + r_base(2);
t26 = cos(qJ(1));
t31 = t26 * pkin(1) + r_base(1);
t30 = t10 * pkin(2) + t32;
t29 = qJ(2) + t33;
t28 = t12 * pkin(2) + t10 * pkin(7) + t31;
t27 = -t12 * pkin(7) + t30;
t13 = qJ(6) + t18;
t11 = cos(t18);
t9 = sin(t18);
t7 = cos(t13);
t6 = sin(t13);
t2 = t21 * pkin(4) + pkin(5) * t9;
t1 = pkin(5) * t11 + t8;
t3 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t26 * rSges(2,1) - t23 * rSges(2,2) + r_base(1)) + g(2) * (t23 * rSges(2,1) + t26 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t33)) - m(3) * (g(1) * (t12 * rSges(3,1) - t10 * rSges(3,2) + t31) + g(2) * (t10 * rSges(3,1) + t12 * rSges(3,2) + t32) + g(3) * (rSges(3,3) + t29)) - m(4) * (g(1) * (t10 * rSges(4,3) + t28) + g(2) * (rSges(4,1) * t40 - t10 * t42 + t30) + g(3) * (t22 * rSges(4,1) + t25 * rSges(4,2) + t29) + (g(1) * (rSges(4,1) * t25 - t42) + g(2) * (-rSges(4,3) - pkin(7))) * t12) - m(5) * (g(1) * (pkin(3) * t38 + (t12 * t36 + t41) * rSges(5,1) + (t10 * t24 - t12 * t37) * rSges(5,2) + t28) + g(2) * (pkin(3) * t40 + (t10 * t36 - t39) * rSges(5,1) + (-t10 * t37 - t12 * t24) * rSges(5,2) + t27) + g(3) * (-t49 * t25 + t29) + (g(3) * (rSges(5,1) * t24 - rSges(5,2) * t21 + pkin(3)) + t46 * t49) * t22) - m(6) * (g(1) * (t8 * t38 + pkin(4) * t41 + (t10 * t9 + t11 * t38) * rSges(6,1) + (t10 * t11 - t9 * t38) * rSges(6,2) + t28) + g(2) * (t8 * t40 - pkin(4) * t39 + (t11 * t40 - t12 * t9) * rSges(6,1) + (-t12 * t11 - t9 * t40) * rSges(6,2) + t27) + g(3) * (-t48 * t25 + t29) + (g(3) * (rSges(6,1) * t11 - rSges(6,2) * t9 + t8) + t46 * t48) * t22) - m(7) * (g(1) * (t1 * t38 + t10 * t2 + (t10 * t6 + t7 * t38) * rSges(7,1) + (t10 * t7 - t6 * t38) * rSges(7,2) + t28) + g(2) * (t1 * t40 - t12 * t2 + (-t12 * t6 + t7 * t40) * rSges(7,1) + (-t12 * t7 - t6 * t40) * rSges(7,2) + t27) + g(3) * (-t47 * t25 + t29) + (g(3) * (rSges(7,1) * t7 - rSges(7,2) * t6 + t1) + t46 * t47) * t22);
U  = t3;
