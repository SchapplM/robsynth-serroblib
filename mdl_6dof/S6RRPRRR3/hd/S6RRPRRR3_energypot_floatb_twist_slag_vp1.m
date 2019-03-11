% Calculate potential energy for
% S6RRPRRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 13:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRR3_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR3_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPRRR3_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR3_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR3_energypot_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR3_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR3_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:21:39
% EndTime: 2019-03-09 13:21:40
% DurationCPUTime: 0.54s
% Computational Cost: add. (257->114), mult. (212->130), div. (0->0), fcn. (196->12), ass. (0->44)
t54 = rSges(5,3) + pkin(8);
t27 = -pkin(9) - pkin(8);
t53 = rSges(6,3) - t27;
t52 = rSges(7,3) + pkin(10) - t27;
t23 = sin(qJ(1));
t26 = cos(qJ(1));
t51 = g(1) * t26 + g(2) * t23;
t48 = rSges(3,3) + pkin(7);
t24 = cos(qJ(4));
t8 = t24 * pkin(4) + pkin(3);
t17 = qJ(2) + pkin(11);
t10 = sin(t17);
t46 = rSges(4,2) * t10;
t11 = cos(t17);
t45 = t11 * t23;
t44 = t11 * t26;
t19 = qJ(4) + qJ(5);
t12 = sin(t19);
t43 = t12 * t23;
t42 = t12 * t26;
t13 = cos(t19);
t41 = t13 * t23;
t40 = t13 * t26;
t21 = sin(qJ(4));
t39 = t21 * t23;
t38 = t21 * t26;
t37 = t23 * t24;
t36 = t24 * t26;
t33 = pkin(6) + r_base(3);
t25 = cos(qJ(2));
t9 = pkin(2) * t25 + pkin(1);
t32 = t26 * t9 + r_base(1);
t20 = -qJ(3) - pkin(7);
t31 = t26 * t20 + t23 * t9 + r_base(2);
t22 = sin(qJ(2));
t30 = t22 * pkin(2) + t33;
t29 = -t23 * t20 + t32;
t28 = rSges(3,1) * t25 - rSges(3,2) * t22 + pkin(1);
t14 = qJ(6) + t19;
t7 = cos(t14);
t6 = sin(t14);
t2 = pkin(4) * t21 + pkin(5) * t12;
t1 = pkin(5) * t13 + t8;
t3 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t26 - rSges(2,2) * t23 + r_base(1)) + g(2) * (rSges(2,1) * t23 + rSges(2,2) * t26 + r_base(2)) + g(3) * (rSges(2,3) + t33)) - m(3) * (g(1) * r_base(1) + g(2) * r_base(2) + g(3) * (rSges(3,1) * t22 + rSges(3,2) * t25 + t33) + (g(1) * t28 - g(2) * t48) * t26 + (g(1) * t48 + g(2) * t28) * t23) - m(4) * (g(1) * (rSges(4,1) * t44 - t26 * t46 + t32) + g(2) * (-t26 * rSges(4,3) + t31) + g(3) * (rSges(4,1) * t10 + rSges(4,2) * t11 + t30) + (g(1) * (rSges(4,3) - t20) + g(2) * (rSges(4,1) * t11 - t46)) * t23) - m(5) * (g(1) * (pkin(3) * t44 + (t11 * t36 + t39) * rSges(5,1) + (-t11 * t38 + t37) * rSges(5,2) + t29) + g(2) * (pkin(3) * t45 + (t11 * t37 - t38) * rSges(5,1) + (-t11 * t39 - t36) * rSges(5,2) + t31) + g(3) * (-t54 * t11 + t30) + (g(3) * (rSges(5,1) * t24 - rSges(5,2) * t21 + pkin(3)) + t51 * t54) * t10) - m(6) * (g(1) * (t8 * t44 + pkin(4) * t39 + (t11 * t40 + t43) * rSges(6,1) + (-t11 * t42 + t41) * rSges(6,2) + t29) + g(2) * (t8 * t45 - pkin(4) * t38 + (t11 * t41 - t42) * rSges(6,1) + (-t11 * t43 - t40) * rSges(6,2) + t31) + g(3) * (-t53 * t11 + t30) + (g(3) * (rSges(6,1) * t13 - rSges(6,2) * t12 + t8) + t51 * t53) * t10) - m(7) * (g(1) * (t1 * t44 + t23 * t2 + (t23 * t6 + t7 * t44) * rSges(7,1) + (t23 * t7 - t6 * t44) * rSges(7,2) + t29) + g(2) * (t1 * t45 - t26 * t2 + (-t26 * t6 + t7 * t45) * rSges(7,1) + (-t26 * t7 - t6 * t45) * rSges(7,2) + t31) + g(3) * (-t52 * t11 + t30) + (g(3) * (rSges(7,1) * t7 - rSges(7,2) * t6 + t1) + t51 * t52) * t10);
U  = t3;
