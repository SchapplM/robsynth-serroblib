% Calculate potential energy for
% S6RRRRRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 03:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRRR3_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR3_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRRRR3_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR3_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR3_energypot_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR3_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRR3_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:39:06
% EndTime: 2019-03-10 03:39:07
% DurationCPUTime: 0.54s
% Computational Cost: add. (257->114), mult. (212->130), div. (0->0), fcn. (196->12), ass. (0->40)
t50 = rSges(5,3) + pkin(9);
t26 = -pkin(10) - pkin(9);
t49 = rSges(6,3) - t26;
t48 = rSges(7,3) + pkin(11) - t26;
t22 = sin(qJ(1));
t25 = cos(qJ(1));
t47 = g(1) * t25 + g(2) * t22;
t44 = rSges(3,3) + pkin(7);
t23 = cos(qJ(4));
t7 = t23 * pkin(4) + pkin(3);
t19 = qJ(2) + qJ(3);
t11 = sin(t19);
t42 = rSges(4,2) * t11;
t13 = cos(t19);
t41 = t13 * t22;
t40 = t13 * t25;
t20 = sin(qJ(4));
t39 = t20 * t22;
t38 = t20 * t25;
t37 = t22 * t23;
t36 = t23 * t25;
t18 = qJ(4) + qJ(5);
t33 = pkin(6) + r_base(3);
t24 = cos(qJ(2));
t8 = pkin(2) * t24 + pkin(1);
t32 = t25 * t8 + r_base(1);
t27 = -pkin(8) - pkin(7);
t31 = t22 * t8 + t25 * t27 + r_base(2);
t21 = sin(qJ(2));
t30 = t21 * pkin(2) + t33;
t29 = -t22 * t27 + t32;
t28 = rSges(3,1) * t24 - rSges(3,2) * t21 + pkin(1);
t14 = qJ(6) + t18;
t12 = cos(t18);
t10 = sin(t18);
t6 = cos(t14);
t5 = sin(t14);
t2 = pkin(4) * t20 + pkin(5) * t10;
t1 = pkin(5) * t12 + t7;
t3 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t25 - rSges(2,2) * t22 + r_base(1)) + g(2) * (rSges(2,1) * t22 + rSges(2,2) * t25 + r_base(2)) + g(3) * (rSges(2,3) + t33)) - m(3) * (g(1) * r_base(1) + g(2) * r_base(2) + g(3) * (rSges(3,1) * t21 + rSges(3,2) * t24 + t33) + (g(1) * t28 - g(2) * t44) * t25 + (g(1) * t44 + g(2) * t28) * t22) - m(4) * (g(1) * (rSges(4,1) * t40 - t25 * t42 + t32) + g(2) * (-t25 * rSges(4,3) + t31) + g(3) * (rSges(4,1) * t11 + rSges(4,2) * t13 + t30) + (g(1) * (rSges(4,3) - t27) + g(2) * (rSges(4,1) * t13 - t42)) * t22) - m(5) * (g(1) * (pkin(3) * t40 + (t13 * t36 + t39) * rSges(5,1) + (-t13 * t38 + t37) * rSges(5,2) + t29) + g(2) * (pkin(3) * t41 + (t13 * t37 - t38) * rSges(5,1) + (-t13 * t39 - t36) * rSges(5,2) + t31) + g(3) * (-t50 * t13 + t30) + (g(3) * (rSges(5,1) * t23 - rSges(5,2) * t20 + pkin(3)) + t47 * t50) * t11) - m(6) * (g(1) * (t7 * t40 + pkin(4) * t39 + (t10 * t22 + t12 * t40) * rSges(6,1) + (-t10 * t40 + t12 * t22) * rSges(6,2) + t29) + g(2) * (t7 * t41 - pkin(4) * t38 + (-t10 * t25 + t12 * t41) * rSges(6,1) + (-t10 * t41 - t12 * t25) * rSges(6,2) + t31) + g(3) * (-t49 * t13 + t30) + (g(3) * (rSges(6,1) * t12 - rSges(6,2) * t10 + t7) + t47 * t49) * t11) - m(7) * (g(1) * (t1 * t40 + t22 * t2 + (t22 * t5 + t6 * t40) * rSges(7,1) + (t22 * t6 - t5 * t40) * rSges(7,2) + t29) + g(2) * (t1 * t41 - t25 * t2 + (-t25 * t5 + t6 * t41) * rSges(7,1) + (-t25 * t6 - t5 * t41) * rSges(7,2) + t31) + g(3) * (-t48 * t13 + t30) + (g(3) * (rSges(7,1) * t6 - rSges(7,2) * t5 + t1) + t47 * t48) * t11);
U  = t3;
