% Calculate potential energy for
% S6RRRPRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 18:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPRR4_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR4_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRPRR4_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR4_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR4_energypot_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR4_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR4_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:15:17
% EndTime: 2019-03-09 18:15:18
% DurationCPUTime: 0.53s
% Computational Cost: add. (257->114), mult. (212->130), div. (0->0), fcn. (196->12), ass. (0->40)
t22 = -pkin(9) - qJ(4);
t50 = rSges(6,3) - t22;
t49 = rSges(7,3) + pkin(10) - t22;
t24 = sin(qJ(1));
t26 = cos(qJ(1));
t48 = g(1) * t26 + g(2) * t24;
t47 = rSges(5,3) + qJ(4);
t44 = rSges(3,3) + pkin(7);
t21 = cos(pkin(11));
t7 = t21 * pkin(4) + pkin(3);
t19 = qJ(2) + qJ(3);
t13 = sin(t19);
t43 = rSges(4,2) * t13;
t14 = cos(t19);
t42 = t14 * t24;
t41 = t14 * t26;
t20 = sin(pkin(11));
t40 = t20 * t24;
t39 = t20 * t26;
t38 = t21 * t24;
t37 = t21 * t26;
t33 = pkin(6) + r_base(3);
t25 = cos(qJ(2));
t8 = pkin(2) * t25 + pkin(1);
t32 = t26 * t8 + r_base(1);
t18 = pkin(11) + qJ(5);
t27 = -pkin(8) - pkin(7);
t31 = t24 * t8 + t26 * t27 + r_base(2);
t23 = sin(qJ(2));
t30 = t23 * pkin(2) + t33;
t29 = -t24 * t27 + t32;
t28 = rSges(3,1) * t25 - rSges(3,2) * t23 + pkin(1);
t12 = qJ(6) + t18;
t11 = cos(t18);
t10 = sin(t18);
t6 = cos(t12);
t5 = sin(t12);
t2 = pkin(4) * t20 + pkin(5) * t10;
t1 = pkin(5) * t11 + t7;
t3 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t26 - rSges(2,2) * t24 + r_base(1)) + g(2) * (rSges(2,1) * t24 + rSges(2,2) * t26 + r_base(2)) + g(3) * (rSges(2,3) + t33)) - m(3) * (g(1) * r_base(1) + g(2) * r_base(2) + g(3) * (rSges(3,1) * t23 + rSges(3,2) * t25 + t33) + (g(1) * t28 - g(2) * t44) * t26 + (g(1) * t44 + g(2) * t28) * t24) - m(4) * (g(1) * (rSges(4,1) * t41 - t26 * t43 + t32) + g(2) * (-t26 * rSges(4,3) + t31) + g(3) * (rSges(4,1) * t13 + rSges(4,2) * t14 + t30) + (g(1) * (rSges(4,3) - t27) + g(2) * (rSges(4,1) * t14 - t43)) * t24) - m(5) * (g(1) * (pkin(3) * t41 + (t14 * t37 + t40) * rSges(5,1) + (-t14 * t39 + t38) * rSges(5,2) + t29) + g(2) * (pkin(3) * t42 + (t14 * t38 - t39) * rSges(5,1) + (-t14 * t40 - t37) * rSges(5,2) + t31) + g(3) * (-t47 * t14 + t30) + (g(3) * (rSges(5,1) * t21 - rSges(5,2) * t20 + pkin(3)) + t48 * t47) * t13) - m(6) * (g(1) * (t7 * t41 + pkin(4) * t40 + (t10 * t24 + t11 * t41) * rSges(6,1) + (-t10 * t41 + t11 * t24) * rSges(6,2) + t29) + g(2) * (t7 * t42 - pkin(4) * t39 + (-t10 * t26 + t11 * t42) * rSges(6,1) + (-t10 * t42 - t11 * t26) * rSges(6,2) + t31) + g(3) * (-t50 * t14 + t30) + (g(3) * (rSges(6,1) * t11 - rSges(6,2) * t10 + t7) + t48 * t50) * t13) - m(7) * (g(1) * (t1 * t41 + t24 * t2 + (t24 * t5 + t6 * t41) * rSges(7,1) + (t24 * t6 - t5 * t41) * rSges(7,2) + t29) + g(2) * (t1 * t42 - t26 * t2 + (-t26 * t5 + t6 * t42) * rSges(7,1) + (-t26 * t6 - t5 * t42) * rSges(7,2) + t31) + g(3) * (-t49 * t14 + t30) + (g(3) * (rSges(7,1) * t6 - rSges(7,2) * t5 + t1) + t48 * t49) * t13);
U  = t3;
