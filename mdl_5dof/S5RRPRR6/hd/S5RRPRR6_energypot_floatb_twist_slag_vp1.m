% Calculate potential energy for
% S5RRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% m [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:18
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRR6_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRPRR6_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR6_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR6_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR6_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:16:44
% EndTime: 2022-01-20 11:16:45
% DurationCPUTime: 0.43s
% Computational Cost: add. (175->84), mult. (141->96), div. (0->0), fcn. (125->10), ass. (0->31)
t42 = rSges(5,3) + pkin(7);
t41 = rSges(6,3) + pkin(8) + pkin(7);
t13 = sin(pkin(9));
t14 = cos(pkin(9));
t40 = rSges(4,1) * t14 - rSges(4,2) * t13;
t12 = qJ(1) + qJ(2);
t6 = sin(t12);
t8 = cos(t12);
t39 = g(1) * t8 + g(2) * t6;
t36 = t6 * t14;
t15 = sin(qJ(4));
t35 = t6 * t15;
t34 = t8 * t14;
t33 = t8 * t15;
t29 = t14 * t15;
t17 = cos(qJ(4));
t28 = t14 * t17;
t26 = pkin(5) + r_base(3);
t16 = sin(qJ(1));
t25 = t16 * pkin(1) + r_base(2);
t18 = cos(qJ(1));
t24 = t18 * pkin(1) + r_base(1);
t23 = pkin(6) + t26;
t22 = t6 * pkin(2) + t25;
t21 = t8 * pkin(2) + t6 * qJ(3) + t24;
t20 = -t8 * qJ(3) + t22;
t11 = qJ(4) + qJ(5);
t7 = cos(t11);
t5 = sin(t11);
t4 = pkin(4) * t17 + pkin(3);
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t18 - rSges(2,2) * t16 + r_base(1)) + g(2) * (rSges(2,1) * t16 + rSges(2,2) * t18 + r_base(2)) + g(3) * (rSges(2,3) + t26)) - m(3) * (g(1) * (rSges(3,1) * t8 - rSges(3,2) * t6 + t24) + g(2) * (rSges(3,1) * t6 + rSges(3,2) * t8 + t25) + g(3) * (rSges(3,3) + t23)) - m(4) * (g(1) * (t6 * rSges(4,3) + t21) + g(2) * (t40 * t6 + t22) + g(3) * (rSges(4,1) * t13 + rSges(4,2) * t14 + t23) + (g(1) * t40 + g(2) * (-rSges(4,3) - qJ(3))) * t8) - m(5) * (g(1) * (pkin(3) * t34 + (t28 * t8 + t35) * rSges(5,1) + (t17 * t6 - t29 * t8) * rSges(5,2) + t21) + g(2) * (pkin(3) * t36 + (t28 * t6 - t33) * rSges(5,1) + (-t17 * t8 - t29 * t6) * rSges(5,2) + t20) + g(3) * (-t42 * t14 + t23) + (g(3) * (rSges(5,1) * t17 - rSges(5,2) * t15 + pkin(3)) + t39 * t42) * t13) - m(6) * (g(1) * (t4 * t34 + pkin(4) * t35 + (t34 * t7 + t5 * t6) * rSges(6,1) + (-t34 * t5 + t6 * t7) * rSges(6,2) + t21) + g(2) * (t4 * t36 - pkin(4) * t33 + (t36 * t7 - t5 * t8) * rSges(6,1) + (-t36 * t5 - t7 * t8) * rSges(6,2) + t20) + g(3) * (-t41 * t14 + t23) + (g(3) * (rSges(6,1) * t7 - rSges(6,2) * t5 + t4) + t39 * t41) * t13);
U = t1;
