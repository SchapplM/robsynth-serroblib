% Calculate potential energy for
% S5RPRRR7
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% m_mdh [6x1]
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
% Datum: 2019-12-31 19:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRR7_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR7_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPRRR7_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR7_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR7_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR7_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR7_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:03:07
% EndTime: 2019-12-31 19:03:08
% DurationCPUTime: 0.42s
% Computational Cost: add. (175->84), mult. (141->96), div. (0->0), fcn. (125->10), ass. (0->31)
t42 = rSges(5,3) + pkin(7);
t41 = rSges(6,3) + pkin(8) + pkin(7);
t14 = sin(qJ(3));
t17 = cos(qJ(3));
t40 = rSges(4,1) * t17 - rSges(4,2) * t14;
t11 = qJ(1) + pkin(9);
t5 = sin(t11);
t6 = cos(t11);
t39 = g(1) * t6 + g(2) * t5;
t13 = sin(qJ(4));
t36 = t5 * t13;
t35 = t5 * t17;
t34 = t6 * t13;
t33 = t6 * t17;
t29 = t13 * t17;
t16 = cos(qJ(4));
t28 = t16 * t17;
t26 = pkin(5) + r_base(3);
t15 = sin(qJ(1));
t25 = t15 * pkin(1) + r_base(2);
t18 = cos(qJ(1));
t24 = t18 * pkin(1) + r_base(1);
t23 = t5 * pkin(2) + t25;
t22 = qJ(2) + t26;
t21 = t6 * pkin(2) + t5 * pkin(6) + t24;
t20 = -t6 * pkin(6) + t23;
t12 = qJ(4) + qJ(5);
t8 = cos(t12);
t7 = sin(t12);
t4 = t16 * pkin(4) + pkin(3);
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t18 * rSges(2,1) - t15 * rSges(2,2) + r_base(1)) + g(2) * (t15 * rSges(2,1) + t18 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t26)) - m(3) * (g(1) * (rSges(3,1) * t6 - rSges(3,2) * t5 + t24) + g(2) * (rSges(3,1) * t5 + rSges(3,2) * t6 + t25) + g(3) * (rSges(3,3) + t22)) - m(4) * (g(1) * (t5 * rSges(4,3) + t21) + g(2) * (t40 * t5 + t23) + g(3) * (rSges(4,1) * t14 + rSges(4,2) * t17 + t22) + (g(1) * t40 + g(2) * (-rSges(4,3) - pkin(6))) * t6) - m(5) * (g(1) * (pkin(3) * t33 + (t28 * t6 + t36) * rSges(5,1) + (t16 * t5 - t29 * t6) * rSges(5,2) + t21) + g(2) * (pkin(3) * t35 + (t28 * t5 - t34) * rSges(5,1) + (-t16 * t6 - t29 * t5) * rSges(5,2) + t20) + g(3) * (-t42 * t17 + t22) + (g(3) * (rSges(5,1) * t16 - rSges(5,2) * t13 + pkin(3)) + t39 * t42) * t14) - m(6) * (g(1) * (t4 * t33 + pkin(4) * t36 + (t33 * t8 + t5 * t7) * rSges(6,1) + (-t33 * t7 + t5 * t8) * rSges(6,2) + t21) + g(2) * (t4 * t35 - pkin(4) * t34 + (t35 * t8 - t6 * t7) * rSges(6,1) + (-t7 * t35 - t6 * t8) * rSges(6,2) + t20) + g(3) * (-t41 * t17 + t22) + (g(3) * (rSges(6,1) * t8 - rSges(6,2) * t7 + t4) + t39 * t41) * t14);
U = t1;
