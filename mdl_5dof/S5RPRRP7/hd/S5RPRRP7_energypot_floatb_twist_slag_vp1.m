% Calculate potential energy for
% S5RPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRP7_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP7_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPRRP7_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP7_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP7_energypot_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP7_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP7_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:44:29
% EndTime: 2019-12-31 18:44:29
% DurationCPUTime: 0.33s
% Computational Cost: add. (173->79), mult. (154->89), div. (0->0), fcn. (142->8), ass. (0->30)
t40 = rSges(6,1) + pkin(4);
t39 = rSges(6,3) + qJ(5);
t18 = qJ(1) + pkin(8);
t13 = sin(t18);
t20 = sin(qJ(3));
t38 = t13 * t20;
t23 = cos(qJ(3));
t37 = t13 * t23;
t14 = cos(t18);
t36 = t14 * t20;
t19 = sin(qJ(4));
t35 = t19 * t23;
t22 = cos(qJ(4));
t34 = t22 * t23;
t33 = pkin(5) + r_base(3);
t21 = sin(qJ(1));
t32 = t21 * pkin(1) + r_base(2);
t24 = cos(qJ(1));
t31 = t24 * pkin(1) + r_base(1);
t30 = qJ(2) + t33;
t29 = t13 * pkin(2) + t32;
t28 = t20 * pkin(3) + t30;
t27 = t14 * pkin(2) + t13 * pkin(6) + t31;
t26 = t14 * t23 * pkin(3) + pkin(7) * t36 + t27;
t25 = pkin(3) * t37 - t14 * pkin(6) + pkin(7) * t38 + t29;
t4 = t13 * t19 + t14 * t34;
t3 = -t13 * t22 + t14 * t35;
t2 = t13 * t34 - t14 * t19;
t1 = t13 * t35 + t14 * t22;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t24 * rSges(2,1) - t21 * rSges(2,2) + r_base(1)) + g(2) * (t21 * rSges(2,1) + t24 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t33)) - m(3) * (g(1) * (t14 * rSges(3,1) - t13 * rSges(3,2) + t31) + g(2) * (t13 * rSges(3,1) + t14 * rSges(3,2) + t32) + g(3) * (rSges(3,3) + t30)) - m(4) * (g(1) * (t13 * rSges(4,3) + t27) + g(2) * (rSges(4,1) * t37 - rSges(4,2) * t38 + t29) + g(3) * (t20 * rSges(4,1) + t23 * rSges(4,2) + t30) + (g(1) * (rSges(4,1) * t23 - rSges(4,2) * t20) + g(2) * (-rSges(4,3) - pkin(6))) * t14) - m(5) * (g(1) * (t4 * rSges(5,1) - t3 * rSges(5,2) + rSges(5,3) * t36 + t26) + g(2) * (t2 * rSges(5,1) - t1 * rSges(5,2) + rSges(5,3) * t38 + t25) + g(3) * ((-rSges(5,3) - pkin(7)) * t23 + (rSges(5,1) * t22 - rSges(5,2) * t19) * t20 + t28)) - m(6) * (g(1) * (t39 * t3 + t40 * t4 + t26) + g(2) * (t39 * t1 + t40 * t2 + t25) + g(3) * (t28 + (-rSges(6,2) - pkin(7)) * t23) + (g(3) * (t39 * t19 + t40 * t22) + (g(1) * t14 + g(2) * t13) * rSges(6,2)) * t20);
U = t5;
