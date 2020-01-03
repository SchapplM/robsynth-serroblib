% Calculate potential energy for
% S5RPRRR13
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
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
% Datum: 2019-12-31 19:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRR13_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR13_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPRRR13_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR13_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR13_energypot_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR13_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR13_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:14:27
% EndTime: 2019-12-31 19:14:28
% DurationCPUTime: 0.39s
% Computational Cost: add. (124->85), mult. (143->95), div. (0->0), fcn. (127->8), ass. (0->28)
t37 = rSges(5,3) + pkin(7);
t36 = rSges(6,3) + pkin(8) + pkin(7);
t12 = sin(qJ(1));
t15 = cos(qJ(1));
t35 = -g(1) * t12 + g(2) * t15;
t14 = cos(qJ(3));
t31 = rSges(4,2) * t14;
t10 = sin(qJ(4));
t30 = t12 * t10;
t11 = sin(qJ(3));
t29 = t12 * t11;
t13 = cos(qJ(4));
t28 = t12 * t13;
t27 = t15 * t10;
t26 = t15 * t11;
t25 = t15 * t13;
t23 = pkin(5) + r_base(3);
t22 = t12 * pkin(1) + r_base(2);
t21 = pkin(2) + t23;
t20 = t15 * pkin(1) + t12 * qJ(2) + r_base(1);
t19 = t12 * pkin(6) + t22;
t18 = t15 * pkin(6) + t20;
t17 = -t15 * qJ(2) + t19;
t9 = qJ(4) + qJ(5);
t3 = cos(t9);
t2 = sin(t9);
t1 = t13 * pkin(4) + pkin(3);
t4 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t15 * rSges(2,1) - t12 * rSges(2,2) + r_base(1)) + g(2) * (t12 * rSges(2,1) + t15 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t23)) - m(3) * (g(1) * (-t15 * rSges(3,2) + t12 * rSges(3,3) + t20) + g(2) * (-t12 * rSges(3,2) + (-rSges(3,3) - qJ(2)) * t15 + t22) + g(3) * (rSges(3,1) + t23)) - m(4) * (g(1) * (rSges(4,1) * t29 + t12 * t31 + t18) + g(2) * (t12 * rSges(4,3) + t19) + g(3) * (t14 * rSges(4,1) - t11 * rSges(4,2) + t21) + (g(1) * rSges(4,3) + g(2) * (-rSges(4,1) * t11 - qJ(2) - t31)) * t15) - m(5) * (g(1) * (pkin(3) * t29 + (t11 * t28 + t27) * rSges(5,1) + (-t10 * t29 + t25) * rSges(5,2) + t18) + g(2) * (-pkin(3) * t26 + (-t11 * t25 + t30) * rSges(5,1) + (t10 * t26 + t28) * rSges(5,2) + t17) + g(3) * (t37 * t11 + t21) + (g(3) * (rSges(5,1) * t13 - rSges(5,2) * t10 + pkin(3)) + t35 * t37) * t14) - m(6) * (g(1) * (t1 * t29 + pkin(4) * t27 + (t15 * t2 + t3 * t29) * rSges(6,1) + (t15 * t3 - t2 * t29) * rSges(6,2) + t18) + g(2) * (-t1 * t26 + pkin(4) * t30 + (t12 * t2 - t3 * t26) * rSges(6,1) + (t12 * t3 + t2 * t26) * rSges(6,2) + t17) + g(3) * (t36 * t11 + t21) + (g(3) * (rSges(6,1) * t3 - rSges(6,2) * t2 + t1) + t35 * t36) * t14);
U = t4;
