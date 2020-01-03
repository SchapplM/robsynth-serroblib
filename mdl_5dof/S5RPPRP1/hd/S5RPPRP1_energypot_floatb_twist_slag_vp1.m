% Calculate potential energy for
% S5RPPRP1
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
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
% Datum: 2020-01-03 11:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPRP1_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPPRP1_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP1_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_energypot_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP1_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP1_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:25:03
% EndTime: 2020-01-03 11:25:04
% DurationCPUTime: 0.41s
% Computational Cost: add. (165->82), mult. (141->91), div. (0->0), fcn. (125->8), ass. (0->30)
t36 = rSges(5,3) + pkin(6);
t35 = rSges(6,3) + qJ(5) + pkin(6);
t10 = qJ(1) + pkin(7);
t7 = sin(t10);
t8 = cos(t10);
t34 = g(2) * t7 - g(3) * t8;
t12 = cos(pkin(8));
t31 = t12 * t7;
t14 = sin(qJ(4));
t30 = t7 * t14;
t28 = t12 * t14;
t16 = cos(qJ(4));
t27 = t12 * t16;
t25 = -rSges(4,3) - qJ(3);
t24 = pkin(5) + r_base(1);
t15 = sin(qJ(1));
t23 = t15 * pkin(1) + r_base(2);
t22 = t7 * pkin(2) + t23;
t21 = qJ(2) + t24;
t17 = cos(qJ(1));
t20 = -t17 * pkin(1) + r_base(3);
t11 = sin(pkin(8));
t19 = rSges(4,1) * t12 - rSges(4,2) * t11;
t18 = -t7 * qJ(3) + t20;
t6 = t16 * pkin(4) + pkin(3);
t4 = -t27 * t8 - t30;
t3 = -t7 * t16 + t28 * t8;
t2 = -t8 * t14 + t27 * t7;
t1 = -t8 * t16 - t28 * t7;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,3) + t24) + g(2) * (t15 * rSges(2,1) + t17 * rSges(2,2) + r_base(2)) + g(3) * (-t17 * rSges(2,1) + t15 * rSges(2,2) + r_base(3))) - m(3) * (g(1) * (rSges(3,3) + t21) + g(2) * (t7 * rSges(3,1) + t8 * rSges(3,2) + t23) + g(3) * (-t8 * rSges(3,1) + t7 * rSges(3,2) + t20)) - m(4) * (g(1) * (t11 * rSges(4,1) + t12 * rSges(4,2) + t21) + g(2) * t22 + g(3) * t20 + (g(2) * t19 + g(3) * t25) * t7 + (g(2) * t25 + g(3) * (-pkin(2) - t19)) * t8) - m(5) * (g(1) * (-t36 * t12 + t21) + g(2) * (t2 * rSges(5,1) + t1 * rSges(5,2) + pkin(3) * t31 - t8 * qJ(3) + t22) + g(3) * (t4 * rSges(5,1) + t3 * rSges(5,2) + t18 + (-pkin(3) * t12 - pkin(2)) * t8) + (g(1) * (rSges(5,1) * t16 - rSges(5,2) * t14 + pkin(3)) + t34 * t36) * t11) - m(6) * (g(1) * (-t35 * t12 + t21) + g(2) * (t2 * rSges(6,1) + t1 * rSges(6,2) + t6 * t31 + t22) + g(3) * (t4 * rSges(6,1) + t3 * rSges(6,2) - pkin(4) * t30 + t18) + (g(2) * (-pkin(4) * t14 - qJ(3)) + g(3) * (-t12 * t6 - pkin(2))) * t8 + (g(1) * (rSges(6,1) * t16 - rSges(6,2) * t14 + t6) + t34 * t35) * t11);
U = t5;
