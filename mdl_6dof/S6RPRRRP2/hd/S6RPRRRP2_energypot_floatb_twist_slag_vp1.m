% Calculate potential energy for
% S6RPRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 06:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRP2_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP2_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRRRP2_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP2_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP2_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP2_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP2_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:58:41
% EndTime: 2019-03-09 05:58:42
% DurationCPUTime: 0.47s
% Computational Cost: add. (256->105), mult. (198->120), div. (0->0), fcn. (182->10), ass. (0->41)
t52 = rSges(5,3) + pkin(8);
t27 = -pkin(9) - pkin(8);
t51 = rSges(6,3) - t27;
t50 = rSges(7,3) + qJ(6) - t27;
t19 = qJ(1) + pkin(10);
t11 = sin(t19);
t12 = cos(t19);
t49 = g(1) * t12 + g(2) * t11;
t24 = cos(qJ(4));
t10 = t24 * pkin(4) + pkin(3);
t22 = sin(qJ(3));
t45 = rSges(4,2) * t22;
t21 = sin(qJ(4));
t44 = t11 * t21;
t25 = cos(qJ(3));
t43 = t11 * t25;
t42 = t12 * t21;
t41 = t12 * t25;
t20 = qJ(4) + qJ(5);
t13 = sin(t20);
t40 = t13 * t25;
t14 = cos(t20);
t39 = t14 * t25;
t38 = t21 * t25;
t37 = t24 * t25;
t34 = pkin(6) + r_base(3);
t23 = sin(qJ(1));
t33 = t23 * pkin(1) + r_base(2);
t26 = cos(qJ(1));
t32 = t26 * pkin(1) + r_base(1);
t31 = t11 * pkin(2) + t33;
t30 = qJ(2) + t34;
t29 = t12 * pkin(2) + t11 * pkin(7) + t32;
t28 = -t12 * pkin(7) + t31;
t6 = t21 * pkin(4) + pkin(5) * t13;
t5 = pkin(5) * t14 + t10;
t4 = t11 * t13 + t12 * t39;
t3 = t11 * t14 - t12 * t40;
t2 = t11 * t39 - t12 * t13;
t1 = -t11 * t40 - t12 * t14;
t7 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t26 * rSges(2,1) - t23 * rSges(2,2) + r_base(1)) + g(2) * (t23 * rSges(2,1) + t26 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t34)) - m(3) * (g(1) * (t12 * rSges(3,1) - t11 * rSges(3,2) + t32) + g(2) * (t11 * rSges(3,1) + t12 * rSges(3,2) + t33) + g(3) * (rSges(3,3) + t30)) - m(4) * (g(1) * (t11 * rSges(4,3) + t29) + g(2) * (rSges(4,1) * t43 - t11 * t45 + t31) + g(3) * (t22 * rSges(4,1) + t25 * rSges(4,2) + t30) + (g(1) * (rSges(4,1) * t25 - t45) + g(2) * (-rSges(4,3) - pkin(7))) * t12) - m(5) * (g(1) * (pkin(3) * t41 + (t12 * t37 + t44) * rSges(5,1) + (t11 * t24 - t12 * t38) * rSges(5,2) + t29) + g(2) * (pkin(3) * t43 + (t11 * t37 - t42) * rSges(5,1) + (-t11 * t38 - t12 * t24) * rSges(5,2) + t28) + g(3) * (-t52 * t25 + t30) + (g(3) * (rSges(5,1) * t24 - rSges(5,2) * t21 + pkin(3)) + t49 * t52) * t22) - m(6) * (g(1) * (t4 * rSges(6,1) + t3 * rSges(6,2) + pkin(4) * t44 + t10 * t41 + t29) + g(2) * (t2 * rSges(6,1) + t1 * rSges(6,2) - pkin(4) * t42 + t10 * t43 + t28) + g(3) * (-t51 * t25 + t30) + (g(3) * (rSges(6,1) * t14 - rSges(6,2) * t13 + t10) + t49 * t51) * t22) - m(7) * (g(1) * (t4 * rSges(7,1) + t3 * rSges(7,2) + t11 * t6 + t5 * t41 + t29) + g(2) * (t2 * rSges(7,1) + t1 * rSges(7,2) - t12 * t6 + t5 * t43 + t28) + g(3) * (-t50 * t25 + t30) + (g(3) * (rSges(7,1) * t14 - rSges(7,2) * t13 + t5) + t49 * t50) * t22);
U  = t7;
