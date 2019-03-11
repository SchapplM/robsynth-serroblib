% Calculate potential energy for
% S6RPRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta5]';
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
% Datum: 2019-03-09 04:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRPP6_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP6_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRRPP6_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP6_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP6_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP6_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPP6_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:46:16
% EndTime: 2019-03-09 04:46:16
% DurationCPUTime: 0.48s
% Computational Cost: add. (192->106), mult. (215->115), div. (0->0), fcn. (203->8), ass. (0->40)
t50 = rSges(7,1) + pkin(5);
t49 = rSges(5,3) + pkin(8);
t20 = -qJ(5) - pkin(8);
t48 = rSges(7,2) - t20;
t47 = rSges(7,3) + qJ(6);
t23 = sin(qJ(1));
t46 = g(1) * t23;
t26 = cos(qJ(1));
t45 = g(2) * t26;
t21 = sin(qJ(4));
t43 = t23 * t21;
t22 = sin(qJ(3));
t42 = t23 * t22;
t24 = cos(qJ(4));
t41 = t23 * t24;
t25 = cos(qJ(3));
t40 = t23 * t25;
t39 = t26 * t21;
t38 = t26 * t22;
t37 = t26 * t24;
t36 = rSges(6,3) - t20;
t35 = pkin(6) + r_base(3);
t34 = t23 * pkin(1) + r_base(2);
t33 = pkin(2) + t35;
t32 = t26 * pkin(1) + t23 * qJ(2) + r_base(1);
t31 = t23 * pkin(7) + t34;
t11 = t24 * pkin(4) + pkin(3);
t30 = t25 * t11 + t33;
t29 = t26 * pkin(7) + t32;
t28 = -t26 * qJ(2) + t31;
t27 = pkin(4) * t39 + t11 * t42 + t20 * t40 + t29;
t19 = qJ(4) + pkin(9);
t13 = cos(t19);
t12 = sin(t19);
t9 = pkin(4) * t43;
t4 = t23 * t12 - t13 * t38;
t3 = t12 * t38 + t23 * t13;
t2 = t26 * t12 + t13 * t42;
t1 = t12 * t42 - t26 * t13;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t26 * rSges(2,1) - t23 * rSges(2,2) + r_base(1)) + g(2) * (t23 * rSges(2,1) + t26 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t35)) - m(3) * (g(1) * (-t26 * rSges(3,2) + t23 * rSges(3,3) + t32) + g(2) * (-t23 * rSges(3,2) + (-rSges(3,3) - qJ(2)) * t26 + t34) + g(3) * (rSges(3,1) + t35)) - m(4) * (g(1) * (rSges(4,1) * t42 + rSges(4,2) * t40 + t29) + g(2) * (t23 * rSges(4,3) + t31) + g(3) * (t25 * rSges(4,1) - t22 * rSges(4,2) + t33) + (g(1) * rSges(4,3) + g(2) * (-rSges(4,1) * t22 - rSges(4,2) * t25 - qJ(2))) * t26) - m(5) * (g(1) * (pkin(3) * t42 + (t22 * t41 + t39) * rSges(5,1) + (-t21 * t42 + t37) * rSges(5,2) + t29) + g(2) * (-pkin(3) * t38 + (-t22 * t37 + t43) * rSges(5,1) + (t21 * t38 + t41) * rSges(5,2) + t28) + g(3) * (t49 * t22 + t33) + (g(3) * (rSges(5,1) * t24 - rSges(5,2) * t21 + pkin(3)) + (-t46 + t45) * t49) * t25) - m(6) * (g(1) * (t2 * rSges(6,1) - t1 * rSges(6,2) - rSges(6,3) * t40 + t27) + g(2) * (t4 * rSges(6,1) + t3 * rSges(6,2) + t31 + t9) + g(3) * ((rSges(6,1) * t13 - rSges(6,2) * t12) * t25 + t36 * t22 + t30) + (-t11 * t22 + t36 * t25 - qJ(2)) * t45) - m(7) * (g(1) * (t47 * t1 + t50 * t2 + t27) + g(2) * (-t11 * t38 - t47 * t3 + t50 * t4 + t28 + t9) + g(3) * (t22 * t48 + t30) + (-rSges(7,2) * t46 + g(3) * (t47 * t12 + t50 * t13) + t48 * t45) * t25);
U  = t5;
