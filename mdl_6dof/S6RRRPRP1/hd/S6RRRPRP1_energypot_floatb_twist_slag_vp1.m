% Calculate potential energy for
% S6RRRPRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-03-09 16:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPRP1_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP1_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRPRP1_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP1_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP1_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP1_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP1_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:30:38
% EndTime: 2019-03-09 16:30:38
% DurationCPUTime: 0.40s
% Computational Cost: add. (246->99), mult. (186->107), div. (0->0), fcn. (166->10), ass. (0->43)
t53 = rSges(6,3) + pkin(9);
t52 = rSges(7,3) + qJ(6) + pkin(9);
t24 = sin(qJ(1));
t27 = cos(qJ(1));
t51 = g(1) * t27 + g(2) * t24;
t28 = -pkin(8) - pkin(7);
t20 = qJ(2) + qJ(3);
t14 = pkin(10) + t20;
t9 = sin(t14);
t50 = rSges(5,2) * t9;
t47 = rSges(3,3) + pkin(7);
t26 = cos(qJ(2));
t13 = t26 * pkin(2) + pkin(1);
t10 = cos(t14);
t45 = t10 * t24;
t44 = t10 * t27;
t22 = sin(qJ(5));
t43 = t24 * t22;
t25 = cos(qJ(5));
t42 = t24 * t25;
t41 = t27 * t22;
t40 = t27 * t25;
t39 = rSges(4,3) - t28;
t37 = pkin(6) + r_base(3);
t16 = cos(t20);
t7 = pkin(3) * t16 + t13;
t36 = t27 * t7 + r_base(1);
t19 = -qJ(4) + t28;
t35 = t27 * t19 + t24 * t7 + r_base(2);
t23 = sin(qJ(2));
t34 = t23 * pkin(2) + t37;
t15 = sin(t20);
t33 = pkin(3) * t15 + t34;
t32 = -t24 * t19 + t36;
t31 = rSges(3,1) * t26 - rSges(3,2) * t23 + pkin(1);
t30 = rSges(4,1) * t16 - rSges(4,2) * t15 + t13;
t29 = g(1) * r_base(1) + g(2) * r_base(2);
t12 = t25 * pkin(5) + pkin(4);
t4 = t10 * t40 + t43;
t3 = -t10 * t41 + t42;
t2 = t10 * t42 - t41;
t1 = -t10 * t43 - t40;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t27 * rSges(2,1) - t24 * rSges(2,2) + r_base(1)) + g(2) * (t24 * rSges(2,1) + t27 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t37)) - m(3) * (g(3) * (t23 * rSges(3,1) + t26 * rSges(3,2) + t37) + (g(1) * t31 - g(2) * t47) * t27 + (g(1) * t47 + g(2) * t31) * t24 + t29) - m(4) * (g(3) * (t15 * rSges(4,1) + t16 * rSges(4,2) + t34) + (g(1) * t30 - g(2) * t39) * t27 + (g(1) * t39 + g(2) * t30) * t24 + t29) - m(5) * (g(1) * (rSges(5,1) * t44 - t27 * t50 + t36) + g(2) * (-t27 * rSges(5,3) + t35) + g(3) * (t9 * rSges(5,1) + t10 * rSges(5,2) + t33) + (g(1) * (rSges(5,3) - t19) + g(2) * (rSges(5,1) * t10 - t50)) * t24) - m(6) * (g(1) * (t4 * rSges(6,1) + t3 * rSges(6,2) + pkin(4) * t44 + t32) + g(2) * (t2 * rSges(6,1) + t1 * rSges(6,2) + pkin(4) * t45 + t35) + g(3) * (-t53 * t10 + t33) + (g(3) * (rSges(6,1) * t25 - rSges(6,2) * t22 + pkin(4)) + t51 * t53) * t9) - m(7) * (g(1) * (t4 * rSges(7,1) + t3 * rSges(7,2) + pkin(5) * t43 + t12 * t44 + t32) + g(2) * (t2 * rSges(7,1) + t1 * rSges(7,2) - pkin(5) * t41 + t12 * t45 + t35) + g(3) * (-t52 * t10 + t33) + (g(3) * (rSges(7,1) * t25 - rSges(7,2) * t22 + t12) + t51 * t52) * t9);
U  = t5;
