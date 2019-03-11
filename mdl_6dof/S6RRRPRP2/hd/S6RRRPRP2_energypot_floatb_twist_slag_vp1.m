% Calculate potential energy for
% S6RRRPRP2
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
% Datum: 2019-03-09 16:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPRP2_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP2_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRPRP2_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP2_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP2_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP2_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP2_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:35:01
% EndTime: 2019-03-09 16:35:02
% DurationCPUTime: 0.35s
% Computational Cost: add. (260->98), mult. (199->107), div. (0->0), fcn. (183->10), ass. (0->43)
t54 = rSges(7,1) + pkin(5);
t53 = rSges(7,3) + qJ(6);
t32 = -pkin(8) - pkin(7);
t52 = rSges(3,3) + pkin(7);
t30 = cos(qJ(2));
t18 = t30 * pkin(2) + pkin(1);
t25 = qJ(2) + qJ(3);
t19 = pkin(10) + t25;
t14 = sin(t19);
t28 = sin(qJ(1));
t51 = t14 * t28;
t31 = cos(qJ(1));
t50 = t14 * t31;
t15 = cos(t19);
t49 = t15 * t31;
t26 = sin(qJ(5));
t48 = t28 * t26;
t29 = cos(qJ(5));
t47 = t28 * t29;
t46 = t31 * t26;
t45 = t31 * t29;
t44 = rSges(4,3) - t32;
t43 = pkin(6) + r_base(3);
t21 = cos(t25);
t7 = pkin(3) * t21 + t18;
t42 = t31 * t7 + r_base(1);
t27 = sin(qJ(2));
t41 = t27 * pkin(2) + t43;
t24 = -qJ(4) + t32;
t40 = t31 * t24 + t28 * t7 + r_base(2);
t20 = sin(t25);
t39 = pkin(3) * t20 + t41;
t38 = t28 * t15 * pkin(4) + pkin(9) * t51 + t40;
t37 = t14 * pkin(4) + t39;
t36 = rSges(3,1) * t30 - rSges(3,2) * t27 + pkin(1);
t35 = rSges(4,1) * t21 - rSges(4,2) * t20 + t18;
t34 = g(1) * r_base(1) + g(2) * r_base(2);
t33 = pkin(4) * t49 + pkin(9) * t50 - t28 * t24 + t42;
t4 = t15 * t45 + t48;
t3 = t15 * t46 - t47;
t2 = t15 * t47 - t46;
t1 = t15 * t48 + t45;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t31 * rSges(2,1) - t28 * rSges(2,2) + r_base(1)) + g(2) * (t28 * rSges(2,1) + t31 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t43)) - m(3) * (g(3) * (t27 * rSges(3,1) + t30 * rSges(3,2) + t43) + (g(1) * t36 - g(2) * t52) * t31 + (g(1) * t52 + g(2) * t36) * t28 + t34) - m(4) * (g(3) * (t20 * rSges(4,1) + t21 * rSges(4,2) + t41) + (g(1) * t35 - g(2) * t44) * t31 + (g(1) * t44 + g(2) * t35) * t28 + t34) - m(5) * (g(1) * (rSges(5,1) * t49 - rSges(5,2) * t50 + t42) + g(2) * (-t31 * rSges(5,3) + t40) + g(3) * (t14 * rSges(5,1) + t15 * rSges(5,2) + t39) + (g(1) * (rSges(5,3) - t24) + g(2) * (rSges(5,1) * t15 - rSges(5,2) * t14)) * t28) - m(6) * (g(1) * (t4 * rSges(6,1) - t3 * rSges(6,2) + rSges(6,3) * t50 + t33) + g(2) * (t2 * rSges(6,1) - t1 * rSges(6,2) + rSges(6,3) * t51 + t38) + g(3) * ((-rSges(6,3) - pkin(9)) * t15 + (rSges(6,1) * t29 - rSges(6,2) * t26) * t14 + t37)) - m(7) * (g(1) * (t53 * t3 + t54 * t4 + t33) + g(2) * (t53 * t1 + t54 * t2 + t38) + g(3) * (t37 + (-rSges(7,2) - pkin(9)) * t15) + (g(3) * (t53 * t26 + t54 * t29) + (g(1) * t31 + g(2) * t28) * rSges(7,2)) * t14);
U  = t5;
