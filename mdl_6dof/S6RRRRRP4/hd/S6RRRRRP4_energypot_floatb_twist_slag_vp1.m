% Calculate potential energy for
% S6RRRRRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 01:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRRP4_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP4_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRRRP4_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP4_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP4_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP4_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP4_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:11:54
% EndTime: 2019-03-10 01:11:55
% DurationCPUTime: 0.52s
% Computational Cost: add. (260->105), mult. (227->117), div. (0->0), fcn. (215->10), ass. (0->46)
t60 = rSges(7,1) + pkin(5);
t59 = rSges(5,3) + pkin(9);
t25 = sin(qJ(1));
t28 = cos(qJ(1));
t58 = g(1) * t28 + g(2) * t25;
t57 = rSges(7,3) + qJ(6);
t54 = rSges(3,3) + pkin(7);
t22 = qJ(2) + qJ(3);
t17 = sin(t22);
t52 = rSges(4,2) * t17;
t19 = cos(t22);
t51 = t19 * t25;
t50 = t19 * t28;
t21 = qJ(4) + qJ(5);
t16 = sin(t21);
t49 = t25 * t16;
t18 = cos(t21);
t48 = t25 * t18;
t23 = sin(qJ(4));
t47 = t25 * t23;
t26 = cos(qJ(4));
t46 = t25 * t26;
t45 = t28 * t16;
t44 = t28 * t18;
t43 = t28 * t23;
t42 = t28 * t26;
t39 = pkin(6) + r_base(3);
t27 = cos(qJ(2));
t14 = t27 * pkin(2) + pkin(1);
t38 = t28 * t14 + r_base(1);
t24 = sin(qJ(2));
t37 = t24 * pkin(2) + t39;
t30 = -pkin(8) - pkin(7);
t36 = t25 * t14 + t28 * t30 + r_base(2);
t35 = -t25 * t30 + t38;
t13 = t26 * pkin(4) + pkin(3);
t29 = -pkin(10) - pkin(9);
t34 = t17 * t13 + t19 * t29 + t37;
t33 = rSges(3,1) * t27 - rSges(3,2) * t24 + pkin(1);
t32 = pkin(4) * t47 + t13 * t50 + t35;
t31 = -pkin(4) * t43 + t13 * t51 + t36;
t4 = t19 * t44 + t49;
t3 = t19 * t45 - t48;
t2 = t19 * t48 - t45;
t1 = t19 * t49 + t44;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t28 * rSges(2,1) - t25 * rSges(2,2) + r_base(1)) + g(2) * (t25 * rSges(2,1) + t28 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t39)) - m(3) * (g(1) * r_base(1) + g(2) * r_base(2) + g(3) * (t24 * rSges(3,1) + t27 * rSges(3,2) + t39) + (g(1) * t33 - g(2) * t54) * t28 + (g(1) * t54 + g(2) * t33) * t25) - m(4) * (g(1) * (rSges(4,1) * t50 - t28 * t52 + t38) + g(2) * (-t28 * rSges(4,3) + t36) + g(3) * (t17 * rSges(4,1) + t19 * rSges(4,2) + t37) + (g(1) * (rSges(4,3) - t30) + g(2) * (rSges(4,1) * t19 - t52)) * t25) - m(5) * (g(1) * (pkin(3) * t50 + (t19 * t42 + t47) * rSges(5,1) + (-t19 * t43 + t46) * rSges(5,2) + t35) + g(2) * (pkin(3) * t51 + (t19 * t46 - t43) * rSges(5,1) + (-t19 * t47 - t42) * rSges(5,2) + t36) + g(3) * (-t59 * t19 + t37) + (g(3) * (rSges(5,1) * t26 - rSges(5,2) * t23 + pkin(3)) + t58 * t59) * t17) - m(6) * (g(1) * (t4 * rSges(6,1) - t3 * rSges(6,2) + t32) + g(2) * (t2 * rSges(6,1) - t1 * rSges(6,2) + t31) + g(3) * (-t19 * rSges(6,3) + t34) + (g(3) * (rSges(6,1) * t18 - rSges(6,2) * t16) + t58 * (rSges(6,3) - t29)) * t17) - m(7) * (g(1) * (t57 * t3 + t60 * t4 + t32) + g(2) * (t57 * t1 + t60 * t2 + t31) + g(3) * (-t19 * rSges(7,2) + t34) + (g(3) * (t57 * t16 + t60 * t18) + t58 * (rSges(7,2) - t29)) * t17);
U  = t5;
