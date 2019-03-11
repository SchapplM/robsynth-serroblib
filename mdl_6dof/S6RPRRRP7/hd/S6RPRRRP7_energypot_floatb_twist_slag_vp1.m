% Calculate potential energy for
% S6RPRRRP7
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
% Datum: 2019-03-09 06:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRP7_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP7_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRRRP7_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP7_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP7_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP7_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP7_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:18:32
% EndTime: 2019-03-09 06:18:33
% DurationCPUTime: 0.53s
% Computational Cost: add. (260->105), mult. (227->117), div. (0->0), fcn. (215->10), ass. (0->46)
t60 = rSges(7,1) + pkin(5);
t59 = rSges(5,3) + pkin(8);
t27 = sin(qJ(1));
t29 = cos(qJ(1));
t58 = g(1) * t29 + g(2) * t27;
t57 = rSges(7,3) + qJ(6);
t21 = pkin(10) + qJ(3);
t16 = sin(t21);
t53 = rSges(4,2) * t16;
t17 = cos(t21);
t52 = t17 * t27;
t51 = t17 * t29;
t22 = qJ(4) + qJ(5);
t18 = sin(t22);
t50 = t27 * t18;
t19 = cos(t22);
t49 = t27 * t19;
t26 = sin(qJ(4));
t48 = t27 * t26;
t28 = cos(qJ(4));
t47 = t27 * t28;
t46 = t29 * t18;
t45 = t29 * t19;
t44 = t29 * t26;
t43 = t29 * t28;
t40 = rSges(3,3) + qJ(2);
t39 = pkin(6) + r_base(3);
t24 = cos(pkin(10));
t13 = t24 * pkin(2) + pkin(1);
t38 = t29 * t13 + r_base(1);
t23 = sin(pkin(10));
t37 = t23 * pkin(2) + t39;
t25 = -pkin(7) - qJ(2);
t36 = t27 * t13 + t29 * t25 + r_base(2);
t35 = -t27 * t25 + t38;
t15 = t28 * pkin(4) + pkin(3);
t30 = -pkin(9) - pkin(8);
t34 = t16 * t15 + t17 * t30 + t37;
t33 = rSges(3,1) * t24 - rSges(3,2) * t23 + pkin(1);
t32 = pkin(4) * t48 + t15 * t51 + t35;
t31 = -pkin(4) * t44 + t15 * t52 + t36;
t4 = t17 * t45 + t50;
t3 = t17 * t46 - t49;
t2 = t17 * t49 - t46;
t1 = t17 * t50 + t45;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t29 * rSges(2,1) - t27 * rSges(2,2) + r_base(1)) + g(2) * (t27 * rSges(2,1) + t29 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t39)) - m(3) * (g(1) * r_base(1) + g(2) * r_base(2) + g(3) * (t23 * rSges(3,1) + t24 * rSges(3,2) + t39) + (g(1) * t33 - g(2) * t40) * t29 + (g(1) * t40 + g(2) * t33) * t27) - m(4) * (g(1) * (rSges(4,1) * t51 - t29 * t53 + t38) + g(2) * (-t29 * rSges(4,3) + t36) + g(3) * (t16 * rSges(4,1) + t17 * rSges(4,2) + t37) + (g(1) * (rSges(4,3) - t25) + g(2) * (rSges(4,1) * t17 - t53)) * t27) - m(5) * (g(1) * (pkin(3) * t51 + (t17 * t43 + t48) * rSges(5,1) + (-t17 * t44 + t47) * rSges(5,2) + t35) + g(2) * (pkin(3) * t52 + (t17 * t47 - t44) * rSges(5,1) + (-t17 * t48 - t43) * rSges(5,2) + t36) + g(3) * (-t59 * t17 + t37) + (g(3) * (rSges(5,1) * t28 - rSges(5,2) * t26 + pkin(3)) + t58 * t59) * t16) - m(6) * (g(1) * (t4 * rSges(6,1) - t3 * rSges(6,2) + t32) + g(2) * (t2 * rSges(6,1) - t1 * rSges(6,2) + t31) + g(3) * (-t17 * rSges(6,3) + t34) + (g(3) * (rSges(6,1) * t19 - rSges(6,2) * t18) + t58 * (rSges(6,3) - t30)) * t16) - m(7) * (g(1) * (t57 * t3 + t60 * t4 + t32) + g(2) * (t57 * t1 + t60 * t2 + t31) + g(3) * (-t17 * rSges(7,2) + t34) + (g(3) * (t57 * t18 + t60 * t19) + t58 * (rSges(7,2) - t30)) * t16);
U  = t5;
