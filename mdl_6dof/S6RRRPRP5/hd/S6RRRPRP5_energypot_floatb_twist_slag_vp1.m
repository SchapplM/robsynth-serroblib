% Calculate potential energy for
% S6RRRPRP5
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
% Datum: 2019-03-09 16:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPRP5_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP5_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRPRP5_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP5_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP5_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP5_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP5_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:47:38
% EndTime: 2019-03-09 16:47:38
% DurationCPUTime: 0.60s
% Computational Cost: add. (268->114), mult. (255->130), div. (0->0), fcn. (247->10), ass. (0->39)
t54 = rSges(7,1) + pkin(5);
t53 = rSges(4,3) + pkin(8);
t25 = -qJ(4) - pkin(8);
t52 = rSges(5,3) - t25;
t28 = sin(qJ(1));
t31 = cos(qJ(1));
t51 = g(1) * t31 + g(2) * t28;
t50 = rSges(7,3) + qJ(6);
t29 = cos(qJ(3));
t15 = t29 * pkin(3) + pkin(2);
t27 = sin(qJ(2));
t46 = rSges(3,2) * t27;
t26 = sin(qJ(3));
t45 = t26 * t31;
t44 = t28 * t26;
t30 = cos(qJ(2));
t43 = t28 * t30;
t42 = t30 * t31;
t38 = pkin(6) + r_base(3);
t24 = qJ(3) + pkin(10);
t37 = t28 * pkin(1) + r_base(2);
t36 = t31 * pkin(1) + t28 * pkin(7) + r_base(1);
t23 = -pkin(9) + t25;
t17 = cos(t24);
t9 = pkin(4) * t17 + t15;
t35 = t30 * t23 + t27 * t9 + t38;
t34 = -t31 * pkin(7) + t37;
t16 = sin(t24);
t10 = t26 * pkin(3) + pkin(4) * t16;
t33 = t28 * t10 + t9 * t42 + t36;
t32 = -t10 * t31 + t9 * t43 + t34;
t18 = qJ(5) + t24;
t14 = cos(t18);
t13 = sin(t18);
t4 = t28 * t13 + t14 * t42;
t3 = t13 * t42 - t28 * t14;
t2 = -t13 * t31 + t14 * t43;
t1 = t13 * t43 + t14 * t31;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t31 - t28 * rSges(2,2) + r_base(1)) + g(2) * (t28 * rSges(2,1) + rSges(2,2) * t31 + r_base(2)) + g(3) * (rSges(2,3) + t38)) - m(3) * (g(1) * (t28 * rSges(3,3) + t36) + g(2) * (rSges(3,1) * t43 - t28 * t46 + t37) + g(3) * (rSges(3,1) * t27 + rSges(3,2) * t30 + t38) + (g(1) * (rSges(3,1) * t30 - t46) + g(2) * (-rSges(3,3) - pkin(7))) * t31) - m(4) * (g(1) * (pkin(2) * t42 + (t29 * t42 + t44) * rSges(4,1) + (-t26 * t42 + t28 * t29) * rSges(4,2) + t36) + g(2) * (pkin(2) * t43 + (t29 * t43 - t45) * rSges(4,1) + (-t26 * t43 - t29 * t31) * rSges(4,2) + t34) + g(3) * (-t53 * t30 + t38) + (g(3) * (rSges(4,1) * t29 - rSges(4,2) * t26 + pkin(2)) + t51 * t53) * t27) - m(5) * (g(1) * (t15 * t42 + pkin(3) * t44 + (t28 * t16 + t17 * t42) * rSges(5,1) + (-t16 * t42 + t28 * t17) * rSges(5,2) + t36) + g(2) * (t15 * t43 - pkin(3) * t45 + (-t16 * t31 + t17 * t43) * rSges(5,1) + (-t16 * t43 - t17 * t31) * rSges(5,2) + t34) + g(3) * (-t52 * t30 + t38) + (g(3) * (rSges(5,1) * t17 - rSges(5,2) * t16 + t15) + t51 * t52) * t27) - m(6) * (g(1) * (t4 * rSges(6,1) - t3 * rSges(6,2) + t33) + g(2) * (t2 * rSges(6,1) - t1 * rSges(6,2) + t32) + g(3) * (-rSges(6,3) * t30 + t35) + (g(3) * (rSges(6,1) * t14 - rSges(6,2) * t13) + t51 * (rSges(6,3) - t23)) * t27) - m(7) * (g(1) * (t50 * t3 + t54 * t4 + t33) + g(2) * (t50 * t1 + t54 * t2 + t32) + g(3) * (-rSges(7,2) * t30 + t35) + (g(3) * (t50 * t13 + t54 * t14) + t51 * (rSges(7,2) - t23)) * t27);
U  = t5;
