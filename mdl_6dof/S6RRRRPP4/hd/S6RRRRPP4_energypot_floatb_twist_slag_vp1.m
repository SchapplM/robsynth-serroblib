% Calculate potential energy for
% S6RRRRPP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
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
% Datum: 2019-03-09 21:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRPP4_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP4_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRRPP4_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP4_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP4_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP4_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPP4_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:58:22
% EndTime: 2019-03-09 20:58:23
% DurationCPUTime: 0.61s
% Computational Cost: add. (268->114), mult. (255->130), div. (0->0), fcn. (247->10), ass. (0->39)
t54 = rSges(7,1) + pkin(5);
t53 = rSges(4,3) + pkin(8);
t31 = -pkin(9) - pkin(8);
t52 = rSges(5,3) - t31;
t27 = sin(qJ(1));
t30 = cos(qJ(1));
t51 = g(1) * t30 + g(2) * t27;
t50 = rSges(7,3) + qJ(6);
t28 = cos(qJ(3));
t15 = t28 * pkin(3) + pkin(2);
t26 = sin(qJ(2));
t46 = rSges(3,2) * t26;
t25 = sin(qJ(3));
t45 = t25 * t27;
t44 = t25 * t30;
t29 = cos(qJ(2));
t43 = t27 * t29;
t42 = t29 * t30;
t24 = qJ(3) + qJ(4);
t38 = pkin(6) + r_base(3);
t37 = t27 * pkin(1) + r_base(2);
t36 = t30 * pkin(1) + t27 * pkin(7) + r_base(1);
t23 = -qJ(5) + t31;
t18 = cos(t24);
t9 = pkin(4) * t18 + t15;
t35 = t29 * t23 + t26 * t9 + t38;
t34 = -t30 * pkin(7) + t37;
t17 = sin(t24);
t10 = t25 * pkin(3) + pkin(4) * t17;
t33 = t27 * t10 + t9 * t42 + t36;
t32 = -t10 * t30 + t9 * t43 + t34;
t16 = pkin(10) + t24;
t14 = cos(t16);
t13 = sin(t16);
t4 = t13 * t27 + t14 * t42;
t3 = t13 * t42 - t27 * t14;
t2 = -t13 * t30 + t14 * t43;
t1 = t13 * t43 + t14 * t30;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t30 - rSges(2,2) * t27 + r_base(1)) + g(2) * (rSges(2,1) * t27 + rSges(2,2) * t30 + r_base(2)) + g(3) * (rSges(2,3) + t38)) - m(3) * (g(1) * (rSges(3,3) * t27 + t36) + g(2) * (rSges(3,1) * t43 - t27 * t46 + t37) + g(3) * (rSges(3,1) * t26 + rSges(3,2) * t29 + t38) + (g(1) * (rSges(3,1) * t29 - t46) + g(2) * (-rSges(3,3) - pkin(7))) * t30) - m(4) * (g(1) * (pkin(2) * t42 + (t28 * t42 + t45) * rSges(4,1) + (-t25 * t42 + t27 * t28) * rSges(4,2) + t36) + g(2) * (pkin(2) * t43 + (t28 * t43 - t44) * rSges(4,1) + (-t25 * t43 - t28 * t30) * rSges(4,2) + t34) + g(3) * (-t53 * t29 + t38) + (g(3) * (rSges(4,1) * t28 - rSges(4,2) * t25 + pkin(2)) + t51 * t53) * t26) - m(5) * (g(1) * (t15 * t42 + pkin(3) * t45 + (t17 * t27 + t18 * t42) * rSges(5,1) + (-t17 * t42 + t18 * t27) * rSges(5,2) + t36) + g(2) * (t15 * t43 - pkin(3) * t44 + (-t17 * t30 + t18 * t43) * rSges(5,1) + (-t17 * t43 - t18 * t30) * rSges(5,2) + t34) + g(3) * (-t52 * t29 + t38) + (g(3) * (rSges(5,1) * t18 - rSges(5,2) * t17 + t15) + t51 * t52) * t26) - m(6) * (g(1) * (rSges(6,1) * t4 - rSges(6,2) * t3 + t33) + g(2) * (rSges(6,1) * t2 - rSges(6,2) * t1 + t32) + g(3) * (-rSges(6,3) * t29 + t35) + (g(3) * (rSges(6,1) * t14 - rSges(6,2) * t13) + t51 * (rSges(6,3) - t23)) * t26) - m(7) * (g(1) * (t50 * t3 + t54 * t4 + t33) + g(2) * (t50 * t1 + t54 * t2 + t32) + g(3) * (-rSges(7,2) * t29 + t35) + (g(3) * (t50 * t13 + t54 * t14) + t51 * (rSges(7,2) - t23)) * t26);
U  = t5;
