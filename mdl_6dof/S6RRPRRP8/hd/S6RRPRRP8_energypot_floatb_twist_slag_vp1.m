% Calculate potential energy for
% S6RRPRRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 12:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRP8_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP8_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPRRP8_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP8_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP8_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP8_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP8_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:21:51
% EndTime: 2019-03-09 12:21:52
% DurationCPUTime: 0.60s
% Computational Cost: add. (268->114), mult. (255->130), div. (0->0), fcn. (247->10), ass. (0->39)
t54 = rSges(7,1) + pkin(5);
t27 = -pkin(8) - qJ(3);
t53 = rSges(5,3) - t27;
t29 = sin(qJ(1));
t31 = cos(qJ(1));
t52 = g(1) * t31 + g(2) * t29;
t51 = rSges(4,3) + qJ(3);
t50 = rSges(7,3) + qJ(6);
t26 = cos(pkin(10));
t15 = t26 * pkin(3) + pkin(2);
t28 = sin(qJ(2));
t47 = rSges(3,2) * t28;
t25 = sin(pkin(10));
t46 = t25 * t31;
t45 = t29 * t25;
t30 = cos(qJ(2));
t44 = t29 * t30;
t43 = t30 * t31;
t38 = pkin(6) + r_base(3);
t24 = pkin(10) + qJ(4);
t37 = t29 * pkin(1) + r_base(2);
t36 = t31 * pkin(1) + t29 * pkin(7) + r_base(1);
t23 = -pkin(9) + t27;
t17 = cos(t24);
t9 = pkin(4) * t17 + t15;
t35 = t30 * t23 + t28 * t9 + t38;
t34 = -t31 * pkin(7) + t37;
t16 = sin(t24);
t10 = t25 * pkin(3) + pkin(4) * t16;
t33 = t29 * t10 + t9 * t43 + t36;
t32 = -t10 * t31 + t9 * t44 + t34;
t18 = qJ(5) + t24;
t14 = cos(t18);
t13 = sin(t18);
t4 = t29 * t13 + t14 * t43;
t3 = t13 * t43 - t29 * t14;
t2 = -t13 * t31 + t14 * t44;
t1 = t13 * t44 + t14 * t31;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t31 - t29 * rSges(2,2) + r_base(1)) + g(2) * (t29 * rSges(2,1) + rSges(2,2) * t31 + r_base(2)) + g(3) * (rSges(2,3) + t38)) - m(3) * (g(1) * (t29 * rSges(3,3) + t36) + g(2) * (rSges(3,1) * t44 - t29 * t47 + t37) + g(3) * (rSges(3,1) * t28 + rSges(3,2) * t30 + t38) + (g(1) * (rSges(3,1) * t30 - t47) + g(2) * (-rSges(3,3) - pkin(7))) * t31) - m(4) * (g(1) * (pkin(2) * t43 + (t26 * t43 + t45) * rSges(4,1) + (-t25 * t43 + t29 * t26) * rSges(4,2) + t36) + g(2) * (pkin(2) * t44 + (t26 * t44 - t46) * rSges(4,1) + (-t25 * t44 - t26 * t31) * rSges(4,2) + t34) + g(3) * (-t51 * t30 + t38) + (g(3) * (rSges(4,1) * t26 - rSges(4,2) * t25 + pkin(2)) + t52 * t51) * t28) - m(5) * (g(1) * (t15 * t43 + pkin(3) * t45 + (t29 * t16 + t17 * t43) * rSges(5,1) + (-t16 * t43 + t29 * t17) * rSges(5,2) + t36) + g(2) * (t15 * t44 - pkin(3) * t46 + (-t16 * t31 + t17 * t44) * rSges(5,1) + (-t16 * t44 - t17 * t31) * rSges(5,2) + t34) + g(3) * (-t53 * t30 + t38) + (g(3) * (rSges(5,1) * t17 - rSges(5,2) * t16 + t15) + t52 * t53) * t28) - m(6) * (g(1) * (t4 * rSges(6,1) - t3 * rSges(6,2) + t33) + g(2) * (t2 * rSges(6,1) - t1 * rSges(6,2) + t32) + g(3) * (-rSges(6,3) * t30 + t35) + (g(3) * (rSges(6,1) * t14 - rSges(6,2) * t13) + t52 * (rSges(6,3) - t23)) * t28) - m(7) * (g(1) * (t50 * t3 + t54 * t4 + t33) + g(2) * (t50 * t1 + t54 * t2 + t32) + g(3) * (-rSges(7,2) * t30 + t35) + (g(3) * (t50 * t13 + t54 * t14) + t52 * (rSges(7,2) - t23)) * t28);
U  = t5;
