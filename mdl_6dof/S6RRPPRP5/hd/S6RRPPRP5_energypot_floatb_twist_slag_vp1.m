% Calculate potential energy for
% S6RRPPRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta4]';
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
% Datum: 2019-03-09 08:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPPRP5_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP5_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPPRP5_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP5_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP5_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP5_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRP5_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:41:19
% EndTime: 2019-03-09 08:41:20
% DurationCPUTime: 0.50s
% Computational Cost: add. (204->109), mult. (252->122), div. (0->0), fcn. (240->8), ass. (0->40)
t56 = rSges(7,1) + pkin(5);
t25 = sin(qJ(1));
t27 = cos(qJ(1));
t55 = g(1) * t27 + g(2) * t25;
t54 = rSges(5,3) + qJ(4);
t53 = rSges(7,3) + qJ(6);
t24 = sin(qJ(2));
t50 = t24 * t25;
t49 = t24 * t27;
t20 = pkin(9) + qJ(5);
t15 = cos(t20);
t48 = t25 * t15;
t21 = sin(pkin(9));
t47 = t25 * t21;
t22 = cos(pkin(9));
t46 = t25 * t22;
t26 = cos(qJ(2));
t45 = t25 * t26;
t42 = qJ(3) * t24;
t40 = pkin(6) + r_base(3);
t39 = t25 * pkin(1) + r_base(2);
t38 = t21 * t49;
t37 = t24 * t47;
t36 = t24 * pkin(2) + t40;
t35 = -pkin(4) * t21 - qJ(3);
t34 = t27 * pkin(1) + t25 * pkin(7) + r_base(1);
t33 = pkin(2) * t45 + t25 * t42 + t39;
t32 = t34 + (pkin(2) * t26 + t42) * t27;
t23 = -pkin(8) - qJ(4);
t31 = -t24 * t23 + t36;
t30 = -t27 * pkin(7) + t33;
t13 = pkin(4) * t22 + pkin(3);
t29 = pkin(4) * t38 + t25 * t13 + t32;
t28 = pkin(4) * t37 - t27 * t13 + t30;
t14 = sin(t20);
t4 = t14 * t50 - t15 * t27;
t3 = t14 * t27 + t24 * t48;
t2 = t14 * t49 + t48;
t1 = t14 * t25 - t15 * t49;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t27 - t25 * rSges(2,2) + r_base(1)) + g(2) * (t25 * rSges(2,1) + rSges(2,2) * t27 + r_base(2)) + g(3) * (rSges(2,3) + t40)) - m(3) * (g(1) * (t25 * rSges(3,3) + t34) + g(2) * (rSges(3,1) * t45 - rSges(3,2) * t50 + t39) + g(3) * (rSges(3,1) * t24 + rSges(3,2) * t26 + t40) + (g(1) * (rSges(3,1) * t26 - rSges(3,2) * t24) + g(2) * (-rSges(3,3) - pkin(7))) * t27) - m(4) * (g(1) * (t25 * rSges(4,1) + t32) + g(2) * (-rSges(4,2) * t45 + rSges(4,3) * t50 + t33) + g(3) * (-t24 * rSges(4,2) + (-rSges(4,3) - qJ(3)) * t26 + t36) + (g(1) * (-rSges(4,2) * t26 + rSges(4,3) * t24) + g(2) * (-rSges(4,1) - pkin(7))) * t27) - m(5) * (g(1) * (t25 * pkin(3) + (t38 + t46) * rSges(5,1) + (t22 * t49 - t47) * rSges(5,2) + t32) + g(2) * (-t27 * pkin(3) + (-t22 * t27 + t37) * rSges(5,1) + (t21 * t27 + t24 * t46) * rSges(5,2) + t30) + g(3) * (t54 * t24 + t36) + (g(3) * (-rSges(5,1) * t21 - rSges(5,2) * t22 - qJ(3)) + t55 * t54) * t26) - m(6) * (g(1) * (t2 * rSges(6,1) - t1 * rSges(6,2) + t29) + g(2) * (t4 * rSges(6,1) + t3 * rSges(6,2) + t28) + g(3) * (t24 * rSges(6,3) + t31) + (g(3) * (-rSges(6,1) * t14 - rSges(6,2) * t15 + t35) + t55 * (rSges(6,3) - t23)) * t26) - m(7) * (g(1) * (t53 * t1 + t2 * t56 + t29) + g(2) * (-t53 * t3 + t4 * t56 + t28) + g(3) * (rSges(7,2) * t24 + t31) + (g(3) * (-t14 * t56 + t53 * t15 + t35) + t55 * (rSges(7,2) - t23)) * t26);
U  = t5;
