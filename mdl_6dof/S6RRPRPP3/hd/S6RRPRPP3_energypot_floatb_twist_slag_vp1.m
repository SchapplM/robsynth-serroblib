% Calculate potential energy for
% S6RRPRPP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
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
% Datum: 2019-03-09 09:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRPP3_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP3_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPRPP3_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP3_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP3_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP3_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPP3_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:54:08
% EndTime: 2019-03-09 09:54:09
% DurationCPUTime: 0.50s
% Computational Cost: add. (245->107), mult. (280->122), div. (0->0), fcn. (278->8), ass. (0->36)
t55 = -rSges(7,1) - pkin(5);
t27 = sin(qJ(1));
t29 = cos(qJ(1));
t54 = g(1) * t29 + g(2) * t27;
t53 = rSges(4,3) + qJ(3);
t52 = rSges(7,3) + qJ(6);
t26 = sin(qJ(2));
t49 = rSges(3,2) * t26;
t23 = sin(pkin(9));
t48 = t23 * t29;
t47 = t27 * t23;
t28 = cos(qJ(2));
t46 = t27 * t28;
t45 = t28 * t29;
t22 = pkin(9) + qJ(4);
t17 = sin(t22);
t44 = t29 * t17;
t40 = pkin(6) + r_base(3);
t39 = t27 * pkin(1) + r_base(2);
t37 = t29 * pkin(1) + t27 * pkin(7) + r_base(1);
t24 = cos(pkin(9));
t15 = pkin(3) * t24 + pkin(2);
t25 = -pkin(8) - qJ(3);
t36 = t26 * t15 + t28 * t25 + t40;
t35 = -t29 * pkin(7) + t39;
t34 = pkin(3) * t47 + t15 * t45 + t37;
t18 = cos(t22);
t33 = t36 + (pkin(4) * t18 + qJ(5) * t17) * t26;
t5 = -t27 * t18 + t28 * t44;
t6 = t27 * t17 + t18 * t45;
t32 = t6 * pkin(4) + t5 * qJ(5) + t34;
t31 = -pkin(3) * t48 + t15 * t46 + t35;
t3 = t17 * t46 + t18 * t29;
t4 = t18 * t46 - t44;
t30 = t4 * pkin(4) + t3 * qJ(5) + t31;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t29 - t27 * rSges(2,2) + r_base(1)) + g(2) * (t27 * rSges(2,1) + rSges(2,2) * t29 + r_base(2)) + g(3) * (rSges(2,3) + t40)) - m(3) * (g(1) * (t27 * rSges(3,3) + t37) + g(2) * (rSges(3,1) * t46 - t27 * t49 + t39) + g(3) * (rSges(3,1) * t26 + rSges(3,2) * t28 + t40) + (g(1) * (rSges(3,1) * t28 - t49) + g(2) * (-rSges(3,3) - pkin(7))) * t29) - m(4) * (g(1) * (pkin(2) * t45 + (t24 * t45 + t47) * rSges(4,1) + (-t23 * t45 + t27 * t24) * rSges(4,2) + t37) + g(2) * (pkin(2) * t46 + (t24 * t46 - t48) * rSges(4,1) + (-t23 * t46 - t24 * t29) * rSges(4,2) + t35) + g(3) * (-t28 * t53 + t40) + (g(3) * (rSges(4,1) * t24 - rSges(4,2) * t23 + pkin(2)) + t54 * t53) * t26) - m(5) * (g(1) * (t6 * rSges(5,1) - t5 * rSges(5,2) + t34) + g(2) * (t4 * rSges(5,1) - t3 * rSges(5,2) + t31) + g(3) * (-t28 * rSges(5,3) + t36) + (g(3) * (rSges(5,1) * t18 - rSges(5,2) * t17) + t54 * (rSges(5,3) - t25)) * t26) - m(6) * (g(1) * (-t6 * rSges(6,2) + t5 * rSges(6,3) + t32) + g(2) * (-t4 * rSges(6,2) + t3 * rSges(6,3) + t30) + g(3) * (-t28 * rSges(6,1) + t33) + (g(3) * (-rSges(6,2) * t18 + rSges(6,3) * t17) + t54 * (rSges(6,1) - t25)) * t26) - m(7) * (g(1) * (t5 * rSges(7,2) + t52 * t6 + t32) + g(2) * (t3 * rSges(7,2) + t4 * t52 + t30) + g(3) * (t55 * t28 + t33) + (g(3) * (rSges(7,2) * t17 + t18 * t52) + t54 * (-t25 - t55)) * t26);
U  = t1;
