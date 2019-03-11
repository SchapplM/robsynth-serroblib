% Calculate potential energy for
% S6RRRRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
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
% Datum: 2019-03-09 22:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRPR8_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR8_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRRPR8_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR8_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR8_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR8_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR8_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:38:44
% EndTime: 2019-03-09 22:38:44
% DurationCPUTime: 0.60s
% Computational Cost: add. (265->115), mult. (306->137), div. (0->0), fcn. (314->10), ass. (0->36)
t54 = rSges(4,3) + pkin(8);
t53 = pkin(10) + rSges(7,3);
t25 = sin(qJ(1));
t29 = cos(qJ(1));
t52 = g(1) * t29 + g(2) * t25;
t24 = sin(qJ(2));
t48 = rSges(3,2) * t24;
t23 = sin(qJ(3));
t47 = t23 * t25;
t46 = t23 * t29;
t28 = cos(qJ(2));
t45 = t25 * t28;
t44 = t28 * t29;
t41 = pkin(6) + r_base(3);
t40 = t25 * pkin(1) + r_base(2);
t38 = t29 * pkin(1) + t25 * pkin(7) + r_base(1);
t27 = cos(qJ(3));
t14 = pkin(3) * t27 + pkin(2);
t30 = -pkin(9) - pkin(8);
t37 = t24 * t14 + t28 * t30 + t41;
t36 = -t29 * pkin(7) + t40;
t35 = pkin(3) * t47 + t14 * t44 + t38;
t21 = qJ(3) + qJ(4);
t16 = sin(t21);
t17 = cos(t21);
t34 = t37 + (pkin(4) * t17 + qJ(5) * t16) * t24;
t5 = t16 * t44 - t25 * t17;
t6 = t16 * t25 + t17 * t44;
t33 = t6 * pkin(4) + t5 * qJ(5) + t35;
t32 = -pkin(3) * t46 + t14 * t45 + t36;
t3 = t16 * t45 + t17 * t29;
t4 = -t16 * t29 + t17 * t45;
t31 = t4 * pkin(4) + t3 * qJ(5) + t32;
t26 = cos(qJ(6));
t22 = sin(qJ(6));
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t29 - rSges(2,2) * t25 + r_base(1)) + g(2) * (rSges(2,1) * t25 + rSges(2,2) * t29 + r_base(2)) + g(3) * (rSges(2,3) + t41)) - m(3) * (g(1) * (rSges(3,3) * t25 + t38) + g(2) * (rSges(3,1) * t45 - t25 * t48 + t40) + g(3) * (rSges(3,1) * t24 + rSges(3,2) * t28 + t41) + (g(1) * (rSges(3,1) * t28 - t48) + g(2) * (-rSges(3,3) - pkin(7))) * t29) - m(4) * (g(1) * (pkin(2) * t44 + (t27 * t44 + t47) * rSges(4,1) + (-t23 * t44 + t25 * t27) * rSges(4,2) + t38) + g(2) * (pkin(2) * t45 + (t27 * t45 - t46) * rSges(4,1) + (-t23 * t45 - t27 * t29) * rSges(4,2) + t36) + g(3) * (-t54 * t28 + t41) + (g(3) * (rSges(4,1) * t27 - rSges(4,2) * t23 + pkin(2)) + t52 * t54) * t24) - m(5) * (g(1) * (t6 * rSges(5,1) - t5 * rSges(5,2) + t35) + g(2) * (t4 * rSges(5,1) - t3 * rSges(5,2) + t32) + g(3) * (-rSges(5,3) * t28 + t37) + (g(3) * (rSges(5,1) * t17 - rSges(5,2) * t16) + t52 * (rSges(5,3) - t30)) * t24) - m(6) * (g(1) * (t6 * rSges(6,1) + t5 * rSges(6,3) + t33) + g(2) * (t4 * rSges(6,1) + t3 * rSges(6,3) + t31) + g(3) * (-rSges(6,2) * t28 + t34) + (g(3) * (rSges(6,1) * t17 + rSges(6,3) * t16) + t52 * (rSges(6,2) - t30)) * t24) - m(7) * (g(1) * (t6 * pkin(5) + (t22 * t5 + t26 * t6) * rSges(7,1) + (-t22 * t6 + t26 * t5) * rSges(7,2) + t33) + g(2) * (t4 * pkin(5) + (t22 * t3 + t26 * t4) * rSges(7,1) + (-t22 * t4 + t26 * t3) * rSges(7,2) + t31) + g(3) * (t53 * t28 + t34) + (g(3) * (t17 * pkin(5) + (t16 * t22 + t17 * t26) * rSges(7,1) + (t16 * t26 - t17 * t22) * rSges(7,2)) + t52 * (-t30 - t53)) * t24);
U  = t1;
