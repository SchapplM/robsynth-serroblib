% Calculate potential energy for
% S6RRPRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-03-09 10:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRPR8_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR8_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPRPR8_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR8_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR8_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR8_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR8_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:49:26
% EndTime: 2019-03-09 10:49:26
% DurationCPUTime: 0.60s
% Computational Cost: add. (265->115), mult. (306->137), div. (0->0), fcn. (314->10), ass. (0->36)
t54 = pkin(9) + rSges(7,3);
t27 = sin(qJ(1));
t30 = cos(qJ(1));
t53 = g(1) * t30 + g(2) * t27;
t52 = rSges(4,3) + qJ(3);
t26 = sin(qJ(2));
t49 = rSges(3,2) * t26;
t22 = sin(pkin(10));
t48 = t22 * t30;
t47 = t27 * t22;
t29 = cos(qJ(2));
t46 = t27 * t29;
t45 = t29 * t30;
t41 = pkin(6) + r_base(3);
t40 = t27 * pkin(1) + r_base(2);
t38 = t30 * pkin(1) + t27 * pkin(7) + r_base(1);
t23 = cos(pkin(10));
t14 = pkin(3) * t23 + pkin(2);
t24 = -pkin(8) - qJ(3);
t37 = t26 * t14 + t29 * t24 + t41;
t36 = -t30 * pkin(7) + t40;
t35 = pkin(3) * t47 + t14 * t45 + t38;
t21 = pkin(10) + qJ(4);
t16 = sin(t21);
t17 = cos(t21);
t34 = t37 + (pkin(4) * t17 + qJ(5) * t16) * t26;
t5 = t16 * t45 - t27 * t17;
t6 = t27 * t16 + t17 * t45;
t33 = t6 * pkin(4) + t5 * qJ(5) + t35;
t32 = -pkin(3) * t48 + t14 * t46 + t36;
t3 = t16 * t46 + t17 * t30;
t4 = -t16 * t30 + t17 * t46;
t31 = t4 * pkin(4) + t3 * qJ(5) + t32;
t28 = cos(qJ(6));
t25 = sin(qJ(6));
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t30 - t27 * rSges(2,2) + r_base(1)) + g(2) * (t27 * rSges(2,1) + rSges(2,2) * t30 + r_base(2)) + g(3) * (rSges(2,3) + t41)) - m(3) * (g(1) * (t27 * rSges(3,3) + t38) + g(2) * (rSges(3,1) * t46 - t27 * t49 + t40) + g(3) * (rSges(3,1) * t26 + rSges(3,2) * t29 + t41) + (g(1) * (rSges(3,1) * t29 - t49) + g(2) * (-rSges(3,3) - pkin(7))) * t30) - m(4) * (g(1) * (pkin(2) * t45 + (t23 * t45 + t47) * rSges(4,1) + (-t22 * t45 + t27 * t23) * rSges(4,2) + t38) + g(2) * (pkin(2) * t46 + (t23 * t46 - t48) * rSges(4,1) + (-t22 * t46 - t23 * t30) * rSges(4,2) + t36) + g(3) * (-t52 * t29 + t41) + (g(3) * (rSges(4,1) * t23 - rSges(4,2) * t22 + pkin(2)) + t53 * t52) * t26) - m(5) * (g(1) * (t6 * rSges(5,1) - t5 * rSges(5,2) + t35) + g(2) * (t4 * rSges(5,1) - t3 * rSges(5,2) + t32) + g(3) * (-rSges(5,3) * t29 + t37) + (g(3) * (rSges(5,1) * t17 - rSges(5,2) * t16) + t53 * (rSges(5,3) - t24)) * t26) - m(6) * (g(1) * (t6 * rSges(6,1) + t5 * rSges(6,3) + t33) + g(2) * (t4 * rSges(6,1) + t3 * rSges(6,3) + t31) + g(3) * (-rSges(6,2) * t29 + t34) + (g(3) * (rSges(6,1) * t17 + rSges(6,3) * t16) + t53 * (rSges(6,2) - t24)) * t26) - m(7) * (g(1) * (t6 * pkin(5) + (t25 * t5 + t28 * t6) * rSges(7,1) + (-t25 * t6 + t28 * t5) * rSges(7,2) + t33) + g(2) * (t4 * pkin(5) + (t25 * t3 + t28 * t4) * rSges(7,1) + (-t25 * t4 + t28 * t3) * rSges(7,2) + t31) + g(3) * (t54 * t29 + t34) + (g(3) * (t17 * pkin(5) + (t16 * t25 + t17 * t28) * rSges(7,1) + (t16 * t28 - t17 * t25) * rSges(7,2)) + t53 * (-t24 - t54)) * t26);
U  = t1;
