% Calculate potential energy for
% S6RPRPRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-03-09 03:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPRP3_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP3_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRPRP3_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP3_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP3_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP3_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP3_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:07:43
% EndTime: 2019-03-09 03:07:43
% DurationCPUTime: 0.46s
% Computational Cost: add. (271->101), mult. (213->113), div. (0->0), fcn. (201->10), ass. (0->40)
t55 = rSges(7,1) + pkin(5);
t22 = qJ(1) + pkin(9);
t16 = sin(t22);
t18 = cos(t22);
t54 = g(1) * t18 + g(2) * t16;
t53 = rSges(5,3) + qJ(4);
t52 = rSges(7,3) + qJ(6);
t26 = sin(qJ(3));
t49 = rSges(4,2) * t26;
t23 = sin(pkin(10));
t48 = t16 * t23;
t28 = cos(qJ(3));
t47 = t16 * t28;
t46 = t18 * t23;
t45 = t18 * t28;
t44 = t23 * t28;
t24 = cos(pkin(10));
t43 = t24 * t28;
t39 = pkin(6) + r_base(3);
t27 = sin(qJ(1));
t38 = t27 * pkin(1) + r_base(2);
t29 = cos(qJ(1));
t37 = t29 * pkin(1) + r_base(1);
t36 = qJ(2) + t39;
t35 = t16 * pkin(2) + t38;
t34 = t18 * pkin(2) + t16 * pkin(7) + t37;
t13 = t24 * pkin(4) + pkin(3);
t25 = -pkin(8) - qJ(4);
t33 = t26 * t13 + t28 * t25 + t36;
t32 = -t18 * pkin(7) + t35;
t31 = pkin(4) * t48 + t13 * t45 + t34;
t30 = -pkin(4) * t46 + t13 * t47 + t32;
t21 = pkin(10) + qJ(5);
t17 = cos(t21);
t15 = sin(t21);
t4 = t16 * t15 + t17 * t45;
t3 = t15 * t45 - t16 * t17;
t2 = -t18 * t15 + t17 * t47;
t1 = t15 * t47 + t18 * t17;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t29 * rSges(2,1) - t27 * rSges(2,2) + r_base(1)) + g(2) * (t27 * rSges(2,1) + t29 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t39)) - m(3) * (g(1) * (t18 * rSges(3,1) - t16 * rSges(3,2) + t37) + g(2) * (t16 * rSges(3,1) + t18 * rSges(3,2) + t38) + g(3) * (rSges(3,3) + t36)) - m(4) * (g(1) * (t16 * rSges(4,3) + t34) + g(2) * (rSges(4,1) * t47 - t16 * t49 + t35) + g(3) * (t26 * rSges(4,1) + t28 * rSges(4,2) + t36) + (g(1) * (rSges(4,1) * t28 - t49) + g(2) * (-rSges(4,3) - pkin(7))) * t18) - m(5) * (g(1) * (pkin(3) * t45 + (t18 * t43 + t48) * rSges(5,1) + (t16 * t24 - t18 * t44) * rSges(5,2) + t34) + g(2) * (pkin(3) * t47 + (t16 * t43 - t46) * rSges(5,1) + (-t16 * t44 - t18 * t24) * rSges(5,2) + t32) + g(3) * (-t53 * t28 + t36) + (g(3) * (rSges(5,1) * t24 - rSges(5,2) * t23 + pkin(3)) + t54 * t53) * t26) - m(6) * (g(1) * (t4 * rSges(6,1) - t3 * rSges(6,2) + t31) + g(2) * (t2 * rSges(6,1) - t1 * rSges(6,2) + t30) + g(3) * (-t28 * rSges(6,3) + t33) + (g(3) * (rSges(6,1) * t17 - rSges(6,2) * t15) + t54 * (rSges(6,3) - t25)) * t26) - m(7) * (g(1) * (t52 * t3 + t55 * t4 + t31) + g(2) * (t52 * t1 + t55 * t2 + t30) + g(3) * (-t28 * rSges(7,2) + t33) + (g(3) * (t52 * t15 + t55 * t17) + t54 * (rSges(7,2) - t25)) * t26);
U  = t5;
