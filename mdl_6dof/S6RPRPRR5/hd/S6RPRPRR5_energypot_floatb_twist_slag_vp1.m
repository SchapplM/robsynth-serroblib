% Calculate potential energy for
% S6RPRPRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
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
% Datum: 2019-03-09 03:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPRR5_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR5_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRPRR5_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR5_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR5_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR5_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR5_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:47:51
% EndTime: 2019-03-09 03:47:51
% DurationCPUTime: 0.42s
% Computational Cost: add. (252->104), mult. (247->116), div. (0->0), fcn. (244->10), ass. (0->40)
t53 = pkin(9) + rSges(7,3);
t28 = -pkin(7) - qJ(2);
t54 = -pkin(8) - t28;
t25 = pkin(10) + qJ(3);
t21 = sin(t25);
t33 = cos(qJ(5));
t52 = t21 * t33;
t34 = cos(qJ(1));
t51 = t21 * t34;
t22 = cos(t25);
t31 = sin(qJ(1));
t50 = t22 * t31;
t49 = t22 * t34;
t48 = qJ(4) * t21;
t47 = rSges(3,3) + qJ(2);
t46 = pkin(6) + r_base(3);
t27 = cos(pkin(10));
t19 = t27 * pkin(2) + pkin(1);
t45 = t34 * t19 + r_base(1);
t26 = sin(pkin(10));
t44 = t26 * pkin(2) + t46;
t43 = t31 * t19 + t34 * t28 + r_base(2);
t42 = t21 * pkin(3) + t44;
t41 = pkin(3) * t49 + t34 * t48 + t45;
t30 = sin(qJ(5));
t5 = t21 * t30 + t22 * t33;
t40 = pkin(3) * t50 + t31 * t48 + t43;
t39 = pkin(4) * t49 + t41;
t38 = rSges(3,1) * t27 - rSges(3,2) * t26 + pkin(1);
t29 = sin(qJ(6));
t32 = cos(qJ(6));
t37 = t32 * rSges(7,1) - t29 * rSges(7,2) + pkin(5);
t36 = pkin(4) * t50 + t34 * pkin(8) + t40;
t35 = t21 * pkin(4) - t22 * qJ(4) + t42;
t6 = -t22 * t30 + t52;
t4 = t5 * t34;
t3 = t30 * t49 - t33 * t51;
t2 = t5 * t31;
t1 = t30 * t50 - t31 * t52;
t7 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t34 * rSges(2,1) - t31 * rSges(2,2) + r_base(1)) + g(2) * (t31 * rSges(2,1) + t34 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t46)) - m(3) * (g(1) * r_base(1) + g(2) * r_base(2) + g(3) * (t26 * rSges(3,1) + t27 * rSges(3,2) + t46) + (g(1) * t38 - g(2) * t47) * t34 + (g(1) * t47 + g(2) * t38) * t31) - m(4) * (g(1) * (rSges(4,1) * t49 - rSges(4,2) * t51 + t45) + g(2) * (-t34 * rSges(4,3) + t43) + g(3) * (t21 * rSges(4,1) + t22 * rSges(4,2) + t44) + (g(1) * (rSges(4,3) - t28) + g(2) * (rSges(4,1) * t22 - rSges(4,2) * t21)) * t31) - m(5) * (g(1) * (rSges(5,1) * t49 + rSges(5,3) * t51 + t41) + g(2) * (-t34 * rSges(5,2) + t40) + g(3) * (t21 * rSges(5,1) + (-rSges(5,3) - qJ(4)) * t22 + t42) + (g(1) * (rSges(5,2) - t28) + g(2) * (rSges(5,1) * t22 + rSges(5,3) * t21)) * t31) - m(6) * (g(2) * (t2 * rSges(6,1) - t1 * rSges(6,2) + t34 * rSges(6,3) + t36) + g(3) * (t6 * rSges(6,1) - t5 * rSges(6,2) + t35) + (t4 * rSges(6,1) - t3 * rSges(6,2) + t39 + (-rSges(6,3) + t54) * t31) * g(1)) - m(7) * (g(1) * (t37 * t4 + (-t29 * rSges(7,1) - t32 * rSges(7,2) + t54) * t31 + t53 * t3 + t39) + g(2) * (t2 * pkin(5) + (t2 * t32 + t34 * t29) * rSges(7,1) + (-t2 * t29 + t34 * t32) * rSges(7,2) + t53 * t1 + t36) + (t37 * t6 + t53 * t5 + t35) * g(3));
U  = t7;
