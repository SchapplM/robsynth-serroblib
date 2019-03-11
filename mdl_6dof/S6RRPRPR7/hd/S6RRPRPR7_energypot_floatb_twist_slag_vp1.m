% Calculate potential energy for
% S6RRPRPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
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
% Datum: 2019-03-09 10:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRPR7_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR7_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPRPR7_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR7_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR7_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR7_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR7_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:45:11
% EndTime: 2019-03-09 10:45:11
% DurationCPUTime: 0.47s
% Computational Cost: add. (229->113), mult. (290->127), div. (0->0), fcn. (293->10), ass. (0->45)
t57 = pkin(9) + rSges(7,3);
t27 = -qJ(5) - pkin(8);
t59 = -pkin(7) - t27;
t58 = -pkin(8) - rSges(5,3);
t26 = qJ(4) + pkin(10);
t21 = cos(t26);
t30 = sin(qJ(2));
t56 = t21 * t30;
t29 = sin(qJ(4));
t55 = t29 * t30;
t34 = cos(qJ(2));
t54 = t29 * t34;
t31 = sin(qJ(1));
t53 = t30 * t31;
t52 = t31 * t34;
t35 = cos(qJ(1));
t51 = t34 * t35;
t50 = qJ(3) * t30;
t49 = pkin(6) + r_base(3);
t48 = pkin(4) * t55;
t47 = t31 * pkin(1) + r_base(2);
t46 = t30 * pkin(2) + t49;
t45 = t35 * pkin(1) + t31 * pkin(7) + r_base(1);
t44 = pkin(2) * t52 + t31 * t50 + t47;
t20 = sin(t26);
t5 = t20 * t30 + t21 * t34;
t33 = cos(qJ(4));
t43 = t30 * t33 - t54;
t42 = t33 * t34 + t55;
t41 = pkin(2) * t51 + t35 * t50 + t45;
t28 = sin(qJ(6));
t32 = cos(qJ(6));
t40 = rSges(7,1) * t32 - rSges(7,2) * t28 + pkin(5);
t39 = -t34 * qJ(3) + t46;
t19 = pkin(4) * t33 + pkin(3);
t38 = t19 * t52 + t31 * t48 + t44;
t37 = t19 * t51 + t31 * t27 + t35 * t48 + t41;
t36 = t42 * rSges(5,1) + t43 * rSges(5,2) + t34 * pkin(3);
t13 = t30 * t19;
t6 = -t20 * t34 + t56;
t4 = t5 * t35;
t3 = t20 * t51 - t35 * t56;
t2 = t5 * t31;
t1 = t20 * t52 - t21 * t53;
t7 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t35 - t31 * rSges(2,2) + r_base(1)) + g(2) * (t31 * rSges(2,1) + rSges(2,2) * t35 + r_base(2)) + g(3) * (rSges(2,3) + t49)) - m(3) * (g(1) * (t31 * rSges(3,3) + t45) + g(2) * (rSges(3,1) * t52 - rSges(3,2) * t53 + t47) + g(3) * (rSges(3,1) * t30 + rSges(3,2) * t34 + t49) + (g(1) * (rSges(3,1) * t34 - rSges(3,2) * t30) + g(2) * (-rSges(3,3) - pkin(7))) * t35) - m(4) * (g(1) * (t31 * rSges(4,2) + t41) + g(2) * (rSges(4,1) * t52 + rSges(4,3) * t53 + t44) + g(3) * (rSges(4,1) * t30 + (-rSges(4,3) - qJ(3)) * t34 + t46) + (g(1) * (rSges(4,1) * t34 + rSges(4,3) * t30) + g(2) * (-rSges(4,2) - pkin(7))) * t35) - m(5) * (g(1) * t41 + g(2) * t44 + g(3) * (t43 * rSges(5,1) - t42 * rSges(5,2) + t30 * pkin(3) + t39) + (g(1) * t58 + g(2) * t36) * t31 + (g(1) * t36 + g(2) * (-pkin(7) - t58)) * t35) - m(6) * (g(1) * (rSges(6,1) * t4 - rSges(6,2) * t3 - rSges(6,3) * t31 + t37) + g(3) * (rSges(6,1) * t6 - rSges(6,2) * t5 + t13 + (-pkin(4) * t29 - qJ(3)) * t34 + t46) + (t2 * rSges(6,1) - t1 * rSges(6,2) + t38 + (rSges(6,3) + t59) * t35) * g(2)) - m(7) * (g(1) * (t4 * pkin(5) + (-t28 * t31 + t32 * t4) * rSges(7,1) + (-t28 * t4 - t31 * t32) * rSges(7,2) + t57 * t3 + t37) + g(2) * (t40 * t2 + (rSges(7,1) * t28 + rSges(7,2) * t32 + t59) * t35 + t57 * t1 + t38) + (-pkin(4) * t54 + t40 * t6 + t57 * t5 + t13 + t39) * g(3));
U  = t7;
