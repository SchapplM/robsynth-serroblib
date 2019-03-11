% Calculate potential energy for
% S6PRRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-03-09 00:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRRRP1_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP1_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRRRRP1_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP1_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP1_energypot_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP1_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRP1_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:55:38
% EndTime: 2019-03-08 23:55:38
% DurationCPUTime: 0.40s
% Computational Cost: add. (364->131), mult. (545->162), div. (0->0), fcn. (633->12), ass. (0->56)
t43 = sin(qJ(3));
t71 = pkin(3) * t43;
t70 = rSges(6,3) + pkin(10);
t69 = pkin(8) + rSges(4,3);
t37 = sin(pkin(11));
t39 = cos(pkin(11));
t44 = sin(qJ(2));
t40 = cos(pkin(6));
t47 = cos(qJ(2));
t59 = t40 * t47;
t17 = t37 * t44 - t39 * t59;
t42 = sin(qJ(5));
t68 = t17 * t42;
t19 = t37 * t59 + t39 * t44;
t67 = t19 * t42;
t38 = sin(pkin(6));
t66 = t37 * t38;
t65 = t38 * t39;
t64 = t38 * t43;
t63 = t38 * t44;
t46 = cos(qJ(3));
t62 = t38 * t46;
t61 = t38 * t47;
t60 = t40 * t44;
t58 = rSges(7,3) + qJ(6) + pkin(10);
t57 = t37 * pkin(1) + r_base(2);
t56 = t37 * t64;
t55 = t42 * t61;
t54 = qJ(1) + r_base(3);
t53 = t39 * pkin(1) + pkin(7) * t66 + r_base(1);
t52 = t40 * pkin(7) + t54;
t20 = -t37 * t60 + t39 * t47;
t30 = t46 * pkin(3) + pkin(2);
t48 = -pkin(9) - pkin(8);
t51 = pkin(3) * t56 - t19 * t48 + t20 * t30 + t53;
t50 = t30 * t63 + t40 * t71 + t48 * t61 + t52;
t18 = t37 * t47 + t39 * t60;
t49 = t18 * t30 + (-pkin(7) - t71) * t65 - t17 * t48 + t57;
t45 = cos(qJ(5));
t36 = qJ(3) + qJ(4);
t32 = cos(t36);
t31 = sin(t36);
t29 = t45 * pkin(5) + pkin(4);
t14 = t40 * t31 + t32 * t63;
t13 = t31 * t63 - t40 * t32;
t10 = t14 * t45 - t55;
t9 = -t14 * t42 - t45 * t61;
t8 = t20 * t32 + t31 * t66;
t7 = t20 * t31 - t32 * t66;
t6 = t18 * t32 - t31 * t65;
t5 = t18 * t31 + t32 * t65;
t4 = t8 * t45 + t67;
t3 = t19 * t45 - t8 * t42;
t2 = t6 * t45 + t68;
t1 = t17 * t45 - t6 * t42;
t11 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t39 * rSges(2,1) - t37 * rSges(2,2) + r_base(1)) + g(2) * (t37 * rSges(2,1) + t39 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t54)) - m(3) * (g(1) * (t20 * rSges(3,1) - t19 * rSges(3,2) + t53) + g(2) * (t18 * rSges(3,1) - t17 * rSges(3,2) + t57) + g(3) * (t40 * rSges(3,3) + t52) + (g(1) * rSges(3,3) * t37 + g(3) * (rSges(3,1) * t44 + rSges(3,2) * t47) + g(2) * (-rSges(3,3) - pkin(7)) * t39) * t38) - m(4) * (g(1) * (t20 * pkin(2) + (t20 * t46 + t56) * rSges(4,1) + (-t20 * t43 + t37 * t62) * rSges(4,2) + t69 * t19 + t53) + g(2) * (t18 * pkin(2) - pkin(7) * t65 + (t18 * t46 - t39 * t64) * rSges(4,1) + (-t18 * t43 - t39 * t62) * rSges(4,2) + t69 * t17 + t57) + g(3) * ((t43 * rSges(4,1) + t46 * rSges(4,2)) * t40 + (-t69 * t47 + (t46 * rSges(4,1) - t43 * rSges(4,2) + pkin(2)) * t44) * t38 + t52)) - m(5) * (g(1) * (t8 * rSges(5,1) - t7 * rSges(5,2) + t19 * rSges(5,3) + t51) + g(2) * (t6 * rSges(5,1) - t5 * rSges(5,2) + t17 * rSges(5,3) + t49) + g(3) * (t14 * rSges(5,1) - t13 * rSges(5,2) - rSges(5,3) * t61 + t50)) - m(6) * (g(1) * (t4 * rSges(6,1) + t3 * rSges(6,2) + t8 * pkin(4) + t70 * t7 + t51) + g(2) * (t2 * rSges(6,1) + t1 * rSges(6,2) + t6 * pkin(4) + t70 * t5 + t49) + g(3) * (t10 * rSges(6,1) + t9 * rSges(6,2) + t14 * pkin(4) + t70 * t13 + t50)) - m(7) * (g(1) * (t4 * rSges(7,1) + t3 * rSges(7,2) + pkin(5) * t67 + t8 * t29 + t58 * t7 + t51) + g(2) * (t2 * rSges(7,1) + t1 * rSges(7,2) + pkin(5) * t68 + t6 * t29 + t58 * t5 + t49) + g(3) * (t10 * rSges(7,1) + t9 * rSges(7,2) - pkin(5) * t55 + t58 * t13 + t14 * t29 + t50));
U  = t11;
