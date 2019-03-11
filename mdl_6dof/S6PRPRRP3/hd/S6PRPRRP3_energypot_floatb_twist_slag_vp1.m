% Calculate potential energy for
% S6PRPRRP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-03-08 20:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRRP3_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRPRRP3_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP3_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_energypot_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP3_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRP3_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:04:35
% EndTime: 2019-03-08 20:04:36
% DurationCPUTime: 0.39s
% Computational Cost: add. (364->131), mult. (545->160), div. (0->0), fcn. (633->12), ass. (0->54)
t37 = sin(pkin(11));
t69 = pkin(3) * t37;
t68 = rSges(6,3) + pkin(9);
t38 = sin(pkin(10));
t41 = cos(pkin(10));
t46 = sin(qJ(2));
t42 = cos(pkin(6));
t48 = cos(qJ(2));
t60 = t42 * t48;
t17 = t38 * t46 - t41 * t60;
t45 = sin(qJ(5));
t67 = t17 * t45;
t19 = t38 * t60 + t41 * t46;
t66 = t19 * t45;
t39 = sin(pkin(6));
t65 = t38 * t39;
t64 = t39 * t41;
t63 = t39 * t46;
t62 = t39 * t48;
t61 = t42 * t46;
t59 = rSges(7,3) + qJ(6) + pkin(9);
t58 = qJ(3) + rSges(4,3);
t57 = t38 * pkin(1) + r_base(2);
t56 = t37 * t65;
t55 = t45 * t62;
t54 = qJ(1) + r_base(3);
t53 = t41 * pkin(1) + pkin(7) * t65 + r_base(1);
t52 = t42 * pkin(7) + t54;
t20 = -t38 * t61 + t41 * t48;
t40 = cos(pkin(11));
t29 = t40 * pkin(3) + pkin(2);
t44 = -pkin(8) - qJ(3);
t51 = pkin(3) * t56 - t19 * t44 + t20 * t29 + t53;
t50 = t29 * t63 + t42 * t69 + t44 * t62 + t52;
t18 = t38 * t48 + t41 * t61;
t49 = t18 * t29 + (-pkin(7) - t69) * t64 - t17 * t44 + t57;
t47 = cos(qJ(5));
t36 = pkin(11) + qJ(4);
t32 = cos(t36);
t31 = sin(t36);
t30 = t47 * pkin(5) + pkin(4);
t14 = t42 * t31 + t32 * t63;
t13 = t31 * t63 - t42 * t32;
t10 = t14 * t47 - t55;
t9 = -t14 * t45 - t47 * t62;
t8 = t20 * t32 + t31 * t65;
t7 = t20 * t31 - t32 * t65;
t6 = t18 * t32 - t31 * t64;
t5 = t18 * t31 + t32 * t64;
t4 = t8 * t47 + t66;
t3 = t19 * t47 - t8 * t45;
t2 = t6 * t47 + t67;
t1 = t17 * t47 - t6 * t45;
t11 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t41 * rSges(2,1) - t38 * rSges(2,2) + r_base(1)) + g(2) * (t38 * rSges(2,1) + t41 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t54)) - m(3) * (g(1) * (t20 * rSges(3,1) - t19 * rSges(3,2) + t53) + g(2) * (t18 * rSges(3,1) - t17 * rSges(3,2) + t57) + g(3) * (t42 * rSges(3,3) + t52) + (g(1) * rSges(3,3) * t38 + g(3) * (rSges(3,1) * t46 + rSges(3,2) * t48) + g(2) * (-rSges(3,3) - pkin(7)) * t41) * t39) - m(4) * (g(1) * (t20 * pkin(2) + (t20 * t40 + t56) * rSges(4,1) + (-t20 * t37 + t40 * t65) * rSges(4,2) + t58 * t19 + t53) + g(2) * (t18 * pkin(2) - pkin(7) * t64 + (t18 * t40 - t37 * t64) * rSges(4,1) + (-t18 * t37 - t40 * t64) * rSges(4,2) + t58 * t17 + t57) + g(3) * ((t37 * rSges(4,1) + t40 * rSges(4,2)) * t42 + (-t58 * t48 + (t40 * rSges(4,1) - t37 * rSges(4,2) + pkin(2)) * t46) * t39 + t52)) - m(5) * (g(1) * (t8 * rSges(5,1) - t7 * rSges(5,2) + t19 * rSges(5,3) + t51) + g(2) * (t6 * rSges(5,1) - t5 * rSges(5,2) + t17 * rSges(5,3) + t49) + g(3) * (t14 * rSges(5,1) - t13 * rSges(5,2) - rSges(5,3) * t62 + t50)) - m(6) * (g(1) * (t4 * rSges(6,1) + t3 * rSges(6,2) + t8 * pkin(4) + t68 * t7 + t51) + g(2) * (t2 * rSges(6,1) + t1 * rSges(6,2) + t6 * pkin(4) + t68 * t5 + t49) + g(3) * (t10 * rSges(6,1) + t9 * rSges(6,2) + t14 * pkin(4) + t68 * t13 + t50)) - m(7) * (g(1) * (t4 * rSges(7,1) + t3 * rSges(7,2) + pkin(5) * t66 + t8 * t30 + t59 * t7 + t51) + g(2) * (t2 * rSges(7,1) + t1 * rSges(7,2) + pkin(5) * t67 + t6 * t30 + t59 * t5 + t49) + g(3) * (t10 * rSges(7,1) + t9 * rSges(7,2) - pkin(5) * t55 + t59 * t13 + t14 * t30 + t50));
U  = t11;
