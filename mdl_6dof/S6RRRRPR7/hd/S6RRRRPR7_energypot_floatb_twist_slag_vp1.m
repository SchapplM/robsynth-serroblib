% Calculate potential energy for
% S6RRRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 22:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRPR7_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR7_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRRPR7_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR7_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR7_energypot_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR7_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR7_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:26:24
% EndTime: 2019-03-09 22:26:25
% DurationCPUTime: 0.66s
% Computational Cost: add. (378->133), mult. (471->159), div. (0->0), fcn. (534->14), ass. (0->55)
t47 = -pkin(10) - pkin(9);
t70 = -t47 + rSges(5,3);
t40 = sin(qJ(3));
t69 = pkin(3) * t40;
t42 = sin(qJ(1));
t68 = g(1) * t42;
t46 = cos(qJ(1));
t67 = g(2) * t46;
t66 = pkin(9) + rSges(4,3);
t65 = pkin(11) + rSges(7,3);
t44 = cos(qJ(3));
t27 = t44 * pkin(3) + pkin(2);
t37 = sin(pkin(6));
t41 = sin(qJ(2));
t64 = t37 * t41;
t63 = t37 * t42;
t45 = cos(qJ(2));
t62 = t37 * t45;
t61 = t37 * t46;
t60 = t41 * t42;
t59 = t41 * t46;
t58 = t42 * t45;
t57 = t45 * t46;
t36 = qJ(3) + qJ(4);
t56 = pkin(7) + r_base(3);
t55 = t42 * pkin(1) + r_base(2);
t38 = cos(pkin(6));
t54 = t38 * pkin(8) + t56;
t53 = t46 * pkin(1) + pkin(8) * t63 + r_base(1);
t29 = sin(t36);
t30 = cos(t36);
t52 = rSges(5,1) * t30 - rSges(5,2) * t29 + t27;
t19 = pkin(4) * t30 + t27;
t20 = pkin(4) * t29 + t69;
t35 = -qJ(5) + t47;
t51 = t19 * t64 + t38 * t20 + t35 * t62 + t54;
t15 = t38 * t58 + t59;
t16 = -t38 * t60 + t57;
t50 = -t15 * t35 + t16 * t19 + t20 * t63 + t53;
t49 = rSges(5,1) * t29 + rSges(5,2) * t30 + t69;
t13 = -t38 * t57 + t60;
t14 = t38 * t59 + t58;
t48 = t14 * t19 + (-pkin(8) - t20) * t61 - t13 * t35 + t55;
t43 = cos(qJ(6));
t39 = sin(qJ(6));
t28 = pkin(12) + t36;
t26 = cos(t28);
t25 = sin(t28);
t8 = t25 * t38 + t26 * t64;
t7 = t25 * t64 - t38 * t26;
t4 = t16 * t26 + t25 * t63;
t3 = t16 * t25 - t26 * t63;
t2 = t14 * t26 - t25 * t61;
t1 = t14 * t25 + t26 * t61;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t46 - rSges(2,2) * t42 + r_base(1)) + g(2) * (rSges(2,1) * t42 + rSges(2,2) * t46 + r_base(2)) + g(3) * (rSges(2,3) + t56)) - m(3) * (g(1) * (rSges(3,1) * t16 - rSges(3,2) * t15 + t53) + g(2) * (rSges(3,1) * t14 - rSges(3,2) * t13 + t55) + g(3) * (rSges(3,3) * t38 + t54) + (rSges(3,3) * t68 + g(3) * (rSges(3,1) * t41 + rSges(3,2) * t45) + (-rSges(3,3) - pkin(8)) * t67) * t37) - m(4) * (g(1) * (t16 * pkin(2) + (t16 * t44 + t40 * t63) * rSges(4,1) + (-t16 * t40 + t44 * t63) * rSges(4,2) + t66 * t15 + t53) + g(2) * (t14 * pkin(2) - pkin(8) * t61 + (t14 * t44 - t40 * t61) * rSges(4,1) + (-t14 * t40 - t44 * t61) * rSges(4,2) + t66 * t13 + t55) + g(3) * ((rSges(4,1) * t40 + rSges(4,2) * t44) * t38 + (-t66 * t45 + (rSges(4,1) * t44 - rSges(4,2) * t40 + pkin(2)) * t41) * t37 + t54)) - m(5) * ((t49 * t68 + (-pkin(8) - t49) * t67) * t37 + (t54 + t49 * t38 + (t52 * t41 - t70 * t45) * t37) * g(3) + (t70 * t13 + t52 * t14 + t55) * g(2) + (t70 * t15 + t52 * t16 + t53) * g(1)) - m(6) * (g(1) * (rSges(6,1) * t4 - rSges(6,2) * t3 + rSges(6,3) * t15 + t50) + g(2) * (rSges(6,1) * t2 - rSges(6,2) * t1 + rSges(6,3) * t13 + t48) + g(3) * (rSges(6,1) * t8 - rSges(6,2) * t7 - rSges(6,3) * t62 + t51)) - m(7) * (g(1) * (t4 * pkin(5) + (t15 * t39 + t4 * t43) * rSges(7,1) + (t15 * t43 - t39 * t4) * rSges(7,2) + t65 * t3 + t50) + g(2) * (t2 * pkin(5) + (t13 * t39 + t2 * t43) * rSges(7,1) + (t13 * t43 - t2 * t39) * rSges(7,2) + t65 * t1 + t48) + g(3) * (t8 * pkin(5) + (-t39 * t62 + t43 * t8) * rSges(7,1) + (-t39 * t8 - t43 * t62) * rSges(7,2) + t65 * t7 + t51));
U  = t5;
