% Calculate potential energy for
% S6RRPRRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 14:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRR9_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR9_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPRRR9_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR9_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR9_energypot_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR9_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR9_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:08:13
% EndTime: 2019-03-09 14:08:13
% DurationCPUTime: 0.66s
% Computational Cost: add. (378->133), mult. (471->159), div. (0->0), fcn. (534->14), ass. (0->55)
t41 = -pkin(9) - qJ(3);
t70 = -t41 + rSges(5,3);
t37 = sin(pkin(12));
t69 = pkin(3) * t37;
t44 = sin(qJ(1));
t68 = g(1) * t44;
t47 = cos(qJ(1));
t67 = g(2) * t47;
t66 = pkin(11) + rSges(7,3);
t39 = cos(pkin(12));
t27 = pkin(3) * t39 + pkin(2);
t38 = sin(pkin(6));
t43 = sin(qJ(2));
t65 = t38 * t43;
t64 = t38 * t44;
t46 = cos(qJ(2));
t63 = t38 * t46;
t62 = t38 * t47;
t61 = t43 * t44;
t60 = t43 * t47;
t59 = t44 * t46;
t58 = t46 * t47;
t57 = qJ(3) + rSges(4,3);
t56 = pkin(7) + r_base(3);
t36 = pkin(12) + qJ(4);
t55 = pkin(1) * t44 + r_base(2);
t40 = cos(pkin(6));
t54 = pkin(8) * t40 + t56;
t53 = pkin(1) * t47 + pkin(8) * t64 + r_base(1);
t28 = sin(t36);
t29 = cos(t36);
t52 = rSges(5,1) * t29 - rSges(5,2) * t28 + t27;
t18 = pkin(4) * t29 + t27;
t20 = pkin(4) * t28 + t69;
t35 = -pkin(10) + t41;
t51 = t18 * t65 + t20 * t40 + t35 * t63 + t54;
t15 = t40 * t59 + t60;
t16 = -t40 * t61 + t58;
t50 = -t15 * t35 + t16 * t18 + t20 * t64 + t53;
t49 = rSges(5,1) * t28 + rSges(5,2) * t29 + t69;
t13 = -t40 * t58 + t61;
t14 = t40 * t60 + t59;
t48 = t14 * t18 + (-pkin(8) - t20) * t62 - t13 * t35 + t55;
t45 = cos(qJ(6));
t42 = sin(qJ(6));
t30 = qJ(5) + t36;
t26 = cos(t30);
t25 = sin(t30);
t8 = t25 * t40 + t26 * t65;
t7 = t25 * t65 - t26 * t40;
t4 = t16 * t26 + t25 * t64;
t3 = t16 * t25 - t26 * t64;
t2 = t14 * t26 - t25 * t62;
t1 = t14 * t25 + t26 * t62;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t47 - rSges(2,2) * t44 + r_base(1)) + g(2) * (rSges(2,1) * t44 + rSges(2,2) * t47 + r_base(2)) + g(3) * (rSges(2,3) + t56)) - m(3) * (g(1) * (rSges(3,1) * t16 - rSges(3,2) * t15 + t53) + g(2) * (t14 * rSges(3,1) - t13 * rSges(3,2) + t55) + g(3) * (rSges(3,3) * t40 + t54) + (rSges(3,3) * t68 + g(3) * (rSges(3,1) * t43 + rSges(3,2) * t46) + (-rSges(3,3) - pkin(8)) * t67) * t38) - m(4) * (g(1) * (t16 * pkin(2) + (t16 * t39 + t37 * t64) * rSges(4,1) + (-t16 * t37 + t39 * t64) * rSges(4,2) + t57 * t15 + t53) + g(2) * (t14 * pkin(2) - pkin(8) * t62 + (t14 * t39 - t37 * t62) * rSges(4,1) + (-t14 * t37 - t39 * t62) * rSges(4,2) + t57 * t13 + t55) + g(3) * ((rSges(4,1) * t37 + rSges(4,2) * t39) * t40 + (-t57 * t46 + (rSges(4,1) * t39 - rSges(4,2) * t37 + pkin(2)) * t43) * t38 + t54)) - m(5) * ((t49 * t68 + (-pkin(8) - t49) * t67) * t38 + (t54 + t49 * t40 + (t52 * t43 - t46 * t70) * t38) * g(3) + (t13 * t70 + t52 * t14 + t55) * g(2) + (t15 * t70 + t52 * t16 + t53) * g(1)) - m(6) * (g(1) * (rSges(6,1) * t4 - rSges(6,2) * t3 + rSges(6,3) * t15 + t50) + g(2) * (t2 * rSges(6,1) - t1 * rSges(6,2) + t13 * rSges(6,3) + t48) + g(3) * (rSges(6,1) * t8 - rSges(6,2) * t7 - rSges(6,3) * t63 + t51)) - m(7) * (g(1) * (t4 * pkin(5) + (t15 * t42 + t4 * t45) * rSges(7,1) + (t15 * t45 - t4 * t42) * rSges(7,2) + t66 * t3 + t50) + g(2) * (t2 * pkin(5) + (t13 * t42 + t2 * t45) * rSges(7,1) + (t13 * t45 - t2 * t42) * rSges(7,2) + t66 * t1 + t48) + g(3) * (t8 * pkin(5) + (-t42 * t63 + t45 * t8) * rSges(7,1) + (-t42 * t8 - t45 * t63) * rSges(7,2) + t66 * t7 + t51));
U  = t5;
