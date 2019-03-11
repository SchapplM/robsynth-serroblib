% Calculate potential energy for
% S6PRRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 21:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRPRR1_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR1_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRRPRR1_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR1_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR1_energypot_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR1_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRR1_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:50:55
% EndTime: 2019-03-08 21:50:56
% DurationCPUTime: 0.66s
% Computational Cost: add. (378->133), mult. (471->163), div. (0->0), fcn. (534->14), ass. (0->55)
t41 = -qJ(4) - pkin(8);
t70 = -t41 + rSges(5,3);
t43 = sin(qJ(3));
t69 = pkin(3) * t43;
t37 = sin(pkin(11));
t68 = g(1) * t37;
t39 = cos(pkin(11));
t67 = g(2) * t39;
t66 = pkin(8) + rSges(4,3);
t65 = pkin(10) + rSges(7,3);
t46 = cos(qJ(3));
t27 = pkin(3) * t46 + pkin(2);
t38 = sin(pkin(6));
t64 = t37 * t38;
t63 = t38 * t39;
t62 = t38 * t43;
t44 = sin(qJ(2));
t61 = t38 * t44;
t60 = t38 * t46;
t47 = cos(qJ(2));
t59 = t38 * t47;
t40 = cos(pkin(6));
t58 = t40 * t44;
t57 = t40 * t47;
t36 = qJ(3) + pkin(12);
t56 = pkin(1) * t37 + r_base(2);
t55 = qJ(1) + r_base(3);
t54 = pkin(1) * t39 + pkin(7) * t64 + r_base(1);
t53 = pkin(7) * t40 + t55;
t28 = sin(t36);
t29 = cos(t36);
t52 = rSges(5,1) * t29 - rSges(5,2) * t28 + t27;
t15 = t37 * t57 + t39 * t44;
t16 = -t37 * t58 + t39 * t47;
t19 = pkin(4) * t29 + t27;
t20 = pkin(4) * t28 + t69;
t35 = -pkin(9) + t41;
t51 = -t15 * t35 + t16 * t19 + t20 * t64 + t54;
t50 = t19 * t61 + t20 * t40 + t35 * t59 + t53;
t49 = rSges(5,1) * t28 + rSges(5,2) * t29 + t69;
t13 = t37 * t44 - t39 * t57;
t14 = t37 * t47 + t39 * t58;
t48 = t14 * t19 + (-pkin(7) - t20) * t63 - t13 * t35 + t56;
t45 = cos(qJ(6));
t42 = sin(qJ(6));
t30 = qJ(5) + t36;
t26 = cos(t30);
t25 = sin(t30);
t8 = t25 * t40 + t26 * t61;
t7 = t25 * t61 - t26 * t40;
t4 = t16 * t26 + t25 * t64;
t3 = t16 * t25 - t26 * t64;
t2 = t14 * t26 - t25 * t63;
t1 = t14 * t25 + t26 * t63;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t39 - rSges(2,2) * t37 + r_base(1)) + g(2) * (rSges(2,1) * t37 + rSges(2,2) * t39 + r_base(2)) + g(3) * (rSges(2,3) + t55)) - m(3) * (g(1) * (rSges(3,1) * t16 - rSges(3,2) * t15 + t54) + g(2) * (rSges(3,1) * t14 - rSges(3,2) * t13 + t56) + g(3) * (t40 * rSges(3,3) + t53) + (rSges(3,3) * t68 + g(3) * (rSges(3,1) * t44 + rSges(3,2) * t47) + (-rSges(3,3) - pkin(7)) * t67) * t38) - m(4) * (g(1) * (t16 * pkin(2) + (t16 * t46 + t37 * t62) * rSges(4,1) + (-t16 * t43 + t37 * t60) * rSges(4,2) + t66 * t15 + t54) + g(2) * (t14 * pkin(2) - pkin(7) * t63 + (t14 * t46 - t39 * t62) * rSges(4,1) + (-t14 * t43 - t39 * t60) * rSges(4,2) + t66 * t13 + t56) + g(3) * ((rSges(4,1) * t43 + rSges(4,2) * t46) * t40 + (-t66 * t47 + (rSges(4,1) * t46 - rSges(4,2) * t43 + pkin(2)) * t44) * t38 + t53)) - m(5) * ((t49 * t68 + (-pkin(7) - t49) * t67) * t38 + (t53 + t49 * t40 + (t52 * t44 - t47 * t70) * t38) * g(3) + (t13 * t70 + t52 * t14 + t56) * g(2) + (t15 * t70 + t52 * t16 + t54) * g(1)) - m(6) * (g(1) * (rSges(6,1) * t4 - rSges(6,2) * t3 + rSges(6,3) * t15 + t51) + g(2) * (rSges(6,1) * t2 - rSges(6,2) * t1 + rSges(6,3) * t13 + t48) + g(3) * (rSges(6,1) * t8 - rSges(6,2) * t7 - rSges(6,3) * t59 + t50)) - m(7) * (g(1) * (t4 * pkin(5) + (t15 * t42 + t4 * t45) * rSges(7,1) + (t15 * t45 - t4 * t42) * rSges(7,2) + t65 * t3 + t51) + g(2) * (t2 * pkin(5) + (t13 * t42 + t2 * t45) * rSges(7,1) + (t13 * t45 - t2 * t42) * rSges(7,2) + t65 * t1 + t48) + g(3) * (t8 * pkin(5) + (-t42 * t59 + t45 * t8) * rSges(7,1) + (-t42 * t8 - t45 * t59) * rSges(7,2) + t65 * t7 + t50));
U  = t5;
