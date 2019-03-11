% Calculate potential energy for
% S6PRPRRP4
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
% Datum: 2019-03-08 20:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRRP4_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP4_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRPRRP4_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP4_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP4_energypot_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP4_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRP4_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:09:08
% EndTime: 2019-03-08 20:09:08
% DurationCPUTime: 0.38s
% Computational Cost: add. (391->126), mult. (591->153), div. (0->0), fcn. (697->12), ass. (0->55)
t41 = sin(pkin(11));
t74 = pkin(3) * t41;
t73 = rSges(7,1) + pkin(5);
t72 = rSges(7,2) + pkin(9);
t71 = rSges(6,3) + pkin(9);
t42 = sin(pkin(10));
t43 = sin(pkin(6));
t70 = t42 * t43;
t45 = cos(pkin(10));
t69 = t43 * t45;
t49 = sin(qJ(2));
t68 = t43 * t49;
t51 = cos(qJ(2));
t67 = t43 * t51;
t46 = cos(pkin(6));
t66 = t46 * t49;
t65 = t46 * t51;
t64 = rSges(7,3) + qJ(6);
t63 = qJ(3) + rSges(4,3);
t62 = t42 * pkin(1) + r_base(2);
t61 = t41 * t70;
t60 = qJ(1) + r_base(3);
t59 = t45 * pkin(1) + pkin(7) * t70 + r_base(1);
t58 = t46 * pkin(7) + t60;
t24 = t42 * t65 + t45 * t49;
t25 = -t42 * t66 + t45 * t51;
t44 = cos(pkin(11));
t34 = t44 * pkin(3) + pkin(2);
t47 = -pkin(8) - qJ(3);
t57 = pkin(3) * t61 - t24 * t47 + t25 * t34 + t59;
t56 = t34 * t68 + t46 * t74 + t47 * t67 + t58;
t40 = pkin(11) + qJ(4);
t35 = sin(t40);
t36 = cos(t40);
t10 = t25 * t36 + t35 * t70;
t55 = t10 * pkin(4) + t57;
t17 = t46 * t35 + t36 * t68;
t54 = t17 * pkin(4) + t56;
t22 = t42 * t49 - t45 * t65;
t23 = t42 * t51 + t45 * t66;
t53 = t23 * t34 + (-pkin(7) - t74) * t69 - t22 * t47 + t62;
t8 = t23 * t36 - t35 * t69;
t52 = t8 * pkin(4) + t53;
t50 = cos(qJ(5));
t48 = sin(qJ(5));
t16 = t35 * t68 - t46 * t36;
t12 = t17 * t50 - t48 * t67;
t11 = t17 * t48 + t50 * t67;
t9 = t25 * t35 - t36 * t70;
t7 = t23 * t35 + t36 * t69;
t4 = t10 * t50 + t24 * t48;
t3 = t10 * t48 - t24 * t50;
t2 = t22 * t48 + t8 * t50;
t1 = -t22 * t50 + t8 * t48;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t45 * rSges(2,1) - t42 * rSges(2,2) + r_base(1)) + g(2) * (t42 * rSges(2,1) + t45 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t60)) - m(3) * (g(1) * (t25 * rSges(3,1) - t24 * rSges(3,2) + t59) + g(2) * (t23 * rSges(3,1) - t22 * rSges(3,2) + t62) + g(3) * (t46 * rSges(3,3) + t58) + (g(1) * rSges(3,3) * t42 + g(3) * (rSges(3,1) * t49 + rSges(3,2) * t51) + g(2) * (-rSges(3,3) - pkin(7)) * t45) * t43) - m(4) * (g(1) * (t25 * pkin(2) + (t25 * t44 + t61) * rSges(4,1) + (-t25 * t41 + t44 * t70) * rSges(4,2) + t63 * t24 + t59) + g(2) * (t23 * pkin(2) - pkin(7) * t69 + (t23 * t44 - t41 * t69) * rSges(4,1) + (-t23 * t41 - t44 * t69) * rSges(4,2) + t63 * t22 + t62) + g(3) * ((t41 * rSges(4,1) + t44 * rSges(4,2)) * t46 + (-t63 * t51 + (t44 * rSges(4,1) - t41 * rSges(4,2) + pkin(2)) * t49) * t43 + t58)) - m(5) * (g(1) * (t10 * rSges(5,1) - t9 * rSges(5,2) + t24 * rSges(5,3) + t57) + g(2) * (t8 * rSges(5,1) - t7 * rSges(5,2) + t22 * rSges(5,3) + t53) + g(3) * (t17 * rSges(5,1) - t16 * rSges(5,2) - rSges(5,3) * t67 + t56)) - m(6) * (g(1) * (t4 * rSges(6,1) - t3 * rSges(6,2) + t71 * t9 + t55) + g(2) * (t2 * rSges(6,1) - t1 * rSges(6,2) + t71 * t7 + t52) + g(3) * (t12 * rSges(6,1) - t11 * rSges(6,2) + t71 * t16 + t54)) - m(7) * (g(1) * (t64 * t3 + t73 * t4 + t72 * t9 + t55) + g(2) * (t64 * t1 + t73 * t2 + t72 * t7 + t52) + g(3) * (t64 * t11 + t73 * t12 + t72 * t16 + t54));
U  = t5;
