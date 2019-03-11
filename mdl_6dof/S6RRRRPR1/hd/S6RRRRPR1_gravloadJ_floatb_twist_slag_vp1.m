% Calculate Gravitation load on the joints for
% S6RRRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPR1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR1_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:52:56
% EndTime: 2019-03-09 21:52:57
% DurationCPUTime: 0.67s
% Computational Cost: add. (631->127), mult. (455->155), div. (0->0), fcn. (391->12), ass. (0->75)
t104 = rSges(7,3) + pkin(10);
t43 = qJ(2) + qJ(3);
t40 = qJ(4) + t43;
t33 = pkin(11) + t40;
t28 = sin(t33);
t29 = cos(t33);
t44 = sin(qJ(6));
t88 = rSges(7,2) * t44;
t103 = t29 * t104 + t28 * t88;
t66 = t29 * rSges(6,1) - rSges(6,2) * t28;
t34 = sin(t40);
t35 = cos(t40);
t73 = t35 * rSges(5,1) - rSges(5,2) * t34;
t37 = sin(t43);
t38 = cos(t43);
t74 = t38 * rSges(4,1) - rSges(4,2) * t37;
t67 = -rSges(5,1) * t34 - rSges(5,2) * t35;
t97 = pkin(3) * t37;
t101 = t67 - t97;
t46 = sin(qJ(1));
t49 = cos(qJ(1));
t100 = g(1) * t49 + g(2) * t46;
t57 = t29 * pkin(5) + t104 * t28;
t47 = cos(qJ(6));
t92 = rSges(7,1) * t47;
t99 = (-pkin(5) - t92) * t28;
t50 = -pkin(8) - pkin(7);
t45 = sin(qJ(2));
t98 = pkin(2) * t45;
t96 = pkin(4) * t34;
t93 = rSges(3,3) + pkin(7);
t48 = cos(qJ(2));
t41 = t48 * pkin(2);
t36 = t41 + pkin(1);
t85 = t44 * t46;
t84 = t44 * t49;
t83 = t46 * t47;
t82 = t47 * t49;
t81 = rSges(4,3) - t50;
t42 = -pkin(9) + t50;
t80 = rSges(5,3) - t42;
t39 = -qJ(5) + t42;
t79 = rSges(6,3) - t39;
t32 = pkin(3) * t38;
t15 = t32 + t36;
t77 = t103 * t46;
t76 = t103 * t49;
t72 = t32 + t73;
t30 = pkin(4) * t35;
t71 = t30 + t66;
t14 = -t96 - t97;
t70 = rSges(3,1) * t48 - rSges(3,2) * t45;
t68 = -rSges(4,1) * t37 - rSges(4,2) * t38;
t65 = -rSges(6,1) * t28 - rSges(6,2) * t29;
t64 = t32 + t71;
t63 = pkin(1) + t70;
t62 = t36 + t74;
t61 = t15 + t73;
t59 = t65 * t46;
t58 = t65 * t49;
t54 = t30 + t57 + (-t88 + t92) * t29;
t53 = t32 + t54;
t52 = t100 * t99;
t11 = t49 * t14;
t10 = t46 * t14;
t9 = t14 - t98;
t8 = t30 + t15;
t7 = t29 * t82 + t85;
t6 = -t29 * t84 + t83;
t5 = -t29 * t83 + t84;
t4 = t29 * t85 + t82;
t3 = t49 * t9;
t2 = t46 * t9;
t1 = t49 * t8;
t12 = [-m(2) * (g(1) * (-rSges(2,1) * t46 - rSges(2,2) * t49) + g(2) * (rSges(2,1) * t49 - rSges(2,2) * t46)) - m(3) * ((g(1) * t93 + g(2) * t63) * t49 + (-g(1) * t63 + g(2) * t93) * t46) - m(4) * ((g(1) * t81 + g(2) * t62) * t49 + (-g(1) * t62 + g(2) * t81) * t46) - m(5) * ((g(1) * t80 + g(2) * t61) * t49 + (-g(1) * t61 + g(2) * t80) * t46) - m(6) * (g(2) * t1 + (g(1) * t79 + g(2) * t66) * t49 + (g(1) * (-t66 - t8) + g(2) * t79) * t46) - m(7) * (g(1) * (rSges(7,1) * t5 + rSges(7,2) * t4) + g(2) * (rSges(7,1) * t7 + rSges(7,2) * t6 + t1) + (-g(1) * t39 + g(2) * t57) * t49 + (g(1) * (-t57 - t8) - g(2) * t39) * t46) -m(3) * (g(3) * t70 + t100 * (-rSges(3,1) * t45 - rSges(3,2) * t48)) - m(4) * (g(3) * (t41 + t74) + t100 * (t68 - t98)) - m(5) * (g(3) * (t41 + t72) + t100 * (-t98 + t101)) - m(6) * (g(1) * (t3 + t58) + g(2) * (t2 + t59) + g(3) * (t41 + t64)) - m(7) * (g(1) * (t3 + t76) + g(2) * (t2 + t77) + g(3) * (t41 + t53) + t52) -m(4) * (g(3) * t74 + t100 * t68) - m(5) * (g(3) * t72 + t100 * t101) - m(6) * (g(1) * (t11 + t58) + g(2) * (t10 + t59) + g(3) * t64) - m(7) * (g(1) * (t11 + t76) + g(2) * (t10 + t77) + g(3) * t53 + t52) -m(7) * (g(1) * t76 + g(2) * t77) + (-m(5) * t73 - m(6) * t71 - m(7) * t54) * g(3) + t100 * (-m(5) * t67 - m(6) * (t65 - t96) - m(7) * (-t96 + t99)) (-m(6) - m(7)) * (g(1) * t46 - g(2) * t49) -m(7) * (g(1) * (rSges(7,1) * t6 - rSges(7,2) * t7) + g(2) * (-rSges(7,1) * t4 + rSges(7,2) * t5) + g(3) * (-rSges(7,1) * t44 - rSges(7,2) * t47) * t28)];
taug  = t12(:);
