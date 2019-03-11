% Calculate potential energy for
% S6PPPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
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
% Datum: 2019-03-08 18:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PPPRRR1_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PPPRRR1_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPPRRR1_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_energypot_floatb_twist_slag_vp1: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPPRRR1_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PPPRRR1_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:39:02
% EndTime: 2019-03-08 18:39:03
% DurationCPUTime: 0.62s
% Computational Cost: add. (972->147), mult. (2565->201), div. (0->0), fcn. (3324->18), ass. (0->74)
t63 = cos(pkin(13));
t66 = cos(pkin(7));
t67 = cos(pkin(6));
t60 = sin(pkin(7));
t61 = sin(pkin(6));
t98 = t60 * t61;
t44 = -t63 * t98 + t66 * t67;
t57 = sin(pkin(13));
t64 = cos(pkin(12));
t58 = sin(pkin(12));
t99 = t58 * t67;
t47 = -t57 * t64 - t63 * t99;
t96 = t61 * t66;
t39 = -t47 * t60 + t58 * t96;
t56 = sin(pkin(14));
t62 = cos(pkin(14));
t95 = t63 * t66;
t97 = t60 * t67;
t36 = t62 * t97 + (-t56 * t57 + t62 * t95) * t61;
t59 = sin(pkin(8));
t65 = cos(pkin(8));
t29 = -t36 * t59 + t44 * t65;
t48 = -t57 * t99 + t63 * t64;
t82 = t47 * t66 + t58 * t98;
t27 = -t48 * t56 + t82 * t62;
t19 = -t27 * t59 + t39 * t65;
t94 = t64 * t67;
t46 = t57 * t94 + t58 * t63;
t45 = -t57 * t58 + t63 * t94;
t83 = t45 * t66 - t64 * t98;
t25 = -t46 * t56 + t83 * t62;
t38 = -t45 * t60 - t64 * t96;
t18 = -t25 * t59 + t38 * t65;
t110 = rSges(6,3) + pkin(10);
t109 = pkin(11) + rSges(7,3);
t108 = cos(qJ(4));
t100 = t57 * t61;
t92 = t61 * qJ(2);
t91 = t58 * pkin(1) + r_base(2);
t88 = qJ(1) + r_base(3);
t87 = t59 * t108;
t86 = t65 * t108;
t85 = t64 * pkin(1) + t58 * t92 + r_base(1);
t84 = t67 * qJ(2) + t88;
t81 = t48 * pkin(2) + t39 * qJ(3) + t85;
t80 = pkin(2) * t100 + t44 * qJ(3) + t84;
t28 = t48 * t62 + t82 * t56;
t79 = t28 * pkin(3) + t19 * pkin(9) + t81;
t70 = sin(qJ(4));
t12 = t28 * t108 + (t27 * t65 + t39 * t59) * t70;
t78 = t12 * pkin(4) + t79;
t37 = t62 * t100 + (t61 * t95 + t97) * t56;
t77 = t37 * pkin(3) + t29 * pkin(9) + t80;
t76 = t46 * pkin(2) + t38 * qJ(3) - t64 * t92 + t91;
t17 = t37 * t108 + (t36 * t65 + t44 * t59) * t70;
t75 = t17 * pkin(4) + t77;
t26 = t46 * t62 + t83 * t56;
t74 = t26 * pkin(3) + t18 * pkin(9) + t76;
t10 = t26 * t108 + (t25 * t65 + t38 * t59) * t70;
t73 = t10 * pkin(4) + t74;
t72 = cos(qJ(5));
t71 = cos(qJ(6));
t69 = sin(qJ(5));
t68 = sin(qJ(6));
t16 = -t36 * t86 + t37 * t70 - t44 * t87;
t11 = -t27 * t86 + t28 * t70 - t39 * t87;
t9 = -t25 * t86 + t26 * t70 - t38 * t87;
t8 = t17 * t72 + t29 * t69;
t7 = t17 * t69 - t29 * t72;
t4 = t12 * t72 + t19 * t69;
t3 = t12 * t69 - t19 * t72;
t2 = t10 * t72 + t18 * t69;
t1 = t10 * t69 - t18 * t72;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t64 - rSges(2,2) * t58 + r_base(1)) + g(2) * (rSges(2,1) * t58 + rSges(2,2) * t64 + r_base(2)) + g(3) * (rSges(2,3) + t88)) - m(3) * (g(1) * (rSges(3,1) * t48 + rSges(3,2) * t47 + t85) + g(2) * (rSges(3,1) * t46 + rSges(3,2) * t45 + t91) + g(3) * (rSges(3,3) * t67 + t84) + (g(1) * rSges(3,3) * t58 + g(3) * (rSges(3,1) * t57 + rSges(3,2) * t63) + g(2) * (-rSges(3,3) - qJ(2)) * t64) * t61) - m(4) * (g(1) * (rSges(4,1) * t28 + rSges(4,2) * t27 + rSges(4,3) * t39 + t81) + g(2) * (rSges(4,1) * t26 + rSges(4,2) * t25 + rSges(4,3) * t38 + t76) + g(3) * (rSges(4,1) * t37 + rSges(4,2) * t36 + rSges(4,3) * t44 + t80)) - m(5) * (g(1) * (rSges(5,1) * t12 - rSges(5,2) * t11 + rSges(5,3) * t19 + t79) + g(2) * (rSges(5,1) * t10 - rSges(5,2) * t9 + rSges(5,3) * t18 + t74) + g(3) * (rSges(5,1) * t17 - rSges(5,2) * t16 + rSges(5,3) * t29 + t77)) - m(6) * (g(1) * (rSges(6,1) * t4 - t3 * rSges(6,2) + t110 * t11 + t78) + g(2) * (rSges(6,1) * t2 - rSges(6,2) * t1 + t110 * t9 + t73) + g(3) * (rSges(6,1) * t8 - rSges(6,2) * t7 + t110 * t16 + t75)) - m(7) * (g(1) * ((t11 * t68 + t4 * t71) * rSges(7,1) + (t11 * t71 - t4 * t68) * rSges(7,2) + t109 * t3 + t11 * pkin(10) + t4 * pkin(5) + t78) + g(2) * (t2 * pkin(5) + t9 * pkin(10) + (t2 * t71 + t68 * t9) * rSges(7,1) + (-t2 * t68 + t71 * t9) * rSges(7,2) + t109 * t1 + t73) + g(3) * (t109 * t7 + (t16 * t68 + t71 * t8) * rSges(7,1) + (t16 * t71 - t68 * t8) * rSges(7,2) + t16 * pkin(10) + t8 * pkin(5) + t75));
U  = t5;
