% Calculate potential energy for
% S6PRPRRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 20:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRRR7_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRPRRR7_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR7_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_energypot_floatb_twist_slag_vp1: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR7_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRR7_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:50:28
% EndTime: 2019-03-08 20:50:29
% DurationCPUTime: 0.62s
% Computational Cost: add. (972->147), mult. (2565->202), div. (0->0), fcn. (3324->18), ass. (0->74)
t59 = sin(pkin(7));
t60 = sin(pkin(6));
t100 = t59 * t60;
t64 = cos(pkin(7));
t65 = cos(pkin(6));
t72 = cos(qJ(2));
t44 = -t72 * t100 + t65 * t64;
t57 = sin(pkin(13));
t62 = cos(pkin(13));
t69 = sin(qJ(2));
t92 = t65 * t72;
t47 = -t57 * t92 - t62 * t69;
t97 = t60 * t64;
t39 = -t47 * t59 + t57 * t97;
t56 = sin(pkin(14));
t61 = cos(pkin(14));
t95 = t64 * t72;
t99 = t59 * t65;
t36 = t61 * t99 + (-t56 * t69 + t61 * t95) * t60;
t58 = sin(pkin(8));
t63 = cos(pkin(8));
t25 = -t36 * t58 + t44 * t63;
t93 = t65 * t69;
t48 = -t57 * t93 + t62 * t72;
t82 = t57 * t100 + t47 * t64;
t28 = -t48 * t56 + t82 * t61;
t19 = -t28 * t58 + t39 * t63;
t46 = t57 * t72 + t62 * t93;
t45 = -t57 * t69 + t62 * t92;
t98 = t60 * t62;
t83 = t45 * t64 - t59 * t98;
t26 = -t46 * t56 + t83 * t61;
t38 = -t45 * t59 - t62 * t97;
t18 = -t26 * t58 + t38 * t63;
t110 = rSges(6,3) + pkin(11);
t109 = pkin(12) + rSges(7,3);
t108 = cos(qJ(4));
t96 = t60 * t69;
t91 = t57 * pkin(1) + r_base(2);
t88 = qJ(1) + r_base(3);
t87 = t58 * t108;
t86 = t63 * t108;
t85 = t57 * t60 * pkin(9) + t62 * pkin(1) + r_base(1);
t84 = t65 * pkin(9) + t88;
t81 = t48 * pkin(2) + t39 * qJ(3) + t85;
t80 = pkin(2) * t96 + t44 * qJ(3) + t84;
t29 = t48 * t61 + t82 * t56;
t79 = t29 * pkin(3) + t19 * pkin(10) + t81;
t68 = sin(qJ(4));
t12 = t29 * t108 + (t28 * t63 + t39 * t58) * t68;
t78 = t12 * pkin(4) + t79;
t77 = t46 * pkin(2) - pkin(9) * t98 + t38 * qJ(3) + t91;
t37 = t61 * t96 + (t60 * t95 + t99) * t56;
t76 = t37 * pkin(3) + t25 * pkin(10) + t80;
t17 = t37 * t108 + (t36 * t63 + t44 * t58) * t68;
t75 = t17 * pkin(4) + t76;
t27 = t46 * t61 + t83 * t56;
t74 = t27 * pkin(3) + t18 * pkin(10) + t77;
t10 = t27 * t108 + (t26 * t63 + t38 * t58) * t68;
t73 = t10 * pkin(4) + t74;
t71 = cos(qJ(5));
t70 = cos(qJ(6));
t67 = sin(qJ(5));
t66 = sin(qJ(6));
t16 = -t36 * t86 + t37 * t68 - t44 * t87;
t11 = -t28 * t86 + t29 * t68 - t39 * t87;
t9 = -t26 * t86 + t27 * t68 - t38 * t87;
t6 = t17 * t71 + t25 * t67;
t5 = t17 * t67 - t25 * t71;
t4 = t12 * t71 + t19 * t67;
t3 = t12 * t67 - t19 * t71;
t2 = t10 * t71 + t18 * t67;
t1 = t10 * t67 - t18 * t71;
t7 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t62 - rSges(2,2) * t57 + r_base(1)) + g(2) * (rSges(2,1) * t57 + rSges(2,2) * t62 + r_base(2)) + g(3) * (rSges(2,3) + t88)) - m(3) * (g(1) * (rSges(3,1) * t48 + rSges(3,2) * t47 + t85) + g(2) * (rSges(3,1) * t46 + rSges(3,2) * t45 + t91) + g(3) * (t65 * rSges(3,3) + t84) + (g(1) * rSges(3,3) * t57 + g(3) * (rSges(3,1) * t69 + rSges(3,2) * t72) + g(2) * (-rSges(3,3) - pkin(9)) * t62) * t60) - m(4) * (g(1) * (rSges(4,1) * t29 + rSges(4,2) * t28 + rSges(4,3) * t39 + t81) + g(2) * (rSges(4,1) * t27 + rSges(4,2) * t26 + rSges(4,3) * t38 + t77) + g(3) * (t37 * rSges(4,1) + t36 * rSges(4,2) + t44 * rSges(4,3) + t80)) - m(5) * (g(1) * (rSges(5,1) * t12 - rSges(5,2) * t11 + rSges(5,3) * t19 + t79) + g(2) * (rSges(5,1) * t10 - rSges(5,2) * t9 + rSges(5,3) * t18 + t74) + g(3) * (t17 * rSges(5,1) - t16 * rSges(5,2) + t25 * rSges(5,3) + t76)) - m(6) * (g(1) * (rSges(6,1) * t4 - rSges(6,2) * t3 + t110 * t11 + t78) + g(2) * (rSges(6,1) * t2 - rSges(6,2) * t1 + t110 * t9 + t73) + g(3) * (t6 * rSges(6,1) - t5 * rSges(6,2) + t110 * t16 + t75)) - m(7) * (g(1) * ((t11 * t66 + t4 * t70) * rSges(7,1) + (t11 * t70 - t4 * t66) * rSges(7,2) + t109 * t3 + t11 * pkin(11) + t4 * pkin(5) + t78) + g(2) * (t2 * pkin(5) + t9 * pkin(11) + (t2 * t70 + t66 * t9) * rSges(7,1) + (-t2 * t66 + t70 * t9) * rSges(7,2) + t109 * t1 + t73) + g(3) * (t109 * t5 + (t16 * t70 - t6 * t66) * rSges(7,2) + (t16 * t66 + t6 * t70) * rSges(7,1) + t16 * pkin(11) + t6 * pkin(5) + t75));
U  = t7;
