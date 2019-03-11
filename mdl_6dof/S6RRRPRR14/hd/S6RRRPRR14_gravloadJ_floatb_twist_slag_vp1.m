% Calculate Gravitation load on the joints for
% S6RRRPRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 20:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRR14_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR14_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR14_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR14_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR14_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR14_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:12:24
% EndTime: 2019-03-09 20:12:27
% DurationCPUTime: 1.39s
% Computational Cost: add. (722->199), mult. (1605->279), div. (0->0), fcn. (1919->12), ass. (0->80)
t54 = sin(qJ(5));
t58 = cos(qJ(5));
t121 = -t54 * rSges(6,1) - t58 * rSges(6,2);
t59 = cos(qJ(3));
t55 = sin(qJ(3));
t86 = qJ(4) * t55;
t118 = -pkin(3) * t59 - t86;
t111 = pkin(4) + pkin(9);
t112 = t58 * rSges(6,1) - t54 * rSges(6,2) + t111;
t100 = pkin(10) + rSges(6,3);
t87 = pkin(11) + pkin(10) + rSges(7,3);
t56 = sin(qJ(2));
t57 = sin(qJ(1));
t60 = cos(qJ(2));
t84 = cos(pkin(6));
t99 = cos(qJ(1));
t75 = t84 * t99;
t33 = t56 * t57 - t60 * t75;
t78 = t57 * t84;
t35 = t56 * t99 + t60 * t78;
t117 = g(1) * t35 + g(2) * t33;
t34 = t56 * t75 + t57 * t60;
t53 = sin(pkin(6));
t80 = t53 * t99;
t17 = t34 * t59 - t55 * t80;
t36 = -t56 * t78 + t60 * t99;
t97 = t53 * t57;
t21 = t36 * t59 + t55 * t97;
t96 = t53 * t59;
t32 = t55 * t84 + t56 * t96;
t116 = g(1) * t21 + g(2) * t17 + g(3) * t32;
t114 = t121 * t55;
t103 = t54 * pkin(5);
t52 = qJ(5) + qJ(6);
t49 = sin(t52);
t50 = cos(t52);
t65 = t49 * rSges(7,1) + t50 * rSges(7,2) + t103;
t113 = t55 * t65 + t59 * t87;
t104 = g(3) * t53;
t102 = rSges(5,1) + pkin(9);
t101 = rSges(4,3) + pkin(9);
t98 = t53 * t56;
t95 = t53 * t60;
t93 = t54 * t60;
t91 = t58 * t60;
t90 = t59 * t60;
t89 = t99 * pkin(1) + pkin(8) * t97;
t88 = pkin(2) * t95 + pkin(9) * t98;
t85 = rSges(5,3) + qJ(4);
t27 = t33 * pkin(2);
t83 = t118 * t33 - t27;
t29 = t35 * pkin(2);
t82 = t118 * t35 - t29;
t81 = t36 * pkin(2) + t89;
t79 = -t57 * pkin(1) + pkin(8) * t80;
t77 = t21 * pkin(3) + t81;
t76 = -t34 * pkin(2) + t79;
t74 = rSges(4,1) * t59 - rSges(4,2) * t55;
t73 = rSges(5,2) * t59 - rSges(5,3) * t55;
t16 = t34 * t55 + t59 * t80;
t72 = t16 * t58 - t33 * t54;
t20 = t36 * t55 - t57 * t96;
t8 = t20 * t58 - t35 * t54;
t71 = g(3) * (t53 * pkin(3) * t90 + t86 * t95 + t88);
t70 = -pkin(3) * t17 + t76;
t31 = t55 * t98 - t59 * t84;
t69 = t31 * t58 + t53 * t93;
t48 = pkin(5) * t58 + pkin(4);
t68 = t50 * rSges(7,1) - t49 * rSges(7,2) + t48;
t67 = qJ(4) - t121;
t66 = pkin(9) + t68;
t64 = qJ(4) + t65;
t6 = t20 * t50 - t35 * t49;
t7 = t20 * t49 + t35 * t50;
t63 = m(7) * (g(1) * (t6 * rSges(7,1) - t7 * rSges(7,2)) + g(2) * ((t16 * t50 - t33 * t49) * rSges(7,1) + (-t16 * t49 - t33 * t50) * rSges(7,2)) + g(3) * ((t31 * t50 + t49 * t95) * rSges(7,1) + (-t31 * t49 + t50 * t95) * rSges(7,2)));
t26 = t31 * pkin(3);
t14 = t20 * pkin(3);
t12 = t16 * pkin(3);
t9 = t20 * t54 + t35 * t58;
t1 = [-m(2) * (g(1) * (-t57 * rSges(2,1) - rSges(2,2) * t99) + g(2) * (rSges(2,1) * t99 - t57 * rSges(2,2))) - m(3) * (g(1) * (-t34 * rSges(3,1) + t33 * rSges(3,2) + rSges(3,3) * t80 + t79) + g(2) * (rSges(3,1) * t36 - rSges(3,2) * t35 + rSges(3,3) * t97 + t89)) - m(4) * (g(1) * (-rSges(4,1) * t17 + rSges(4,2) * t16 - t101 * t33 + t76) + g(2) * (rSges(4,1) * t21 - rSges(4,2) * t20 + t101 * t35 + t81)) - m(5) * (g(1) * (rSges(5,2) * t17 - t102 * t33 - t16 * t85 + t70) + g(2) * (-rSges(5,2) * t21 + t102 * t35 + t20 * t85 + t77)) - m(6) * (g(1) * (-t100 * t17 - t112 * t33 - t67 * t16 + t70) + g(2) * (rSges(6,1) * t9 + rSges(6,2) * t8 + t20 * qJ(4) + t100 * t21 + t111 * t35 + t77)) - m(7) * (g(1) * (-t16 * t64 - t17 * t87 - t33 * t66 + t70) + g(2) * (t7 * rSges(7,1) + t6 * rSges(7,2) + (pkin(9) + t48) * t35 + t87 * t21 + (qJ(4) + t103) * t20 + t77)) -m(3) * (g(1) * (-rSges(3,1) * t35 - rSges(3,2) * t36) + g(2) * (-rSges(3,1) * t33 - rSges(3,2) * t34) + (rSges(3,1) * t60 - rSges(3,2) * t56) * t104) - m(4) * (g(1) * (t101 * t36 - t35 * t74 - t29) + g(2) * (t101 * t34 - t33 * t74 - t27) + g(3) * t88 + (rSges(4,3) * t56 + t60 * t74) * t104) - m(5) * (g(1) * (t102 * t36 + t35 * t73 + t82) + g(2) * (t102 * t34 + t33 * t73 + t83) + t71 + (rSges(5,1) * t56 - t60 * t73) * t104) - m(6) * (g(1) * (t112 * t36 + t114 * t35 + t82) + g(2) * (t112 * t34 + t114 * t33 + t83) + t71 + (t56 * pkin(4) + (t55 * t93 + t56 * t58) * rSges(6,1) + (-t54 * t56 + t55 * t91) * rSges(6,2)) * t104 + (t90 * t104 - t117 * t59) * t100) - m(7) * (t71 + (t113 * t60 + t68 * t56) * t104 - t117 * t113 + (t66 * t34 + t83) * g(2) + (t66 * t36 + t82) * g(1)) -m(4) * (g(1) * (-rSges(4,1) * t20 - rSges(4,2) * t21) + g(2) * (-rSges(4,1) * t16 - rSges(4,2) * t17) + g(3) * (-rSges(4,1) * t31 - rSges(4,2) * t32)) - m(5) * (g(1) * (rSges(5,2) * t20 + t21 * t85 - t14) + g(2) * (rSges(5,2) * t16 + t17 * t85 - t12) + g(3) * (rSges(5,2) * t31 + t32 * t85 - t26)) + (-g(1) * (-t87 * t20 - t14) - g(2) * (-t87 * t16 - t12) - g(3) * (-t87 * t31 - t26) - t116 * t64) * m(7) + (-g(1) * (-t100 * t20 - t14) - g(2) * (-t100 * t16 - t12) - g(3) * (-t100 * t31 - t26) - t116 * t67) * m(6) (-m(5) - m(6) - m(7)) * (g(1) * t20 + g(2) * t16 + g(3) * t31) -m(6) * (g(1) * (rSges(6,1) * t8 - rSges(6,2) * t9) + g(2) * (t72 * rSges(6,1) + (-t16 * t54 - t33 * t58) * rSges(6,2)) + g(3) * (t69 * rSges(6,1) + (-t31 * t54 + t53 * t91) * rSges(6,2))) - t63 - m(7) * (g(1) * t8 + g(2) * t72 + g(3) * t69) * pkin(5), -t63];
taug  = t1(:);
