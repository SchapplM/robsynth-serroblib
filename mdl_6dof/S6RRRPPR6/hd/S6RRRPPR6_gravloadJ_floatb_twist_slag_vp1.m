% Calculate Gravitation load on the joints for
% S6RRRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-03-09 15:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPPR6_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR6_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR6_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR6_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR6_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPR6_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:47:57
% EndTime: 2019-03-09 15:48:00
% DurationCPUTime: 1.29s
% Computational Cost: add. (728->197), mult. (1259->275), div. (0->0), fcn. (1466->12), ass. (0->77)
t60 = sin(qJ(2));
t61 = sin(qJ(1));
t64 = cos(qJ(2));
t109 = cos(qJ(1));
t93 = cos(pkin(6));
t78 = t93 * t109;
t33 = t60 * t78 + t61 * t64;
t59 = sin(qJ(3));
t63 = cos(qJ(3));
t56 = sin(pkin(6));
t89 = t56 * t109;
t66 = t33 * t59 + t63 * t89;
t65 = t66 * pkin(3);
t55 = qJ(3) + pkin(11);
t52 = sin(t55);
t53 = cos(t55);
t9 = t33 * t52 + t53 * t89;
t123 = -t9 * pkin(4) - t65;
t102 = t56 * t60;
t122 = -t59 * t102 + t93 * t63;
t121 = pkin(4) * t53 + qJ(5) * t52;
t100 = t56 * t63;
t86 = t61 * t93;
t35 = t109 * t64 - t60 * t86;
t15 = t61 * t100 - t35 * t59;
t32 = t60 * t61 - t64 * t78;
t58 = sin(qJ(6));
t62 = cos(qJ(6));
t120 = -t32 * t62 - t58 * t9;
t119 = -t32 * t58 + t62 * t9;
t110 = pkin(10) + rSges(7,3);
t118 = t110 * t53;
t113 = g(2) * t32;
t34 = t109 * t60 + t64 * t86;
t117 = -g(1) * t34 - t113;
t116 = -m(6) - m(7);
t112 = g(3) * t56;
t111 = -pkin(9) - rSges(4,3);
t104 = t52 * t58;
t103 = t52 * t62;
t101 = t56 * t61;
t99 = t56 * t64;
t51 = pkin(3) * t63 + pkin(2);
t57 = -qJ(4) - pkin(9);
t98 = -t32 * t51 - t33 * t57;
t97 = -t34 * t51 - t35 * t57;
t96 = t109 * pkin(1) + pkin(8) * t101;
t94 = rSges(6,3) + qJ(5);
t92 = t59 * t101;
t88 = -t61 * pkin(1) + pkin(8) * t89;
t10 = t33 * t53 - t52 * t89;
t45 = t59 * t89;
t87 = -t33 * t63 + t45;
t84 = -t121 * t32 + t98;
t83 = -t121 * t34 + t97;
t38 = t51 * t99;
t82 = g(3) * (t121 * t99 + t38);
t80 = t15 * pkin(3);
t79 = pkin(3) * t92 - t34 * t57 + t35 * t51 + t96;
t77 = rSges(5,1) * t53 - rSges(5,2) * t52;
t76 = rSges(7,1) * t58 + rSges(7,2) * t62;
t75 = rSges(6,2) * t53 - rSges(6,3) * t52;
t13 = -t53 * t101 + t35 * t52;
t74 = -t13 * pkin(4) + t80;
t73 = t122 * pkin(3);
t14 = t52 * t101 + t35 * t53;
t72 = t14 * pkin(4) + t79;
t71 = rSges(4,1) * t63 - rSges(4,2) * t59 + pkin(2);
t26 = t52 * t102 - t93 * t53;
t69 = -t26 * pkin(4) + t73;
t68 = pkin(3) * t45 + t32 * t57 - t33 * t51 + t88;
t67 = -pkin(4) * t10 + t68;
t27 = t53 * t102 + t93 * t52;
t16 = t35 * t63 + t92;
t3 = t13 * t58 + t34 * t62;
t2 = t13 * t62 - t34 * t58;
t1 = [-m(2) * (g(1) * (-t61 * rSges(2,1) - t109 * rSges(2,2)) + g(2) * (t109 * rSges(2,1) - t61 * rSges(2,2))) - m(3) * (g(1) * (-t33 * rSges(3,1) + t32 * rSges(3,2) + rSges(3,3) * t89 + t88) + g(2) * (rSges(3,1) * t35 - rSges(3,2) * t34 + rSges(3,3) * t101 + t96)) - m(4) * (g(1) * (t87 * rSges(4,1) + t66 * rSges(4,2) - t33 * pkin(2) + t111 * t32 + t88) + g(2) * (rSges(4,1) * t16 + rSges(4,2) * t15 + pkin(2) * t35 - t111 * t34 + t96)) - m(5) * (g(1) * (-rSges(5,1) * t10 + rSges(5,2) * t9 - rSges(5,3) * t32 + t68) + g(2) * (rSges(5,1) * t14 - rSges(5,2) * t13 + rSges(5,3) * t34 + t79)) - m(6) * (g(1) * (-rSges(6,1) * t32 + rSges(6,2) * t10 - t9 * t94 + t67) + g(2) * (rSges(6,1) * t34 - rSges(6,2) * t14 + t94 * t13 + t72)) - m(7) * (g(1) * (rSges(7,1) * t120 - rSges(7,2) * t119 - t32 * pkin(5) - t9 * qJ(5) - t110 * t10 + t67) + g(2) * (rSges(7,1) * t3 + rSges(7,2) * t2 + pkin(5) * t34 + qJ(5) * t13 + t110 * t14 + t72)) -m(3) * (g(1) * (-rSges(3,1) * t34 - rSges(3,2) * t35) + g(2) * (-rSges(3,1) * t32 - rSges(3,2) * t33) + (rSges(3,1) * t64 - rSges(3,2) * t60) * t112) - m(4) * (g(1) * (-t111 * t35 - t71 * t34) - g(2) * t111 * t33 - t71 * t113 + (-t111 * t60 + t71 * t64) * t112) - m(5) * (g(1) * (rSges(5,3) * t35 - t77 * t34 + t97) + g(2) * (rSges(5,3) * t33 - t77 * t32 + t98) + g(3) * t38 + (t77 * t64 + (rSges(5,3) - t57) * t60) * t112) - m(6) * (g(1) * (rSges(6,1) * t35 + t75 * t34 + t83) + g(2) * (rSges(6,1) * t33 + t75 * t32 + t84) + t82 + (-t75 * t64 + (rSges(6,1) - t57) * t60) * t112) - m(7) * (g(1) * (t35 * pkin(5) + (-t34 * t104 + t35 * t62) * rSges(7,1) + (-t34 * t103 - t35 * t58) * rSges(7,2) + t83) + g(2) * (t33 * pkin(5) + (-t32 * t104 + t33 * t62) * rSges(7,1) + (-t32 * t103 - t33 * t58) * rSges(7,2) + t84) + t82 + t117 * t118 + ((t62 * rSges(7,1) - t58 * rSges(7,2) + pkin(5) - t57) * t60 + (t76 * t52 + t118) * t64) * t112) -m(4) * (g(1) * (rSges(4,1) * t15 - rSges(4,2) * t16) + g(2) * (-t66 * rSges(4,1) + t87 * rSges(4,2)) + g(3) * (t122 * rSges(4,1) + (-t60 * t100 - t93 * t59) * rSges(4,2))) - m(5) * (g(1) * (-rSges(5,1) * t13 - rSges(5,2) * t14 + t80) + g(2) * (-t9 * rSges(5,1) - t10 * rSges(5,2) - t65) + g(3) * (-rSges(5,1) * t26 - rSges(5,2) * t27 + t73)) - m(6) * (g(1) * (rSges(6,2) * t13 + t94 * t14 + t74) + g(2) * (t9 * rSges(6,2) + t94 * t10 + t123) + g(3) * (rSges(6,2) * t26 + t94 * t27 + t69)) + (-g(1) * (-t110 * t13 + t74) - g(2) * (-t110 * t9 + t123) - g(3) * (-t110 * t26 + t69) - (g(1) * t14 + g(2) * t10 + g(3) * t27) * (qJ(5) + t76)) * m(7) (-m(5) + t116) * (-g(3) * t99 - t117) t116 * (g(1) * t13 + g(2) * t9 + g(3) * t26) -m(7) * (g(1) * (rSges(7,1) * t2 - rSges(7,2) * t3) + g(2) * (rSges(7,1) * t119 + rSges(7,2) * t120) + g(3) * ((t26 * t62 + t58 * t99) * rSges(7,1) + (-t26 * t58 + t62 * t99) * rSges(7,2)))];
taug  = t1(:);
