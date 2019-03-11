% Calculate Gravitation load on the joints for
% S6RRPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-03-09 10:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRPR6_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR6_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR6_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR6_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR6_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR6_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:38:11
% EndTime: 2019-03-09 10:38:15
% DurationCPUTime: 1.25s
% Computational Cost: add. (799->176), mult. (1935->250), div. (0->0), fcn. (2422->12), ass. (0->76)
t60 = sin(qJ(4));
t64 = cos(qJ(4));
t118 = pkin(4) * t64 + qJ(5) * t60;
t62 = sin(qJ(1));
t65 = cos(qJ(1));
t106 = cos(qJ(2));
t61 = sin(qJ(2));
t91 = sin(pkin(11));
t92 = cos(pkin(11));
t44 = -t106 * t92 + t61 * t91;
t58 = cos(pkin(6));
t66 = t58 * t44;
t68 = t106 * t91 + t61 * t92;
t26 = t62 * t66 - t65 * t68;
t89 = t62 * t106;
t97 = t65 * t61;
t42 = -t58 * t89 - t97;
t69 = t42 * pkin(2);
t67 = t26 * pkin(3) + t69;
t123 = t118 * t26 + t67;
t59 = sin(qJ(6));
t63 = cos(qJ(6));
t122 = t59 * rSges(7,1) + t63 * rSges(7,2);
t88 = t65 * t106;
t99 = t62 * t61;
t119 = t58 * t88 - t99;
t112 = -pkin(5) - pkin(9);
t114 = t63 * rSges(7,1) - t59 * rSges(7,2) - t112;
t107 = pkin(10) + rSges(7,3);
t96 = t68 * t58;
t27 = -t65 * t44 - t62 * t96;
t22 = t62 * t44 - t65 * t96;
t115 = t122 * t60;
t113 = -m(6) - m(7);
t110 = rSges(6,1) + pkin(9);
t108 = rSges(5,3) + pkin(9);
t57 = sin(pkin(6));
t105 = t57 * t62;
t104 = t57 * t65;
t90 = t106 * pkin(2);
t56 = t90 + pkin(1);
t100 = t62 * t56;
t37 = t44 * t57;
t53 = t57 * t90;
t95 = -t37 * pkin(3) + t53;
t93 = rSges(6,3) + qJ(5);
t39 = pkin(2) * t58 * t61 + (-pkin(8) - qJ(3)) * t57;
t87 = rSges(4,3) * t57 - t39;
t9 = -t60 * t104 - t22 * t64;
t84 = -t118 * t37 + t95;
t81 = t119 * pkin(2);
t50 = t65 * t56;
t80 = t27 * pkin(3) - t39 * t62 + t50;
t77 = rSges(5,1) * t64 - rSges(5,2) * t60;
t76 = rSges(6,2) * t64 - rSges(6,3) * t60;
t23 = -t62 * t68 - t65 * t66;
t75 = t23 * pkin(3) + t81;
t13 = t105 * t60 + t27 * t64;
t74 = t13 * pkin(4) + t80;
t8 = t104 * t64 - t22 * t60;
t73 = pkin(3) * t22 - t65 * t39 - t100;
t72 = qJ(5) + t122;
t71 = -pkin(4) * t9 + t73;
t70 = t118 * t23 + t75;
t43 = -t58 * t99 + t88;
t41 = -t58 * t97 - t89;
t38 = t68 * t57;
t30 = t38 * t64 + t58 * t60;
t29 = t38 * t60 - t58 * t64;
t28 = t29 * pkin(4);
t12 = -t105 * t64 + t27 * t60;
t6 = t12 * pkin(4);
t4 = t8 * pkin(4);
t3 = t12 * t59 - t26 * t63;
t2 = t12 * t63 + t26 * t59;
t1 = [-m(2) * (g(1) * (-t62 * rSges(2,1) - rSges(2,2) * t65) + g(2) * (rSges(2,1) * t65 - t62 * rSges(2,2))) - m(3) * (g(1) * (t41 * rSges(3,1) - rSges(3,2) * t119 - t62 * pkin(1)) + g(2) * (t43 * rSges(3,1) + t42 * rSges(3,2) + pkin(1) * t65) + (g(1) * t65 + g(2) * t62) * t57 * (rSges(3,3) + pkin(8))) - m(4) * (g(1) * (rSges(4,1) * t22 - t23 * rSges(4,2) + t65 * t87 - t100) + g(2) * (rSges(4,1) * t27 + rSges(4,2) * t26 + t62 * t87 + t50)) - m(5) * (g(1) * (-rSges(5,1) * t9 + rSges(5,2) * t8 + t108 * t23 + t73) + g(2) * (rSges(5,1) * t13 - rSges(5,2) * t12 - t108 * t26 + t80)) - m(6) * (g(1) * (rSges(6,2) * t9 + t110 * t23 - t8 * t93 + t71) + g(2) * (-rSges(6,2) * t13 - t110 * t26 + t12 * t93 + t74)) - m(7) * (g(1) * (-t107 * t9 + t114 * t23 - t72 * t8 + t71) + g(2) * (rSges(7,1) * t3 + rSges(7,2) * t2 + qJ(5) * t12 + t107 * t13 + t112 * t26 + t74)) -m(3) * (g(1) * (rSges(3,1) * t42 - rSges(3,2) * t43) + g(2) * (rSges(3,1) * t119 + rSges(3,2) * t41) + g(3) * (rSges(3,1) * t106 - rSges(3,2) * t61) * t57) - m(4) * (g(1) * (t26 * rSges(4,1) - rSges(4,2) * t27 + t69) + g(2) * (rSges(4,1) * t23 + rSges(4,2) * t22 + t81) + g(3) * (-rSges(4,1) * t37 - rSges(4,2) * t38 + t53)) - m(5) * (g(1) * (t108 * t27 + t26 * t77 + t67) + g(2) * (-t108 * t22 + t23 * t77 + t75) + g(3) * (t108 * t38 - t37 * t77 + t95)) - m(6) * (g(1) * (t110 * t27 - t26 * t76 + t123) + g(2) * (-t110 * t22 - t23 * t76 + t70) + g(3) * (t110 * t38 + t37 * t76 + t84)) - m(7) * (g(1) * (t114 * t27 + t115 * t26 + t123) + g(2) * (-t114 * t22 + t115 * t23 + t70) + g(3) * (t114 * t38 - t115 * t37 + t84) + (g(1) * t26 + g(2) * t23 - g(3) * t37) * t64 * t107) (-m(4) - m(5) + t113) * (g(3) * t58 + (g(1) * t62 - g(2) * t65) * t57) -m(5) * (g(1) * (-rSges(5,1) * t12 - rSges(5,2) * t13) + g(2) * (-rSges(5,1) * t8 - rSges(5,2) * t9) + g(3) * (-rSges(5,1) * t29 - rSges(5,2) * t30)) - m(6) * (g(1) * (rSges(6,2) * t12 + t13 * t93 - t6) + g(2) * (rSges(6,2) * t8 + t9 * t93 - t4) + g(3) * (rSges(6,2) * t29 + t30 * t93 - t28)) + (-g(1) * (-t107 * t12 - t6) - g(2) * (-t107 * t8 - t4) - g(3) * (-t107 * t29 - t28) - (g(1) * t13 + g(2) * t9 + g(3) * t30) * t72) * m(7), t113 * (g(1) * t12 + g(2) * t8 + g(3) * t29) -m(7) * (g(1) * (rSges(7,1) * t2 - rSges(7,2) * t3) + g(2) * ((t23 * t59 + t63 * t8) * rSges(7,1) + (t23 * t63 - t59 * t8) * rSges(7,2)) + g(3) * ((t29 * t63 - t37 * t59) * rSges(7,1) + (-t29 * t59 - t37 * t63) * rSges(7,2)))];
taug  = t1(:);
