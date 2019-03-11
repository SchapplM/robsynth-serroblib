% Calculate Gravitation load on the joints for
% S6RRRPRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
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
% Datum: 2019-03-09 17:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRP9_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP9_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP9_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP9_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP9_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP9_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:21:53
% EndTime: 2019-03-09 17:21:57
% DurationCPUTime: 1.21s
% Computational Cost: add. (474->189), mult. (1099->253), div. (0->0), fcn. (1232->8), ass. (0->80)
t58 = sin(qJ(1));
t60 = cos(qJ(3));
t61 = cos(qJ(2));
t56 = sin(qJ(3));
t62 = cos(qJ(1));
t92 = t62 * t56;
t33 = -t58 * t60 + t61 * t92;
t93 = t61 * t62;
t34 = t58 * t56 + t60 * t93;
t55 = sin(qJ(5));
t59 = cos(qJ(5));
t11 = t33 * t55 + t34 * t59;
t73 = -t33 * t59 + t34 * t55;
t120 = -rSges(6,1) * t73 - rSges(6,2) * t11;
t105 = rSges(7,1) + pkin(5);
t87 = rSges(7,3) + qJ(6);
t119 = -t105 * t73 + t11 * t87;
t96 = t58 * t61;
t31 = t56 * t96 + t60 * t62;
t94 = t60 * t61;
t32 = t58 * t94 - t92;
t114 = -t31 * t59 + t32 * t55;
t5 = t31 * t55 + t32 * t59;
t118 = -rSges(6,1) * t114 - rSges(6,2) * t5;
t117 = -t105 * t114 + t5 * t87;
t108 = g(2) * t58;
t113 = g(1) * t62 + t108;
t57 = sin(qJ(2));
t101 = t56 * t57;
t95 = t59 * t60;
t19 = -t101 * t55 - t57 * t95;
t100 = t56 * t59;
t98 = t57 * t60;
t20 = -t100 * t57 + t55 * t98;
t112 = -t105 * t20 - t19 * t87;
t111 = -pkin(3) - pkin(4);
t110 = g(1) * t58;
t107 = g(3) * t57;
t52 = t61 * pkin(2);
t106 = -rSges(5,1) - pkin(3);
t104 = -rSges(7,2) - pkin(9);
t103 = -rSges(6,3) - pkin(9);
t99 = t56 * t61;
t97 = t57 * t62;
t91 = t57 * pkin(8) + t52;
t90 = t62 * pkin(1) + t58 * pkin(7);
t89 = qJ(4) * t56;
t88 = rSges(5,3) + qJ(4);
t86 = -pkin(1) - t52;
t85 = t62 * t104;
t84 = t103 * t62;
t83 = pkin(3) * t94 + t61 * t89 + t91;
t82 = pkin(2) * t93 + pkin(8) * t97 + t90;
t24 = t31 * pkin(3);
t81 = -t31 * pkin(4) + qJ(4) * t32 - t24;
t53 = t62 * pkin(7);
t80 = -t32 * pkin(3) - qJ(4) * t31 + t53;
t28 = t33 * pkin(3);
t79 = -t33 * pkin(4) + qJ(4) * t34 - t28;
t78 = t34 * pkin(3) + t82;
t77 = pkin(4) * t94 + t83;
t76 = rSges(3,1) * t61 - rSges(3,2) * t57;
t74 = -rSges(6,1) * t20 + rSges(6,2) * t19;
t72 = t55 * t56 + t95;
t39 = qJ(4) * t98;
t71 = t101 * t111 + t39;
t69 = t58 * t57 * pkin(9) - t32 * pkin(4) + t80;
t68 = t57 * (t55 * t60 - t100);
t67 = t57 * t72;
t65 = t34 * pkin(4) + t33 * qJ(4) + t78;
t63 = t113 * (t111 * t60 - pkin(2) - t89);
t46 = pkin(8) * t93;
t42 = pkin(8) * t96;
t22 = t72 * t61;
t21 = t55 * t94 - t59 * t99;
t15 = t62 * t67;
t14 = t62 * t68;
t13 = t58 * t67;
t12 = t58 * t68;
t1 = [-m(2) * (g(1) * (-t58 * rSges(2,1) - rSges(2,2) * t62) + g(2) * (rSges(2,1) * t62 - t58 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) * t62 + t53) + g(2) * (rSges(3,1) * t93 - rSges(3,2) * t97 + t90) + (g(1) * (-pkin(1) - t76) + g(2) * rSges(3,3)) * t58) - m(4) * (g(1) * (-rSges(4,1) * t32 + rSges(4,2) * t31 + t53) + g(2) * (t34 * rSges(4,1) - t33 * rSges(4,2) + rSges(4,3) * t97 + t82) + ((-rSges(4,3) - pkin(8)) * t57 + t86) * t110) - m(5) * (g(1) * (-rSges(5,1) * t32 - rSges(5,3) * t31 + t80) + g(2) * (t34 * rSges(5,1) + rSges(5,2) * t97 + t88 * t33 + t78) + ((-rSges(5,2) - pkin(8)) * t57 + t86) * t110) - m(6) * (g(1) * (-rSges(6,1) * t5 + rSges(6,2) * t114 + t69) + g(2) * (t11 * rSges(6,1) - rSges(6,2) * t73 + t57 * t84 + t65) + ((rSges(6,3) - pkin(8)) * t57 + t86) * t110) - m(7) * (g(1) * (-t105 * t5 - t114 * t87 + t69) + g(2) * (t105 * t11 + t57 * t85 + t73 * t87 + t65) + ((rSges(7,2) - pkin(8)) * t57 + t86) * t110) -m(3) * (g(3) * t76 + t113 * (-rSges(3,1) * t57 - rSges(3,2) * t61)) - m(4) * (g(1) * (rSges(4,3) * t93 + t46) + g(2) * (rSges(4,3) * t96 + t42) + g(3) * (rSges(4,1) * t94 - rSges(4,2) * t99 + t91) + (g(3) * rSges(4,3) + t113 * (-rSges(4,1) * t60 + rSges(4,2) * t56 - pkin(2))) * t57) - m(5) * (g(1) * (rSges(5,2) * t93 + t46) + g(2) * (rSges(5,2) * t96 + t42) + g(3) * (rSges(5,1) * t94 + rSges(5,3) * t99 + t83) + (g(3) * rSges(5,2) + t113 * (t106 * t60 - t56 * t88 - pkin(2))) * t57) - m(6) * (g(1) * (-t15 * rSges(6,1) + t14 * rSges(6,2) + t46) + g(2) * (-rSges(6,1) * t13 + rSges(6,2) * t12 + t42) + g(3) * (rSges(6,1) * t22 - rSges(6,2) * t21 + t77) + (g(1) * t84 + t103 * t108) * t61 + (g(3) * t103 + t63) * t57) - m(7) * (g(1) * (-t105 * t15 - t87 * t14 + t46) + g(2) * (-t105 * t13 - t87 * t12 + t42) + g(3) * (t105 * t22 + t87 * t21 + t77) + (g(1) * t85 + t104 * t108) * t61 + (g(3) * t104 + t63) * t57) -m(4) * (g(1) * (-rSges(4,1) * t33 - rSges(4,2) * t34) + g(2) * (-rSges(4,1) * t31 - rSges(4,2) * t32) + (-rSges(4,1) * t56 - rSges(4,2) * t60) * t107) - m(5) * (g(1) * (-rSges(5,1) * t33 + t34 * t88 - t28) + g(2) * (-rSges(5,1) * t31 + t32 * t88 - t24) + g(3) * t39 + (rSges(5,3) * t60 + t106 * t56) * t107) - m(6) * (g(1) * (-t120 + t79) + g(2) * (-t118 + t81) + g(3) * (t71 - t74)) - m(7) * (g(1) * (-t119 + t79) + g(2) * (-t117 + t81) + g(3) * (-t112 + t71)) (-m(5) - m(6) - m(7)) * (g(1) * t33 + g(2) * t31 + g(3) * t101) -m(6) * (g(1) * t120 + g(2) * t118 + g(3) * t74) - m(7) * (g(1) * t119 + g(2) * t117 + g(3) * t112) -m(7) * (g(1) * t73 + g(2) * t114 + g(3) * t20)];
taug  = t1(:);
