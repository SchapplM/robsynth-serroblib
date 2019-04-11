% Calculate Gravitation load on the joints for
% S6RRRRRR10V2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,d6]';
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
% Datum: 2019-04-11 14:56
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRR10V2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(6,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10V2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR10V2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S6RRRRRR10V2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR10V2_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRR10V2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-11 14:41:38
% EndTime: 2019-04-11 14:41:44
% DurationCPUTime: 1.13s
% Computational Cost: add. (697->181), mult. (990->271), div. (0->0), fcn. (1076->12), ass. (0->90)
t123 = rSges(7,3) + pkin(6);
t60 = qJ(2) + qJ(3);
t57 = sin(t60);
t65 = sin(qJ(1));
t118 = t57 * t65;
t63 = sin(qJ(4));
t70 = cos(qJ(1));
t103 = t70 * t63;
t68 = cos(qJ(4));
t107 = t65 * t68;
t58 = cos(t60);
t33 = t107 * t58 - t103;
t62 = sin(qJ(5));
t67 = cos(qJ(5));
t10 = t118 * t62 + t33 * t67;
t104 = t68 * t70;
t108 = t65 * t63;
t32 = t108 * t58 + t104;
t61 = sin(qJ(6));
t66 = cos(qJ(6));
t137 = t10 * t61 - t32 * t66;
t136 = -t10 * t66 - t32 * t61;
t84 = rSges(7,1) * t66 - rSges(7,2) * t61;
t73 = -t123 * t62 - t67 * t84;
t110 = t62 * t68;
t78 = t110 * t57 + t58 * t67;
t21 = t78 * t65;
t105 = t67 * t68;
t29 = t105 * t57 - t58 * t62;
t22 = t29 * t65;
t134 = -rSges(6,1) * t22 + rSges(6,2) * t21;
t23 = t78 * t70;
t24 = t29 * t70;
t133 = -rSges(6,1) * t24 + rSges(6,2) * t23;
t101 = t57 * t108;
t114 = t58 * t65;
t132 = rSges(5,2) * t101 + rSges(5,3) * t114;
t113 = t58 * t70;
t99 = t57 * t103;
t131 = rSges(5,2) * t99 + rSges(5,3) * t113;
t130 = rSges(4,1) * t58 - rSges(4,2) * t57;
t129 = g(1) * t70 + g(2) * t65;
t128 = t57 * t129;
t64 = sin(qJ(2));
t127 = pkin(2) * t64;
t126 = g(1) * t65;
t52 = t57 * pkin(5);
t53 = t58 * pkin(3);
t122 = rSges(5,1) * t68;
t117 = t57 * t67;
t116 = t57 * t70;
t115 = t58 * t63;
t111 = t61 * t63;
t109 = t63 * t66;
t102 = t52 + t53;
t100 = t57 * t109;
t69 = cos(qJ(2));
t59 = t69 * pkin(2);
t56 = t59 + pkin(1);
t39 = t70 * t56;
t98 = pkin(3) * t113 + pkin(5) * t116 + t39;
t95 = -t56 - t53;
t44 = pkin(5) * t114;
t93 = -t127 * t65 + t44;
t46 = pkin(5) * t113;
t92 = -t127 * t70 + t46;
t30 = t110 * t58 - t117;
t31 = t105 * t58 + t57 * t62;
t91 = rSges(6,1) * t31 - rSges(6,2) * t30 + rSges(6,3) * t115 + t102;
t90 = -pkin(3) * t57 - t127;
t89 = rSges(3,1) * t69 - rSges(3,2) * t64;
t86 = -rSges(4,1) * t57 - rSges(4,2) * t58;
t85 = -rSges(6,1) * t67 + rSges(6,2) * t62;
t83 = rSges(7,1) * t61 + rSges(7,2) * t66;
t82 = pkin(1) + t89;
t81 = (-t100 * t65 + t22 * t61) * rSges(7,2) + (-t101 * t61 - t22 * t66) * rSges(7,1) + t44 - t123 * t21;
t80 = t46 + (t24 * t61 - t66 * t99) * rSges(7,2) + (-t24 * t66 - t61 * t99) * rSges(7,1) - t123 * t23;
t79 = t117 * t65 - t33 * t62;
t77 = t102 + (t109 * t58 - t31 * t61) * rSges(7,2) + (t111 * t58 + t31 * t66) * rSges(7,1) + t123 * t30;
t76 = -rSges(5,2) * t115 + rSges(5,3) * t57 + t122 * t58 + t102;
t74 = (t95 - t52) * t126;
t72 = (-pkin(3) - t122) * t128;
t71 = (-rSges(6,3) * t63 - pkin(3)) * t128;
t35 = t104 * t58 + t108;
t34 = t103 * t58 - t107;
t14 = t116 * t62 + t35 * t67;
t13 = -t116 * t67 + t35 * t62;
t2 = t14 * t66 + t34 * t61;
t1 = -t14 * t61 + t34 * t66;
t3 = [-m(2) * (g(1) * (-rSges(2,1) * t65 - rSges(2,2) * t70) + g(2) * (rSges(2,1) * t70 - rSges(2,2) * t65)) - m(3) * ((g(1) * rSges(3,3) + g(2) * t82) * t70 + (g(2) * rSges(3,3) - g(1) * t82) * t65) - m(4) * (g(2) * t39 + (g(1) * rSges(4,3) + g(2) * t130) * t70 + (g(1) * (-t56 - t130) + g(2) * rSges(4,3)) * t65) - m(5) * (g(1) * (-rSges(5,1) * t33 + rSges(5,2) * t32) + g(2) * (rSges(5,1) * t35 - rSges(5,2) * t34 + rSges(5,3) * t116 + t98) + ((-rSges(5,3) - pkin(5)) * t57 + t95) * t126) - m(6) * (g(1) * (-rSges(6,1) * t10 - rSges(6,2) * t79 - t32 * rSges(6,3)) + g(2) * (rSges(6,1) * t14 - rSges(6,2) * t13 + rSges(6,3) * t34 + t98) + t74) - m(7) * (g(1) * (rSges(7,1) * t136 + rSges(7,2) * t137 + t123 * t79) + g(2) * (rSges(7,1) * t2 + rSges(7,2) * t1 + t123 * t13 + t98) + t74) -m(3) * (g(3) * t89 + t129 * (-rSges(3,1) * t64 - rSges(3,2) * t69)) - m(4) * (g(3) * (t59 + t130) + t129 * (t86 - t127)) - m(5) * (g(1) * (t92 + t131) + g(2) * (t93 + t132) + g(3) * (t59 + t76) + t72) - m(6) * (g(1) * (t92 + t133) + g(2) * (t93 + t134) + g(3) * (t59 + t91) + t71) - m(7) * (g(1) * (t70 * t90 + t80) + g(2) * (t65 * t90 + t81) + g(3) * (t59 + t77)) -m(4) * (g(3) * t130 + t129 * t86) - m(5) * (g(1) * (t46 + t131) + g(2) * (t44 + t132) + g(3) * t76 + t72) - m(6) * (g(1) * (t46 + t133) + g(2) * (t44 + t134) + g(3) * t91 + t71) - m(7) * (g(1) * (-pkin(3) * t116 + t80) + g(2) * (-pkin(3) * t118 + t81) + g(3) * t77) -m(5) * (g(1) * (-rSges(5,1) * t34 - rSges(5,2) * t35) + g(2) * (-rSges(5,1) * t32 - rSges(5,2) * t33)) - m(6) * (g(1) * (rSges(6,3) * t35 + t34 * t85) + g(2) * (rSges(6,3) * t33 + t32 * t85)) - m(7) * (g(1) * (t34 * t73 + t35 * t83) + g(2) * (t32 * t73 + t33 * t83)) + ((m(5) * rSges(5,2) - m(6) * rSges(6,3) - m(7) * t83) * t68 + (m(5) * rSges(5,1) - m(6) * t85 - m(7) * t73) * t63) * g(3) * t57, -m(6) * (g(1) * (-rSges(6,1) * t13 - rSges(6,2) * t14) + g(2) * (rSges(6,1) * t79 - rSges(6,2) * t10) + g(3) * (-rSges(6,1) * t78 - rSges(6,2) * t29)) - m(7) * (g(1) * (t123 * t14 - t13 * t84) + g(2) * (t10 * t123 + t79 * t84) + g(3) * (t123 * t29 - t78 * t84)) -m(7) * (g(1) * (rSges(7,1) * t1 - rSges(7,2) * t2) + g(2) * (-rSges(7,1) * t137 + rSges(7,2) * t136) + g(3) * ((-t29 * t61 + t100) * rSges(7,1) + (-t111 * t57 - t29 * t66) * rSges(7,2)))];
taug  = t3(:);
