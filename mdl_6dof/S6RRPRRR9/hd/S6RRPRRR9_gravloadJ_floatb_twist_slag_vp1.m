% Calculate Gravitation load on the joints for
% S6RRPRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 14:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRR9_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR9_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR9_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR9_gravloadJ_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR9_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR9_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:08:32
% EndTime: 2019-03-09 14:08:34
% DurationCPUTime: 1.11s
% Computational Cost: add. (869->181), mult. (1172->253), div. (0->0), fcn. (1354->14), ass. (0->79)
t119 = pkin(11) + rSges(7,3);
t100 = cos(pkin(6));
t67 = sin(pkin(6));
t71 = sin(qJ(2));
t109 = t67 * t71;
t65 = pkin(12) + qJ(4);
t59 = sin(t65);
t60 = cos(t65);
t133 = t100 * t60 - t59 * t109;
t72 = sin(qJ(1));
t108 = t67 * t72;
t117 = cos(qJ(1));
t74 = cos(qJ(2));
t92 = t72 * t100;
t43 = t117 * t74 - t71 * t92;
t20 = t60 * t108 - t43 * t59;
t70 = sin(qJ(6));
t73 = cos(qJ(6));
t132 = rSges(7,1) * t73 - rSges(7,2) * t70 + pkin(5);
t87 = t100 * t117;
t41 = t71 * t87 + t72 * t74;
t61 = qJ(5) + t65;
t56 = sin(t61);
t57 = cos(t61);
t96 = t67 * t117;
t15 = t41 * t57 - t56 * t96;
t40 = t71 * t72 - t74 * t87;
t131 = t15 * t70 - t40 * t73;
t130 = -t15 * t73 - t40 * t70;
t107 = t67 * t74;
t129 = g(3) * t107;
t128 = rSges(7,1) * t70 + rSges(7,2) * t73;
t123 = g(2) * t40;
t42 = t117 * t71 + t74 * t92;
t127 = g(1) * t42 + t123;
t126 = t119 * t56 + t132 * t57;
t94 = -t41 * t56 - t57 * t96;
t125 = rSges(6,1) * t94 - t15 * rSges(6,2);
t122 = g(2) * t41;
t68 = cos(pkin(12));
t58 = t68 * pkin(3) + pkin(2);
t44 = pkin(4) * t60 + t58;
t121 = t44 * t129;
t120 = g(3) * t67;
t18 = -t57 * t108 + t43 * t56;
t19 = t56 * t108 + t43 * t57;
t118 = -t18 * rSges(6,1) - t19 * rSges(6,2);
t69 = -pkin(9) - qJ(3);
t64 = -pkin(10) + t69;
t106 = -t40 * t44 - t41 * t64;
t105 = -t42 * t44 - t43 * t64;
t31 = t100 * t57 - t56 * t109;
t32 = t100 * t56 + t57 * t109;
t104 = t31 * rSges(6,1) - t32 * rSges(6,2);
t103 = t117 * pkin(1) + pkin(8) * t108;
t102 = t69 - rSges(5,3);
t101 = qJ(3) + rSges(4,3);
t66 = sin(pkin(12));
t98 = t66 * t108;
t95 = -t72 * pkin(1) + pkin(8) * t96;
t93 = -t41 * t60 + t59 * t96;
t90 = t66 * t96;
t89 = t20 * pkin(4);
t49 = pkin(3) * t66 + pkin(4) * t59;
t88 = t49 * t108 - t42 * t64 + t43 * t44 + t103;
t86 = rSges(6,1) * t57 - rSges(6,2) * t56;
t85 = t133 * pkin(4);
t84 = rSges(4,1) * t68 - rSges(4,2) * t66 + pkin(2);
t82 = rSges(5,1) * t60 - rSges(5,2) * t59 + t58;
t81 = t40 * t64 - t41 * t44 + t49 * t96 + t95;
t80 = t41 * t59 + t60 * t96;
t79 = t119 * t15 + t132 * t94;
t78 = t119 * t19 - t132 * t18;
t77 = t119 * t32 + t132 * t31;
t76 = t80 * pkin(4);
t21 = t59 * t108 + t43 * t60;
t2 = t19 * t73 + t42 * t70;
t1 = -t19 * t70 + t42 * t73;
t3 = [-m(2) * (g(1) * (-t72 * rSges(2,1) - t117 * rSges(2,2)) + g(2) * (t117 * rSges(2,1) - t72 * rSges(2,2))) - m(3) * (g(1) * (-t41 * rSges(3,1) + t40 * rSges(3,2) + rSges(3,3) * t96 + t95) + g(2) * (rSges(3,1) * t43 - rSges(3,2) * t42 + rSges(3,3) * t108 + t103)) - m(4) * (g(1) * (-t41 * pkin(2) + (-t41 * t68 + t90) * rSges(4,1) + (t41 * t66 + t68 * t96) * rSges(4,2) - t101 * t40 + t95) + g(2) * (t43 * pkin(2) + (t43 * t68 + t98) * rSges(4,1) + (t68 * t108 - t43 * t66) * rSges(4,2) + t101 * t42 + t103)) - m(5) * (g(1) * (t93 * rSges(5,1) + t80 * rSges(5,2) + pkin(3) * t90 + t102 * t40 - t41 * t58 + t95) + g(2) * (rSges(5,1) * t21 + rSges(5,2) * t20 + pkin(3) * t98 - t102 * t42 + t43 * t58 + t103)) - m(6) * (g(1) * (-rSges(6,1) * t15 - rSges(6,2) * t94 - rSges(6,3) * t40 + t81) + g(2) * (rSges(6,1) * t19 - rSges(6,2) * t18 + rSges(6,3) * t42 + t88)) - m(7) * (g(1) * (rSges(7,1) * t130 + rSges(7,2) * t131 - t15 * pkin(5) + t119 * t94 + t81) + g(2) * (rSges(7,1) * t2 + rSges(7,2) * t1 + pkin(5) * t19 + t119 * t18 + t88)) -m(3) * (g(1) * (-rSges(3,1) * t42 - rSges(3,2) * t43) + g(2) * (-rSges(3,1) * t40 - rSges(3,2) * t41) + (rSges(3,1) * t74 - rSges(3,2) * t71) * t120) - m(4) * (g(1) * (t101 * t43 - t84 * t42) + t101 * t122 - t84 * t123 + (t101 * t71 + t84 * t74) * t120) - m(5) * (g(1) * (-t102 * t43 - t82 * t42) - t102 * t122 - t82 * t123 + (-t102 * t71 + t82 * t74) * t120) - m(6) * (g(1) * (rSges(6,3) * t43 - t86 * t42 + t105) + g(2) * (rSges(6,3) * t41 - t86 * t40 + t106) + t121 + (t86 * t74 + (rSges(6,3) - t64) * t71) * t120) - m(7) * (g(1) * (t128 * t43 + t105) + g(2) * (t128 * t41 + t106) + t121 + ((-t64 + t128) * t71 + t126 * t74) * t120 - t127 * t126) (-m(4) - m(5) - m(6) - m(7)) * (t127 - t129) -m(5) * (g(1) * (rSges(5,1) * t20 - rSges(5,2) * t21) + g(2) * (-t80 * rSges(5,1) + t93 * rSges(5,2)) + g(3) * (t133 * rSges(5,1) + (-t100 * t59 - t60 * t109) * rSges(5,2))) - m(6) * (g(1) * (t89 + t118) + g(2) * (-t76 + t125) + g(3) * (t85 + t104)) - m(7) * (g(1) * (t78 + t89) + g(2) * (-t76 + t79) + g(3) * (t77 + t85)) -m(6) * (g(1) * t118 + g(2) * t125 + g(3) * t104) - m(7) * (g(1) * t78 + g(2) * t79 + g(3) * t77) -m(7) * (g(1) * (rSges(7,1) * t1 - rSges(7,2) * t2) + g(2) * (-rSges(7,1) * t131 + rSges(7,2) * t130) + g(3) * ((-t73 * t107 - t32 * t70) * rSges(7,1) + (t70 * t107 - t32 * t73) * rSges(7,2)))];
taug  = t3(:);
