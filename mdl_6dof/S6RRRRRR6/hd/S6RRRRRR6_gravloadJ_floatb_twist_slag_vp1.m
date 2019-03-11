% Calculate Gravitation load on the joints for
% S6RRRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 04:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRR6_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR6_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR6_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR6_gravloadJ_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR6_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRR6_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 04:11:47
% EndTime: 2019-03-10 04:11:50
% DurationCPUTime: 1.43s
% Computational Cost: add. (1035->204), mult. (1644->288), div. (0->0), fcn. (1941->14), ass. (0->90)
t149 = pkin(11) + rSges(6,3);
t125 = pkin(12) + pkin(11) + rSges(7,3);
t124 = cos(pkin(6));
t82 = sin(pkin(6));
t85 = sin(qJ(2));
t135 = t82 * t85;
t84 = sin(qJ(3));
t88 = cos(qJ(3));
t166 = t124 * t88 - t84 * t135;
t133 = t82 * t88;
t86 = sin(qJ(1));
t116 = t86 * t124;
t148 = cos(qJ(1));
t89 = cos(qJ(2));
t59 = -t116 * t85 + t148 * t89;
t35 = t86 * t133 - t59 * t84;
t83 = sin(qJ(5));
t87 = cos(qJ(5));
t165 = rSges(6,1) * t87 - rSges(6,2) * t83 + pkin(4);
t73 = pkin(5) * t87 + pkin(4);
t80 = qJ(5) + qJ(6);
t75 = sin(t80);
t77 = cos(t80);
t164 = -rSges(7,1) * t77 + rSges(7,2) * t75 - t73;
t120 = t82 * t148;
t112 = t124 * t148;
t57 = t112 * t85 + t86 * t89;
t81 = qJ(3) + qJ(4);
t76 = sin(t81);
t78 = cos(t81);
t30 = -t76 * t120 + t57 * t78;
t56 = -t112 * t89 + t85 * t86;
t163 = t30 * t83 - t56 * t87;
t139 = t56 * t83;
t162 = -t30 * t87 - t139;
t161 = t30 * t75 - t56 * t77;
t160 = -t30 * t77 - t56 * t75;
t159 = rSges(6,1) * t83 + rSges(6,2) * t87;
t154 = g(2) * t56;
t58 = t116 * t89 + t148 * t85;
t158 = g(1) * t58 + t154;
t157 = t125 * t76 - t164 * t78;
t156 = t149 * t76 + t165 * t78;
t153 = g(2) * t57;
t132 = t82 * t89;
t74 = pkin(3) * t88 + pkin(2);
t152 = g(3) * t74 * t132;
t151 = g(3) * t82;
t150 = -pkin(9) - rSges(4,3);
t137 = t58 * t83;
t134 = t82 * t86;
t118 = -t78 * t120 - t57 * t76;
t131 = rSges(5,1) * t118 - t30 * rSges(5,2);
t33 = -t134 * t78 + t59 * t76;
t34 = t134 * t76 + t59 * t78;
t130 = -t33 * rSges(5,1) - t34 * rSges(5,2);
t91 = -pkin(10) - pkin(9);
t129 = -t56 * t74 - t57 * t91;
t128 = -t58 * t74 - t59 * t91;
t50 = t124 * t78 - t135 * t76;
t51 = t124 * t76 + t135 * t78;
t127 = t50 * rSges(5,1) - t51 * rSges(5,2);
t126 = t148 * pkin(1) + pkin(8) * t134;
t123 = t84 * t134;
t119 = -t86 * pkin(1) + pkin(8) * t120;
t68 = t84 * t120;
t117 = -t57 * t88 + t68;
t114 = t35 * pkin(3);
t113 = pkin(3) * t123 - t58 * t91 + t59 * t74 + t126;
t111 = rSges(5,1) * t78 - rSges(5,2) * t76;
t7 = -t34 * t83 + t58 * t87;
t109 = t166 * pkin(3);
t108 = rSges(4,1) * t88 - rSges(4,2) * t84 + pkin(2);
t106 = -t132 * t87 - t51 * t83;
t104 = pkin(3) * t68 + t56 * t91 - t57 * t74 + t119;
t103 = -t164 * t118 + t125 * t30;
t102 = t125 * t34 + t164 * t33;
t101 = t125 * t51 - t164 * t50;
t100 = t120 * t88 + t57 * t84;
t99 = rSges(7,1) * t75 + rSges(7,2) * t77 + pkin(5) * t83;
t98 = t165 * t118 + t149 * t30;
t97 = t149 * t34 - t165 * t33;
t96 = t149 * t51 + t165 * t50;
t95 = t100 * pkin(3);
t5 = -t34 * t75 + t58 * t77;
t6 = t34 * t77 + t58 * t75;
t92 = m(7) * (g(1) * (t5 * rSges(7,1) - t6 * rSges(7,2)) + g(2) * (-t161 * rSges(7,1) + t160 * rSges(7,2)) + g(3) * ((-t132 * t77 - t51 * t75) * rSges(7,1) + (t132 * t75 - t51 * t77) * rSges(7,2)));
t36 = t59 * t88 + t123;
t8 = t34 * t87 + t137;
t1 = [-m(2) * (g(1) * (-t86 * rSges(2,1) - rSges(2,2) * t148) + g(2) * (rSges(2,1) * t148 - t86 * rSges(2,2))) - m(3) * (g(1) * (-t57 * rSges(3,1) + t56 * rSges(3,2) + rSges(3,3) * t120 + t119) + g(2) * (rSges(3,1) * t59 - rSges(3,2) * t58 + rSges(3,3) * t134 + t126)) - m(4) * (g(1) * (rSges(4,1) * t117 + rSges(4,2) * t100 - t57 * pkin(2) + t150 * t56 + t119) + g(2) * (rSges(4,1) * t36 + rSges(4,2) * t35 + pkin(2) * t59 - t150 * t58 + t126)) - m(5) * (g(1) * (-rSges(5,1) * t30 - rSges(5,2) * t118 - rSges(5,3) * t56 + t104) + g(2) * (rSges(5,1) * t34 - rSges(5,2) * t33 + rSges(5,3) * t58 + t113)) - m(6) * (g(1) * (t162 * rSges(6,1) + t163 * rSges(6,2) - t30 * pkin(4) + t149 * t118 + t104) + g(2) * (rSges(6,1) * t8 + rSges(6,2) * t7 + pkin(4) * t34 + t149 * t33 + t113)) - m(7) * (g(1) * (t160 * rSges(7,1) + t161 * rSges(7,2) - pkin(5) * t139 + t125 * t118 - t30 * t73 + t104) + g(2) * (rSges(7,1) * t6 + rSges(7,2) * t5 + pkin(5) * t137 + t125 * t33 + t34 * t73 + t113)) -m(3) * (g(1) * (-rSges(3,1) * t58 - rSges(3,2) * t59) + g(2) * (-rSges(3,1) * t56 - rSges(3,2) * t57) + (rSges(3,1) * t89 - rSges(3,2) * t85) * t151) - m(4) * (g(1) * (-t108 * t58 - t150 * t59) - t150 * t153 - t108 * t154 + (t108 * t89 - t150 * t85) * t151) - m(5) * (g(1) * (rSges(5,3) * t59 - t111 * t58 + t128) + g(2) * (rSges(5,3) * t57 - t111 * t56 + t129) + t152 + (t111 * t89 + (rSges(5,3) - t91) * t85) * t151) - m(6) * (g(1) * (t159 * t59 + t128) + g(2) * (t159 * t57 + t129) + t152 + ((-t91 + t159) * t85 + t156 * t89) * t151 - t158 * t156) - m(7) * (g(2) * t129 + t152 + t99 * t153 + ((-t91 + t99) * t85 + t157 * t89) * t151 - t158 * t157 + (t99 * t59 + t128) * g(1)) -m(4) * (g(1) * (rSges(4,1) * t35 - rSges(4,2) * t36) + g(2) * (-rSges(4,1) * t100 + rSges(4,2) * t117) + g(3) * (t166 * rSges(4,1) + (-t124 * t84 - t133 * t85) * rSges(4,2))) - m(5) * (g(1) * (t114 + t130) + g(2) * (-t95 + t131) + g(3) * (t109 + t127)) - m(6) * (g(1) * (t114 + t97) + g(2) * (-t95 + t98) + g(3) * (t109 + t96)) - m(7) * (g(1) * (t102 + t114) + g(2) * (t103 - t95) + g(3) * (t101 + t109)) -m(5) * (g(1) * t130 + g(2) * t131 + g(3) * t127) - m(6) * (g(1) * t97 + g(2) * t98 + g(3) * t96) - m(7) * (g(1) * t102 + g(2) * t103 + g(3) * t101) -m(6) * (g(1) * (rSges(6,1) * t7 - rSges(6,2) * t8) + g(2) * (-rSges(6,1) * t163 + t162 * rSges(6,2)) + g(3) * (t106 * rSges(6,1) + (t132 * t83 - t51 * t87) * rSges(6,2))) - t92 - m(7) * (g(1) * t7 - g(2) * t163 + g(3) * t106) * pkin(5), -t92];
taug  = t1(:);
