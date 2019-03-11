% Calculate Gravitation load on the joints for
% S6RRRRRR5
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
% Datum: 2019-03-10 04:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRR5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR5_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRR5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:56:55
% EndTime: 2019-03-10 03:56:58
% DurationCPUTime: 1.29s
% Computational Cost: add. (1081->204), mult. (1423->285), div. (0->0), fcn. (1646->14), ass. (0->89)
t140 = pkin(12) + rSges(7,3);
t119 = cos(pkin(6));
t81 = sin(pkin(6));
t84 = sin(qJ(2));
t133 = t81 * t84;
t80 = qJ(3) + qJ(4);
t74 = sin(t80);
t75 = cos(t80);
t152 = t119 * t75 - t74 * t133;
t85 = sin(qJ(1));
t132 = t81 * t85;
t111 = t85 * t119;
t88 = cos(qJ(2));
t89 = cos(qJ(1));
t54 = -t84 * t111 + t88 * t89;
t23 = t75 * t132 - t54 * t74;
t82 = sin(qJ(6));
t86 = cos(qJ(6));
t151 = rSges(7,1) * t86 - rSges(7,2) * t82 + pkin(5);
t129 = t81 * t89;
t110 = t89 * t119;
t52 = t84 * t110 + t85 * t88;
t76 = qJ(5) + t80;
t71 = sin(t76);
t72 = cos(t76);
t14 = -t71 * t129 + t52 * t72;
t51 = -t88 * t110 + t84 * t85;
t150 = t14 * t82 - t51 * t86;
t149 = -t14 * t86 - t51 * t82;
t148 = t82 * rSges(7,1) + t86 * rSges(7,2);
t147 = t140 * t71 + t151 * t72;
t90 = -pkin(10) - pkin(9);
t114 = -t72 * t129 - t52 * t71;
t146 = rSges(6,1) * t114 - t14 * rSges(6,2);
t145 = g(2) * t51;
t144 = g(2) * t52;
t130 = t81 * t88;
t87 = cos(qJ(3));
t77 = t87 * pkin(3);
t62 = pkin(4) * t75 + t77;
t59 = pkin(2) + t62;
t143 = g(3) * t59 * t130;
t142 = g(3) * t81;
t141 = -pkin(9) - rSges(4,3);
t17 = -t72 * t132 + t54 * t71;
t18 = t71 * t132 + t54 * t72;
t139 = -t17 * rSges(6,1) - t18 * rSges(6,2);
t131 = t81 * t87;
t79 = -pkin(11) + t90;
t126 = -t51 * t59 - t52 * t79;
t53 = t88 * t111 + t89 * t84;
t125 = -t53 * t59 - t54 * t79;
t38 = t119 * t72 - t71 * t133;
t39 = t119 * t71 + t72 * t133;
t124 = t38 * rSges(6,1) - t39 * rSges(6,2);
t83 = sin(qJ(3));
t61 = pkin(3) * t83 + pkin(4) * t74;
t123 = t62 * t132 - t54 * t61;
t122 = t119 * t62 - t61 * t133;
t121 = t89 * pkin(1) + pkin(8) * t132;
t120 = t90 - rSges(5,3);
t117 = t83 * t132;
t66 = t83 * t129;
t115 = -t85 * pkin(1) + pkin(8) * t129;
t113 = t74 * t129 - t52 * t75;
t112 = -t52 * t87 + t66;
t108 = t23 * pkin(4);
t107 = -t62 * t129 - t52 * t61;
t106 = t61 * t132 - t53 * t79 + t54 * t59 + t121;
t105 = rSges(6,1) * t72 - rSges(6,2) * t71;
t104 = t152 * pkin(4);
t103 = rSges(4,1) * t87 - rSges(4,2) * t83 + pkin(2);
t101 = t75 * t129 + t52 * t74;
t100 = t87 * t129 + t52 * t83;
t25 = t85 * t131 - t54 * t83;
t73 = t77 + pkin(2);
t99 = rSges(5,1) * t75 - rSges(5,2) * t74 + t73;
t98 = t61 * t129 + t51 * t79 - t52 * t59 + t115;
t97 = t101 * pkin(4);
t96 = t119 * t87 - t83 * t133;
t95 = t151 * t114 + t14 * t140;
t94 = t140 * t18 - t151 * t17;
t93 = t140 * t39 + t151 * t38;
t24 = t74 * t132 + t54 * t75;
t92 = g(1) * (t23 * rSges(5,1) - t24 * rSges(5,2)) + g(2) * (-t101 * rSges(5,1) + t113 * rSges(5,2)) + g(3) * (t152 * rSges(5,1) + (-t119 * t74 - t75 * t133) * rSges(5,2));
t26 = t54 * t87 + t117;
t2 = t18 * t86 + t53 * t82;
t1 = -t18 * t82 + t53 * t86;
t3 = [-m(2) * (g(1) * (-rSges(2,1) * t85 - rSges(2,2) * t89) + g(2) * (rSges(2,1) * t89 - rSges(2,2) * t85)) - m(3) * (g(1) * (-rSges(3,1) * t52 + rSges(3,2) * t51 + rSges(3,3) * t129 + t115) + g(2) * (rSges(3,1) * t54 - rSges(3,2) * t53 + rSges(3,3) * t132 + t121)) - m(4) * (g(1) * (t112 * rSges(4,1) + t100 * rSges(4,2) - t52 * pkin(2) + t141 * t51 + t115) + g(2) * (rSges(4,1) * t26 + rSges(4,2) * t25 + pkin(2) * t54 - t141 * t53 + t121)) - m(5) * (g(1) * (t113 * rSges(5,1) + t101 * rSges(5,2) + pkin(3) * t66 + t120 * t51 - t52 * t73 + t115) + g(2) * (t24 * rSges(5,1) + t23 * rSges(5,2) + pkin(3) * t117 - t120 * t53 + t54 * t73 + t121)) - m(6) * (g(1) * (-rSges(6,1) * t14 - rSges(6,2) * t114 - rSges(6,3) * t51 + t98) + g(2) * (rSges(6,1) * t18 - rSges(6,2) * t17 + rSges(6,3) * t53 + t106)) - m(7) * (g(1) * (rSges(7,1) * t149 + rSges(7,2) * t150 - t14 * pkin(5) + t140 * t114 + t98) + g(2) * (rSges(7,1) * t2 + rSges(7,2) * t1 + pkin(5) * t18 + t140 * t17 + t106)) -m(3) * (g(1) * (-rSges(3,1) * t53 - rSges(3,2) * t54) + g(2) * (-rSges(3,1) * t51 - rSges(3,2) * t52) + (rSges(3,1) * t88 - rSges(3,2) * t84) * t142) - m(4) * (g(1) * (-t103 * t53 - t141 * t54) - t141 * t144 - t103 * t145 + (t103 * t88 - t141 * t84) * t142) - m(5) * (g(1) * (-t120 * t54 - t99 * t53) - t120 * t144 - t99 * t145 + (-t120 * t84 + t99 * t88) * t142) - m(6) * (g(1) * (rSges(6,3) * t54 - t105 * t53 + t125) + g(2) * (rSges(6,3) * t52 - t105 * t51 + t126) + t143 + (t105 * t88 + (rSges(6,3) - t79) * t84) * t142) - m(7) * (g(2) * (t148 * t52 + t126) + t143 - t147 * t145 + ((-t79 + t148) * t84 + t147 * t88) * t142 + (-t147 * t53 + t148 * t54 + t125) * g(1)) -m(4) * (g(1) * (rSges(4,1) * t25 - rSges(4,2) * t26) + g(2) * (-t100 * rSges(4,1) + t112 * rSges(4,2)) + g(3) * (t96 * rSges(4,1) + (-t119 * t83 - t84 * t131) * rSges(4,2))) - m(5) * ((g(1) * t25 - g(2) * t100 + g(3) * t96) * pkin(3) + t92) - m(6) * (g(1) * (t123 + t139) + g(2) * (t107 + t146) + g(3) * (t122 + t124)) - m(7) * (g(1) * (t94 + t123) + g(2) * (t107 + t95) + g(3) * (t93 + t122)) -m(5) * t92 - m(6) * (g(1) * (t108 + t139) + g(2) * (-t97 + t146) + g(3) * (t104 + t124)) - m(7) * (g(1) * (t108 + t94) + g(2) * (-t97 + t95) + g(3) * (t104 + t93)) -m(6) * (g(1) * t139 + g(2) * t146 + g(3) * t124) - m(7) * (g(1) * t94 + g(2) * t95 + g(3) * t93) -m(7) * (g(1) * (rSges(7,1) * t1 - rSges(7,2) * t2) + g(2) * (-t150 * rSges(7,1) + rSges(7,2) * t149) + g(3) * ((-t86 * t130 - t39 * t82) * rSges(7,1) + (t82 * t130 - t39 * t86) * rSges(7,2)))];
taug  = t3(:);
