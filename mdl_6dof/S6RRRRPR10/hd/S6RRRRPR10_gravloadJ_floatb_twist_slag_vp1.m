% Calculate Gravitation load on the joints for
% S6RRRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
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
% Datum: 2019-03-09 23:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPR10_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR10_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR10_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR10_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR10_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR10_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:00:59
% EndTime: 2019-03-09 23:01:01
% DurationCPUTime: 1.29s
% Computational Cost: add. (881->208), mult. (1460->287), div. (0->0), fcn. (1712->12), ass. (0->89)
t144 = pkin(11) + rSges(7,3);
t77 = sin(qJ(6));
t81 = cos(qJ(6));
t154 = t77 * rSges(7,1) + t81 * rSges(7,2);
t123 = cos(pkin(6));
t76 = sin(pkin(6));
t79 = sin(qJ(2));
t136 = t76 * t79;
t78 = sin(qJ(3));
t82 = cos(qJ(3));
t156 = t123 * t82 - t78 * t136;
t134 = t76 * t82;
t110 = t79 * t123;
t142 = cos(qJ(2));
t143 = cos(qJ(1));
t80 = sin(qJ(1));
t55 = -t80 * t110 + t143 * t142;
t28 = t80 * t134 - t55 * t78;
t75 = qJ(3) + qJ(4);
t72 = sin(t75);
t125 = qJ(5) * t72;
t73 = cos(t75);
t155 = -pkin(4) * t73 - t125;
t116 = t76 * t143;
t53 = t143 * t110 + t80 * t142;
t22 = t73 * t116 + t53 * t72;
t100 = t123 * t142;
t52 = -t143 * t100 + t79 * t80;
t153 = -t22 * t77 - t52 * t81;
t152 = t22 * t81 - t52 * t77;
t23 = -t72 * t116 + t53 * t73;
t151 = t22 * rSges(6,2) + t23 * rSges(6,3);
t150 = rSges(4,1) * t82 - rSges(4,2) * t78 + pkin(2);
t149 = -t144 * t22 + t154 * t23;
t147 = g(3) * t76;
t145 = -pkin(9) - rSges(4,3);
t135 = t76 * t80;
t131 = -t22 * rSges(5,1) - t23 * rSges(5,2);
t26 = -t73 * t135 + t55 * t72;
t27 = t72 * t135 + t55 * t73;
t130 = -t26 * rSges(5,1) - t27 * rSges(5,2);
t71 = pkin(3) * t82 + pkin(2);
t83 = -pkin(10) - pkin(9);
t129 = -t52 * t71 - t53 * t83;
t54 = t80 * t100 + t143 * t79;
t128 = -t54 * t71 - t55 * t83;
t46 = -t123 * t73 + t72 * t136;
t47 = t123 * t72 + t73 * t136;
t127 = -t46 * rSges(5,1) - t47 * rSges(5,2);
t126 = t143 * pkin(1) + pkin(8) * t135;
t124 = rSges(6,3) + qJ(5);
t122 = t78 * t135;
t115 = t76 * t142;
t103 = t73 * t115;
t57 = t71 * t115;
t119 = pkin(4) * t103 + t115 * t125 + t57;
t118 = t72 * t142;
t117 = t73 * t142;
t114 = t77 * t142;
t113 = t81 * t142;
t112 = -t80 * pkin(1) + pkin(8) * t116;
t65 = t78 * t116;
t111 = -t53 * t82 + t65;
t108 = -t22 * pkin(4) + t23 * qJ(5);
t107 = -t26 * pkin(4) + qJ(5) * t27;
t106 = -t46 * pkin(4) + qJ(5) * t47;
t105 = t155 * t52 + t129;
t104 = t155 * t54 + t128;
t102 = t28 * pkin(3);
t101 = pkin(3) * t122 - t54 * t83 + t55 * t71 + t126;
t99 = -rSges(5,1) * t73 + rSges(5,2) * t72;
t98 = rSges(6,2) * t73 - rSges(6,3) * t72;
t97 = t156 * pkin(3);
t96 = t26 * rSges(6,2) + t27 * rSges(6,3) + t107;
t95 = t46 * rSges(6,2) + t47 * rSges(6,3) + t106;
t94 = t27 * pkin(4) + t101;
t92 = t81 * rSges(7,1) - t77 * rSges(7,2) + pkin(5);
t91 = pkin(3) * t65 + t52 * t83 - t53 * t71 + t112;
t90 = -pkin(4) * t23 + t91;
t89 = t82 * t116 + t53 * t78;
t88 = -t144 * t26 + t154 * t27 + t107;
t87 = t89 * pkin(3);
t86 = -t144 * t46 + t154 * t47 + t106;
t85 = -t144 * t73 - t154 * t72;
t84 = t108 - t87;
t29 = t55 * t82 + t122;
t3 = t26 * t77 + t54 * t81;
t2 = t26 * t81 - t54 * t77;
t1 = [-m(2) * (g(1) * (-t80 * rSges(2,1) - t143 * rSges(2,2)) + g(2) * (t143 * rSges(2,1) - t80 * rSges(2,2))) - m(3) * (g(1) * (-t53 * rSges(3,1) + t52 * rSges(3,2) + rSges(3,3) * t116 + t112) + g(2) * (rSges(3,1) * t55 - rSges(3,2) * t54 + rSges(3,3) * t135 + t126)) - m(4) * (g(1) * (t111 * rSges(4,1) + t89 * rSges(4,2) - t53 * pkin(2) + t145 * t52 + t112) + g(2) * (rSges(4,1) * t29 + rSges(4,2) * t28 + pkin(2) * t55 - t145 * t54 + t126)) - m(5) * (g(1) * (-rSges(5,1) * t23 + rSges(5,2) * t22 - rSges(5,3) * t52 + t91) + g(2) * (rSges(5,1) * t27 - rSges(5,2) * t26 + rSges(5,3) * t54 + t101)) - m(6) * (g(1) * (-rSges(6,1) * t52 + rSges(6,2) * t23 - t124 * t22 + t90) + g(2) * (rSges(6,1) * t54 - rSges(6,2) * t27 + t124 * t26 + t94)) - m(7) * (g(1) * (rSges(7,1) * t153 - rSges(7,2) * t152 - t52 * pkin(5) - t22 * qJ(5) - t144 * t23 + t90) + g(2) * (rSges(7,1) * t3 + rSges(7,2) * t2 + pkin(5) * t54 + qJ(5) * t26 + t144 * t27 + t94)) -m(3) * (g(1) * (-rSges(3,1) * t54 - rSges(3,2) * t55) + g(2) * (-rSges(3,1) * t52 - rSges(3,2) * t53) + (t142 * rSges(3,1) - rSges(3,2) * t79) * t147) - m(4) * (g(1) * (-t145 * t55 - t150 * t54) + (-t145 * t53 - t150 * t52) * g(2) + (t150 * t142 - t145 * t79) * t147) - m(5) * (g(1) * (rSges(5,3) * t55 + t99 * t54 + t128) + g(2) * (rSges(5,3) * t53 + t99 * t52 + t129) + g(3) * t57 + (rSges(5,1) * t117 - rSges(5,2) * t118 + (rSges(5,3) - t83) * t79) * t147) - m(6) * (g(1) * (rSges(6,1) * t55 + t98 * t54 + t104) + g(2) * (rSges(6,1) * t53 + t98 * t52 + t105) + g(3) * t119 + (-rSges(6,2) * t117 + rSges(6,3) * t118 + (rSges(6,1) - t83) * t79) * t147) - m(7) * (g(1) * (t85 * t54 + t92 * t55 + t104) + g(2) * (t85 * t52 + t92 * t53 + t105) + g(3) * (t119 + ((t72 * t114 + t79 * t81) * rSges(7,1) + (t72 * t113 - t77 * t79) * rSges(7,2)) * t76 + (pkin(5) - t83) * t136 + t144 * t103)) -m(4) * (g(1) * (rSges(4,1) * t28 - rSges(4,2) * t29) + g(2) * (-t89 * rSges(4,1) + t111 * rSges(4,2)) + g(3) * (t156 * rSges(4,1) + (-t123 * t78 - t79 * t134) * rSges(4,2))) - m(5) * (g(1) * (t102 + t130) + g(2) * (-t87 + t131) + g(3) * (t97 + t127)) - m(6) * (g(1) * (t102 + t96) + g(2) * (t84 + t151) + g(3) * (t95 + t97)) - m(7) * (g(1) * (t102 + t88) + g(2) * (t84 + t149) + g(3) * (t86 + t97)) -m(5) * (g(1) * t130 + g(2) * t131 + g(3) * t127) - m(6) * (g(1) * t96 + g(2) * (t108 + t151) + g(3) * t95) - m(7) * (g(1) * t88 + g(2) * (t108 + t149) + g(3) * t86) (-m(6) - m(7)) * (g(1) * t26 + g(2) * t22 + g(3) * t46) -m(7) * (g(1) * (rSges(7,1) * t2 - rSges(7,2) * t3) + g(2) * (rSges(7,1) * t152 + rSges(7,2) * t153) + g(3) * ((t76 * t114 + t46 * t81) * rSges(7,1) + (t76 * t113 - t46 * t77) * rSges(7,2)))];
taug  = t1(:);
