% Calculate Gravitation load on the joints for
% S6RRRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 22:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPR7_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR7_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR7_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR7_gravloadJ_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR7_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR7_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:26:34
% EndTime: 2019-03-09 22:26:38
% DurationCPUTime: 1.21s
% Computational Cost: add. (951->201), mult. (1305->279), div. (0->0), fcn. (1503->14), ass. (0->91)
t141 = pkin(11) + rSges(7,3);
t120 = cos(pkin(6));
t82 = sin(pkin(6));
t85 = sin(qJ(2));
t134 = t82 * t85;
t81 = qJ(3) + qJ(4);
t76 = sin(t81);
t77 = cos(t81);
t156 = t120 * t77 - t76 * t134;
t86 = sin(qJ(1));
t133 = t82 * t86;
t112 = t86 * t120;
t89 = cos(qJ(2));
t90 = cos(qJ(1));
t55 = -t85 * t112 + t89 * t90;
t24 = t77 * t133 - t55 * t76;
t83 = sin(qJ(6));
t87 = cos(qJ(6));
t155 = rSges(7,1) * t87 - rSges(7,2) * t83 + pkin(5);
t130 = t82 * t90;
t111 = t90 * t120;
t53 = t85 * t111 + t86 * t89;
t75 = pkin(12) + t81;
t71 = sin(t75);
t72 = cos(t75);
t15 = -t71 * t130 + t53 * t72;
t52 = -t89 * t111 + t85 * t86;
t154 = t15 * t83 - t52 * t87;
t153 = -t15 * t87 - t52 * t83;
t131 = t82 * t89;
t152 = g(3) * t131;
t151 = t83 * rSges(7,1) + t87 * rSges(7,2);
t146 = g(2) * t52;
t54 = t89 * t112 + t90 * t85;
t150 = g(1) * t54 + t146;
t149 = t141 * t71 + t155 * t72;
t91 = -pkin(10) - pkin(9);
t115 = -t72 * t130 - t53 * t71;
t148 = rSges(6,1) * t115 - t15 * rSges(6,2);
t145 = g(2) * t53;
t88 = cos(qJ(3));
t78 = t88 * pkin(3);
t63 = pkin(4) * t77 + t78;
t60 = pkin(2) + t63;
t144 = t60 * t152;
t143 = g(3) * t82;
t142 = -pkin(9) - rSges(4,3);
t18 = -t72 * t133 + t55 * t71;
t19 = t71 * t133 + t55 * t72;
t140 = -t18 * rSges(6,1) - t19 * rSges(6,2);
t132 = t82 * t88;
t80 = -qJ(5) + t91;
t127 = -t52 * t60 - t53 * t80;
t126 = -t54 * t60 - t55 * t80;
t39 = t120 * t72 - t71 * t134;
t40 = t120 * t71 + t72 * t134;
t125 = t39 * rSges(6,1) - t40 * rSges(6,2);
t84 = sin(qJ(3));
t62 = pkin(3) * t84 + pkin(4) * t76;
t124 = t63 * t133 - t55 * t62;
t123 = t120 * t63 - t62 * t134;
t122 = t90 * pkin(1) + pkin(8) * t133;
t121 = t91 - rSges(5,3);
t118 = t84 * t133;
t67 = t84 * t130;
t116 = -t86 * pkin(1) + pkin(8) * t130;
t114 = t76 * t130 - t53 * t77;
t113 = -t53 * t88 + t67;
t109 = t24 * pkin(4);
t108 = -t63 * t130 - t53 * t62;
t107 = t62 * t133 - t54 * t80 + t55 * t60 + t122;
t106 = rSges(6,1) * t72 - rSges(6,2) * t71;
t105 = t156 * pkin(4);
t104 = rSges(4,1) * t88 - rSges(4,2) * t84 + pkin(2);
t102 = t77 * t130 + t53 * t76;
t101 = t88 * t130 + t53 * t84;
t26 = t86 * t132 - t55 * t84;
t74 = t78 + pkin(2);
t100 = rSges(5,1) * t77 - rSges(5,2) * t76 + t74;
t99 = t62 * t130 + t52 * t80 - t53 * t60 + t116;
t98 = t102 * pkin(4);
t97 = t120 * t88 - t84 * t134;
t96 = t155 * t115 + t141 * t15;
t95 = t141 * t19 - t155 * t18;
t94 = t141 * t40 + t155 * t39;
t25 = t76 * t133 + t55 * t77;
t93 = g(1) * (t24 * rSges(5,1) - t25 * rSges(5,2)) + g(2) * (-t102 * rSges(5,1) + t114 * rSges(5,2)) + g(3) * (t156 * rSges(5,1) + (-t120 * t76 - t77 * t134) * rSges(5,2));
t27 = t55 * t88 + t118;
t2 = t19 * t87 + t54 * t83;
t1 = -t19 * t83 + t54 * t87;
t3 = [-m(2) * (g(1) * (-rSges(2,1) * t86 - rSges(2,2) * t90) + g(2) * (rSges(2,1) * t90 - rSges(2,2) * t86)) - m(3) * (g(1) * (-rSges(3,1) * t53 + rSges(3,2) * t52 + rSges(3,3) * t130 + t116) + g(2) * (rSges(3,1) * t55 - rSges(3,2) * t54 + rSges(3,3) * t133 + t122)) - m(4) * (g(1) * (t113 * rSges(4,1) + t101 * rSges(4,2) - t53 * pkin(2) + t142 * t52 + t116) + g(2) * (rSges(4,1) * t27 + rSges(4,2) * t26 + pkin(2) * t55 - t142 * t54 + t122)) - m(5) * (g(1) * (t114 * rSges(5,1) + t102 * rSges(5,2) + pkin(3) * t67 + t121 * t52 - t53 * t74 + t116) + g(2) * (t25 * rSges(5,1) + t24 * rSges(5,2) + pkin(3) * t118 - t121 * t54 + t55 * t74 + t122)) - m(6) * (g(1) * (-rSges(6,1) * t15 - rSges(6,2) * t115 - rSges(6,3) * t52 + t99) + g(2) * (rSges(6,1) * t19 - rSges(6,2) * t18 + rSges(6,3) * t54 + t107)) - m(7) * (g(1) * (rSges(7,1) * t153 + t154 * rSges(7,2) - t15 * pkin(5) + t141 * t115 + t99) + g(2) * (rSges(7,1) * t2 + rSges(7,2) * t1 + pkin(5) * t19 + t141 * t18 + t107)) -m(3) * (g(1) * (-rSges(3,1) * t54 - rSges(3,2) * t55) + g(2) * (-rSges(3,1) * t52 - rSges(3,2) * t53) + (rSges(3,1) * t89 - rSges(3,2) * t85) * t143) - m(4) * (g(1) * (-t104 * t54 - t142 * t55) - t142 * t145 - t104 * t146 + (t104 * t89 - t142 * t85) * t143) - m(5) * (g(1) * (-t100 * t54 - t121 * t55) - t121 * t145 - t100 * t146 + (t100 * t89 - t121 * t85) * t143) - m(6) * (g(1) * (rSges(6,3) * t55 - t106 * t54 + t126) + g(2) * (rSges(6,3) * t53 - t106 * t52 + t127) + t144 + (t106 * t89 + (rSges(6,3) - t80) * t85) * t143) - m(7) * (g(1) * (t151 * t55 + t126) + g(2) * (t151 * t53 + t127) + t144 + ((-t80 + t151) * t85 + t149 * t89) * t143 - t150 * t149) -m(4) * (g(1) * (rSges(4,1) * t26 - rSges(4,2) * t27) + g(2) * (-t101 * rSges(4,1) + t113 * rSges(4,2)) + g(3) * (t97 * rSges(4,1) + (-t120 * t84 - t85 * t132) * rSges(4,2))) - m(5) * ((g(1) * t26 - g(2) * t101 + g(3) * t97) * pkin(3) + t93) - m(6) * (g(1) * (t124 + t140) + g(2) * (t108 + t148) + g(3) * (t123 + t125)) - m(7) * (g(1) * (t95 + t124) + g(2) * (t108 + t96) + g(3) * (t94 + t123)) -m(5) * t93 - m(6) * (g(1) * (t109 + t140) + g(2) * (-t98 + t148) + g(3) * (t105 + t125)) - m(7) * (g(1) * (t109 + t95) + g(2) * (-t98 + t96) + g(3) * (t105 + t94)) (-m(6) - m(7)) * (t150 - t152) -m(7) * (g(1) * (rSges(7,1) * t1 - rSges(7,2) * t2) + g(2) * (-t154 * rSges(7,1) + rSges(7,2) * t153) + g(3) * ((-t87 * t131 - t40 * t83) * rSges(7,1) + (t83 * t131 - t40 * t87) * rSges(7,2)))];
taug  = t3(:);
