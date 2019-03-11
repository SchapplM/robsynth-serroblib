% Calculate Gravitation load on the joints for
% S6PRRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-03-09 01:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRRR4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR4_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRR4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:54:52
% EndTime: 2019-03-09 00:54:55
% DurationCPUTime: 1.37s
% Computational Cost: add. (1223->222), mult. (2870->343), div. (0->0), fcn. (3624->16), ass. (0->102)
t146 = pkin(12) + rSges(7,3);
t120 = sin(pkin(6));
t143 = sin(qJ(2));
t107 = t120 * t143;
t144 = cos(qJ(3));
t122 = cos(pkin(7));
t103 = t122 * t120;
t123 = cos(pkin(6));
t145 = cos(qJ(2));
t79 = sin(pkin(7));
t159 = t145 * t103 + t123 * t79;
t82 = sin(qJ(3));
t49 = t144 * t107 + t159 * t82;
t108 = t145 * t120;
t63 = -t108 * t79 + t122 * t123;
t81 = sin(qJ(4));
t84 = cos(qJ(4));
t165 = -t49 * t81 + t63 * t84;
t119 = sin(pkin(13));
t101 = t120 * t119;
t104 = t123 * t119;
t121 = cos(pkin(13));
t89 = t104 * t145 + t121 * t143;
t160 = -t79 * t101 + t89 * t122;
t65 = -t104 * t143 + t121 * t145;
t33 = t65 * t144 - t160 * t82;
t51 = t101 * t122 + t79 * t89;
t164 = -t33 * t81 + t51 * t84;
t102 = t121 * t120;
t105 = t123 * t121;
t88 = -t105 * t145 + t119 * t143;
t161 = t79 * t102 + t88 * t122;
t64 = t105 * t143 + t119 * t145;
t31 = t144 * t64 - t161 * t82;
t50 = -t102 * t122 + t79 * t88;
t163 = -t31 * t81 + t50 * t84;
t80 = sin(qJ(6));
t83 = cos(qJ(6));
t162 = t83 * rSges(7,1) - t80 * rSges(7,2) + pkin(5);
t158 = g(1) * t65 + g(2) * t64;
t30 = t144 * t161 + t64 * t82;
t32 = t144 * t160 + t65 * t82;
t48 = t107 * t82 - t144 * t159;
t157 = g(1) * t32 + g(2) * t30 + g(3) * t48;
t78 = qJ(4) + qJ(5);
t76 = sin(t78);
t77 = cos(t78);
t11 = -t31 * t76 + t50 * t77;
t12 = t31 * t77 + t50 * t76;
t156 = t11 * rSges(6,1) - t12 * rSges(6,2);
t13 = -t33 * t76 + t51 * t77;
t14 = t33 * t77 + t51 * t76;
t155 = t13 * rSges(6,1) - t14 * rSges(6,2);
t150 = t77 * pkin(5);
t149 = t79 * pkin(9);
t147 = pkin(10) + rSges(5,3);
t136 = t76 * t79;
t135 = t77 * t79;
t134 = t77 * t80;
t133 = t77 * t83;
t132 = t79 * t81;
t131 = t79 * t84;
t75 = pkin(4) * t84 + pkin(3);
t85 = -pkin(11) - pkin(10);
t128 = -t30 * t75 - t31 * t85;
t127 = -t32 * t75 - t33 * t85;
t26 = -t49 * t76 + t63 * t77;
t27 = t49 * t77 + t63 * t76;
t126 = t26 * rSges(6,1) - t27 * rSges(6,2);
t125 = -t48 * t75 - t49 * t85;
t99 = t79 * t107;
t124 = pkin(2) * t108 + pkin(9) * t99;
t109 = t122 * t144;
t38 = t109 * t64 - t82 * t88;
t115 = t82 * t122;
t39 = -t115 * t64 - t144 * t88;
t61 = t88 * pkin(2);
t118 = -t38 * t85 + t39 * t75 - t61;
t40 = t109 * t65 - t82 * t89;
t41 = -t115 * t65 - t144 * t89;
t62 = t89 * pkin(2);
t117 = -t40 * t85 + t41 * t75 - t62;
t113 = t163 * pkin(4);
t112 = t164 * pkin(4);
t111 = t165 * pkin(4);
t94 = t143 * t103;
t58 = t108 * t82 + t144 * t94;
t59 = t108 * t144 - t82 * t94;
t96 = t81 * t99;
t110 = pkin(4) * t96 - t58 * t85 + t59 * t75 + t124;
t106 = -rSges(6,1) * t77 + rSges(6,2) * t76;
t93 = t162 * t11 + t12 * t146;
t92 = t162 * t13 + t14 * t146;
t91 = t146 * t27 + t162 * t26;
t90 = t158 * (pkin(4) * t81 + pkin(9)) * t79;
t43 = t59 * t77 + t76 * t99;
t42 = t59 * t76 - t77 * t99;
t18 = t136 * t65 + t41 * t77;
t17 = -t135 * t65 + t41 * t76;
t16 = t136 * t64 + t39 * t77;
t15 = -t135 * t64 + t39 * t76;
t1 = [(-m(2) - m(3) - m(4) - m(5) - m(6) - m(7)) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t89 - t65 * rSges(3,2)) + g(2) * (-rSges(3,1) * t88 - t64 * rSges(3,2)) + g(3) * (rSges(3,1) * t108 - rSges(3,2) * t107)) - m(4) * (g(1) * (rSges(4,1) * t41 - rSges(4,2) * t40 - t62) + g(2) * (rSges(4,1) * t39 - rSges(4,2) * t38 - t61) + g(3) * (t59 * rSges(4,1) - t58 * rSges(4,2) + t124) + (rSges(4,3) * g(3) * t107 + t158 * (rSges(4,3) + pkin(9))) * t79) - m(5) * (g(1) * (t41 * pkin(3) - t62 + t65 * t149 + (t132 * t65 + t41 * t84) * rSges(5,1) + (t131 * t65 - t41 * t81) * rSges(5,2) + t147 * t40) + g(2) * (t39 * pkin(3) - t61 + t64 * t149 + (t132 * t64 + t39 * t84) * rSges(5,1) + (t131 * t64 - t39 * t81) * rSges(5,2) + t147 * t38) + g(3) * (t59 * pkin(3) + (t59 * t84 + t96) * rSges(5,1) + (-t59 * t81 + t84 * t99) * rSges(5,2) + t147 * t58 + t124)) - m(6) * (g(1) * (rSges(6,1) * t18 - rSges(6,2) * t17 + rSges(6,3) * t40 + t117) + g(2) * (rSges(6,1) * t16 - rSges(6,2) * t15 + rSges(6,3) * t38 + t118) + g(3) * (rSges(6,1) * t43 - t42 * rSges(6,2) + rSges(6,3) * t58 + t110) + t90) - m(7) * (g(1) * (t18 * pkin(5) + (t18 * t83 + t40 * t80) * rSges(7,1) + (-t18 * t80 + t40 * t83) * rSges(7,2) + t117 + t146 * t17) + g(2) * (t16 * pkin(5) + (t16 * t83 + t38 * t80) * rSges(7,1) + (-t16 * t80 + t38 * t83) * rSges(7,2) + t118 + t146 * t15) + g(3) * (t43 * pkin(5) + (t43 * t83 + t58 * t80) * rSges(7,1) + (-t43 * t80 + t58 * t83) * rSges(7,2) + t146 * t42 + t110) + t90) -m(4) * (g(1) * (-rSges(4,1) * t32 - rSges(4,2) * t33) + g(2) * (-rSges(4,1) * t30 - rSges(4,2) * t31) + g(3) * (-rSges(4,1) * t48 - rSges(4,2) * t49)) - m(5) * ((g(1) * t33 + g(2) * t31 + g(3) * t49) * t147 + t157 * (-t84 * rSges(5,1) + t81 * rSges(5,2) - pkin(3))) - m(6) * (g(1) * (rSges(6,3) * t33 + t106 * t32 + t127) + g(2) * (rSges(6,3) * t31 + t106 * t30 + t128) + g(3) * (rSges(6,3) * t49 + t106 * t48 + t125)) + (-g(1) * (-t32 * t150 + (-t133 * t32 + t33 * t80) * rSges(7,1) + (t134 * t32 + t33 * t83) * rSges(7,2) + t127) - g(2) * (-t30 * t150 + (-t133 * t30 + t31 * t80) * rSges(7,1) + (t134 * t30 + t31 * t83) * rSges(7,2) + t128) - g(3) * (-t48 * t150 + (-t133 * t48 + t49 * t80) * rSges(7,1) + (t134 * t48 + t49 * t83) * rSges(7,2) + t125) + t157 * t76 * t146) * m(7), -m(5) * (g(1) * (t164 * rSges(5,1) + (-t33 * t84 - t51 * t81) * rSges(5,2)) + g(2) * (t163 * rSges(5,1) + (-t31 * t84 - t50 * t81) * rSges(5,2)) + g(3) * (t165 * rSges(5,1) + (-t49 * t84 - t63 * t81) * rSges(5,2))) - m(6) * (g(1) * (t112 + t155) + g(2) * (t113 + t156) + g(3) * (t111 + t126)) - m(7) * (g(1) * (t112 + t92) + g(2) * (t113 + t93) + g(3) * (t111 + t91)) -m(6) * (g(1) * t155 + g(2) * t156 + g(3) * t126) - m(7) * (g(1) * t92 + g(2) * t93 + g(3) * t91) -m(7) * (g(1) * ((-t14 * t80 + t32 * t83) * rSges(7,1) + (-t14 * t83 - t32 * t80) * rSges(7,2)) + g(2) * ((-t12 * t80 + t30 * t83) * rSges(7,1) + (-t12 * t83 - t30 * t80) * rSges(7,2)) + g(3) * ((-t27 * t80 + t48 * t83) * rSges(7,1) + (-t27 * t83 - t48 * t80) * rSges(7,2)))];
taug  = t1(:);
