% Calculate Gravitation load on the joints for
% S5RRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d2,d3,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRR12_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR12_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR12_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR12_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR12_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR12_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:46:36
% EndTime: 2019-12-31 22:46:44
% DurationCPUTime: 1.67s
% Computational Cost: add. (927->209), mult. (2503->329), div. (0->0), fcn. (3167->14), ass. (0->89)
t71 = sin(pkin(5));
t77 = cos(qJ(1));
t113 = t71 * t77;
t70 = sin(pkin(6));
t106 = t70 * t113;
t119 = cos(qJ(3));
t117 = sin(qJ(2));
t118 = sin(qJ(1));
t108 = cos(pkin(5));
t120 = cos(qJ(2));
t94 = t108 * t120;
t53 = t117 * t118 - t77 * t94;
t93 = t108 * t117;
t54 = t118 * t120 + t77 * t93;
t74 = sin(qJ(3));
t107 = cos(pkin(6));
t92 = t107 * t119;
t19 = t106 * t119 + t53 * t92 + t54 * t74;
t99 = t74 * t107;
t20 = -t74 * t106 + t119 * t54 - t53 * t99;
t100 = t71 * t107;
t39 = -t77 * t100 + t53 * t70;
t73 = sin(qJ(4));
t76 = cos(qJ(4));
t4 = t20 * t76 + t39 * t73;
t72 = sin(qJ(5));
t75 = cos(qJ(5));
t131 = -t19 * t75 + t4 * t72;
t130 = -t19 * t72 - t4 * t75;
t83 = t117 * t77 + t118 * t94;
t78 = t118 * t100 + t83 * t70;
t127 = -t20 * t73 + t39 * t76;
t103 = t71 * t118;
t126 = -t70 * t103 + t83 * t107;
t125 = pkin(9) * t70;
t124 = t76 * pkin(4);
t122 = rSges(5,3) + pkin(10);
t121 = rSges(6,3) + pkin(11);
t115 = t70 * t73;
t114 = t70 * t76;
t112 = t72 * t76;
t111 = t75 * t76;
t104 = t71 * t120;
t102 = t71 * t117;
t98 = t70 * t102;
t110 = pkin(2) * t104 + pkin(9) * t98;
t109 = t77 * pkin(1) + pkin(8) * t103;
t91 = t107 * t117;
t47 = (t119 * t120 - t74 * t91) * t71;
t105 = t47 * pkin(3) + t110;
t101 = t70 * t108;
t95 = -pkin(1) * t118 + pkin(8) * t113;
t90 = -rSges(5,1) * t76 + rSges(5,2) * t73;
t28 = -t119 * t53 - t54 * t99;
t48 = t53 * pkin(2);
t89 = t28 * pkin(3) + t125 * t54 - t48;
t55 = -t118 * t93 + t120 * t77;
t30 = -t119 * t83 - t55 * t99;
t50 = t83 * pkin(2);
t88 = t30 * pkin(3) + t125 * t55 - t50;
t87 = t75 * rSges(6,1) - t72 * rSges(6,2) + pkin(4);
t85 = -t54 * pkin(2) - t39 * pkin(9) + t95;
t84 = -pkin(3) * t20 + t85;
t80 = t55 * pkin(2) + t78 * pkin(9) + t109;
t24 = t55 * t119 - t126 * t74;
t79 = t24 * pkin(3) + t80;
t52 = -t104 * t70 + t107 * t108;
t46 = (t119 * t91 + t120 * t74) * t71;
t37 = t74 * t101 + (t117 * t119 + t120 * t99) * t71;
t36 = -t101 * t119 + t102 * t74 - t104 * t92;
t35 = t36 * pkin(3);
t32 = t47 * t76 + t73 * t98;
t31 = t47 * t73 - t76 * t98;
t29 = t55 * t92 - t74 * t83;
t27 = -t53 * t74 + t54 * t92;
t23 = t119 * t126 + t55 * t74;
t18 = t37 * t76 + t52 * t73;
t17 = -t37 * t73 + t52 * t76;
t15 = t23 * pkin(3);
t13 = t19 * pkin(3);
t12 = t115 * t55 + t30 * t76;
t11 = -t114 * t55 + t30 * t73;
t10 = t115 * t54 + t28 * t76;
t9 = -t114 * t54 + t28 * t73;
t8 = t24 * t76 + t73 * t78;
t7 = t24 * t73 - t76 * t78;
t2 = t23 * t72 + t75 * t8;
t1 = t23 * t75 - t72 * t8;
t3 = [-m(2) * (g(1) * (-rSges(2,1) * t118 - t77 * rSges(2,2)) + g(2) * (t77 * rSges(2,1) - rSges(2,2) * t118)) - m(3) * (g(1) * (-t54 * rSges(3,1) + t53 * rSges(3,2) + rSges(3,3) * t113 + t95) + g(2) * (t55 * rSges(3,1) - rSges(3,2) * t83 + rSges(3,3) * t103 + t109)) - m(4) * (g(1) * (-rSges(4,1) * t20 + rSges(4,2) * t19 - rSges(4,3) * t39 + t85) + g(2) * (t24 * rSges(4,1) - t23 * rSges(4,2) + rSges(4,3) * t78 + t80)) - m(5) * (g(1) * (-rSges(5,1) * t4 - rSges(5,2) * t127 - t122 * t19 + t84) + g(2) * (t8 * rSges(5,1) - t7 * rSges(5,2) + t122 * t23 + t79)) - m(6) * (g(1) * (t130 * rSges(6,1) + t131 * rSges(6,2) - t4 * pkin(4) - t19 * pkin(10) + t121 * t127 + t84) + g(2) * (t2 * rSges(6,1) + t1 * rSges(6,2) + t8 * pkin(4) + t23 * pkin(10) + t121 * t7 + t79)), -m(3) * (g(1) * (-rSges(3,1) * t83 - t55 * rSges(3,2)) + g(2) * (-rSges(3,1) * t53 - rSges(3,2) * t54) + g(3) * (rSges(3,1) * t120 - rSges(3,2) * t117) * t71) - m(4) * (g(1) * (t30 * rSges(4,1) - t29 * rSges(4,2) - t50) + g(2) * (t28 * rSges(4,1) - t27 * rSges(4,2) - t48) + g(3) * (t47 * rSges(4,1) - t46 * rSges(4,2) + t110) + (rSges(4,3) * g(3) * t102 + (g(1) * t55 + g(2) * t54) * (rSges(4,3) + pkin(9))) * t70) - m(5) * (g(1) * (t12 * rSges(5,1) - t11 * rSges(5,2) + t122 * t29 + t88) + g(2) * (t10 * rSges(5,1) - t9 * rSges(5,2) + t122 * t27 + t89) + g(3) * (t32 * rSges(5,1) - t31 * rSges(5,2) + t122 * t46 + t105)) - m(6) * (g(1) * (t12 * pkin(4) + t29 * pkin(10) + (t12 * t75 + t29 * t72) * rSges(6,1) + (-t12 * t72 + t29 * t75) * rSges(6,2) + t121 * t11 + t88) + g(2) * (t10 * pkin(4) + t27 * pkin(10) + (t10 * t75 + t27 * t72) * rSges(6,1) + (-t10 * t72 + t27 * t75) * rSges(6,2) + t121 * t9 + t89) + g(3) * (t32 * pkin(4) + t46 * pkin(10) + (t32 * t75 + t46 * t72) * rSges(6,1) + (-t32 * t72 + t46 * t75) * rSges(6,2) + t121 * t31 + t105)), -m(4) * (g(1) * (-rSges(4,1) * t23 - rSges(4,2) * t24) + g(2) * (-rSges(4,1) * t19 - rSges(4,2) * t20) + g(3) * (-rSges(4,1) * t36 - rSges(4,2) * t37)) - m(5) * (g(1) * (t122 * t24 + t23 * t90 - t15) + g(2) * (t122 * t20 + t19 * t90 - t13) + g(3) * (t122 * t37 + t36 * t90 - t35)) + (-g(1) * (-t23 * t124 - t15 + t24 * pkin(10) + (-t111 * t23 + t24 * t72) * rSges(6,1) + (t112 * t23 + t24 * t75) * rSges(6,2)) - g(2) * (-t19 * t124 - t13 + t20 * pkin(10) + (-t111 * t19 + t20 * t72) * rSges(6,1) + (t112 * t19 + t20 * t75) * rSges(6,2)) - g(3) * (-t36 * t124 - t35 + t37 * pkin(10) + (-t111 * t36 + t37 * t72) * rSges(6,1) + (t112 * t36 + t37 * t75) * rSges(6,2)) - (-g(1) * t23 - g(2) * t19 - g(3) * t36) * t73 * t121) * m(6), -m(5) * (g(1) * (-rSges(5,1) * t7 - rSges(5,2) * t8) + g(2) * (rSges(5,1) * t127 - rSges(5,2) * t4) + g(3) * (rSges(5,1) * t17 - rSges(5,2) * t18)) - m(6) * (g(1) * (t121 * t8 - t7 * t87) + (t121 * t18 + t87 * t17) * g(3) + (t121 * t4 + t127 * t87) * g(2)), -m(6) * (g(1) * (rSges(6,1) * t1 - rSges(6,2) * t2) + g(2) * (-t131 * rSges(6,1) + t130 * rSges(6,2)) + g(3) * ((-t18 * t72 + t36 * t75) * rSges(6,1) + (-t18 * t75 - t36 * t72) * rSges(6,2)))];
taug = t3(:);
