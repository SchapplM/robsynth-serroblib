% Calculate Gravitation load on the joints for
% S6RRRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
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
% Datum: 2019-03-09 22:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPR5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR5_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:12:07
% EndTime: 2019-03-09 22:12:10
% DurationCPUTime: 1.10s
% Computational Cost: add. (639->182), mult. (844->226), div. (0->0), fcn. (864->10), ass. (0->90)
t50 = qJ(2) + qJ(3);
t47 = sin(t50);
t51 = sin(qJ(6));
t52 = sin(qJ(4));
t55 = cos(qJ(6));
t56 = cos(qJ(4));
t76 = t51 * t56 - t52 * t55;
t10 = t76 * t47;
t75 = t51 * t52 + t55 * t56;
t125 = t75 * t47;
t130 = -rSges(7,1) * t125 + rSges(7,2) * t10;
t54 = sin(qJ(1));
t101 = t52 * t54;
t48 = cos(t50);
t58 = cos(qJ(1));
t99 = t56 * t58;
t17 = t48 * t101 + t99;
t100 = t54 * t56;
t98 = t58 * t52;
t18 = t48 * t100 - t98;
t77 = t17 * t51 + t18 * t55;
t87 = -t17 * t55 + t18 * t51;
t129 = -t87 * rSges(7,1) - t77 * rSges(7,2);
t128 = t130 * t54;
t127 = t130 * t58;
t126 = g(3) * t47;
t110 = -rSges(7,3) - pkin(10);
t104 = t48 * t54;
t107 = t47 * t52;
t92 = rSges(5,2) * t107;
t122 = rSges(5,3) * t104 + t54 * t92;
t102 = t48 * t58;
t121 = rSges(5,3) * t102 + t58 * t92;
t120 = t48 * rSges(4,1) - rSges(4,2) * t47;
t119 = g(1) * t58 + g(2) * t54;
t118 = t119 * t47;
t117 = -pkin(4) - pkin(5);
t53 = sin(qJ(2));
t116 = pkin(2) * t53;
t59 = -pkin(8) - pkin(7);
t113 = g(2) * t59;
t42 = t47 * pkin(9);
t43 = t48 * pkin(3);
t112 = -rSges(6,1) - pkin(4);
t111 = rSges(3,3) + pkin(7);
t40 = t47 * rSges(6,2);
t39 = t47 * rSges(5,3);
t106 = t47 * t58;
t105 = t48 * t52;
t103 = t48 * t56;
t97 = t58 * t59;
t96 = rSges(4,3) - t59;
t95 = t42 + t43;
t94 = qJ(5) * t52;
t93 = rSges(6,3) + qJ(5);
t57 = cos(qJ(2));
t49 = t57 * pkin(2);
t46 = t49 + pkin(1);
t25 = t58 * t46;
t91 = pkin(3) * t102 + pkin(9) * t106 + t25;
t89 = -t46 - t43;
t86 = pkin(4) * t103 + t48 * t94 + t95;
t33 = pkin(9) * t104;
t85 = -t54 * t116 + t33;
t37 = pkin(9) * t102;
t84 = -t58 * t116 + t37;
t19 = t48 * t98 - t100;
t20 = t48 * t99 + t101;
t2 = t19 * t55 - t20 * t51;
t3 = t19 * t51 + t20 * t55;
t83 = rSges(7,1) * t2 - rSges(7,2) * t3;
t82 = rSges(3,1) * t57 - rSges(3,2) * t53;
t79 = -rSges(4,1) * t47 - rSges(4,2) * t48;
t78 = -rSges(7,1) * t10 - rSges(7,2) * t125;
t74 = t89 - t42;
t73 = pkin(1) + t82;
t72 = -t18 * pkin(4) - t17 * qJ(5) - t97;
t71 = pkin(5) * t103 + t86 + (rSges(7,1) * t75 - rSges(7,2) * t76) * t48;
t70 = t20 * pkin(4) + t19 * qJ(5) + t91;
t69 = rSges(6,1) * t103 + rSges(6,3) * t105 + t40 + t86;
t67 = rSges(5,1) * t103 - rSges(5,2) * t105 + t39 + t95;
t62 = (-rSges(5,1) * t56 - pkin(3)) * t118;
t61 = (t112 * t56 - t93 * t52 - pkin(3)) * t118;
t60 = (t117 * t56 - pkin(3) - t94) * t118 + (t119 * t48 + t126) * t110;
t32 = rSges(6,2) * t102;
t27 = rSges(6,2) * t104;
t23 = t47 * t56 * qJ(5);
t15 = t19 * pkin(4);
t13 = t17 * pkin(4);
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t54 - rSges(2,2) * t58) + g(2) * (rSges(2,1) * t58 - rSges(2,2) * t54)) - m(3) * ((g(1) * t111 + g(2) * t73) * t58 + (-g(1) * t73 + g(2) * t111) * t54) - m(4) * (g(2) * t25 + (g(1) * t96 + g(2) * t120) * t58 + (g(1) * (-t46 - t120) + g(2) * t96) * t54) - m(5) * (g(1) * (-t18 * rSges(5,1) + t17 * rSges(5,2) - t97) + g(2) * (t20 * rSges(5,1) - t19 * rSges(5,2) + rSges(5,3) * t106 + t91) + (g(1) * (t74 - t39) - t113) * t54) - m(6) * (g(1) * (-t18 * rSges(6,1) - t17 * rSges(6,3) + t72) + g(2) * (t20 * rSges(6,1) + rSges(6,2) * t106 + t19 * rSges(6,3) + t70) + (g(1) * (t74 - t40) - t113) * t54) - m(7) * (g(1) * (-t77 * rSges(7,1) + t87 * rSges(7,2) - t18 * pkin(5) + t72) + g(2) * (t3 * rSges(7,1) + t2 * rSges(7,2) + t20 * pkin(5) + t106 * t110 + t70) + (-t113 + (t89 + (-pkin(9) - t110) * t47) * g(1)) * t54) -m(3) * (g(3) * t82 + t119 * (-rSges(3,1) * t53 - rSges(3,2) * t57)) - m(4) * (g(3) * (t49 + t120) + t119 * (t79 - t116)) - m(5) * (g(1) * (t84 + t121) + g(2) * (t85 + t122) + g(3) * (t49 + t67) + t62) - m(6) * (g(1) * (t32 + t84) + g(2) * (t27 + t85) + g(3) * (t49 + t69) + t61) - m(7) * (g(1) * (t84 + t127) + g(2) * (t85 + t128) + g(3) * (t49 + t71) + t60) -m(4) * (g(3) * t120 + t119 * t79) - m(5) * (g(1) * (t37 + t121) + g(2) * (t33 + t122) + g(3) * t67 + t62) - m(6) * (g(1) * (t32 + t37) + g(2) * (t27 + t33) + g(3) * t69 + t61) - m(7) * (g(1) * (t37 + t127) + g(2) * (t33 + t128) + g(3) * t71 + t60) -m(5) * (g(1) * (-rSges(5,1) * t19 - rSges(5,2) * t20) + g(2) * (-rSges(5,1) * t17 - rSges(5,2) * t18)) - m(6) * (g(1) * (-rSges(6,1) * t19 + t93 * t20 - t15) + g(2) * (-rSges(6,1) * t17 + t93 * t18 - t13) + g(3) * t23) - m(7) * (g(1) * (-t19 * pkin(5) + t20 * qJ(5) - t15 - t83) + g(2) * (-t17 * pkin(5) + t18 * qJ(5) - t129 - t13) + g(3) * (t23 - t78)) + ((m(5) * rSges(5,2) - m(6) * rSges(6,3)) * t56 + (m(5) * rSges(5,1) - m(6) * t112 - m(7) * t117) * t52) * t126 (-m(6) - m(7)) * (g(1) * t19 + g(2) * t17 + g(3) * t107) -m(7) * (g(1) * t83 + g(2) * t129 + g(3) * t78)];
taug  = t1(:);
