% Calculate Gravitation load on the joints for
% S6RRRPPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-03-09 16:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPPR9_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR9_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR9_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR9_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR9_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPR9_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:11:39
% EndTime: 2019-03-09 16:11:43
% DurationCPUTime: 1.66s
% Computational Cost: add. (814->217), mult. (2024->314), div. (0->0), fcn. (2494->12), ass. (0->92)
t130 = cos(qJ(1));
t73 = sin(pkin(6));
t105 = t73 * t130;
t129 = cos(qJ(3));
t128 = sin(qJ(1));
t77 = sin(qJ(2));
t79 = cos(qJ(2));
t114 = cos(pkin(6));
t94 = t114 * t130;
t56 = t128 * t79 + t77 * t94;
t76 = sin(qJ(3));
t30 = -t76 * t105 + t129 * t56;
t55 = t128 * t77 - t79 * t94;
t72 = sin(pkin(11));
t74 = cos(pkin(11));
t7 = t30 * t72 - t55 * t74;
t75 = sin(qJ(6));
t78 = cos(qJ(6));
t8 = t30 * t74 + t55 * t72;
t143 = t7 * t78 - t75 * t8;
t142 = -t7 * t75 - t78 * t8;
t141 = -pkin(3) * t129 - qJ(4) * t76;
t140 = -pkin(4) * t74 - qJ(5) * t72;
t104 = t73 * t129;
t29 = t104 * t130 + t56 * t76;
t103 = t73 * t128;
t93 = t114 * t128;
t58 = t130 * t79 - t77 * t93;
t33 = -t103 * t129 + t58 * t76;
t123 = t73 * t77;
t53 = -t114 * t129 + t123 * t76;
t139 = g(1) * t33 + g(2) * t29 + g(3) * t53;
t138 = -m(6) - m(7);
t133 = g(3) * t73;
t132 = rSges(4,3) + pkin(9);
t131 = rSges(7,3) + pkin(10);
t125 = t55 * t76;
t57 = t130 * t77 + t79 * t93;
t124 = t57 * t76;
t122 = t73 * t79;
t121 = t130 * pkin(1) + pkin(8) * t103;
t120 = pkin(2) * t122 + pkin(9) * t123;
t117 = rSges(6,2) + qJ(4);
t116 = rSges(5,3) + qJ(4);
t115 = rSges(6,3) + qJ(5);
t113 = t76 * t122;
t112 = qJ(4) - t131;
t23 = t29 * pkin(3);
t111 = t140 * t29 - t23;
t25 = t33 * pkin(3);
t110 = t140 * t33 - t25;
t48 = t53 * pkin(3);
t109 = t140 * t53 - t48;
t108 = t58 * pkin(2) + t121;
t106 = t72 * t129;
t102 = t74 * t129;
t101 = t79 * t129;
t97 = t73 * t101;
t100 = pkin(3) * t97 + qJ(4) * t113 + t120;
t99 = g(1) * t112;
t98 = g(2) * t112;
t96 = -pkin(1) * t128 + pkin(8) * t105;
t37 = (t101 * t74 + t72 * t77) * t73;
t95 = t37 * pkin(4) + t100;
t92 = -rSges(5,1) * t74 + rSges(5,2) * t72;
t91 = -rSges(6,1) * t74 - rSges(6,3) * t72;
t49 = t55 * pkin(2);
t90 = t56 * pkin(9) + t141 * t55 - t49;
t51 = t57 * pkin(2);
t89 = t58 * pkin(9) + t141 * t57 - t51;
t88 = -t56 * pkin(2) + t96;
t16 = -t102 * t55 + t56 * t72;
t87 = t16 * pkin(4) + t90;
t18 = -t102 * t57 + t58 * t72;
t86 = t18 * pkin(4) + t89;
t34 = t103 * t76 + t129 * t58;
t85 = t34 * pkin(3) + pkin(9) * t57 + t108;
t84 = rSges(4,1) * t129 - rSges(4,2) * t76;
t12 = t34 * t74 + t57 * t72;
t83 = t12 * pkin(4) + t85;
t82 = -pkin(3) * t30 - t55 * pkin(9) + t88;
t81 = -pkin(4) * t8 + t82;
t54 = t104 * t77 + t114 * t76;
t36 = -t123 * t74 + t72 * t97;
t28 = -t122 * t72 + t54 * t74;
t27 = t122 * t74 + t54 * t72;
t17 = -t106 * t57 - t58 * t74;
t15 = -t106 * t55 - t56 * t74;
t11 = t34 * t72 - t57 * t74;
t3 = t11 * t75 + t12 * t78;
t2 = t11 * t78 - t12 * t75;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t128 - rSges(2,2) * t130) + g(2) * (rSges(2,1) * t130 - rSges(2,2) * t128)) - m(3) * (g(1) * (-t56 * rSges(3,1) + t55 * rSges(3,2) + rSges(3,3) * t105 + t96) + g(2) * (t58 * rSges(3,1) - t57 * rSges(3,2) + rSges(3,3) * t103 + t121)) - m(4) * (g(1) * (-rSges(4,1) * t30 + rSges(4,2) * t29 - t132 * t55 + t88) + g(2) * (rSges(4,1) * t34 - rSges(4,2) * t33 + t132 * t57 + t108)) - m(5) * (g(1) * (-rSges(5,1) * t8 + rSges(5,2) * t7 - t116 * t29 + t82) + g(2) * (rSges(5,1) * t12 - rSges(5,2) * t11 + t116 * t33 + t85)) - m(6) * (g(1) * (-rSges(6,1) * t8 - t115 * t7 - t117 * t29 + t81) + g(2) * (rSges(6,1) * t12 + t11 * t115 + t117 * t33 + t83)) - m(7) * (g(1) * (t142 * rSges(7,1) - t143 * rSges(7,2) - t8 * pkin(5) - t7 * qJ(5) + t81) + g(2) * (rSges(7,1) * t3 + rSges(7,2) * t2 + pkin(5) * t12 + qJ(5) * t11 + t83) + t33 * t98 - t29 * t99) -m(3) * (g(1) * (-rSges(3,1) * t57 - rSges(3,2) * t58) + g(2) * (-rSges(3,1) * t55 - rSges(3,2) * t56) + (rSges(3,1) * t79 - rSges(3,2) * t77) * t133) - m(4) * (g(1) * (t132 * t58 - t84 * t57 - t51) + g(2) * (t132 * t56 - t55 * t84 - t49) + g(3) * t120 + (rSges(4,3) * t77 + t84 * t79) * t133) - m(5) * (g(1) * (rSges(5,1) * t18 - rSges(5,2) * t17 - rSges(5,3) * t124 + t89) + g(2) * (rSges(5,1) * t16 - rSges(5,2) * t15 - rSges(5,3) * t125 + t90) + g(3) * (t37 * rSges(5,1) - t36 * rSges(5,2) + rSges(5,3) * t113 + t100)) - m(6) * (g(1) * (rSges(6,1) * t18 - rSges(6,2) * t124 + t115 * t17 + t86) + g(2) * (rSges(6,1) * t16 - rSges(6,2) * t125 + t115 * t15 + t87) + g(3) * (t37 * rSges(6,1) + rSges(6,2) * t113 + t115 * t36 + t95)) - m(7) * (g(1) * (t18 * pkin(5) + t17 * qJ(5) + (t17 * t75 + t18 * t78) * rSges(7,1) + (t17 * t78 - t18 * t75) * rSges(7,2) + t86) + g(2) * (t16 * pkin(5) + t15 * qJ(5) + (t15 * t75 + t16 * t78) * rSges(7,1) + (t15 * t78 - t16 * t75) * rSges(7,2) + t87) + g(3) * (t37 * pkin(5) + t36 * qJ(5) + (t36 * t75 + t37 * t78) * rSges(7,1) + (t36 * t78 - t37 * t75) * rSges(7,2) + t95) + (g(1) * t57 + g(2) * t55 - g(3) * t122) * t76 * t131) -m(4) * (g(1) * (-rSges(4,1) * t33 - rSges(4,2) * t34) + g(2) * (-rSges(4,1) * t29 - rSges(4,2) * t30) + g(3) * (-rSges(4,1) * t53 - rSges(4,2) * t54)) - m(5) * (g(1) * (t116 * t34 + t33 * t92 - t25) + g(2) * (t116 * t30 + t29 * t92 - t23) + g(3) * (t116 * t54 + t53 * t92 - t48)) - m(6) * (g(1) * (t117 * t34 + t33 * t91 + t110) + g(2) * (t117 * t30 + t29 * t91 + t111) + g(3) * (t117 * t54 + t53 * t91 + t109)) + (-g(1) * t110 - g(2) * t111 - t30 * t98 - t34 * t99 - t139 * (-t74 * pkin(5) + (-t72 * t75 - t74 * t78) * rSges(7,1) + (-t72 * t78 + t74 * t75) * rSges(7,2)) + (-t112 * t54 - t109) * g(3)) * m(7) (-m(5) + t138) * t139, t138 * (g(1) * t11 + g(2) * t7 + g(3) * t27) -m(7) * (g(1) * (rSges(7,1) * t2 - rSges(7,2) * t3) + g(2) * (t143 * rSges(7,1) + t142 * rSges(7,2)) + g(3) * ((t27 * t78 - t28 * t75) * rSges(7,1) + (-t27 * t75 - t28 * t78) * rSges(7,2)))];
taug  = t1(:);
