% Calculate Gravitation load on the joints for
% S6RPRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 07:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRR11_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR11_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR11_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR11_gravloadJ_floatb_twist_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR11_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR11_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:38:28
% EndTime: 2019-03-09 07:38:33
% DurationCPUTime: 2.02s
% Computational Cost: add. (1316->198), mult. (3392->295), div. (0->0), fcn. (4343->16), ass. (0->86)
t120 = cos(qJ(3));
t107 = cos(pkin(7));
t121 = cos(qJ(1));
t103 = sin(pkin(13));
t119 = sin(qJ(1));
t106 = cos(pkin(13));
t108 = cos(pkin(6));
t91 = t108 * t106;
t78 = t119 * t103 - t121 * t91;
t104 = sin(pkin(7));
t105 = sin(pkin(6));
t88 = t105 * t104;
t141 = t78 * t107 + t121 * t88;
t90 = t108 * t103;
t42 = t106 * t119 + t121 * t90;
t60 = sin(qJ(3));
t26 = -t120 * t42 + t141 * t60;
t59 = sin(qJ(4));
t62 = cos(qJ(4));
t89 = t107 * t105;
t68 = -t78 * t104 + t121 * t89;
t14 = t26 * t62 + t59 * t68;
t23 = t141 * t120 + t42 * t60;
t58 = sin(qJ(5));
t61 = cos(qJ(5));
t93 = t14 * t58 + t23 * t61;
t144 = t14 * t61 - t23 * t58;
t140 = t26 * t59 - t62 * t68;
t72 = t103 * t121 + t119 * t91;
t64 = t72 * t104 + t119 * t89;
t137 = t72 * t107 - t119 * t88;
t136 = t104 * t108 + t106 * t89;
t43 = t106 * t121 - t119 * t90;
t28 = t43 * t120 - t137 * t60;
t16 = t28 * t62 + t59 * t64;
t87 = t105 * t103;
t33 = t120 * t87 + t136 * t60;
t41 = -t106 * t88 + t107 * t108;
t22 = t33 * t62 + t41 * t59;
t134 = g(1) * t16 - g(2) * t14 + g(3) * t22;
t15 = t28 * t59 - t62 * t64;
t21 = -t33 * t59 + t41 * t62;
t133 = -g(1) * t15 + g(2) * t140 + g(3) * t21;
t27 = t120 * t137 + t43 * t60;
t32 = -t120 * t136 + t60 * t87;
t132 = (-g(1) * t27 - g(2) * t23 - g(3) * t32) * t59;
t124 = t62 * pkin(4);
t123 = pkin(10) + rSges(5,3);
t122 = pkin(11) + rSges(6,3);
t118 = t26 * t58;
t117 = t28 * t58;
t116 = t33 * t58;
t57 = qJ(5) + qJ(6);
t54 = sin(t57);
t115 = t54 * t62;
t55 = cos(t57);
t114 = t55 * t62;
t113 = t58 * t62;
t112 = t61 * t62;
t53 = pkin(5) * t61 + pkin(4);
t111 = t62 * t53;
t95 = t105 * t119;
t110 = t121 * pkin(1) + qJ(2) * t95;
t109 = pkin(12) + pkin(11) + rSges(7,3);
t102 = pkin(5) * t58 + pkin(10);
t17 = t23 * pkin(3);
t101 = -pkin(10) * t26 - t17;
t19 = t27 * pkin(3);
t100 = t28 * pkin(10) - t19;
t31 = t32 * pkin(3);
t99 = t33 * pkin(10) - t31;
t96 = t121 * t105;
t98 = -pkin(1) * t119 + qJ(2) * t96;
t94 = -rSges(5,1) * t62 + rSges(5,2) * t59;
t7 = -t16 * t58 + t27 * t61;
t92 = -t22 * t58 + t32 * t61;
t84 = t55 * rSges(7,1) - t54 * rSges(7,2) + t53;
t5 = -t16 * t54 + t27 * t55;
t6 = t16 * t55 + t27 * t54;
t77 = m(7) * (g(1) * (t5 * rSges(7,1) - t6 * rSges(7,2)) + g(2) * ((t14 * t54 + t23 * t55) * rSges(7,1) + (t14 * t55 - t23 * t54) * rSges(7,2)) + g(3) * ((-t22 * t54 + t32 * t55) * rSges(7,1) + (-t22 * t55 - t32 * t54) * rSges(7,2)));
t69 = -t42 * pkin(2) + t68 * pkin(9) + t98;
t67 = t26 * pkin(3) + t69;
t66 = t43 * pkin(2) + t64 * pkin(9) + t110;
t65 = t28 * pkin(3) + t66;
t8 = t16 * t61 + t27 * t58;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t119 - rSges(2,2) * t121) + g(2) * (rSges(2,1) * t121 - rSges(2,2) * t119)) - m(3) * (g(1) * (-t42 * rSges(3,1) + rSges(3,2) * t78 + rSges(3,3) * t96 + t98) + g(2) * (t43 * rSges(3,1) - rSges(3,2) * t72 + rSges(3,3) * t95 + t110)) - m(4) * (g(1) * (t26 * rSges(4,1) + rSges(4,2) * t23 + rSges(4,3) * t68 + t69) + g(2) * (t28 * rSges(4,1) - t27 * rSges(4,2) + rSges(4,3) * t64 + t66)) - m(5) * (g(1) * (t14 * rSges(5,1) - rSges(5,2) * t140 - t123 * t23 + t67) + g(2) * (t16 * rSges(5,1) - t15 * rSges(5,2) + t123 * t27 + t65)) - m(6) * (g(1) * (t144 * rSges(6,1) - t93 * rSges(6,2) + t14 * pkin(4) - t23 * pkin(10) + t122 * t140 + t67) + g(2) * (t8 * rSges(6,1) + t7 * rSges(6,2) + t16 * pkin(4) + t27 * pkin(10) + t122 * t15 + t65)) - m(7) * (g(1) * (t84 * t14 - (rSges(7,1) * t54 + rSges(7,2) * t55 + t102) * t23 + t109 * t140 + t67) + g(2) * (t6 * rSges(7,1) + t5 * rSges(7,2) + t102 * t27 + t109 * t15 + t16 * t53 + t65)) (-m(3) - m(4) - m(5) - m(6) - m(7)) * (g(1) * t95 - g(2) * t96 + g(3) * t108) -m(4) * (g(1) * (-rSges(4,1) * t27 - rSges(4,2) * t28) + g(2) * (-rSges(4,1) * t23 + rSges(4,2) * t26) + g(3) * (-rSges(4,1) * t32 - rSges(4,2) * t33)) - m(5) * (g(1) * (t123 * t28 + t27 * t94 - t19) + g(2) * (-t123 * t26 + t23 * t94 - t17) + g(3) * (t123 * t33 + t32 * t94 - t31)) + (-g(1) * (-t27 * t111 + pkin(5) * t117 + (-t114 * t27 + t28 * t54) * rSges(7,1) + (t115 * t27 + t28 * t55) * rSges(7,2) + t100) - g(2) * (-t23 * t111 - pkin(5) * t118 + (-t114 * t23 - t26 * t54) * rSges(7,1) + (t115 * t23 - t26 * t55) * rSges(7,2) + t101) - g(3) * (-t32 * t111 + pkin(5) * t116 + (-t114 * t32 + t33 * t54) * rSges(7,1) + (t115 * t32 + t33 * t55) * rSges(7,2) + t99) - t109 * t132) * m(7) + (-g(1) * (-t27 * t124 + (-t112 * t27 + t117) * rSges(6,1) + (t113 * t27 + t28 * t61) * rSges(6,2) + t100) - g(2) * (-t23 * t124 + (-t112 * t23 - t118) * rSges(6,1) + (t113 * t23 - t26 * t61) * rSges(6,2) + t101) - g(3) * (-t32 * t124 + (-t112 * t32 + t116) * rSges(6,1) + (t113 * t32 + t33 * t61) * rSges(6,2) + t99) - t122 * t132) * m(6), -m(5) * (g(1) * (-rSges(5,1) * t15 - rSges(5,2) * t16) + g(2) * (rSges(5,1) * t140 + rSges(5,2) * t14) + g(3) * (rSges(5,1) * t21 - rSges(5,2) * t22)) - m(6) * (t133 * (t61 * rSges(6,1) - t58 * rSges(6,2) + pkin(4)) + t134 * t122) - m(7) * (t109 * t134 + t133 * t84) -m(6) * (g(1) * (rSges(6,1) * t7 - rSges(6,2) * t8) + g(2) * (t93 * rSges(6,1) + t144 * rSges(6,2)) + g(3) * (t92 * rSges(6,1) + (-t22 * t61 - t32 * t58) * rSges(6,2))) - t77 - m(7) * (g(1) * t7 + g(2) * t93 + g(3) * t92) * pkin(5), -t77];
taug  = t1(:);
