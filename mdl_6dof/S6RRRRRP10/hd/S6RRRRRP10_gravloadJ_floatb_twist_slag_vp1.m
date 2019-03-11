% Calculate Gravitation load on the joints for
% S6RRRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 02:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRP10_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP10_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP10_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP10_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP10_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP10_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:16:20
% EndTime: 2019-03-10 02:16:25
% DurationCPUTime: 1.68s
% Computational Cost: add. (980->215), mult. (1969->305), div. (0->0), fcn. (2391->12), ass. (0->99)
t154 = cos(qJ(1));
t89 = sin(pkin(6));
t124 = t89 * t154;
t129 = cos(pkin(6));
t117 = t129 * t154;
t153 = sin(qJ(1));
t92 = sin(qJ(2));
t95 = cos(qJ(2));
t68 = t92 * t117 + t153 * t95;
t91 = sin(qJ(3));
t94 = cos(qJ(3));
t39 = -t91 * t124 + t68 * t94;
t67 = -t95 * t117 + t153 * t92;
t90 = sin(qJ(4));
t93 = cos(qJ(4));
t171 = -t39 * t90 + t67 * t93;
t157 = rSges(7,1) + pkin(5);
t163 = pkin(4) * t90;
t172 = pkin(9) + t163;
t130 = rSges(7,3) + qJ(6);
t123 = t89 * t153;
t116 = t129 * t153;
t70 = -t92 * t116 + t154 * t95;
t43 = t91 * t123 + t70 * t94;
t69 = t95 * t116 + t154 * t92;
t17 = -t43 * t90 + t69 * t93;
t148 = t67 * t90;
t170 = -t39 * t93 - t148;
t88 = qJ(4) + qJ(5);
t85 = sin(t88);
t86 = cos(t88);
t11 = t39 * t85 - t67 * t86;
t12 = t39 * t86 + t67 * t85;
t122 = -t94 * t124 - t68 * t91;
t42 = -t94 * t123 + t70 * t91;
t141 = t89 * t92;
t65 = t129 * t94 - t91 * t141;
t169 = -g(1) * t42 + g(2) * t122 + g(3) * t65;
t168 = g(1) * t69 + g(2) * t67;
t140 = t89 * t95;
t167 = (g(3) * t140 - t168) * t91;
t105 = t93 * rSges(5,1) - t90 * rSges(5,2) + pkin(3);
t155 = pkin(10) + rSges(5,3);
t166 = t105 * t94 + t155 * t91;
t165 = -t11 * rSges(6,1) - t12 * rSges(6,2);
t15 = t43 * t85 - t69 * t86;
t16 = t43 * t86 + t69 * t85;
t164 = -t15 * rSges(6,1) - t16 * rSges(6,2);
t158 = g(3) * t89;
t156 = rSges(4,3) + pkin(9);
t146 = t69 * t90;
t84 = pkin(4) * t93 + pkin(3);
t144 = t84 * t94;
t143 = t85 * t94;
t142 = t86 * t94;
t139 = t94 * t95;
t96 = -pkin(11) - pkin(10);
t136 = t122 * t84 - t39 * t96;
t135 = -t42 * t84 - t43 * t96;
t66 = t129 * t91 + t94 * t141;
t32 = t86 * t140 + t66 * t85;
t127 = t85 * t140;
t33 = t66 * t86 - t127;
t134 = -t32 * rSges(6,1) - t33 * rSges(6,2);
t133 = t65 * t84 - t66 * t96;
t132 = t154 * pkin(1) + pkin(8) * t123;
t131 = pkin(2) * t140 + pkin(9) * t141;
t126 = t70 * pkin(2) + t132;
t125 = g(3) * t131;
t121 = t89 * t84 * t139 + t141 * t163 + t131;
t120 = t171 * pkin(4);
t119 = t17 * pkin(4);
t118 = -t153 * pkin(1) + pkin(8) * t124;
t115 = rSges(4,1) * t94 - rSges(4,2) * t91;
t114 = t90 * rSges(5,1) + t93 * rSges(5,2);
t113 = rSges(6,1) * t86 - rSges(6,2) * t85;
t112 = -t157 * t11 + t130 * t12;
t61 = t67 * pkin(2);
t111 = -t67 * t144 + t172 * t68 - t61;
t63 = t69 * pkin(2);
t110 = -t69 * t144 + t172 * t70 - t63;
t109 = pkin(9) * t69 + t126;
t108 = t130 * t16 - t157 * t15;
t107 = t130 * t33 - t157 * t32;
t106 = -t68 * pkin(2) + t118;
t104 = pkin(9) + t114;
t103 = -t93 * t140 - t66 * t90;
t102 = t103 * pkin(4);
t101 = -t67 * pkin(9) + t106;
t100 = pkin(4) * t146 - t42 * t96 + t43 * t84 + t109;
t98 = -pkin(4) * t148 - t122 * t96 - t39 * t84 + t101;
t45 = (t86 * t139 + t85 * t92) * t89;
t44 = t94 * t127 - t86 * t141;
t22 = -t69 * t142 + t70 * t85;
t21 = -t69 * t143 - t70 * t86;
t20 = -t67 * t142 + t68 * t85;
t19 = -t67 * t143 - t68 * t86;
t18 = t43 * t93 + t146;
t1 = [-m(2) * (g(1) * (-t153 * rSges(2,1) - t154 * rSges(2,2)) + g(2) * (t154 * rSges(2,1) - t153 * rSges(2,2))) - m(3) * (g(1) * (-t68 * rSges(3,1) + t67 * rSges(3,2) + rSges(3,3) * t124 + t118) + g(2) * (t70 * rSges(3,1) - t69 * rSges(3,2) + rSges(3,3) * t123 + t132)) - m(4) * (g(1) * (-rSges(4,1) * t39 - rSges(4,2) * t122 - t156 * t67 + t106) + g(2) * (rSges(4,1) * t43 - rSges(4,2) * t42 + t156 * t69 + t126)) - m(5) * (g(1) * (t170 * rSges(5,1) - rSges(5,2) * t171 - t39 * pkin(3) + t155 * t122 + t101) + g(2) * (rSges(5,1) * t18 + rSges(5,2) * t17 + pkin(3) * t43 + t155 * t42 + t109)) - m(6) * (g(1) * (-rSges(6,1) * t12 + rSges(6,2) * t11 + rSges(6,3) * t122 + t98) + g(2) * (rSges(6,1) * t16 - rSges(6,2) * t15 + rSges(6,3) * t42 + t100)) - m(7) * (g(1) * (rSges(7,2) * t122 - t11 * t130 - t12 * t157 + t98) + g(2) * (rSges(7,2) * t42 + t130 * t15 + t157 * t16 + t100)) -m(3) * (g(1) * (-rSges(3,1) * t69 - rSges(3,2) * t70) + g(2) * (-rSges(3,1) * t67 - rSges(3,2) * t68) + (rSges(3,1) * t95 - rSges(3,2) * t92) * t158) - m(4) * (g(1) * (-t115 * t69 + t156 * t70 - t63) + g(2) * (-t115 * t67 + t156 * t68 - t61) + t125 + (rSges(4,3) * t92 + t115 * t95) * t158) - m(5) * (t125 + (t114 * t92 + t166 * t95) * t158 - t168 * t166 + (t104 * t68 - t61) * g(2) + (t104 * t70 - t63) * g(1)) - m(6) * (g(1) * (t22 * rSges(6,1) - t21 * rSges(6,2) + t110) + g(2) * (t20 * rSges(6,1) - t19 * rSges(6,2) + t111) + g(3) * (t45 * rSges(6,1) - t44 * rSges(6,2) + t121) + (rSges(6,3) - t96) * t167) - m(7) * (g(1) * (t130 * t21 + t157 * t22 + t110) + g(2) * (t130 * t19 + t157 * t20 + t111) + g(3) * (t130 * t44 + t157 * t45 + t121) + (rSges(7,2) - t96) * t167) -m(4) * (g(1) * (-rSges(4,1) * t42 - rSges(4,2) * t43) + g(2) * (rSges(4,1) * t122 - rSges(4,2) * t39) + g(3) * (rSges(4,1) * t65 - rSges(4,2) * t66)) - m(5) * ((g(1) * t43 + g(2) * t39 + g(3) * t66) * t155 + t169 * t105) - m(6) * (g(1) * (rSges(6,3) * t43 - t113 * t42 + t135) + g(2) * (rSges(6,3) * t39 + t113 * t122 + t136) + g(3) * (rSges(6,3) * t66 + t113 * t65 + t133)) - m(7) * (g(1) * (rSges(7,2) * t43 + t135) + g(2) * (rSges(7,2) * t39 + t136) + g(3) * (rSges(7,2) * t66 + t133) + t169 * (t130 * t85 + t157 * t86)) -m(5) * (g(1) * (rSges(5,1) * t17 - rSges(5,2) * t18) + g(2) * (t171 * rSges(5,1) + t170 * rSges(5,2)) + g(3) * (t103 * rSges(5,1) + (t90 * t140 - t66 * t93) * rSges(5,2))) - m(6) * (g(1) * (t119 + t164) + g(2) * (t120 + t165) + g(3) * (t102 + t134)) - m(7) * (g(1) * (t108 + t119) + g(2) * (t112 + t120) + g(3) * (t102 + t107)) -m(6) * (g(1) * t164 + g(2) * t165 + g(3) * t134) - m(7) * (g(1) * t108 + g(2) * t112 + g(3) * t107) -m(7) * (g(1) * t15 + g(2) * t11 + g(3) * t32)];
taug  = t1(:);
