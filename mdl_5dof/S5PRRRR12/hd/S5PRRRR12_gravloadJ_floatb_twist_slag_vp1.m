% Calculate Gravitation load on the joints for
% S5PRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d2,d3,d4,d5,theta1]';
% m [6x1]
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
% Datum: 2024-09-28 18:09
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRR12_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR12_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR12_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR12_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR12_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR12_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-28 18:07:20
% EndTime: 2024-09-28 18:07:22
% DurationCPUTime: 0.97s
% Computational Cost: add. (681->167), mult. (908->270), div. (0->0), fcn. (925->26), ass. (0->105)
t82 = sin(pkin(6));
t141 = pkin(10) * t82;
t81 = sin(pkin(11));
t112 = t81 * t141;
t84 = cos(pkin(11));
t86 = cos(pkin(5));
t48 = t84 * pkin(4) + t86 * t112;
t111 = t84 * t141;
t130 = t81 * t86;
t50 = pkin(4) * t130 - t111;
t88 = sin(qJ(4));
t92 = cos(qJ(4));
t13 = -t48 * t92 + t50 * t88;
t10 = -t84 * pkin(3) + t13;
t142 = t48 * t88 + t50 * t92;
t7 = -pkin(3) * t130 - t142;
t89 = sin(qJ(3));
t93 = cos(qJ(3));
t146 = t10 * t89 + t7 * t93;
t49 = -t81 * pkin(4) + t86 * t111;
t125 = t84 * t86;
t51 = pkin(4) * t125 + t112;
t14 = t49 * t92 - t51 * t88;
t11 = -pkin(3) * t81 + t14;
t104 = t49 * t88 + t51 * t92;
t8 = pkin(3) * t125 + t104;
t145 = t11 * t89 + t8 * t93;
t144 = t10 * t93 - t7 * t89;
t143 = t11 * t93 - t8 * t89;
t138 = rSges(6,3) * t82;
t58 = -t88 * pkin(4) + t92 * t141;
t133 = t58 * t89;
t80 = qJ(2) + qJ(3);
t79 = qJ(4) + t80;
t72 = cos(t79);
t132 = t81 * t72;
t76 = sin(t80);
t131 = t81 * t76;
t83 = sin(pkin(5));
t129 = t82 * t83;
t94 = cos(qJ(2));
t128 = t83 * t94;
t127 = t84 * t72;
t126 = t84 * t76;
t85 = cos(pkin(6));
t87 = sin(qJ(5));
t124 = t85 * t87;
t91 = cos(qJ(5));
t123 = t85 * t91;
t122 = t86 * t87;
t90 = sin(qJ(2));
t121 = t86 * t90;
t120 = t86 * t91;
t119 = t86 * t94;
t118 = t89 * t90;
t75 = pkin(5) - t80;
t70 = -qJ(4) + t75;
t60 = sin(t70) / 0.2e1;
t74 = pkin(5) + t80;
t69 = qJ(4) + t74;
t65 = sin(t69);
t44 = t60 - t65 / 0.2e1;
t61 = cos(t69) / 0.2e1;
t66 = cos(t70);
t45 = t66 / 0.2e1 + t61;
t71 = sin(t79);
t117 = (t84 * t45 - t81 * t71) * rSges(5,1) + (t84 * t44 - t132) * rSges(5,2);
t116 = (-t81 * t45 - t84 * t71) * rSges(5,1) + (-t81 * t44 - t127) * rSges(5,2);
t115 = (t65 / 0.2e1 + t60) * rSges(5,1) + (t61 - t66 / 0.2e1) * rSges(5,2);
t108 = t85 * t122;
t30 = t81 * t108 - t84 * t91;
t107 = t85 * t120;
t31 = -t81 * t107 - t84 * t87;
t36 = t81 * t120 + t84 * t124;
t37 = t81 * t122 - t84 * t123;
t114 = (t30 * t71 - t36 * t72) * rSges(6,1) + (-t31 * t71 + t37 * t72) * rSges(6,2) + (-t71 * t130 + t127) * t138;
t32 = t84 * t108 + t81 * t91;
t33 = t84 * t107 - t81 * t87;
t34 = -t84 * t120 + t81 * t124;
t35 = -t84 * t122 - t81 * t123;
t113 = (t71 * t125 + t132) * t138 + (-t32 * t71 - t34 * t72) * rSges(6,1) + (-t33 * t71 + t35 * t72) * rSges(6,2);
t110 = t87 * t129;
t109 = t91 * t129;
t106 = t71 * rSges(6,3) * t129 + ((-t71 * t124 + t72 * t91) * rSges(6,1) + (-t71 * t123 - t72 * t87) * rSges(6,2)) * t83;
t52 = t58 * t93;
t57 = -pkin(4) * t92 - t88 * t141;
t56 = pkin(3) - t57;
t103 = (-t89 * t56 + t52) * t83 * t90 + t106;
t100 = -t56 * t93 - t133;
t99 = -pkin(3) * t118 + (t93 * pkin(3) + pkin(2)) * t94;
t98 = t84 * t119 - t81 * t90;
t97 = -t81 * t119 - t84 * t90;
t96 = pkin(3) * (t93 * t94 - t118);
t62 = sin(t75) / 0.2e1;
t67 = sin(t74);
t53 = t62 - t67 / 0.2e1;
t63 = cos(t74) / 0.2e1;
t68 = cos(t75);
t54 = t68 / 0.2e1 + t63;
t77 = cos(t80);
t95 = g(1) * ((-t81 * t54 - t126) * rSges(4,1) + (-t81 * t53 - t84 * t77) * rSges(4,2)) + g(2) * ((t84 * t54 - t131) * rSges(4,1) + (t84 * t53 - t81 * t77) * rSges(4,2)) + g(3) * ((t67 / 0.2e1 + t62) * rSges(4,1) + (t63 - t68 / 0.2e1) * rSges(4,2));
t59 = -t90 * pkin(2) - pkin(3) * t76;
t29 = t86 * t96;
t28 = t99 * t86;
t1 = [(-m(2) - m(3) - m(4) - m(5) - m(6)) * g(3), -m(3) * (g(1) * (t97 * rSges(3,1) + (t81 * t121 - t84 * t94) * rSges(3,2)) + g(2) * (t98 * rSges(3,1) + (-t84 * t121 - t81 * t94) * rSges(3,2)) + g(3) * (rSges(3,1) * t94 - rSges(3,2) * t90) * t83) - m(4) * ((g(1) * t97 + g(2) * t98 + g(3) * t128) * pkin(2) + t95) - m(5) * (g(1) * (-t81 * t28 + t84 * t59 + t116) + g(2) * (t84 * t28 + t81 * t59 + t117) + g(3) * (t99 * t83 + t115)) - m(6) * (g(1) * (-(t84 * pkin(2) - t144) * t90 + (-pkin(2) * t130 + t146) * t94 + t114) + g(2) * (-(t81 * pkin(2) - t143) * t90 + (pkin(2) * t125 + t145) * t94 + t113) + g(3) * ((pkin(2) - t100) * t128 + t103)), -m(4) * t95 - m(5) * (g(1) * (-pkin(3) * t126 - t81 * t29 + t116) + g(2) * (-pkin(3) * t131 + t84 * t29 + t117) + g(3) * (t83 * t96 + t115)) - m(6) * (g(1) * (t144 * t90 + t146 * t94 + t114) + g(2) * (t143 * t90 + t145 * t94 + t113) + g(3) * (-t100 * t128 + t103)), -m(5) * (g(1) * t116 + g(2) * t117 + g(3) * t115) - m(6) * (g(1) * ((t13 * t89 - t142 * t93) * t94 + (t13 * t93 + t142 * t89) * t90 + t114) + g(2) * ((t104 * t93 + t14 * t89) * t94 + (-t104 * t89 + t14 * t93) * t90 + t113) + g(3) * (((t57 * t89 + t52) * t90 - (t57 * t93 - t133) * t94) * t83 + t106)), -m(6) * (g(1) * ((t81 * t109 + t31 * t72 + t37 * t71) * rSges(6,1) + (-t81 * t110 + t30 * t72 + t36 * t71) * rSges(6,2)) + g(2) * ((-t84 * t109 + t33 * t72 + t35 * t71) * rSges(6,1) + (t84 * t110 - t32 * t72 + t34 * t71) * rSges(6,2)) + g(3) * ((rSges(6,1) * t91 - rSges(6,2) * t87) * t86 * t82 + ((t72 * t123 - t71 * t87) * rSges(6,1) + (-t72 * t124 - t71 * t91) * rSges(6,2)) * t83))];
taug = t1(:);
