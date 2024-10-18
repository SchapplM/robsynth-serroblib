% Calculate Gravitation load on the joints for
% S5RRRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,d3,d4,d5]';
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
% Datum: 2024-09-27 18:44
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRR14_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR14_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR14_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR14_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR14_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR14_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 18:42:23
% EndTime: 2024-09-27 18:42:24
% DurationCPUTime: 0.41s
% Computational Cost: add. (724->125), mult. (499->178), div. (0->0), fcn. (450->22), ass. (0->80)
t121 = rSges(4,3) + pkin(8);
t120 = pkin(8) + pkin(9);
t76 = qJ(3) + qJ(4);
t68 = pkin(5) + t76;
t63 = qJ(5) + t68;
t55 = cos(t63) / 0.2e1;
t104 = pkin(5) - t76;
t103 = -qJ(5) + t104;
t59 = cos(t103);
t36 = t59 / 0.2e1 + t55;
t73 = qJ(5) + t76;
t64 = sin(t73);
t77 = qJ(1) + qJ(2);
t70 = sin(t77);
t72 = cos(t77);
t9 = -t72 * t36 + t70 * t64;
t107 = sin(t63) / 0.2e1;
t92 = sin(t103);
t35 = t107 - t92 / 0.2e1;
t65 = cos(t73);
t99 = -t72 * t35 - t70 * t65;
t119 = -t9 * rSges(6,1) + t99 * rSges(6,2);
t10 = -t70 * t36 - t72 * t64;
t98 = t70 * t35 - t72 * t65;
t118 = t10 * rSges(6,1) + t98 * rSges(6,2);
t78 = sin(pkin(5));
t117 = g(3) * t78;
t82 = sin(qJ(1));
t116 = t82 * pkin(1);
t84 = cos(qJ(3));
t67 = t84 * pkin(3) + pkin(2);
t69 = sin(t76);
t115 = t70 * t69;
t114 = t70 * t78;
t113 = t72 * t69;
t112 = t72 * t78;
t79 = cos(pkin(5));
t81 = sin(qJ(3));
t111 = t79 * t81;
t110 = t79 * t84;
t80 = sin(qJ(4));
t109 = t80 * t81;
t108 = (t107 + t92 / 0.2e1) * rSges(6,1) + (t55 - t59 / 0.2e1) * rSges(6,2);
t106 = sin(t68) / 0.2e1;
t105 = t72 * rSges(3,1) - t70 * rSges(3,2);
t30 = -t70 * t110 - t72 * t81;
t31 = -t70 * t111 + t72 * t84;
t102 = t31 * rSges(4,1) + t30 * rSges(4,2) + t72 * pkin(2) + t121 * t114;
t101 = sin(t104);
t100 = -t70 * rSges(3,1) - t72 * rSges(3,2);
t39 = t106 - t101 / 0.2e1;
t71 = cos(t76);
t97 = -t72 * t39 - t70 * t71;
t96 = t70 * t39 - t72 * t71;
t83 = cos(qJ(4));
t66 = t83 * pkin(4) + pkin(3);
t95 = -pkin(4) * t109 + t66 * t84;
t94 = t72 * t110 - t70 * t81;
t21 = -t78 * (pkin(10) + t120) + (t84 * t80 * pkin(4) + t81 * t66) * t79;
t42 = pkin(4) * t71 + t67;
t93 = -t98 * rSges(6,1) + t10 * rSges(6,2) + rSges(6,3) * t114 - t70 * t21 + t72 * t42;
t29 = -t72 * t111 - t70 * t84;
t91 = t29 * rSges(4,1) - t94 * rSges(4,2) - t70 * pkin(2) + t121 * t112;
t58 = cos(t68) / 0.2e1;
t61 = cos(t104);
t40 = t61 / 0.2e1 + t58;
t20 = -t70 * t40 - t113;
t41 = pkin(3) * t111 - t78 * t120;
t90 = -t96 * rSges(5,1) + t20 * rSges(5,2) + rSges(5,3) * t114 - t70 * t41 + t72 * t67;
t89 = pkin(4) * (t83 * t84 - t109);
t88 = t99 * rSges(6,1) + t9 * rSges(6,2) + rSges(6,3) * t112 - t72 * t21 - t70 * t42;
t19 = -t72 * t40 + t115;
t87 = t97 * rSges(5,1) + t19 * rSges(5,2) + rSges(5,3) * t112 - t72 * t41 - t70 * t67;
t86 = g(1) * (t20 * rSges(5,1) + t96 * rSges(5,2)) + g(2) * (-t19 * rSges(5,1) + t97 * rSges(5,2)) + g(3) * ((t106 + t101 / 0.2e1) * rSges(5,1) + (t58 - t61 / 0.2e1) * rSges(5,2));
t85 = cos(qJ(1));
t75 = t85 * pkin(1);
t43 = -t81 * pkin(3) - pkin(4) * t69;
t27 = t79 * t89;
t22 = t95 * t79;
t1 = [-m(2) * (g(1) * (-t82 * rSges(2,1) - t85 * rSges(2,2)) + g(2) * (t85 * rSges(2,1) - t82 * rSges(2,2))) - m(3) * (g(1) * (t100 - t116) + g(2) * (t105 + t75)) - m(4) * (g(1) * (t91 - t116) + g(2) * (t75 + t102)) - m(5) * (g(1) * (t87 - t116) + g(2) * (t75 + t90)) - m(6) * (g(1) * (t88 - t116) + g(2) * (t75 + t93)), -m(3) * (g(1) * t100 + g(2) * t105) - m(4) * (g(1) * t91 + g(2) * t102) - m(5) * (g(1) * t87 + g(2) * t90) - m(6) * (g(1) * t88 + g(2) * t93), -m(4) * (g(1) * (t30 * rSges(4,1) - t31 * rSges(4,2)) + g(2) * (rSges(4,1) * t94 + t29 * rSges(4,2)) + (rSges(4,1) * t84 - rSges(4,2) * t81) * t117) - m(5) * ((g(1) * t30 + g(2) * t94 + t84 * t117) * pkin(3) + t86) - m(6) * (g(1) * (-t70 * t22 + t72 * t43 + t118) + g(2) * (t72 * t22 + t70 * t43 + t119) + g(3) * (t95 * t78 + t108)), -m(5) * t86 - m(6) * (g(1) * (-pkin(4) * t113 - t70 * t27 + t118) + g(2) * (-pkin(4) * t115 + t72 * t27 + t119) + g(3) * (t78 * t89 + t108)), -m(6) * (g(1) * t118 + g(2) * t119 + g(3) * t108)];
taug = t1(:);
