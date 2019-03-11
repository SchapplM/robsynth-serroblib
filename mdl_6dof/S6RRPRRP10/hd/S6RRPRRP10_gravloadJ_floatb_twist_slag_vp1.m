% Calculate Gravitation load on the joints for
% S6RRPRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 12:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRP10_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP10_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP10_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP10_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP10_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP10_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:36:34
% EndTime: 2019-03-09 12:36:38
% DurationCPUTime: 1.26s
% Computational Cost: add. (842->191), mult. (1447->274), div. (0->0), fcn. (1730->12), ass. (0->81)
t74 = sin(qJ(2));
t75 = sin(qJ(1));
t77 = cos(qJ(2));
t115 = cos(qJ(1));
t97 = cos(pkin(6));
t85 = t97 * t115;
t47 = t74 * t85 + t75 * t77;
t68 = pkin(11) + qJ(4);
t65 = sin(t68);
t66 = cos(t68);
t70 = sin(pkin(6));
t93 = t70 * t115;
t19 = t47 * t66 - t65 * t93;
t46 = t74 * t75 - t77 * t85;
t73 = sin(qJ(5));
t76 = cos(qJ(5));
t1 = t19 * t73 - t46 * t76;
t2 = t19 * t76 + t46 * t73;
t117 = rSges(7,2) + pkin(10);
t118 = rSges(7,1) + pkin(5);
t99 = rSges(7,3) + qJ(6);
t78 = t118 * t76 + t99 * t73;
t121 = pkin(4) * t66;
t120 = g(2) * t46;
t119 = g(3) * t70;
t116 = rSges(6,3) + pkin(10);
t114 = t46 * t65;
t90 = t75 * t97;
t48 = t115 * t74 + t77 * t90;
t111 = t48 * t65;
t110 = t65 * t77;
t109 = t66 * t73;
t108 = t66 * t76;
t107 = t70 * t74;
t106 = t70 * t75;
t105 = t70 * t77;
t72 = -pkin(9) - qJ(3);
t104 = t72 * t74;
t103 = t76 * t77;
t71 = cos(pkin(11));
t64 = pkin(3) * t71 + pkin(2);
t102 = -t46 * t64 - t47 * t72;
t49 = t115 * t77 - t74 * t90;
t101 = -t48 * t64 - t49 * t72;
t100 = t115 * pkin(1) + pkin(8) * t106;
t98 = qJ(3) + rSges(4,3);
t69 = sin(pkin(11));
t96 = t69 * t106;
t95 = t73 * t105;
t50 = t64 * t105;
t94 = t50 + (pkin(10) * t65 + t121) * t105;
t92 = -t75 * pkin(1) + pkin(8) * t93;
t91 = -t47 * t65 - t66 * t93;
t89 = -pkin(10) * t114 - t46 * t121 + t102;
t88 = -pkin(10) * t111 - t48 * t121 + t101;
t87 = t69 * t93;
t86 = pkin(3) * t96 - t48 * t72 + t49 * t64 + t100;
t84 = rSges(5,1) * t66 - rSges(5,2) * t65;
t83 = rSges(6,1) * t76 - rSges(6,2) * t73;
t23 = t65 * t106 + t49 * t66;
t82 = t23 * pkin(4) + t86;
t81 = t71 * rSges(4,1) - t69 * rSges(4,2) + pkin(2);
t80 = pkin(3) * t87 + t46 * t72 - t47 * t64 + t92;
t79 = -pkin(4) * t19 + t80;
t36 = t66 * t107 + t97 * t65;
t35 = -t65 * t107 + t97 * t66;
t34 = t35 * pkin(4);
t25 = (t66 * t103 + t73 * t74) * t70;
t24 = -t76 * t107 + t66 * t95;
t22 = -t66 * t106 + t49 * t65;
t17 = t36 * t76 - t95;
t16 = t70 * t103 + t36 * t73;
t14 = t22 * pkin(4);
t12 = t91 * pkin(4);
t10 = -t48 * t108 + t49 * t73;
t9 = -t48 * t109 - t49 * t76;
t8 = -t46 * t108 + t47 * t73;
t7 = -t46 * t109 - t47 * t76;
t6 = t23 * t76 + t48 * t73;
t5 = t23 * t73 - t48 * t76;
t3 = [-m(2) * (g(1) * (-t75 * rSges(2,1) - t115 * rSges(2,2)) + g(2) * (t115 * rSges(2,1) - t75 * rSges(2,2))) - m(3) * (g(1) * (-t47 * rSges(3,1) + t46 * rSges(3,2) + rSges(3,3) * t93 + t92) + g(2) * (rSges(3,1) * t49 - rSges(3,2) * t48 + rSges(3,3) * t106 + t100)) - m(4) * (g(1) * (-t47 * pkin(2) + (-t47 * t71 + t87) * rSges(4,1) + (t47 * t69 + t71 * t93) * rSges(4,2) - t98 * t46 + t92) + g(2) * (t49 * pkin(2) + (t49 * t71 + t96) * rSges(4,1) + (t71 * t106 - t49 * t69) * rSges(4,2) + t98 * t48 + t100)) - m(5) * (g(1) * (-rSges(5,1) * t19 - rSges(5,2) * t91 - rSges(5,3) * t46 + t80) + g(2) * (rSges(5,1) * t23 - t22 * rSges(5,2) + t48 * rSges(5,3) + t86)) - m(6) * (g(1) * (-rSges(6,1) * t2 + rSges(6,2) * t1 + t116 * t91 + t79) + g(2) * (rSges(6,1) * t6 - rSges(6,2) * t5 + t116 * t22 + t82)) - m(7) * (g(1) * (-t1 * t99 + t117 * t91 - t118 * t2 + t79) + g(2) * (t117 * t22 + t118 * t6 + t99 * t5 + t82)) -m(3) * (g(1) * (-rSges(3,1) * t48 - rSges(3,2) * t49) + g(2) * (-rSges(3,1) * t46 - rSges(3,2) * t47) + (rSges(3,1) * t77 - rSges(3,2) * t74) * t119) - m(4) * (g(1) * (-t81 * t48 + t98 * t49) + g(2) * t98 * t47 - t81 * t120 + (t98 * t74 + t81 * t77) * t119) - m(5) * (g(1) * (rSges(5,3) * t49 - t84 * t48 + t101) + g(2) * (rSges(5,3) * t47 - t84 * t46 + t102) + g(3) * t50 + (t84 * t77 + (rSges(5,3) - t72) * t74) * t119) - m(6) * (g(1) * (rSges(6,1) * t10 - rSges(6,2) * t9 - rSges(6,3) * t111 + t88) + g(2) * (rSges(6,1) * t8 - rSges(6,2) * t7 - rSges(6,3) * t114 + t89) + g(3) * (t25 * rSges(6,1) - t24 * rSges(6,2) + (rSges(6,3) * t110 - t104) * t70 + t94)) - m(7) * (g(1) * (-rSges(7,2) * t111 + t118 * t10 + t99 * t9 + t88) + g(2) * (-rSges(7,2) * t114 + t118 * t8 + t99 * t7 + t89) + g(3) * ((rSges(7,2) * t110 - t104) * t70 + t118 * t25 + t99 * t24 + t94)) (-m(4) - m(5) - m(6) - m(7)) * (g(1) * t48 - g(3) * t105 + t120) -m(5) * (g(1) * (-rSges(5,1) * t22 - rSges(5,2) * t23) + g(2) * (rSges(5,1) * t91 - rSges(5,2) * t19) + g(3) * (rSges(5,1) * t35 - rSges(5,2) * t36)) - m(6) * (g(1) * (t116 * t23 - t83 * t22 - t14) + g(2) * (t116 * t19 + t83 * t91 + t12) + g(3) * (t116 * t36 + t83 * t35 + t34)) - m(7) * ((t117 * t36 + t78 * t35 + t34) * g(3) + (t117 * t19 + t78 * t91 + t12) * g(2) + (t117 * t23 - t78 * t22 - t14) * g(1)) -m(6) * (g(1) * (-rSges(6,1) * t5 - rSges(6,2) * t6) + g(2) * (-rSges(6,1) * t1 - rSges(6,2) * t2) + g(3) * (-rSges(6,1) * t16 - rSges(6,2) * t17)) - m(7) * (g(1) * (-t118 * t5 + t99 * t6) + g(2) * (-t118 * t1 + t99 * t2) + g(3) * (-t118 * t16 + t99 * t17)) -m(7) * (g(1) * t5 + g(2) * t1 + g(3) * t16)];
taug  = t3(:);
