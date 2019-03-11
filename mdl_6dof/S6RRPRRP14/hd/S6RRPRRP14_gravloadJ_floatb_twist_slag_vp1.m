% Calculate Gravitation load on the joints for
% S6RRPRRP14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
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
% Datum: 2019-03-09 13:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRP14_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP14_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP14_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP14_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP14_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP14_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:03:47
% EndTime: 2019-03-09 13:03:50
% DurationCPUTime: 1.05s
% Computational Cost: add. (604->181), mult. (1448->257), div. (0->0), fcn. (1735->10), ass. (0->76)
t63 = sin(qJ(4));
t117 = pkin(4) * t63 + qJ(3);
t110 = rSges(7,1) + pkin(5);
t109 = rSges(7,2) + pkin(10);
t94 = rSges(7,3) + qJ(6);
t64 = sin(qJ(2));
t65 = sin(qJ(1));
t68 = cos(qJ(2));
t69 = cos(qJ(1));
t93 = cos(pkin(6));
t84 = t69 * t93;
t44 = t64 * t84 + t65 * t68;
t85 = t65 * t93;
t46 = -t64 * t85 + t68 * t69;
t116 = g(1) * t46 + g(2) * t44;
t62 = sin(qJ(5));
t66 = cos(qJ(5));
t61 = sin(pkin(6));
t101 = t61 * t69;
t43 = t64 * t65 - t68 * t84;
t67 = cos(qJ(4));
t79 = t101 * t67 - t43 * t63;
t3 = t44 * t66 + t62 * t79;
t4 = -t44 * t62 + t66 * t79;
t104 = t61 * t64;
t115 = (-g(3) * t104 - t116) * t67;
t72 = t110 * t66 + t62 * t94;
t111 = g(3) * t61;
t108 = rSges(5,3) + pkin(9);
t107 = rSges(6,3) + pkin(10);
t103 = t61 * t65;
t102 = t61 * t68;
t100 = t62 * t63;
t99 = t63 * t66;
t98 = t64 * t66;
t97 = pkin(2) * t102 + qJ(3) * t104;
t96 = t69 * pkin(1) + pkin(8) * t103;
t95 = rSges(4,3) + qJ(3);
t91 = t63 * t104;
t90 = t46 * pkin(2) + t96;
t89 = pkin(9) * t102 + t97;
t88 = -t65 * pkin(1) + pkin(8) * t101;
t37 = t43 * pkin(2);
t87 = -pkin(9) * t43 - t37;
t45 = t69 * t64 + t68 * t85;
t39 = t45 * pkin(2);
t86 = -pkin(9) * t45 - t39;
t83 = pkin(4) * t91 + t89;
t82 = -t44 * pkin(2) + t88;
t81 = rSges(5,1) * t63 + rSges(5,2) * t67;
t80 = rSges(6,1) * t66 - rSges(6,2) * t62;
t78 = t101 * t63 + t43 * t67;
t76 = pkin(3) * t103 + qJ(3) * t45 + t90;
t75 = t117 * t44 + t87;
t74 = t117 * t46 + t86;
t73 = pkin(3) * t101 - t43 * qJ(3) + t82;
t19 = t103 * t67 + t45 * t63;
t71 = t19 * pkin(4) + pkin(9) * t46 + t76;
t70 = pkin(4) * t79 - pkin(9) * t44 + t73;
t42 = -t102 * t63 + t67 * t93;
t41 = -t102 * t67 - t63 * t93;
t36 = t41 * pkin(4);
t25 = (t62 * t68 + t63 * t98) * t61;
t24 = -t102 * t66 + t62 * t91;
t18 = t103 * t63 - t45 * t67;
t17 = t104 * t62 + t42 * t66;
t16 = t42 * t62 - t61 * t98;
t14 = t78 * pkin(4);
t12 = t18 * pkin(4);
t11 = -t45 * t62 + t46 * t99;
t10 = t100 * t46 + t45 * t66;
t9 = -t43 * t62 + t44 * t99;
t8 = t100 * t44 + t43 * t66;
t2 = t19 * t66 + t46 * t62;
t1 = t19 * t62 - t46 * t66;
t5 = [-m(2) * (g(1) * (-t65 * rSges(2,1) - rSges(2,2) * t69) + g(2) * (rSges(2,1) * t69 - t65 * rSges(2,2))) - m(3) * (g(1) * (-t44 * rSges(3,1) + t43 * rSges(3,2) + rSges(3,3) * t101 + t88) + g(2) * (rSges(3,1) * t46 - rSges(3,2) * t45 + rSges(3,3) * t103 + t96)) - m(4) * (g(1) * (rSges(4,1) * t101 + t44 * rSges(4,2) - t43 * t95 + t82) + g(2) * (rSges(4,1) * t103 - rSges(4,2) * t46 + t45 * t95 + t90)) - m(5) * (g(1) * (rSges(5,1) * t79 - rSges(5,2) * t78 - t108 * t44 + t73) + g(2) * (rSges(5,1) * t19 - rSges(5,2) * t18 + t108 * t46 + t76)) - m(6) * (g(1) * (rSges(6,1) * t4 - rSges(6,2) * t3 + t107 * t78 + t70) + g(2) * (rSges(6,1) * t2 - rSges(6,2) * t1 + t107 * t18 + t71)) - m(7) * (g(1) * (t109 * t78 + t110 * t4 + t3 * t94 + t70) + g(2) * (t1 * t94 + t109 * t18 + t110 * t2 + t71)) -m(3) * (g(1) * (-rSges(3,1) * t45 - rSges(3,2) * t46) + g(2) * (-rSges(3,1) * t43 - rSges(3,2) * t44) + (rSges(3,1) * t68 - rSges(3,2) * t64) * t111) - m(4) * (g(1) * (rSges(4,2) * t45 + t46 * t95 - t39) + g(2) * (rSges(4,2) * t43 + t44 * t95 - t37) + g(3) * ((-rSges(4,2) * t68 + rSges(4,3) * t64) * t61 + t97)) - m(5) * (g(1) * (-rSges(5,3) * t45 + t86) + g(2) * (-rSges(5,3) * t43 + t87) + g(3) * t89 + (rSges(5,3) * t68 + t64 * t81) * t111 + t116 * (qJ(3) + t81)) - m(6) * (g(1) * (rSges(6,1) * t11 - rSges(6,2) * t10 + t74) + g(2) * (rSges(6,1) * t9 - rSges(6,2) * t8 + t75) + g(3) * (rSges(6,1) * t25 - rSges(6,2) * t24 + t83) + t107 * t115) - m(7) * (g(1) * (t94 * t10 + t110 * t11 + t74) + g(2) * (t110 * t9 + t94 * t8 + t75) + g(3) * (t110 * t25 + t94 * t24 + t83) + t109 * t115) (-m(4) - m(5) - m(6) - m(7)) * (g(1) * t45 + g(2) * t43 - g(3) * t102) -m(5) * (g(1) * (-rSges(5,1) * t18 - rSges(5,2) * t19) + g(2) * (rSges(5,1) * t78 + rSges(5,2) * t79) + g(3) * (rSges(5,1) * t41 - rSges(5,2) * t42)) - m(6) * (g(1) * (t107 * t19 - t18 * t80 - t12) + g(2) * (-t107 * t79 + t78 * t80 + t14) + g(3) * (t107 * t42 + t41 * t80 + t36)) - m(7) * ((t109 * t42 + t72 * t41 + t36) * g(3) + (-t109 * t79 + t72 * t78 + t14) * g(2) + (t109 * t19 - t72 * t18 - t12) * g(1)) -m(6) * (g(1) * (-rSges(6,1) * t1 - rSges(6,2) * t2) + g(2) * (rSges(6,1) * t3 + rSges(6,2) * t4) + g(3) * (-rSges(6,1) * t16 - rSges(6,2) * t17)) - m(7) * (g(1) * (-t1 * t110 + t2 * t94) + g(2) * (t110 * t3 - t4 * t94) + g(3) * (-t110 * t16 + t17 * t94)) -m(7) * (g(1) * t1 - g(2) * t3 + g(3) * t16)];
taug  = t5(:);
