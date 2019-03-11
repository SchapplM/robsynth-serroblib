% Calculate Gravitation load on the joints for
% S6PRRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
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
% Datum: 2019-03-08 23:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRPR6_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR6_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR6_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR6_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR6_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRPR6_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:31:12
% EndTime: 2019-03-08 23:31:16
% DurationCPUTime: 1.08s
% Computational Cost: add. (693->190), mult. (1766->276), div. (0->0), fcn. (2182->12), ass. (0->83)
t61 = sin(qJ(4));
t65 = cos(qJ(4));
t117 = pkin(4) * t65 + qJ(5) * t61;
t59 = sin(pkin(6));
t67 = cos(qJ(2));
t101 = t59 * t67;
t66 = cos(qJ(3));
t102 = t59 * t66;
t62 = sin(qJ(3));
t63 = sin(qJ(2));
t94 = cos(pkin(6));
t48 = t63 * t102 + t94 * t62;
t25 = t65 * t101 + t48 * t61;
t90 = t61 * t101;
t26 = t48 * t65 - t90;
t60 = sin(qJ(6));
t64 = cos(qJ(6));
t116 = (t25 * t64 - t26 * t60) * rSges(7,1) - (t25 * t60 + t26 * t64) * rSges(7,2);
t58 = sin(pkin(11));
t86 = t58 * t94;
t93 = cos(pkin(11));
t46 = -t63 * t86 + t93 * t67;
t24 = t58 * t59 * t62 + t46 * t66;
t45 = t93 * t63 + t67 * t86;
t6 = t24 * t61 - t45 * t65;
t7 = t24 * t65 + t45 * t61;
t115 = (t6 * t64 - t60 * t7) * rSges(7,1) - (t6 * t60 + t64 * t7) * rSges(7,2);
t71 = t94 * t93;
t44 = t58 * t67 + t63 * t71;
t85 = t59 * t93;
t22 = t44 * t66 - t62 * t85;
t43 = t58 * t63 - t67 * t71;
t4 = t22 * t61 - t43 * t65;
t5 = t22 * t65 + t43 * t61;
t114 = (t4 * t64 - t5 * t60) * rSges(7,1) - (t4 * t60 + t5 * t64) * rSges(7,2);
t113 = -m(6) - m(7);
t112 = pkin(3) * t66;
t110 = g(3) * t59;
t109 = rSges(6,2) + pkin(9);
t108 = rSges(4,3) + pkin(8);
t107 = rSges(5,3) + pkin(9);
t106 = rSges(7,3) + pkin(10);
t105 = t43 * t62;
t104 = t45 * t62;
t103 = t59 * t63;
t100 = t61 * t66;
t99 = t65 * t66;
t98 = t66 * t67;
t97 = pkin(2) * t101 + pkin(8) * t103;
t95 = rSges(6,3) + qJ(5);
t92 = pkin(9) - t106;
t91 = t62 * t101;
t21 = -t44 * t62 - t66 * t85;
t18 = t21 * pkin(3);
t89 = t117 * t21 + t18;
t23 = t58 * t102 - t46 * t62;
t19 = t23 * pkin(3);
t88 = t117 * t23 + t19;
t47 = -t62 * t103 + t94 * t66;
t42 = t47 * pkin(3);
t87 = t117 * t47 + t42;
t84 = t59 * pkin(3) * t98 + pkin(9) * t91 + t97;
t29 = (t61 * t63 + t65 * t98) * t59;
t83 = t29 * pkin(4) + t84;
t78 = rSges(4,1) * t66 - rSges(4,2) * t62;
t77 = rSges(5,1) * t65 - rSges(5,2) * t61;
t76 = rSges(6,1) * t65 + rSges(6,3) * t61;
t40 = t43 * pkin(2);
t73 = t44 * pkin(8) - pkin(9) * t105 - t43 * t112 - t40;
t41 = t45 * pkin(2);
t72 = t46 * pkin(8) - pkin(9) * t104 - t45 * t112 - t41;
t11 = -t43 * t99 + t44 * t61;
t70 = t11 * pkin(4) + t73;
t13 = -t45 * t99 + t46 * t61;
t69 = t13 * pkin(4) + t72;
t68 = t65 * pkin(5) + (t60 * t61 + t64 * t65) * rSges(7,1) + (-t60 * t65 + t61 * t64) * rSges(7,2);
t28 = -t65 * t103 + t66 * t90;
t20 = t25 * pkin(4);
t12 = -t45 * t100 - t46 * t65;
t10 = -t43 * t100 - t44 * t65;
t3 = t6 * pkin(4);
t2 = t4 * pkin(4);
t1 = [(-m(2) - m(3) - m(4) - m(5) + t113) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t45 - rSges(3,2) * t46) + g(2) * (-rSges(3,1) * t43 - rSges(3,2) * t44) + (rSges(3,1) * t67 - rSges(3,2) * t63) * t110) - m(4) * (g(1) * (t108 * t46 - t78 * t45 - t41) + g(2) * (t108 * t44 - t78 * t43 - t40) + g(3) * t97 + (rSges(4,3) * t63 + t78 * t67) * t110) - m(5) * (g(1) * (rSges(5,1) * t13 - rSges(5,2) * t12 - rSges(5,3) * t104 + t72) + g(2) * (rSges(5,1) * t11 - rSges(5,2) * t10 - rSges(5,3) * t105 + t73) + g(3) * (t29 * rSges(5,1) - t28 * rSges(5,2) + rSges(5,3) * t91 + t84)) - m(6) * (g(1) * (rSges(6,1) * t13 - rSges(6,2) * t104 + t95 * t12 + t69) + g(2) * (rSges(6,1) * t11 - rSges(6,2) * t105 + t95 * t10 + t70) + g(3) * (t29 * rSges(6,1) + rSges(6,2) * t91 + t95 * t28 + t83)) - m(7) * (g(1) * (t13 * pkin(5) + t12 * qJ(5) + (t12 * t60 + t13 * t64) * rSges(7,1) + (t12 * t64 - t13 * t60) * rSges(7,2) + t69) + g(2) * (t11 * pkin(5) + t10 * qJ(5) + (t10 * t60 + t11 * t64) * rSges(7,1) + (t10 * t64 - t11 * t60) * rSges(7,2) + t70) + g(3) * (t29 * pkin(5) + t28 * qJ(5) + (t28 * t60 + t29 * t64) * rSges(7,1) + (t28 * t64 - t29 * t60) * rSges(7,2) + t83) + (g(1) * t45 + g(2) * t43 - g(3) * t101) * t62 * t106) (-m(4) * (rSges(4,1) * t47 - rSges(4,2) * t48) - m(5) * (t107 * t48 + t77 * t47 + t42) - m(6) * (t109 * t48 + t76 * t47 + t87) - m(7) * (t47 * t68 + t48 * t92 + t87)) * g(3) + (-m(4) * (rSges(4,1) * t21 - rSges(4,2) * t22) - m(5) * (t107 * t22 + t77 * t21 + t18) - m(6) * (t109 * t22 + t76 * t21 + t89) - m(7) * (t21 * t68 + t22 * t92 + t89)) * g(2) + (-m(4) * (rSges(4,1) * t23 - rSges(4,2) * t24) - m(5) * (t107 * t24 + t77 * t23 + t19) - m(6) * (t109 * t24 + t76 * t23 + t88) - m(7) * (t23 * t68 + t24 * t92 + t88)) * g(1), -m(5) * (g(1) * (-rSges(5,1) * t6 - rSges(5,2) * t7) + g(2) * (-rSges(5,1) * t4 - rSges(5,2) * t5) + g(3) * (-rSges(5,1) * t25 - rSges(5,2) * t26)) - m(6) * (g(1) * (-rSges(6,1) * t6 + t95 * t7 - t3) + g(2) * (-rSges(6,1) * t4 + t95 * t5 - t2) + g(3) * (-rSges(6,1) * t25 + t95 * t26 - t20)) - m(7) * (g(1) * (-t6 * pkin(5) + t7 * qJ(5) - t115 - t3) + g(2) * (-t4 * pkin(5) + t5 * qJ(5) - t114 - t2) + g(3) * (-t25 * pkin(5) + t26 * qJ(5) - t116 - t20)) t113 * (g(1) * t6 + g(2) * t4 + g(3) * t25) -m(7) * (g(1) * t115 + g(2) * t114 + g(3) * t116)];
taug  = t1(:);
