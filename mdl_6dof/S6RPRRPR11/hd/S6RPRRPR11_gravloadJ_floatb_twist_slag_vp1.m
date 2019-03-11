% Calculate Gravitation load on the joints for
% S6RPRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 05:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRPR11_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR11_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR11_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR11_gravloadJ_floatb_twist_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR11_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR11_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:39:15
% EndTime: 2019-03-09 05:39:20
% DurationCPUTime: 1.73s
% Computational Cost: add. (1176->187), mult. (3033->277), div. (0->0), fcn. (3874->16), ass. (0->81)
t114 = cos(qJ(1));
t112 = sin(qJ(1));
t100 = cos(pkin(6));
t98 = cos(pkin(12));
t85 = t100 * t98;
t95 = sin(pkin(12));
t72 = t112 * t95 - t114 * t85;
t96 = sin(pkin(7));
t97 = sin(pkin(6));
t82 = t97 * t96;
t99 = cos(pkin(7));
t131 = t114 * t82 + t72 * t99;
t113 = cos(qJ(3));
t84 = t100 * t95;
t35 = t112 * t98 + t114 * t84;
t55 = sin(qJ(3));
t19 = -t113 * t35 + t131 * t55;
t54 = sin(qJ(4));
t56 = cos(qJ(4));
t83 = t99 * t97;
t62 = -t114 * t83 + t72 * t96;
t7 = t19 * t56 - t54 * t62;
t6 = t19 * t54 + t56 * t62;
t67 = t112 * t85 + t114 * t95;
t57 = t112 * t83 + t67 * t96;
t16 = t113 * t131 + t35 * t55;
t129 = -t112 * t82 + t67 * t99;
t128 = t96 * t100 + t98 * t83;
t81 = t97 * t95;
t27 = t113 * t81 + t128 * t55;
t66 = t100 * t99 - t82 * t98;
t15 = t27 * t56 + t54 * t66;
t36 = -t112 * t84 + t114 * t98;
t21 = t36 * t113 - t129 * t55;
t9 = t21 * t56 + t54 * t57;
t127 = g(1) * t9 - g(2) * t7 + g(3) * t15;
t14 = t27 * t54 - t56 * t66;
t8 = t21 * t54 - t56 * t57;
t126 = g(1) * t8 - g(2) * t6 + g(3) * t14;
t20 = t113 * t129 + t36 * t55;
t26 = -t113 * t128 + t55 * t81;
t125 = (-g(1) * t20 - g(2) * t16 - g(3) * t26) * t54;
t122 = -m(6) - m(7);
t116 = t56 * pkin(4);
t115 = pkin(10) + rSges(5,3);
t51 = sin(pkin(13));
t111 = t19 * t51;
t110 = t21 * t51;
t109 = t27 * t51;
t50 = pkin(13) + qJ(6);
t47 = sin(t50);
t108 = t47 * t56;
t48 = cos(t50);
t107 = t48 * t56;
t106 = t51 * t56;
t52 = cos(pkin(13));
t105 = t52 * t56;
t46 = pkin(5) * t52 + pkin(4);
t104 = t56 * t46;
t87 = t97 * t112;
t103 = t114 * pkin(1) + qJ(2) * t87;
t102 = pkin(11) + qJ(5) + rSges(7,3);
t101 = qJ(5) + rSges(6,3);
t94 = t51 * pkin(5) + pkin(10);
t10 = t16 * pkin(3);
t93 = -pkin(10) * t19 - t10;
t12 = t20 * pkin(3);
t92 = t21 * pkin(10) - t12;
t25 = t26 * pkin(3);
t91 = t27 * pkin(10) - t25;
t88 = t114 * t97;
t90 = -pkin(1) * t112 + qJ(2) * t88;
t86 = -rSges(5,1) * t56 + rSges(5,2) * t54;
t78 = t48 * rSges(7,1) - t47 * rSges(7,2) + t46;
t63 = -t35 * pkin(2) - t62 * pkin(9) + t90;
t60 = t19 * pkin(3) + t63;
t59 = t36 * pkin(2) + t57 * pkin(9) + t103;
t58 = t21 * pkin(3) + t59;
t3 = t20 * t47 + t48 * t9;
t2 = t20 * t48 - t47 * t9;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t112 - rSges(2,2) * t114) + g(2) * (rSges(2,1) * t114 - rSges(2,2) * t112)) - m(3) * (g(1) * (-t35 * rSges(3,1) + rSges(3,2) * t72 + rSges(3,3) * t88 + t90) + g(2) * (t36 * rSges(3,1) - rSges(3,2) * t67 + rSges(3,3) * t87 + t103)) - m(4) * (g(1) * (t19 * rSges(4,1) + rSges(4,2) * t16 - rSges(4,3) * t62 + t63) + g(2) * (t21 * rSges(4,1) - t20 * rSges(4,2) + rSges(4,3) * t57 + t59)) - m(5) * (g(1) * (t7 * rSges(5,1) - t6 * rSges(5,2) - t115 * t16 + t60) + g(2) * (t9 * rSges(5,1) - t8 * rSges(5,2) + t115 * t20 + t58)) - m(6) * (g(1) * (t7 * pkin(4) - t16 * pkin(10) + (-t16 * t51 + t52 * t7) * rSges(6,1) + (-t16 * t52 - t51 * t7) * rSges(6,2) + t101 * t6 + t60) + g(2) * (t20 * pkin(10) + t9 * pkin(4) + (t20 * t51 + t52 * t9) * rSges(6,1) + (t20 * t52 - t51 * t9) * rSges(6,2) + t101 * t8 + t58)) - m(7) * (g(1) * (t102 * t6 + t78 * t7 - (t47 * rSges(7,1) + t48 * rSges(7,2) + t94) * t16 + t60) + g(2) * (t3 * rSges(7,1) + t2 * rSges(7,2) + t102 * t8 + t20 * t94 + t9 * t46 + t58)) (-m(3) - m(4) - m(5) + t122) * (g(1) * t87 - g(2) * t88 + g(3) * t100) -m(4) * (g(1) * (-t20 * rSges(4,1) - t21 * rSges(4,2)) + g(2) * (-t16 * rSges(4,1) + rSges(4,2) * t19) + g(3) * (-t26 * rSges(4,1) - t27 * rSges(4,2))) - m(5) * (g(1) * (t115 * t21 + t20 * t86 - t12) + g(2) * (-t115 * t19 + t16 * t86 - t10) + g(3) * (t115 * t27 + t26 * t86 - t25)) + (-g(1) * (-t20 * t104 + pkin(5) * t110 + (-t107 * t20 + t21 * t47) * rSges(7,1) + (t108 * t20 + t21 * t48) * rSges(7,2) + t92) - g(2) * (-t16 * t104 - pkin(5) * t111 + (-t107 * t16 - t19 * t47) * rSges(7,1) + (t108 * t16 - t19 * t48) * rSges(7,2) + t93) - g(3) * (-t26 * t104 + pkin(5) * t109 + (-t107 * t26 + t27 * t47) * rSges(7,1) + (t108 * t26 + t27 * t48) * rSges(7,2) + t91) - t102 * t125) * m(7) + (-g(1) * (-t20 * t116 + (-t105 * t20 + t110) * rSges(6,1) + (t106 * t20 + t21 * t52) * rSges(6,2) + t92) - g(2) * (-t16 * t116 + (-t105 * t16 - t111) * rSges(6,1) + (t106 * t16 - t19 * t52) * rSges(6,2) + t93) - g(3) * (-t26 * t116 + (-t105 * t26 + t109) * rSges(6,1) + (t106 * t26 + t27 * t52) * rSges(6,2) + t91) - t101 * t125) * m(6), -m(5) * (g(1) * (-rSges(5,1) * t8 - rSges(5,2) * t9) + g(2) * (rSges(5,1) * t6 + rSges(5,2) * t7) + g(3) * (-rSges(5,1) * t14 - rSges(5,2) * t15)) - m(6) * (t126 * (-t52 * rSges(6,1) + t51 * rSges(6,2) - pkin(4)) + t127 * t101) - m(7) * (t102 * t127 - t126 * t78) t122 * t126, -m(7) * (g(1) * (rSges(7,1) * t2 - rSges(7,2) * t3) + g(2) * ((t16 * t48 + t47 * t7) * rSges(7,1) + (-t16 * t47 + t48 * t7) * rSges(7,2)) + g(3) * ((-t15 * t47 + t26 * t48) * rSges(7,1) + (-t15 * t48 - t26 * t47) * rSges(7,2)))];
taug  = t1(:);
