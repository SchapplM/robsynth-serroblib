% Calculate Gravitation load on the joints for
% S6RPRPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-03-09 04:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPRR11_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR11_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR11_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR11_gravloadJ_floatb_twist_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR11_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR11_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:11:28
% EndTime: 2019-03-09 04:11:32
% DurationCPUTime: 1.49s
% Computational Cost: add. (1016->165), mult. (2403->248), div. (0->0), fcn. (3041->16), ass. (0->73)
t114 = cos(qJ(3));
t103 = cos(pkin(7));
t65 = cos(qJ(1));
t113 = sin(qJ(1));
t102 = cos(pkin(12));
t104 = cos(pkin(6));
t92 = t104 * t102;
t99 = sin(pkin(12));
t80 = t113 * t99 - t65 * t92;
t100 = sin(pkin(7));
t101 = sin(pkin(6));
t89 = t101 * t100;
t128 = t80 * t103 + t65 * t89;
t91 = t104 * t99;
t43 = t102 * t113 + t65 * t91;
t63 = sin(qJ(3));
t20 = t128 * t114 + t43 * t63;
t62 = sin(qJ(6));
t64 = cos(qJ(6));
t23 = -t114 * t43 + t128 * t63;
t90 = t103 * t101;
t34 = -t80 * t100 + t65 * t90;
t58 = pkin(13) + qJ(5);
t55 = sin(t58);
t56 = cos(t58);
t7 = t23 * t56 + t34 * t55;
t132 = t20 * t64 + t62 * t7;
t131 = -t20 * t62 + t64 * t7;
t127 = t23 * t55 - t34 * t56;
t74 = t113 * t92 + t65 * t99;
t67 = t74 * t100 + t113 * t90;
t124 = t74 * t103 - t113 * t89;
t44 = t102 * t65 - t113 * t91;
t24 = t114 * t124 + t44 * t63;
t123 = t100 * t104 + t102 * t90;
t88 = t101 * t99;
t32 = -t114 * t123 + t63 * t88;
t120 = g(1) * t24 + g(2) * t20 + g(3) * t32;
t116 = t56 * pkin(5);
t115 = rSges(7,3) + pkin(11);
t59 = sin(pkin(13));
t112 = t34 * t59;
t111 = t56 * t62;
t110 = t56 * t64;
t60 = cos(pkin(13));
t54 = pkin(4) * t60 + pkin(3);
t61 = -pkin(10) - qJ(4);
t109 = -t20 * t54 + t23 * t61;
t25 = t44 * t114 - t124 * t63;
t108 = -t24 * t54 - t25 * t61;
t33 = t114 * t88 + t123 * t63;
t107 = -t32 * t54 - t33 * t61;
t94 = t101 * t113;
t106 = t65 * pkin(1) + qJ(2) * t94;
t105 = qJ(4) + rSges(5,3);
t98 = -m(5) - m(6) - m(7);
t97 = t65 * t101;
t96 = -pkin(1) * t113 + qJ(2) * t97;
t93 = -rSges(6,1) * t56 + rSges(6,2) * t55;
t85 = t64 * rSges(7,1) - t62 * rSges(7,2) + pkin(5);
t71 = -t43 * pkin(2) + t34 * pkin(9) + t96;
t70 = t44 * pkin(2) + t67 * pkin(9) + t106;
t69 = pkin(4) * t112 + t20 * t61 + t23 * t54 + t71;
t66 = t67 * t59;
t68 = pkin(4) * t66 - t24 * t61 + t25 * t54 + t70;
t42 = -t102 * t89 + t103 * t104;
t15 = t33 * t56 + t42 * t55;
t14 = -t33 * t55 + t42 * t56;
t9 = t25 * t56 + t55 * t67;
t8 = t25 * t55 - t56 * t67;
t2 = t24 * t62 + t64 * t9;
t1 = t24 * t64 - t62 * t9;
t3 = [-m(2) * (g(1) * (-rSges(2,1) * t113 - t65 * rSges(2,2)) + g(2) * (t65 * rSges(2,1) - rSges(2,2) * t113)) - m(3) * (g(1) * (-t43 * rSges(3,1) + rSges(3,2) * t80 + rSges(3,3) * t97 + t96) + g(2) * (t44 * rSges(3,1) - rSges(3,2) * t74 + rSges(3,3) * t94 + t106)) - m(4) * (g(1) * (t23 * rSges(4,1) + rSges(4,2) * t20 + t34 * rSges(4,3) + t71) + g(2) * (t25 * rSges(4,1) - t24 * rSges(4,2) + rSges(4,3) * t67 + t70)) - m(5) * (g(1) * (t23 * pkin(3) + (t23 * t60 + t112) * rSges(5,1) + (-t23 * t59 + t34 * t60) * rSges(5,2) - t105 * t20 + t71) + g(2) * (t25 * pkin(3) + (t25 * t60 + t66) * rSges(5,1) + (-t25 * t59 + t60 * t67) * rSges(5,2) + t105 * t24 + t70)) - m(6) * (g(1) * (t7 * rSges(6,1) - rSges(6,2) * t127 - rSges(6,3) * t20 + t69) + g(2) * (t9 * rSges(6,1) - t8 * rSges(6,2) + t24 * rSges(6,3) + t68)) - m(7) * (g(1) * (t131 * rSges(7,1) - t132 * rSges(7,2) + t7 * pkin(5) + t115 * t127 + t69) + g(2) * (t2 * rSges(7,1) + t1 * rSges(7,2) + t9 * pkin(5) + t115 * t8 + t68)) (-m(3) - m(4) + t98) * (g(1) * t94 - g(2) * t97 + g(3) * t104) -m(4) * (g(1) * (-rSges(4,1) * t24 - rSges(4,2) * t25) + g(2) * (-rSges(4,1) * t20 + rSges(4,2) * t23) + g(3) * (-rSges(4,1) * t32 - rSges(4,2) * t33)) - m(5) * (t120 * (-t60 * rSges(5,1) + t59 * rSges(5,2) - pkin(3)) + (g(1) * t25 - g(2) * t23 + g(3) * t33) * t105) - m(6) * (g(1) * (rSges(6,3) * t25 + t93 * t24 + t108) + g(2) * (-rSges(6,3) * t23 + t93 * t20 + t109) + g(3) * (rSges(6,3) * t33 + t93 * t32 + t107)) + (-g(1) * (-t24 * t116 + (-t24 * t110 + t25 * t62) * rSges(7,1) + (t24 * t111 + t25 * t64) * rSges(7,2) + t108) - g(2) * (-t20 * t116 + (-t20 * t110 - t23 * t62) * rSges(7,1) + (t20 * t111 - t23 * t64) * rSges(7,2) + t109) - g(3) * (-t32 * t116 + (-t32 * t110 + t33 * t62) * rSges(7,1) + (t32 * t111 + t33 * t64) * rSges(7,2) + t107) + t120 * t55 * t115) * m(7), t98 * t120, -m(6) * (g(1) * (-rSges(6,1) * t8 - rSges(6,2) * t9) + g(2) * (rSges(6,1) * t127 + rSges(6,2) * t7) + g(3) * (rSges(6,1) * t14 - rSges(6,2) * t15)) - m(7) * (g(1) * (t115 * t9 - t8 * t85) + (t115 * t15 + t85 * t14) * g(3) + (-t115 * t7 + t127 * t85) * g(2)) -m(7) * (g(1) * (rSges(7,1) * t1 - rSges(7,2) * t2) + g(2) * (t132 * rSges(7,1) + t131 * rSges(7,2)) + g(3) * ((-t15 * t62 + t32 * t64) * rSges(7,1) + (-t15 * t64 - t32 * t62) * rSges(7,2)))];
taug  = t3(:);
