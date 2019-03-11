% Calculate Gravitation load on the joints for
% S6RPRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2]';
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
% Datum: 2019-03-09 05:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRPR12_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR12_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR12_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRPR12_gravloadJ_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR12_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR12_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:48:09
% EndTime: 2019-03-09 05:48:13
% DurationCPUTime: 1.39s
% Computational Cost: add. (1099->170), mult. (2944->243), div. (0->0), fcn. (3757->14), ass. (0->76)
t106 = cos(pkin(7));
t115 = cos(qJ(1));
t102 = sin(pkin(12));
t113 = sin(qJ(1));
t105 = cos(pkin(12));
t107 = cos(pkin(6));
t92 = t107 * t105;
t80 = t102 * t113 - t115 * t92;
t103 = sin(pkin(7));
t104 = sin(pkin(6));
t89 = t104 * t103;
t134 = t106 * t80 + t115 * t89;
t114 = cos(qJ(3));
t91 = t107 * t102;
t46 = t105 * t113 + t115 * t91;
t60 = sin(qJ(3));
t28 = -t46 * t114 + t134 * t60;
t59 = sin(qJ(4));
t62 = cos(qJ(4));
t90 = t106 * t104;
t69 = -t103 * t80 + t115 * t90;
t10 = t28 * t59 - t62 * t69;
t11 = t28 * t62 + t59 * t69;
t58 = sin(qJ(6));
t61 = cos(qJ(6));
t131 = -rSges(7,1) * t58 - rSges(7,2) * t61;
t75 = t102 * t115 + t113 * t92;
t63 = t103 * t75 + t113 * t90;
t128 = -pkin(4) * t62 - qJ(5) * t59;
t120 = pkin(10) + pkin(5);
t122 = rSges(7,1) * t61 - rSges(7,2) * t58 + t120;
t25 = t114 * t134 + t46 * t60;
t116 = pkin(11) + rSges(7,3);
t127 = t106 * t75 - t113 * t89;
t126 = t103 * t107 + t105 * t90;
t123 = t131 * t59;
t121 = -m(6) - m(7);
t118 = pkin(10) + rSges(6,1);
t117 = pkin(10) + rSges(5,3);
t95 = t104 * t113;
t110 = pkin(1) * t115 + qJ(2) * t95;
t108 = qJ(5) + rSges(6,3);
t19 = t25 * pkin(3);
t101 = t128 * t25 - t19;
t47 = t105 * t115 - t113 * t91;
t29 = t114 * t127 + t47 * t60;
t21 = t29 * pkin(3);
t100 = t128 * t29 - t21;
t88 = t104 * t102;
t37 = -t114 * t126 + t60 * t88;
t36 = t37 * pkin(3);
t99 = t128 * t37 - t36;
t96 = t115 * t104;
t98 = -pkin(1) * t113 + qJ(2) * t96;
t94 = -rSges(5,1) * t62 + rSges(5,2) * t59;
t93 = rSges(6,2) * t62 - rSges(6,3) * t59;
t86 = qJ(5) - t131;
t74 = -t105 * t89 + t106 * t107;
t71 = -t46 * pkin(2) + pkin(9) * t69 + t98;
t68 = pkin(3) * t28 + t71;
t67 = pkin(4) * t11 + t68;
t66 = t47 * pkin(2) + pkin(9) * t63 + t110;
t30 = t114 * t47 - t127 * t60;
t65 = pkin(3) * t30 + t66;
t13 = t30 * t62 + t59 * t63;
t64 = pkin(4) * t13 + t65;
t38 = t114 * t88 + t126 * t60;
t24 = t38 * t62 + t59 * t74;
t23 = t38 * t59 - t62 * t74;
t18 = t23 * pkin(4);
t12 = t30 * t59 - t62 * t63;
t6 = t12 * pkin(4);
t4 = t10 * pkin(4);
t3 = t12 * t58 + t29 * t61;
t2 = t12 * t61 - t29 * t58;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t113 - rSges(2,2) * t115) + g(2) * (rSges(2,1) * t115 - rSges(2,2) * t113)) - m(3) * (g(1) * (-t46 * rSges(3,1) + rSges(3,2) * t80 + rSges(3,3) * t96 + t98) + g(2) * (t47 * rSges(3,1) - rSges(3,2) * t75 + rSges(3,3) * t95 + t110)) - m(4) * (g(1) * (t28 * rSges(4,1) + rSges(4,2) * t25 + rSges(4,3) * t69 + t71) + g(2) * (t30 * rSges(4,1) - t29 * rSges(4,2) + rSges(4,3) * t63 + t66)) - m(5) * (g(1) * (t11 * rSges(5,1) - t10 * rSges(5,2) - t117 * t25 + t68) + g(2) * (t13 * rSges(5,1) - t12 * rSges(5,2) + t117 * t29 + t65)) - m(6) * (g(1) * (-t11 * rSges(6,2) + t10 * t108 - t118 * t25 + t67) + g(2) * (-t13 * rSges(6,2) + t108 * t12 + t118 * t29 + t64)) - m(7) * (g(1) * (t86 * t10 + t116 * t11 - t122 * t25 + t67) + g(2) * (t3 * rSges(7,1) + t2 * rSges(7,2) + t12 * qJ(5) + t116 * t13 + t120 * t29 + t64)) (-m(3) - m(4) - m(5) + t121) * (g(1) * t95 - g(2) * t96 + g(3) * t107) -m(4) * (g(1) * (-rSges(4,1) * t29 - rSges(4,2) * t30) + g(2) * (-rSges(4,1) * t25 + rSges(4,2) * t28) + g(3) * (-rSges(4,1) * t37 - rSges(4,2) * t38)) - m(5) * (g(1) * (t117 * t30 + t29 * t94 - t21) + g(2) * (-t117 * t28 + t25 * t94 - t19) + g(3) * (t117 * t38 + t37 * t94 - t36)) - m(6) * (g(1) * (t118 * t30 + t29 * t93 + t100) + g(2) * (-t118 * t28 + t25 * t93 + t101) + g(3) * (t118 * t38 + t37 * t93 + t99)) + (-g(1) * (t122 * t30 + t123 * t29 + t100) - g(2) * (-t122 * t28 + t123 * t25 + t101) - g(3) * (t122 * t38 + t123 * t37 + t99) - (-g(1) * t29 - g(2) * t25 - g(3) * t37) * t62 * t116) * m(7), -m(5) * (g(1) * (-rSges(5,1) * t12 - rSges(5,2) * t13) + g(2) * (rSges(5,1) * t10 + rSges(5,2) * t11) + g(3) * (-rSges(5,1) * t23 - rSges(5,2) * t24)) - m(6) * (g(1) * (rSges(6,2) * t12 + t108 * t13 - t6) + g(2) * (-rSges(6,2) * t10 - t108 * t11 + t4) + g(3) * (rSges(6,2) * t23 + t108 * t24 - t18)) + (-g(1) * (-t116 * t12 - t6) - g(2) * (t10 * t116 + t4) - g(3) * (-t116 * t23 - t18) - (g(1) * t13 - g(2) * t11 + g(3) * t24) * t86) * m(7), t121 * (g(1) * t12 - g(2) * t10 + g(3) * t23) -m(7) * (g(1) * (rSges(7,1) * t2 - rSges(7,2) * t3) + g(2) * ((-t10 * t61 - t25 * t58) * rSges(7,1) + (t10 * t58 - t25 * t61) * rSges(7,2)) + g(3) * ((t23 * t61 - t37 * t58) * rSges(7,1) + (-t23 * t58 - t37 * t61) * rSges(7,2)))];
taug  = t1(:);
