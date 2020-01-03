% Calculate Gravitation load on the joints for
% S5RRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d2,d3,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRR12_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR12_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR12_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR12_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR12_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR12_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:46:36
% EndTime: 2019-12-31 22:46:44
% DurationCPUTime: 1.54s
% Computational Cost: add. (926->143), mult. (2514->220), div. (0->0), fcn. (3167->14), ass. (0->77)
t140 = m(5) + m(6);
t145 = t140 * pkin(10);
t72 = sin(qJ(5));
t75 = cos(qJ(5));
t139 = m(6) * pkin(4) + t75 * mrSges(6,1) - t72 * mrSges(6,2) + mrSges(5,1);
t136 = -m(6) * pkin(11) + mrSges(5,2) - mrSges(6,3);
t115 = cos(pkin(6));
t71 = sin(pkin(5));
t108 = t71 * t115;
t125 = sin(qJ(1));
t70 = sin(pkin(6));
t124 = sin(qJ(2));
t77 = cos(qJ(1));
t116 = cos(pkin(5));
t127 = cos(qJ(2));
t99 = t116 * t127;
t84 = t124 * t77 + t125 * t99;
t144 = t125 * t108 + t84 * t70;
t138 = mrSges(4,2) - mrSges(5,3);
t94 = -t72 * mrSges(6,1) - t75 * mrSges(6,2);
t130 = t94 + t138;
t53 = t124 * t125 - t77 * t99;
t143 = t77 * t108 - t53 * t70;
t73 = sin(qJ(4));
t76 = cos(qJ(4));
t142 = pkin(3) * t140 - t136 * t73 + t139 * t76 + mrSges(4,1);
t74 = sin(qJ(3));
t107 = t74 * t115;
t119 = t71 * t77;
t114 = t70 * t119;
t126 = cos(qJ(3));
t98 = t116 * t124;
t54 = t125 * t127 + t77 * t98;
t20 = -t107 * t53 - t74 * t114 + t126 * t54;
t4 = -t143 * t73 + t20 * t76;
t141 = -t143 * t76 - t20 * t73;
t137 = -t70 * mrSges(4,3) + mrSges(3,2);
t111 = t71 * t125;
t135 = t70 * t111 - t84 * t115;
t131 = t130 - t145;
t128 = pkin(9) * t70;
t121 = t70 * t73;
t120 = t70 * t76;
t110 = t71 * t124;
t105 = t70 * t110;
t112 = t71 * t127;
t118 = pkin(2) * t112 + pkin(9) * t105;
t117 = t77 * pkin(1) + pkin(8) * t111;
t109 = t70 * t116;
t103 = -t53 * pkin(2) + t128 * t54;
t55 = -t125 * t98 + t127 * t77;
t102 = -t84 * pkin(2) + t128 * t55;
t100 = -pkin(1) * t125 + pkin(8) * t119;
t97 = t115 * t126;
t96 = t115 * t124;
t92 = t138 - t145;
t86 = -t54 * pkin(2) + pkin(9) * t143 + t100;
t19 = t114 * t126 + t53 * t97 + t54 * t74;
t80 = t55 * pkin(2) + pkin(9) * t144 + t117;
t24 = t55 * t126 + t135 * t74;
t79 = t24 * pkin(3) + t80;
t52 = -t112 * t70 + t115 * t116;
t47 = (t126 * t127 - t74 * t96) * t71;
t46 = (t126 * t96 + t127 * t74) * t71;
t37 = t74 * t109 + (t107 * t127 + t124 * t126) * t71;
t36 = -t109 * t126 + t110 * t74 - t112 * t97;
t30 = -t107 * t55 - t126 * t84;
t29 = t55 * t97 - t74 * t84;
t28 = -t107 * t54 - t126 * t53;
t27 = -t53 * t74 + t54 * t97;
t23 = -t135 * t126 + t55 * t74;
t18 = t37 * t76 + t52 * t73;
t8 = t144 * t73 + t24 * t76;
t7 = -t144 * t76 + t24 * t73;
t2 = t23 * t72 + t75 * t8;
t1 = t23 * t75 - t72 * t8;
t3 = [(-t77 * mrSges(2,1) + t125 * mrSges(2,2) - m(3) * t117 - t55 * mrSges(3,1) + t84 * mrSges(3,2) - mrSges(3,3) * t111 - m(4) * t80 - t24 * mrSges(4,1) - t144 * mrSges(4,3) - m(5) * t79 - t8 * mrSges(5,1) - m(6) * (t8 * pkin(4) + t79) - t2 * mrSges(6,1) - t1 * mrSges(6,2) + t136 * t7 + t92 * t23) * g(2) + (t125 * mrSges(2,1) - m(3) * t100 + t54 * mrSges(3,1) - t53 * mrSges(3,2) - m(4) * t86 + t20 * mrSges(4,1) - t143 * mrSges(4,3) + (-t71 * mrSges(3,3) + mrSges(2,2)) * t77 + t136 * t141 + t139 * t4 - (t92 + t94) * t19 + t140 * (pkin(3) * t20 - t86)) * g(1), (-(mrSges(3,1) * t127 - mrSges(3,2) * t124) * t71 - m(4) * t118 - t47 * mrSges(4,1) - mrSges(4,3) * t105 - t140 * (t47 * pkin(3) + pkin(10) * t46 + t118) - t139 * (t105 * t73 + t47 * t76) + t130 * t46 + t136 * (-t105 * t76 + t47 * t73)) * g(3) + (-m(4) * t103 + mrSges(3,1) * t53 - t28 * mrSges(4,1) - t140 * (t28 * pkin(3) + pkin(10) * t27 + t103) + t136 * (-t120 * t54 + t28 * t73) + t137 * t54 - t139 * (t121 * t54 + t28 * t76) + t130 * t27) * g(2) + (-m(4) * t102 + t84 * mrSges(3,1) - t30 * mrSges(4,1) - t140 * (t30 * pkin(3) + pkin(10) * t29 + t102) + t137 * t55 - t139 * (t121 * t55 + t30 * t76) + t130 * t29 + t136 * (-t120 * t55 + t30 * t73)) * g(1), (t131 * t37 + t142 * t36) * g(3) + (t131 * t20 + t142 * t19) * g(2) + (t131 * t24 + t142 * t23) * g(1), (t136 * t18 - t139 * (-t37 * t73 + t52 * t76)) * g(3) + (t136 * t4 - t139 * t141) * g(2) + (t136 * t8 + t139 * t7) * g(1), -g(1) * (mrSges(6,1) * t1 - mrSges(6,2) * t2) - g(2) * ((t19 * t75 - t4 * t72) * mrSges(6,1) + (-t19 * t72 - t4 * t75) * mrSges(6,2)) - g(3) * ((-t18 * t72 + t36 * t75) * mrSges(6,1) + (-t18 * t75 - t36 * t72) * mrSges(6,2))];
taug = t3(:);
