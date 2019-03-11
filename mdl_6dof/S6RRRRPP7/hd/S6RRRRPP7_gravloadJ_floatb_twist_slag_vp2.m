% Calculate Gravitation load on the joints for
% S6RRRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,theta5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPP7_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP7_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP7_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPP7_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP7_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP7_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:17:20
% EndTime: 2019-03-09 21:17:23
% DurationCPUTime: 1.57s
% Computational Cost: add. (881->145), mult. (1838->195), div. (0->0), fcn. (2191->12), ass. (0->73)
t143 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t159 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t73 = qJ(4) + pkin(11);
t70 = sin(t73);
t71 = cos(t73);
t164 = t143 * t71 + t159 * t70;
t76 = sin(qJ(4));
t79 = cos(qJ(4));
t148 = -m(5) * pkin(3) - t79 * mrSges(5,1) + t76 * mrSges(5,2) - mrSges(4,1);
t150 = m(6) + m(7);
t69 = pkin(4) * t79 + pkin(3);
t162 = t150 * t69;
t161 = -m(5) * pkin(10) + mrSges(4,2) - mrSges(5,3);
t160 = mrSges(6,3) + mrSges(7,2);
t158 = -m(4) - t150;
t77 = sin(qJ(3));
t80 = cos(qJ(3));
t157 = t148 * t80 + t161 * t77 - mrSges(3,1);
t121 = mrSges(3,2) - mrSges(4,3);
t97 = mrSges(5,1) * t76 + mrSges(5,2) * t79;
t156 = -t97 + t121;
t139 = sin(qJ(1));
t74 = sin(pkin(6));
t108 = t74 * t139;
t140 = cos(qJ(1));
t78 = sin(qJ(2));
t81 = cos(qJ(2));
t115 = cos(pkin(6));
t99 = t115 * t139;
t55 = t140 * t81 - t78 * t99;
t28 = t108 * t77 + t55 * t80;
t54 = t140 * t78 + t81 * t99;
t7 = -t28 * t76 + t54 * t79;
t109 = t74 * t140;
t100 = t115 * t140;
t53 = t100 * t78 + t139 * t81;
t24 = -t109 * t77 + t53 * t80;
t52 = -t100 * t81 + t139 * t78;
t155 = -t24 * t76 + t52 * t79;
t1 = t24 * t70 - t52 * t71;
t154 = t24 * t71 + t52 * t70;
t141 = pkin(4) * t76;
t153 = -m(5) * pkin(9) - t141 * t150 - t143 * t70 + t159 * t71 + t156;
t75 = -qJ(5) - pkin(10);
t152 = -t157 + (-t150 * t75 + t160) * t77 + (t162 + t164) * t80;
t151 = m(4) + m(5);
t145 = -pkin(9) * (t150 + t151) + t121;
t144 = -t148 + t164;
t92 = -t160 + t161;
t134 = t52 * t76;
t132 = t54 * t76;
t127 = t74 * t78;
t126 = t74 * t81;
t122 = t80 * t81;
t117 = pkin(2) * t126 + pkin(9) * t127;
t116 = pkin(1) * t140 + pkin(8) * t108;
t114 = t77 * t126;
t113 = t70 * t126;
t112 = pkin(2) * t55 + t116;
t101 = -pkin(1) * t139 + pkin(8) * t109;
t95 = pkin(2) * t53 - t101;
t51 = t115 * t77 + t127 * t80;
t93 = -t126 * t79 - t51 * t76;
t23 = t109 * t80 + t53 * t77;
t50 = -t115 * t80 + t127 * t77;
t48 = t54 * pkin(2);
t46 = t52 * pkin(2);
t27 = -t108 * t80 + t55 * t77;
t17 = t126 * t71 + t51 * t70;
t8 = t28 * t79 + t132;
t6 = t28 * t71 + t54 * t70;
t5 = t28 * t70 - t54 * t71;
t2 = [(-t140 * mrSges(2,1) + t139 * mrSges(2,2) - m(3) * t116 - t55 * mrSges(3,1) - mrSges(3,3) * t108 - m(4) * t112 - t28 * mrSges(4,1) - m(5) * (pkin(3) * t28 + t112) - t8 * mrSges(5,1) - t7 * mrSges(5,2) - t143 * t6 + t145 * t54 - t159 * t5 + t92 * t27 - t150 * (pkin(4) * t132 - t27 * t75 + t28 * t69 + t112)) * g(2) + (t139 * mrSges(2,1) + t140 * mrSges(2,2) - m(3) * t101 + t53 * mrSges(3,1) - mrSges(3,3) * t109 - t148 * t24 + (-t145 + t97) * t52 + t143 * t154 + t159 * t1 - t92 * t23 + t151 * t95 + t150 * (pkin(4) * t134 - t23 * t75 + t24 * t69 + t95)) * g(1) (-t150 * (-t114 * t75 + t127 * t141 + t117) - t159 * (t113 * t80 - t127 * t71) - t151 * t117 - t160 * t114 + (-t122 * t162 - t143 * (t122 * t71 + t70 * t78) + t157 * t81 + t156 * t78) * t74) * g(3) + (m(5) * t46 + t158 * (pkin(9) * t53 - t46) + t153 * t53 + t152 * t52) * g(2) + (m(5) * t48 + t158 * (pkin(9) * t55 - t48) + t153 * t55 + t152 * t54) * g(1) (-t150 * (-t50 * t69 - t51 * t75) + t92 * t51 + t144 * t50) * g(3) + (-t150 * (-t23 * t69 - t24 * t75) + t92 * t24 + t144 * t23) * g(2) + (-t150 * (-t27 * t69 - t28 * t75) + t92 * t28 + t144 * t27) * g(1) (-t93 * mrSges(5,1) - (t126 * t76 - t51 * t79) * mrSges(5,2) - t159 * (t51 * t71 - t113) + t143 * t17) * g(3) + (-t155 * mrSges(5,1) - (-t24 * t79 - t134) * mrSges(5,2) - t159 * t154 + t143 * t1) * g(2) + (-t7 * mrSges(5,1) + t8 * mrSges(5,2) + t143 * t5 - t159 * t6) * g(1) + (-g(1) * t7 - g(2) * t155 - g(3) * t93) * t150 * pkin(4), t150 * (-g(1) * t27 - g(2) * t23 - g(3) * t50) (-g(1) * t5 - g(2) * t1 - g(3) * t17) * m(7)];
taug  = t2(:);
