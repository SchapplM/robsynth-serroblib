% Calculate Gravitation load on the joints for
% S6RRRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 01:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRP6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP6_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:26:35
% EndTime: 2019-03-10 01:26:39
% DurationCPUTime: 1.18s
% Computational Cost: add. (722->154), mult. (770->175), div. (0->0), fcn. (743->10), ass. (0->82)
t142 = mrSges(6,1) + mrSges(7,1);
t141 = mrSges(6,2) - mrSges(7,3);
t144 = m(6) + m(7);
t143 = mrSges(4,3) + mrSges(5,3);
t140 = -mrSges(6,3) - mrSges(7,2);
t62 = -pkin(9) - pkin(8);
t139 = -m(4) * pkin(8) + m(5) * t62 - t143;
t55 = qJ(3) + qJ(4);
t49 = qJ(5) + t55;
t44 = sin(t49);
t45 = cos(t49);
t59 = cos(qJ(3));
t51 = t59 * pkin(3);
t46 = t51 + pkin(2);
t47 = sin(t55);
t48 = cos(t55);
t56 = sin(qJ(3));
t138 = m(4) * pkin(2) + m(5) * t46 + t59 * mrSges(4,1) + t48 * mrSges(5,1) - t56 * mrSges(4,2) - t47 * mrSges(5,2) - t141 * t44 + t142 * t45;
t58 = sin(qJ(1));
t61 = cos(qJ(1));
t137 = g(1) * t61 + g(2) * t58;
t136 = m(7) * qJ(6) + mrSges(7,3);
t57 = sin(qJ(2));
t135 = t140 * t57;
t60 = cos(qJ(2));
t96 = t60 * t61;
t21 = -t47 * t96 + t48 * t58;
t120 = m(5) * pkin(3);
t132 = mrSges(2,2) - mrSges(3,3);
t109 = mrSges(6,2) * t45;
t131 = mrSges(5,1) * t47 + mrSges(6,1) * t44 + mrSges(5,2) * t48 + t109;
t130 = mrSges(4,1) + t120;
t95 = t61 * t44;
t13 = -t58 * t45 + t60 * t95;
t14 = t44 * t58 + t45 * t96;
t129 = t142 * t13 + t141 * t14;
t97 = t58 * t60;
t11 = t44 * t97 + t45 * t61;
t12 = t45 * t97 - t95;
t128 = t142 * t11 + t141 * t12;
t127 = -m(3) - m(4) - m(5);
t22 = t47 * t58 + t48 * t96;
t126 = -t21 * mrSges(5,1) + t22 * mrSges(5,2) + t129;
t19 = t47 * t97 + t48 * t61;
t20 = t47 * t61 - t48 * t97;
t125 = t19 * mrSges(5,1) - t20 * mrSges(5,2) + t128;
t79 = t60 * mrSges(3,1) - t57 * mrSges(3,2);
t124 = t143 * t57 + mrSges(2,1) + t79;
t123 = m(7) * pkin(5) + t142;
t122 = -mrSges(6,2) + t136;
t115 = pkin(3) * t56;
t114 = pkin(4) * t47;
t113 = pkin(5) * t44;
t110 = g(3) * t57;
t108 = t44 * mrSges(7,1);
t54 = -pkin(10) + t62;
t105 = t54 * t57;
t104 = t56 * t58;
t103 = t56 * t61;
t98 = t57 * t61;
t36 = pkin(4) * t48 + t51;
t34 = pkin(2) + t36;
t30 = t60 * t34;
t35 = t114 + t115;
t33 = t61 * t35;
t92 = -t60 * t33 + t58 * t36;
t91 = t61 * pkin(1) + t58 * pkin(7);
t86 = t136 * t45 * t57;
t84 = -t11 * pkin(5) + qJ(6) * t12;
t83 = -t35 * t97 - t36 * t61;
t81 = -t13 * pkin(5) + qJ(6) * t14;
t80 = pkin(2) * t60 + pkin(8) * t57;
t73 = pkin(5) * t45 + qJ(6) * t44;
t72 = t60 * t46 - t57 * t62;
t71 = t21 * pkin(4);
t28 = -t56 * t96 + t58 * t59;
t26 = t56 * t97 + t59 * t61;
t68 = t19 * pkin(4);
t52 = t61 * pkin(7);
t29 = t59 * t96 + t104;
t27 = -t59 * t97 + t103;
t1 = [(-t104 * t120 - t29 * mrSges(4,1) - t22 * mrSges(5,1) - t28 * mrSges(4,2) - t21 * mrSges(5,2) + t140 * t98 + t127 * t91 - t144 * (t34 * t96 + t58 * t35 - t54 * t98 + t91) + t132 * t58 - t123 * t14 - t122 * t13 + (-m(4) * t80 - m(5) * t72 - t124) * t61) * g(2) + (-t103 * t120 - t27 * mrSges(4,1) - t20 * mrSges(5,1) - t26 * mrSges(4,2) - t19 * mrSges(5,2) - t144 * (t58 * t105 + t33 + t52) + t132 * t61 + t127 * t52 + t123 * t12 + t122 * t11 + (m(3) * pkin(1) - m(4) * (-pkin(1) - t80) - m(5) * (-pkin(1) - t72) - t144 * (-pkin(1) - t30) + t124 - t135) * t58) * g(1) (-t79 - t144 * (t30 - t105) + t135) * g(3) + ((-m(7) * t73 - t138) * g(3) + t137 * (t144 * t54 + mrSges(3,2) + t139 + t140)) * t60 + (t139 * g(3) + t137 * (mrSges(3,1) + m(6) * t34 - m(7) * (-t34 - t73) + t138)) * t57, -g(3) * ((m(7) * (-t35 - t113) - t108) * t57 + t86) + (m(5) * t115 + m(6) * t35 + mrSges(4,1) * t56 + mrSges(4,2) * t59 + t131) * t110 + (-t27 * mrSges(4,2) - m(6) * t83 - m(7) * (t83 + t84) + t130 * t26 + t125) * g(2) + (t29 * mrSges(4,2) - m(6) * t92 - m(7) * (t81 + t92) - t130 * t28 + t126) * g(1), -g(3) * ((m(7) * (-t113 - t114) - t108) * t57 + t86) + (m(6) * t114 + t131) * t110 + (m(6) * t68 - m(7) * (-t68 + t84) + t125) * g(2) + (-m(6) * t71 - m(7) * (t71 + t81) + t126) * g(1) ((t123 * t44 + t109) * t57 - t86) * g(3) + (-m(7) * t84 + t128) * g(2) + (-m(7) * t81 + t129) * g(1) (-g(1) * t13 - g(2) * t11 - t44 * t110) * m(7)];
taug  = t1(:);
