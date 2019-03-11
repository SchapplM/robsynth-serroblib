% Calculate Gravitation load on the joints for
% S6RPRPRR9
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
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPRR9_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR9_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_gravloadJ_floatb_twist_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR9_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR9_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:01:17
% EndTime: 2019-03-09 04:01:21
% DurationCPUTime: 1.37s
% Computational Cost: add. (1198->154), mult. (3200->228), div. (0->0), fcn. (4079->16), ass. (0->86)
t74 = sin(qJ(6));
t78 = cos(qJ(6));
t134 = m(6) + m(7);
t98 = pkin(10) * t134 - mrSges(5,2) + mrSges(6,3);
t135 = mrSges(7,1) * t74 + mrSges(7,2) * t78 + t98;
t138 = m(7) * pkin(5) + mrSges(7,1) * t78 - mrSges(7,2) * t74 + mrSges(6,1);
t105 = -m(7) * pkin(11) + mrSges(6,2) - mrSges(7,3);
t143 = m(5) + t134;
t72 = cos(pkin(7));
t80 = cos(qJ(3));
t120 = t72 * t80;
t67 = sin(pkin(12));
t71 = cos(pkin(12));
t76 = sin(qJ(3));
t142 = t120 * t71 - t67 * t76;
t77 = sin(qJ(1));
t118 = t77 * t71;
t73 = cos(pkin(6));
t81 = cos(qJ(1));
t119 = t73 * t81;
t51 = t119 * t67 + t118;
t129 = t51 * t76;
t125 = t67 * t77;
t50 = -t119 * t71 + t125;
t141 = t120 * t50 + t129;
t69 = sin(pkin(6));
t122 = t69 * t81;
t68 = sin(pkin(7));
t130 = t50 * t68;
t33 = t122 * t72 - t130;
t75 = sin(qJ(5));
t79 = cos(qJ(5));
t66 = sin(pkin(13));
t70 = cos(pkin(13));
t100 = t66 * t80 + t70 * t76;
t44 = t100 * t68;
t46 = t100 * t72;
t54 = t66 * t76 - t70 * t80;
t85 = t122 * t44 + t46 * t50 + t51 * t54;
t140 = -t33 * t79 + t75 * t85;
t139 = t33 * t75 + t79 * t85;
t136 = t105 * t75 - t138 * t79 - mrSges(5,1);
t26 = (t46 * t71 - t54 * t67) * t69 + t44 * t73;
t132 = pkin(3) * t80;
t52 = -t118 * t73 - t67 * t81;
t128 = t52 * t68;
t53 = -t125 * t73 + t71 * t81;
t127 = t53 * t76;
t124 = t68 * t80;
t123 = t69 * t77;
t121 = t72 * t76;
t117 = -mrSges(4,3) - mrSges(5,3);
t116 = pkin(9) + qJ(4);
t114 = qJ(2) * t69;
t115 = pkin(1) * t81 + t114 * t77;
t112 = t68 * t123;
t111 = t68 * t122;
t108 = -t77 * pkin(1) + t114 * t81;
t106 = -m(4) * pkin(9) * t72 - mrSges(3,3);
t47 = pkin(3) * t68 * t76 + t116 * t72;
t48 = pkin(3) * t121 - t116 * t68;
t64 = pkin(2) + t132;
t104 = t123 * t47 + t48 * t52 + t53 * t64 + t115;
t99 = t112 * t132 + (t120 * t52 - t127) * pkin(3);
t22 = t123 * t44 + t46 * t52 - t53 * t54;
t97 = pkin(4) * t22 + t104;
t96 = t52 * t72 + t112;
t95 = (t124 * t73 + t142 * t69) * pkin(3);
t93 = t122 * t47 + t48 * t50 - t51 * t64 + t108;
t92 = t111 * t76 + t121 * t50 - t51 * t80;
t43 = t54 * t68;
t45 = t54 * t72;
t86 = t100 * t51 - t122 * t43 - t45 * t50;
t84 = (-t111 * t80 - t141) * pkin(3);
t49 = -t68 * t69 * t71 + t72 * t73;
t35 = t123 * t72 - t128;
t28 = t53 * t80 + t76 * t96;
t27 = t80 * t96 - t127;
t25 = -t43 * t73 + (-t100 * t67 - t45 * t71) * t69;
t21 = -t100 * t53 - t123 * t43 - t45 * t52;
t14 = t26 * t79 + t49 * t75;
t8 = t22 * t79 + t35 * t75;
t7 = t22 * t75 - t35 * t79;
t2 = -t21 * t74 + t78 * t8;
t1 = -t21 * t78 - t74 * t8;
t3 = [(-t81 * mrSges(2,1) - m(3) * t115 - t53 * mrSges(3,1) - t52 * mrSges(3,2) - m(4) * (pkin(2) * t53 - pkin(9) * t128 + t115) - t28 * mrSges(4,1) - t27 * mrSges(4,2) - m(5) * t104 - t22 * mrSges(5,1) - m(6) * t97 - t8 * mrSges(6,1) - m(7) * (pkin(5) * t8 + t97) - t2 * mrSges(7,1) - t1 * mrSges(7,2) + (t106 * t69 + mrSges(2,2)) * t77 + t105 * t7 + t117 * t35 + t98 * t21) * g(2) + (t77 * mrSges(2,1) - m(3) * t108 + t51 * mrSges(3,1) - t50 * mrSges(3,2) - m(4) * (-t51 * pkin(2) - pkin(9) * t130 + t108) - t92 * mrSges(4,1) - t141 * mrSges(4,2) - m(5) * t93 - t85 * mrSges(5,1) + (mrSges(2,2) + (-mrSges(4,2) * t124 + t106) * t69) * t81 + t105 * t140 + t117 * t33 - t138 * t139 + t135 * t86 + t134 * (-pkin(4) * t85 - t93)) * g(1) (-t73 * g(3) + (-g(1) * t77 + g(2) * t81) * t69) * (m(3) + m(4) + t143) (-(mrSges(4,1) * t80 - mrSges(4,2) * t76) * t73 * t68 - (t142 * mrSges(4,1) + (-t121 * t71 - t67 * t80) * mrSges(4,2)) * t69 - m(5) * t95 - t134 * (pkin(4) * t25 + t95) + t136 * t25 - t135 * t26) * g(3) + (-(-t129 + (-t50 * t72 - t111) * t80) * mrSges(4,1) - t92 * mrSges(4,2) - m(5) * t84 - t134 * (-pkin(4) * t86 + t84) - t136 * t86 + t135 * t85) * g(2) + (-m(5) * t99 - mrSges(4,1) * t27 + mrSges(4,2) * t28 - t134 * (pkin(4) * t21 + t99) + t136 * t21 - t135 * t22) * g(1), t143 * (-g(1) * t35 + g(2) * t33 - g(3) * t49) (t105 * t14 - t138 * (-t26 * t75 + t49 * t79)) * g(3) + (-t105 * t139 - t138 * t140) * g(2) + (t105 * t8 + t138 * t7) * g(1), -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((t139 * t74 + t78 * t86) * mrSges(7,1) + (t139 * t78 - t74 * t86) * mrSges(7,2)) - g(3) * ((-t14 * t74 - t25 * t78) * mrSges(7,1) + (-t14 * t78 + t25 * t74) * mrSges(7,2))];
taug  = t3(:);
