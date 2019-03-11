% Calculate Gravitation load on the joints for
% S6RRRPRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
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
% Datum: 2019-03-09 18:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRP12_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP12_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP12_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP12_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP12_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP12_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:51:43
% EndTime: 2019-03-09 17:51:46
% DurationCPUTime: 1.24s
% Computational Cost: add. (695->134), mult. (1708->188), div. (0->0), fcn. (2032->10), ass. (0->70)
t134 = mrSges(4,2) - mrSges(5,3);
t135 = m(6) + m(7);
t140 = m(5) + t135;
t80 = -qJ(4) * t140 + t134;
t139 = -mrSges(4,1) + mrSges(5,2);
t138 = mrSges(6,3) + mrSges(7,2);
t92 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t90 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t126 = cos(qJ(1));
t66 = sin(pkin(6));
t102 = t66 * t126;
t69 = sin(qJ(2));
t70 = sin(qJ(1));
t73 = cos(qJ(2));
t109 = cos(pkin(6));
t88 = t109 * t126;
t48 = t69 * t88 + t70 * t73;
t68 = sin(qJ(3));
t72 = cos(qJ(3));
t20 = t102 * t72 + t48 * t68;
t47 = t69 * t70 - t73 * t88;
t67 = sin(qJ(5));
t71 = cos(qJ(5));
t137 = t20 * t67 + t47 * t71;
t136 = t20 * t71 - t47 * t67;
t132 = t134 * t68 + t139 * t72 - mrSges(3,1);
t105 = -mrSges(5,1) + mrSges(3,2) - mrSges(4,3);
t131 = pkin(9) * (m(4) + t140) - t105;
t130 = t138 - t139;
t129 = t138 * t72 - t132;
t128 = -t92 * t67 + t90 * t71 + t80;
t123 = t47 * t72;
t95 = t70 * t109;
t49 = t126 * t69 + t73 * t95;
t122 = t49 * t72;
t121 = t66 * t69;
t120 = t66 * t70;
t119 = t66 * t72;
t118 = t66 * t73;
t117 = t67 * t68;
t116 = t67 * t73;
t115 = t68 * t71;
t112 = pkin(2) * t118 + pkin(9) * t121;
t111 = t126 * pkin(1) + pkin(8) * t120;
t110 = qJ(4) * t68;
t107 = t71 * t118;
t106 = t72 * t118;
t50 = t126 * t73 - t69 * t95;
t104 = t50 * pkin(2) + t111;
t101 = -pkin(1) * t70 + pkin(8) * t102;
t100 = -t47 * pkin(2) + pkin(9) * t48;
t99 = -t49 * pkin(2) + pkin(9) * t50;
t21 = -t68 * t102 + t48 * t72;
t25 = t120 * t68 + t50 * t72;
t94 = t25 * pkin(3) + t104;
t93 = pkin(3) * t106 + t110 * t118 + t112;
t91 = -t48 * pkin(2) + t101;
t85 = -pkin(3) * t21 + t91;
t84 = -pkin(3) * t123 - t47 * t110 + t100;
t83 = -pkin(3) * t122 - t49 * t110 + t99;
t76 = -t135 * pkin(10) - t130;
t45 = -t109 * t72 + t121 * t68;
t36 = t45 * pkin(3);
t24 = -t119 * t70 + t50 * t68;
t18 = t116 * t66 + t45 * t71;
t16 = t24 * pkin(3);
t14 = t20 * pkin(3);
t6 = t24 * t67 + t49 * t71;
t5 = -t24 * t71 + t49 * t67;
t1 = [(-t126 * mrSges(2,1) - m(3) * t111 - t50 * mrSges(3,1) - m(4) * t104 - m(5) * t94 + (-mrSges(3,3) * t66 + mrSges(2,2)) * t70 - t92 * t6 - t90 * t5 + t80 * t24 - t131 * t49 + t76 * t25 - t135 * (t49 * pkin(4) + t94)) * g(2) + (t70 * mrSges(2,1) + t126 * mrSges(2,2) - m(3) * t101 + t48 * mrSges(3,1) - mrSges(3,3) * t102 - m(4) * t91 - m(5) * t85 + t92 * t137 - t90 * t136 - t80 * t20 + t131 * t47 - t76 * t21 + t135 * (t47 * pkin(4) - t85)) * g(1) (-m(4) * t112 - m(5) * t93 - t135 * (pkin(4) * t121 + pkin(10) * t106 + t93) - t90 * (-t107 * t68 + t121 * t67) - t138 * t106 + (-t92 * (t116 * t68 + t69 * t71) + t105 * t69 + t132 * t73) * t66) * g(3) + (-m(4) * t100 - m(5) * t84 - t135 * (t48 * pkin(4) - pkin(10) * t123 + t84) - t92 * (-t117 * t47 + t48 * t71) - t90 * (t115 * t47 + t48 * t67) + t105 * t48 + t129 * t47) * g(2) + (-m(4) * t99 - m(5) * t83 - t135 * (t50 * pkin(4) - pkin(10) * t122 + t83) - t90 * (t115 * t49 + t50 * t67) - t92 * (-t117 * t49 + t50 * t71) + t105 * t50 + t129 * t49) * g(1) (m(5) * t36 - t135 * (-pkin(10) * t45 - t36) + t128 * (t109 * t68 + t119 * t69) + t130 * t45) * g(3) + (m(5) * t14 - t135 * (-pkin(10) * t20 - t14) + t128 * t21 + t130 * t20) * g(2) + (m(5) * t16 - t135 * (-pkin(10) * t24 - t16) + t128 * t25 + t130 * t24) * g(1), t140 * (-g(1) * t24 - g(2) * t20 - g(3) * t45) (t90 * (-t45 * t67 + t107) - t92 * t18) * g(3) + (-t136 * t92 - t90 * t137) * g(2) + (t5 * t92 - t6 * t90) * g(1) (-g(1) * t5 + g(2) * t136 + g(3) * t18) * m(7)];
taug  = t1(:);
