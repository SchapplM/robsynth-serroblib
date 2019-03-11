% Calculate potential energy for
% S6PRRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRRPP2_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP2_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPP2_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP2_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP2_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPP2_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:48:39
% EndTime: 2019-03-08 22:48:40
% DurationCPUTime: 0.54s
% Computational Cost: add. (335->95), mult. (748->111), div. (0->0), fcn. (903->10), ass. (0->54)
t146 = -mrSges(4,3) + mrSges(3,2);
t145 = mrSges(5,2) - mrSges(6,3) - mrSges(7,2);
t144 = mrSges(4,2) - mrSges(5,3) - mrSges(6,2) - m(7) * (pkin(9) - qJ(6)) + mrSges(7,3);
t143 = -m(7) * pkin(5) - mrSges(5,1) - mrSges(6,1) - mrSges(7,1);
t111 = cos(pkin(10));
t113 = sin(qJ(3));
t110 = sin(pkin(6));
t139 = cos(qJ(3));
t130 = t110 * t139;
t109 = sin(pkin(10));
t116 = cos(qJ(2));
t114 = sin(qJ(2));
t138 = cos(pkin(6));
t129 = t114 * t138;
t96 = t109 * t116 + t111 * t129;
t84 = t111 * t130 + t96 * t113;
t142 = pkin(9) * t84;
t98 = -t109 * t129 + t111 * t116;
t86 = -t109 * t130 + t113 * t98;
t141 = pkin(9) * t86;
t134 = t110 * t114;
t99 = t113 * t134 - t138 * t139;
t140 = t99 * pkin(9);
t137 = t109 * t110;
t136 = t110 * t111;
t135 = t110 * t113;
t133 = t110 * t116;
t132 = t111 * pkin(1) + pkin(7) * t137;
t131 = t138 * pkin(7) + qJ(1);
t128 = t116 * t138;
t127 = t109 * pkin(1) - pkin(7) * t136;
t97 = t109 * t128 + t111 * t114;
t125 = t98 * pkin(2) + pkin(8) * t97 + t132;
t87 = t109 * t135 + t98 * t139;
t124 = t87 * pkin(3) + t125;
t123 = pkin(2) * t134 - pkin(8) * t133 + t131;
t100 = t138 * t113 + t114 * t130;
t122 = t100 * pkin(3) + t123;
t95 = t109 * t114 - t111 * t128;
t121 = t96 * pkin(2) + pkin(8) * t95 + t127;
t85 = -t111 * t135 + t96 * t139;
t120 = t85 * pkin(3) + t121;
t112 = sin(qJ(4));
t115 = cos(qJ(4));
t79 = t112 * t87 - t97 * t115;
t80 = t112 * t97 + t115 * t87;
t119 = t80 * pkin(4) + t79 * qJ(5) + t124;
t88 = t100 * t112 + t115 * t133;
t89 = t100 * t115 - t112 * t133;
t118 = t89 * pkin(4) + t88 * qJ(5) + t122;
t77 = t112 * t85 - t95 * t115;
t78 = t112 * t95 + t115 * t85;
t117 = t78 * pkin(4) + t77 * qJ(5) + t120;
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - m(3) * t131 - t138 * mrSges(3,3) - (t114 * mrSges(3,1) + t116 * mrSges(3,2)) * t110 - m(4) * t123 - t100 * mrSges(4,1) + mrSges(4,3) * t133 - m(5) * (t122 + t140) - m(6) * (t118 + t140) - m(7) * t118 + t143 * t89 + t145 * t88 + t144 * t99) * g(3) + (-mrSges(1,2) - t109 * mrSges(2,1) - t111 * mrSges(2,2) - m(3) * t127 - t96 * mrSges(3,1) + mrSges(3,3) * t136 - m(4) * t121 - t85 * mrSges(4,1) - m(5) * (t120 + t142) - m(6) * (t117 + t142) - m(7) * t117 + t146 * t95 + t143 * t78 + t145 * t77 + t144 * t84) * g(2) + (-mrSges(1,1) - t111 * mrSges(2,1) + t109 * mrSges(2,2) - m(3) * t132 - t98 * mrSges(3,1) - mrSges(3,3) * t137 - m(4) * t125 - t87 * mrSges(4,1) - m(5) * (t124 + t141) - m(6) * (t119 + t141) - m(7) * t119 + t146 * t97 + t143 * t80 + t145 * t79 + t144 * t86) * g(1);
U  = t1;
