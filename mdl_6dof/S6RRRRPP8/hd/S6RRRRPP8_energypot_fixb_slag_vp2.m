% Calculate potential energy for
% S6RRRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
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
% Datum: 2019-03-09 21:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRPP8_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP8_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP8_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP8_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP8_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP8_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:30:12
% EndTime: 2019-03-09 21:30:13
% DurationCPUTime: 0.53s
% Computational Cost: add. (335->95), mult. (748->110), div. (0->0), fcn. (903->10), ass. (0->53)
t144 = -mrSges(4,3) + mrSges(3,2);
t143 = mrSges(5,2) - mrSges(6,3) - mrSges(7,2);
t142 = mrSges(4,2) - mrSges(5,3) - mrSges(6,2) - m(7) * (pkin(10) - qJ(6)) + mrSges(7,3);
t141 = -m(7) * pkin(5) - mrSges(5,1) - mrSges(6,1) - mrSges(7,1);
t110 = sin(qJ(3));
t112 = sin(qJ(1));
t108 = sin(pkin(6));
t137 = cos(qJ(3));
t129 = t108 * t137;
t111 = sin(qJ(2));
t114 = cos(qJ(2));
t115 = cos(qJ(1));
t135 = cos(pkin(6));
t128 = t112 * t135;
t99 = -t111 * t128 + t115 * t114;
t87 = t110 * t99 - t112 * t129;
t140 = pkin(10) * t87;
t134 = t108 * t111;
t94 = t110 * t134 - t135 * t137;
t139 = pkin(10) * t94;
t127 = t115 * t135;
t97 = t111 * t127 + t112 * t114;
t85 = t97 * t110 + t115 * t129;
t138 = t85 * pkin(10);
t136 = t135 * pkin(8) + pkin(7);
t133 = t108 * t112;
t132 = t108 * t114;
t131 = t108 * t115;
t130 = t115 * pkin(1) + pkin(8) * t133;
t126 = t112 * pkin(1) - pkin(8) * t131;
t98 = t115 * t111 + t114 * t128;
t124 = t99 * pkin(2) + pkin(9) * t98 + t130;
t88 = t110 * t133 + t137 * t99;
t123 = t88 * pkin(3) + t124;
t122 = pkin(2) * t134 - pkin(9) * t132 + t136;
t95 = t110 * t135 + t111 * t129;
t121 = t95 * pkin(3) + t122;
t96 = t111 * t112 - t114 * t127;
t120 = t97 * pkin(2) + t96 * pkin(9) + t126;
t86 = -t110 * t131 + t137 * t97;
t119 = t86 * pkin(3) + t120;
t109 = sin(qJ(4));
t113 = cos(qJ(4));
t78 = t109 * t88 - t98 * t113;
t79 = t109 * t98 + t113 * t88;
t118 = t79 * pkin(4) + t78 * qJ(5) + t123;
t83 = t109 * t95 + t113 * t132;
t84 = -t109 * t132 + t113 * t95;
t117 = t84 * pkin(4) + qJ(5) * t83 + t121;
t76 = t109 * t86 - t96 * t113;
t77 = t109 * t96 + t113 * t86;
t116 = t77 * pkin(4) + t76 * qJ(5) + t119;
t1 = (-mrSges(1,3) - m(2) * pkin(7) - mrSges(2,3) - m(3) * t136 - t135 * mrSges(3,3) - (t111 * mrSges(3,1) + t114 * mrSges(3,2)) * t108 - m(4) * t122 - t95 * mrSges(4,1) + mrSges(4,3) * t132 - m(5) * (t121 + t139) - m(6) * (t117 + t139) - m(7) * t117 + t141 * t84 + t143 * t83 + t142 * t94) * g(3) + (-mrSges(1,2) - t112 * mrSges(2,1) - t115 * mrSges(2,2) - m(3) * t126 - t97 * mrSges(3,1) + mrSges(3,3) * t131 - m(4) * t120 - t86 * mrSges(4,1) - m(5) * (t119 + t138) - m(6) * (t116 + t138) - m(7) * t116 + t144 * t96 + t141 * t77 + t143 * t76 + t142 * t85) * g(2) + (-mrSges(1,1) - t115 * mrSges(2,1) + t112 * mrSges(2,2) - m(3) * t130 - t99 * mrSges(3,1) - mrSges(3,3) * t133 - m(4) * t124 - t88 * mrSges(4,1) - m(5) * (t123 + t140) - m(6) * (t118 + t140) - m(7) * t118 + t144 * t98 + t141 * t79 + t143 * t78 + t142 * t87) * g(1);
U  = t1;
