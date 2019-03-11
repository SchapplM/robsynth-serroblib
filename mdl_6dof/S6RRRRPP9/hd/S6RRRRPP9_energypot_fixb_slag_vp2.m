% Calculate potential energy for
% S6RRRRPP9
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
% Datum: 2019-03-09 21:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRPP9_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP9_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP9_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP9_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP9_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP9_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:41:41
% EndTime: 2019-03-09 21:41:42
% DurationCPUTime: 0.52s
% Computational Cost: add. (335->95), mult. (748->110), div. (0->0), fcn. (903->10), ass. (0->53)
t146 = -mrSges(4,3) + mrSges(3,2);
t145 = mrSges(5,2) - mrSges(6,3) - mrSges(7,2);
t144 = mrSges(4,2) - mrSges(5,3) - mrSges(6,1) - m(7) * (pkin(5) + pkin(10)) - mrSges(7,1);
t143 = -m(7) * qJ(6) - mrSges(5,1) + mrSges(6,2) - mrSges(7,3);
t113 = sin(qJ(2));
t116 = cos(qJ(2));
t117 = cos(qJ(1));
t114 = sin(qJ(1));
t137 = cos(pkin(6));
t130 = t114 * t137;
t100 = -t113 * t130 + t117 * t116;
t112 = sin(qJ(3));
t110 = sin(pkin(6));
t139 = cos(qJ(3));
t131 = t110 * t139;
t88 = t100 * t112 - t114 * t131;
t142 = pkin(10) * t88;
t136 = t110 * t113;
t95 = t112 * t136 - t137 * t139;
t141 = pkin(10) * t95;
t129 = t117 * t137;
t98 = t113 * t129 + t114 * t116;
t86 = t98 * t112 + t117 * t131;
t140 = t86 * pkin(10);
t138 = t137 * pkin(8) + pkin(7);
t135 = t110 * t114;
t134 = t110 * t116;
t133 = t110 * t117;
t132 = t117 * pkin(1) + pkin(8) * t135;
t127 = t114 * pkin(1) - pkin(8) * t133;
t99 = t117 * t113 + t116 * t130;
t126 = t100 * pkin(2) + pkin(9) * t99 + t132;
t89 = t100 * t139 + t112 * t135;
t125 = t89 * pkin(3) + t126;
t124 = pkin(2) * t136 - pkin(9) * t134 + t138;
t96 = t137 * t112 + t113 * t131;
t123 = t96 * pkin(3) + t124;
t97 = t113 * t114 - t116 * t129;
t122 = t98 * pkin(2) + t97 * pkin(9) + t127;
t87 = -t112 * t133 + t98 * t139;
t121 = t87 * pkin(3) + t122;
t111 = sin(qJ(4));
t115 = cos(qJ(4));
t79 = t111 * t89 - t99 * t115;
t80 = t111 * t99 + t115 * t89;
t120 = t80 * pkin(4) + t79 * qJ(5) + t125;
t84 = t111 * t96 + t115 * t134;
t85 = -t111 * t134 + t115 * t96;
t119 = t85 * pkin(4) + qJ(5) * t84 + t123;
t77 = t111 * t87 - t97 * t115;
t78 = t111 * t97 + t115 * t87;
t118 = t78 * pkin(4) + t77 * qJ(5) + t121;
t1 = (-mrSges(1,3) - m(2) * pkin(7) - mrSges(2,3) - m(3) * t138 - t137 * mrSges(3,3) - (t113 * mrSges(3,1) + t116 * mrSges(3,2)) * t110 - m(4) * t124 - t96 * mrSges(4,1) + mrSges(4,3) * t134 - m(5) * (t123 + t141) - m(6) * (t119 + t141) - m(7) * t119 + t143 * t85 + t145 * t84 + t144 * t95) * g(3) + (-mrSges(1,2) - t114 * mrSges(2,1) - t117 * mrSges(2,2) - m(3) * t127 - t98 * mrSges(3,1) + mrSges(3,3) * t133 - m(4) * t122 - t87 * mrSges(4,1) - m(5) * (t121 + t140) - m(6) * (t118 + t140) - m(7) * t118 + t146 * t97 + t143 * t78 + t145 * t77 + t144 * t86) * g(2) + (-mrSges(1,1) - t117 * mrSges(2,1) + t114 * mrSges(2,2) - m(3) * t132 - t100 * mrSges(3,1) - mrSges(3,3) * t135 - m(4) * t126 - t89 * mrSges(4,1) - m(5) * (t125 + t142) - m(6) * (t120 + t142) - m(7) * t120 + t146 * t99 + t143 * t80 + t145 * t79 + t144 * t88) * g(1);
U  = t1;
