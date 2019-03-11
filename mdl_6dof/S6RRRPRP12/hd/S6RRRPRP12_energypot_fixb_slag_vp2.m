% Calculate potential energy for
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
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPRP12_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP12_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP12_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP12_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP12_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP12_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:51:42
% EndTime: 2019-03-09 17:51:42
% DurationCPUTime: 0.53s
% Computational Cost: add. (301->90), mult. (659->107), div. (0->0), fcn. (780->10), ass. (0->45)
t144 = -m(6) - m(7);
t143 = mrSges(4,2) - mrSges(5,3);
t142 = mrSges(4,3) + mrSges(5,1);
t141 = mrSges(3,2) - t142;
t140 = -mrSges(4,1) + mrSges(5,2) - mrSges(6,3) - mrSges(7,2);
t139 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t138 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t137 = cos(qJ(3));
t135 = cos(pkin(6));
t136 = t135 * pkin(8) + pkin(7);
t108 = sin(pkin(6));
t111 = sin(qJ(2));
t134 = t108 * t111;
t114 = cos(qJ(2));
t133 = t108 * t114;
t112 = sin(qJ(1));
t132 = t112 * t108;
t115 = cos(qJ(1));
t131 = t115 * t108;
t130 = t115 * pkin(1) + pkin(8) * t132;
t129 = pkin(9) * t133;
t128 = pkin(2) * t134 + t136;
t127 = t108 * t137;
t126 = t112 * t135;
t125 = t115 * t135;
t124 = t112 * pkin(1) - pkin(8) * t131;
t97 = t115 * t111 + t114 * t126;
t98 = -t111 * t126 + t115 * t114;
t123 = t98 * pkin(2) + t97 * pkin(9) + t130;
t110 = sin(qJ(3));
t93 = t110 * t134 - t135 * t137;
t94 = t135 * t110 + t111 * t127;
t122 = t94 * pkin(3) + t93 * qJ(4) + t128;
t95 = t112 * t111 - t114 * t125;
t96 = t111 * t125 + t112 * t114;
t121 = t96 * pkin(2) + t95 * pkin(9) + t124;
t86 = t98 * t110 - t112 * t127;
t87 = t110 * t132 + t98 * t137;
t120 = t87 * pkin(3) + t86 * qJ(4) + t123;
t84 = t96 * t110 + t115 * t127;
t85 = -t110 * t131 + t96 * t137;
t119 = t85 * pkin(3) + t84 * qJ(4) + t121;
t113 = cos(qJ(5));
t109 = sin(qJ(5));
t1 = (-mrSges(1,3) - m(2) * pkin(7) - mrSges(2,3) - m(3) * t136 - t135 * mrSges(3,3) - (t111 * mrSges(3,1) + t114 * mrSges(3,2)) * t108 - m(4) * (t128 - t129) - m(5) * (t122 - t129) + t143 * t93 + t144 * (t94 * pkin(10) + (-pkin(4) - pkin(9)) * t133 + t122) + t139 * (t93 * t109 - t113 * t133) + t138 * (t109 * t133 + t93 * t113) + t142 * t133 + t140 * t94) * g(3) + (-m(3) * t124 - m(4) * t121 - m(5) * t119 - t112 * mrSges(2,1) - t96 * mrSges(3,1) - t115 * mrSges(2,2) + mrSges(3,3) * t131 - mrSges(1,2) + t143 * t84 + t144 * (t95 * pkin(4) + t85 * pkin(10) + t119) + t139 * (t84 * t109 + t95 * t113) - t138 * (t95 * t109 - t84 * t113) + t141 * t95 + t140 * t85) * g(2) + (-m(3) * t130 - m(4) * t123 - m(5) * t120 - t115 * mrSges(2,1) - t98 * mrSges(3,1) + t112 * mrSges(2,2) - mrSges(3,3) * t132 - mrSges(1,1) + t143 * t86 + t144 * (t97 * pkin(4) + t87 * pkin(10) + t120) + t139 * (t86 * t109 + t97 * t113) - t138 * (t97 * t109 - t86 * t113) + t141 * t97 + t140 * t87) * g(1);
U  = t1;
