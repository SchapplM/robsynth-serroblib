% Calculate potential energy for
% S6RRRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
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
% Datum: 2019-03-09 17:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPRP8_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP8_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP8_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP8_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP8_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP8_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:14:10
% EndTime: 2019-03-09 17:14:11
% DurationCPUTime: 0.63s
% Computational Cost: add. (194->87), mult. (371->92), div. (0->0), fcn. (388->8), ass. (0->43)
t132 = -m(5) - m(6);
t131 = -mrSges(4,1) - mrSges(5,1);
t130 = -mrSges(6,1) - mrSges(7,1);
t129 = mrSges(2,2) - mrSges(3,3);
t128 = mrSges(4,2) - mrSges(5,3);
t127 = -mrSges(6,2) - mrSges(7,2);
t95 = sin(qJ(2));
t98 = cos(qJ(3));
t116 = t95 * t98;
t94 = sin(qJ(3));
t118 = t95 * t94;
t126 = pkin(3) * t116 + qJ(4) * t118;
t99 = cos(qJ(2));
t125 = -t99 * mrSges(3,1) + t95 * mrSges(3,2) - mrSges(2,1);
t124 = mrSges(5,2) + mrSges(4,3) - mrSges(6,3) - mrSges(7,3);
t93 = sin(qJ(5));
t120 = pkin(5) * t93;
t123 = -m(7) * (qJ(4) + t120) + t128;
t97 = cos(qJ(5));
t87 = pkin(5) * t97 + pkin(4);
t122 = -m(6) * pkin(4) - m(7) * t87 + t131;
t92 = -qJ(6) - pkin(9);
t121 = m(6) * pkin(9) - m(7) * t92 - t124;
t119 = t95 * pkin(2) + pkin(6);
t96 = sin(qJ(1));
t117 = t95 * t96;
t115 = t96 * t99;
t100 = cos(qJ(1));
t114 = t100 * pkin(1) + t96 * pkin(7);
t113 = t100 * t95;
t112 = t100 * t99;
t110 = t96 * pkin(1) - pkin(7) * t100;
t109 = t119 + t126;
t108 = pkin(2) * t112 + pkin(8) * t113 + t114;
t107 = -pkin(8) * t99 + t119;
t79 = t112 * t98 + t96 * t94;
t106 = t79 * pkin(3) + t108;
t104 = pkin(2) * t115 + pkin(8) * t117 + t110;
t77 = -t100 * t94 + t115 * t98;
t103 = t77 * pkin(3) + t104;
t78 = t112 * t94 - t96 * t98;
t76 = t100 * t98 + t115 * t94;
t1 = (-m(3) * t110 - m(4) * t104 - m(7) * t103 - mrSges(1,2) + t125 * t96 + t132 * (t76 * qJ(4) + t103) + t122 * t77 + t123 * t76 + t130 * (t76 * t93 + t77 * t97) + t127 * (t76 * t97 - t77 * t93) - t129 * t100 + t121 * t117) * g(2) + (-m(3) * t114 - m(4) * t108 - m(7) * t106 - mrSges(1,1) + t129 * t96 + t132 * (t78 * qJ(4) + t106) + t122 * t79 + t123 * t78 + t130 * (t78 * t93 + t79 * t97) + t127 * (t78 * t97 - t79 * t93) + t125 * t100 + t121 * t113) * g(1) + (-mrSges(1,3) - mrSges(2,3) - m(4) * t107 - m(5) * (t107 + t126) - m(6) * (pkin(4) * t116 + t109) - m(7) * (t87 * t116 + t118 * t120 + t109) + (-m(2) - m(3)) * pkin(6) + (-mrSges(3,2) - m(6) * (-pkin(8) + pkin(9)) - m(7) * (-pkin(8) - t92) + t124) * t99 + (t130 * (t93 * t94 + t97 * t98) + t127 * (-t93 * t98 + t94 * t97) + t128 * t94 + t131 * t98 - mrSges(3,1)) * t95) * g(3);
U  = t1;
