% Calculate potential energy for
% S6RRRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
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
% Datum: 2019-03-09 22:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRPR8_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR8_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR8_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR8_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR8_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR8_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:38:44
% EndTime: 2019-03-09 22:38:45
% DurationCPUTime: 0.66s
% Computational Cost: add. (244->82), mult. (315->75), div. (0->0), fcn. (314->10), ass. (0->39)
t132 = -m(4) * pkin(8) + mrSges(3,2) - mrSges(6,2) - mrSges(4,3) - mrSges(5,3) + mrSges(7,3);
t91 = sin(qJ(3));
t95 = cos(qJ(3));
t131 = -m(4) * pkin(2) - mrSges(4,1) * t95 + mrSges(4,2) * t91 - mrSges(3,1);
t129 = -m(6) - m(7);
t92 = sin(qJ(2));
t96 = cos(qJ(2));
t98 = -pkin(9) - pkin(8);
t128 = t131 * t96 - mrSges(2,1) + (-m(7) * (-pkin(10) - t98) + t132) * t92;
t127 = -m(3) - m(4);
t93 = sin(qJ(1));
t114 = t93 * t96;
t89 = qJ(3) + qJ(4);
t84 = sin(t89);
t85 = cos(t89);
t97 = cos(qJ(1));
t71 = t84 * t114 + t85 * t97;
t72 = t85 * t114 - t84 * t97;
t125 = t72 * pkin(4) + t71 * qJ(5);
t124 = -t91 * mrSges(4,1) - t95 * mrSges(4,2) + mrSges(2,2) - mrSges(3,3);
t90 = sin(qJ(6));
t94 = cos(qJ(6));
t120 = -t90 * mrSges(7,1) - t94 * mrSges(7,2) + mrSges(5,2) - mrSges(6,3);
t119 = -m(7) * pkin(5) - t94 * mrSges(7,1) + t90 * mrSges(7,2) - mrSges(5,1) - mrSges(6,1);
t118 = pkin(3) * t91;
t115 = t92 * t98;
t113 = t97 * t96;
t112 = t97 * pkin(1) + t93 * pkin(7);
t109 = t97 * t115;
t87 = t93 * pkin(1);
t108 = -t97 * pkin(7) + t87;
t82 = pkin(3) * t95 + pkin(2);
t107 = t82 * t113 + t93 * t118 + t112;
t73 = t84 * t113 - t93 * t85;
t74 = t85 * t113 + t84 * t93;
t101 = t74 * pkin(4) + t73 * qJ(5) + t107;
t75 = t82 * t114;
t100 = -t93 * t115 + t75 + t87 + (-pkin(7) - t118) * t97;
t1 = (-mrSges(1,3) - mrSges(2,3) + (-m(2) + t127) * pkin(6) + (-m(7) * pkin(10) - t132) * t96 + (-m(5) + t129) * (t92 * t82 + t96 * t98 + pkin(6)) + (t129 * (pkin(4) * t85 + qJ(5) * t84) + t119 * t85 + t120 * t84 + t131) * t92) * g(3) + (-mrSges(1,2) - m(3) * t108 - m(4) * t87 - m(5) * t100 - m(6) * (t100 + t125) - m(7) * (t108 + t75 + t125) + t119 * t72 + t120 * t71 + (m(4) * pkin(7) + m(7) * t118 - t124) * t97 + t128 * t93) * g(2) + (-mrSges(1,1) - m(5) * (t107 - t109) - m(6) * (t101 - t109) - m(7) * t101 + t119 * t74 + t120 * t73 + t127 * t112 + t124 * t93 + t128 * t97) * g(1);
U  = t1;
