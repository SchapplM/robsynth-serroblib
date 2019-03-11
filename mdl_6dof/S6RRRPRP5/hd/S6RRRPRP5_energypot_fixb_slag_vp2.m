% Calculate potential energy for
% S6RRRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-03-09 16:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPRP5_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP5_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP5_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP5_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP5_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP5_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:47:38
% EndTime: 2019-03-09 16:47:39
% DurationCPUTime: 0.53s
% Computational Cost: add. (247->72), mult. (264->68), div. (0->0), fcn. (247->10), ass. (0->34)
t90 = cos(qJ(3));
t76 = t90 * pkin(3) + pkin(2);
t85 = qJ(3) + pkin(10);
t77 = sin(t85);
t78 = cos(t85);
t87 = sin(qJ(3));
t119 = -m(4) * pkin(2) - m(5) * t76 - t90 * mrSges(4,1) - t78 * mrSges(5,1) + t87 * mrSges(4,2) + t77 * mrSges(5,2) - mrSges(3,1);
t86 = -qJ(4) - pkin(8);
t118 = -m(4) * pkin(8) + m(5) * t86 + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t117 = -m(6) - m(7);
t116 = -mrSges(6,3) - mrSges(7,2);
t115 = m(3) + m(4) + m(5);
t88 = sin(qJ(2));
t91 = cos(qJ(2));
t112 = t118 * t88 + t119 * t91 - mrSges(2,1);
t111 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t110 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t108 = pkin(3) * t87;
t109 = m(5) * t108 + t87 * mrSges(4,1) + t77 * mrSges(5,1) + t90 * mrSges(4,2) + t78 * mrSges(5,2) - mrSges(2,2) + mrSges(3,3);
t89 = sin(qJ(1));
t107 = t88 * t89;
t92 = cos(qJ(1));
t106 = t88 * t92;
t105 = t89 * t91;
t104 = t91 * t92;
t103 = t92 * pkin(1) + t89 * pkin(7);
t84 = -pkin(9) + t86;
t81 = t89 * pkin(1);
t79 = qJ(5) + t85;
t75 = cos(t79);
t74 = sin(t79);
t71 = pkin(4) * t77 + t108;
t70 = pkin(4) * t78 + t76;
t1 = (-mrSges(1,3) - mrSges(2,3) + t117 * (t88 * t70 + t91 * t84 + pkin(6)) + (-m(2) - t115) * pkin(6) + (-t116 - t118) * t91 + (t110 * t74 + t111 * t75 + t119) * t88) * g(3) + (-mrSges(1,2) + t117 * (-t84 * t107 + t70 * t105 + t81 + (-pkin(7) - t71) * t92) + t111 * (t75 * t105 - t74 * t92) + t110 * (t74 * t105 + t75 * t92) + t116 * t107 - t115 * t81 + (t115 * pkin(7) + t109) * t92 + t112 * t89) * g(2) + (-mrSges(1,1) + t117 * (t70 * t104 - t84 * t106 + t89 * t71 + t103) + t111 * (t75 * t104 + t89 * t74) + t110 * (t74 * t104 - t89 * t75) + t116 * t106 - t115 * t103 + t112 * t92 - t109 * t89) * g(1);
U  = t1;
