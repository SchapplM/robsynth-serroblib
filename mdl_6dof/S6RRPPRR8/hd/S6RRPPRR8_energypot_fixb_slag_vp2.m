% Calculate potential energy for
% S6RRPPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
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
% Datum: 2019-03-09 09:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPPRR8_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR8_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR8_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR8_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR8_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR8_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:22:46
% EndTime: 2019-03-09 09:22:47
% DurationCPUTime: 0.70s
% Computational Cost: add. (206->78), mult. (371->76), div. (0->0), fcn. (388->10), ass. (0->36)
t124 = mrSges(5,2) + mrSges(4,3) - mrSges(6,3) - mrSges(7,3) - mrSges(3,2);
t88 = sin(qJ(2));
t91 = cos(qJ(2));
t93 = -pkin(9) - pkin(8);
t123 = -t91 * mrSges(3,1) - mrSges(2,1) + (m(6) * pkin(8) - m(7) * t93 - t124) * t88;
t122 = -m(5) - m(6);
t121 = mrSges(2,2) - mrSges(3,3);
t85 = sin(pkin(10));
t86 = cos(pkin(10));
t120 = (pkin(3) * t86 + qJ(4) * t85) * t88;
t84 = qJ(5) + qJ(6);
t78 = sin(t84);
t79 = cos(t84);
t90 = cos(qJ(5));
t116 = -t78 * mrSges(7,1) - t90 * mrSges(6,2) - t79 * mrSges(7,2) + mrSges(4,2) - mrSges(5,3);
t87 = sin(qJ(5));
t115 = -m(7) * (pkin(5) * t87 + qJ(4)) - t87 * mrSges(6,1) + t116;
t114 = -m(6) * pkin(4) - m(7) * (pkin(5) * t90 + pkin(4)) - t90 * mrSges(6,1) - t79 * mrSges(7,1) + t87 * mrSges(6,2) + t78 * mrSges(7,2) - mrSges(4,1) - mrSges(5,1);
t113 = t88 * pkin(2) + pkin(6);
t89 = sin(qJ(1));
t108 = t89 * t91;
t92 = cos(qJ(1));
t107 = t91 * t92;
t106 = t92 * pkin(1) + t89 * pkin(7);
t105 = qJ(3) * t88;
t104 = t89 * pkin(1) - t92 * pkin(7);
t101 = pkin(2) * t107 + t92 * t105 + t106;
t100 = -t91 * qJ(3) + t113;
t69 = t107 * t86 + t85 * t89;
t99 = t69 * pkin(3) + t101;
t97 = pkin(2) * t108 + t89 * t105 + t104;
t67 = t108 * t86 - t85 * t92;
t96 = t67 * pkin(3) + t97;
t68 = t107 * t85 - t89 * t86;
t66 = t108 * t85 + t86 * t92;
t1 = (-mrSges(1,3) - mrSges(2,3) - m(4) * t100 - m(5) * (t100 + t120) + (-m(6) - m(7)) * (t113 + t120) + (-m(2) - m(3)) * pkin(6) + (-m(6) * (pkin(8) - qJ(3)) - m(7) * (-qJ(3) - t93) + t124) * t91 + (t114 * t86 - mrSges(3,1) + ((-m(7) * pkin(5) - mrSges(6,1)) * t87 + t116) * t85) * t88) * g(3) + (-m(3) * t104 - m(4) * t97 - m(7) * t96 - mrSges(1,2) + t122 * (t66 * qJ(4) + t96) - t121 * t92 + t114 * t67 + t115 * t66 + t123 * t89) * g(2) + (-m(3) * t106 - m(4) * t101 - m(7) * t99 - mrSges(1,1) + t122 * (t68 * qJ(4) + t99) + t121 * t89 + t114 * t69 + t115 * t68 + t123 * t92) * g(1);
U  = t1;
