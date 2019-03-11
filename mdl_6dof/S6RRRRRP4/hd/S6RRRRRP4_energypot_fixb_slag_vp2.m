% Calculate potential energy for
% S6RRRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 01:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRRP4_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP4_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP4_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP4_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP4_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP4_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:11:55
% EndTime: 2019-03-10 01:11:55
% DurationCPUTime: 0.51s
% Computational Cost: add. (239->70), mult. (236->65), div. (0->0), fcn. (215->10), ass. (0->33)
t127 = m(5) * pkin(9) - mrSges(4,2) + mrSges(7,2) + mrSges(5,3) + mrSges(6,3);
t88 = sin(qJ(4));
t91 = cos(qJ(4));
t126 = -m(5) * pkin(3) - t91 * mrSges(5,1) + t88 * mrSges(5,2) - mrSges(4,1);
t122 = -m(6) - m(7);
t87 = qJ(2) + qJ(3);
t82 = sin(t87);
t84 = cos(t87);
t89 = sin(qJ(2));
t92 = cos(qJ(2));
t94 = -pkin(10) - pkin(9);
t124 = -m(3) * pkin(1) - t92 * mrSges(3,1) + t89 * mrSges(3,2) + t126 * t84 - mrSges(2,1) + (-t122 * t94 - t127) * t82;
t123 = -m(4) - m(5);
t117 = m(3) * pkin(7) + t88 * mrSges(5,1) + t91 * mrSges(5,2) - mrSges(2,2) + mrSges(3,3) + mrSges(4,3);
t116 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t115 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t114 = pkin(4) * t88;
t113 = pkin(2) * t89 + pkin(6);
t90 = sin(qJ(1));
t109 = t84 * t90;
t93 = cos(qJ(1));
t108 = t84 * t93;
t86 = qJ(4) + qJ(5);
t83 = cos(t86);
t107 = t90 * t83;
t79 = pkin(2) * t92 + pkin(1);
t95 = -pkin(8) - pkin(7);
t106 = t79 * t90 + t93 * t95;
t74 = t93 * t79;
t104 = -t90 * t95 + t74;
t81 = sin(t86);
t78 = pkin(4) * t91 + pkin(3);
t1 = (-mrSges(3,1) * t89 - mrSges(3,2) * t92 - mrSges(1,3) - mrSges(2,3) + t123 * t113 + t122 * (t78 * t82 + t84 * t94 + t113) + (-m(2) - m(3)) * pkin(6) + t127 * t84 + (t115 * t81 + t116 * t83 + t126) * t82) * g(3) + (-mrSges(1,2) + t122 * (t109 * t78 - t114 * t93 + t106) + t116 * (t107 * t84 - t81 * t93) + t115 * (t109 * t81 + t83 * t93) + t123 * t106 + t117 * t93 + t124 * t90) * g(2) + (-m(4) * t104 - m(5) * t74 - mrSges(1,1) + t122 * (t108 * t78 + t114 * t90 + t104) + t116 * (t108 * t83 + t81 * t90) + t115 * (t108 * t81 - t107) + (m(5) * t95 - t117) * t90 + t124 * t93) * g(1);
U  = t1;
