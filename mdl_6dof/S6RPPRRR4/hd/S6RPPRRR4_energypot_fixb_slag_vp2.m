% Calculate potential energy for
% S6RPPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 02:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPRRR4_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR4_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR4_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR4_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR4_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR4_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:25:18
% EndTime: 2019-03-09 02:25:18
% DurationCPUTime: 0.40s
% Computational Cost: add. (183->59), mult. (282->50), div. (0->0), fcn. (305->10), ass. (0->25)
t102 = m(5) + m(6) + m(7);
t68 = qJ(5) + qJ(6);
t61 = sin(t68);
t62 = cos(t68);
t70 = sin(qJ(5));
t72 = cos(qJ(5));
t101 = mrSges(5,1) + m(6) * pkin(4) + t72 * mrSges(6,1) - t70 * mrSges(6,2) + m(7) * (pkin(5) * t72 + pkin(4)) + t62 * mrSges(7,1) - t61 * mrSges(7,2);
t100 = mrSges(5,2) + m(7) * (-pkin(9) - pkin(8)) - mrSges(7,3) - m(6) * pkin(8) - mrSges(6,3);
t99 = -mrSges(2,1) - mrSges(3,1);
t98 = mrSges(2,2) - mrSges(3,3);
t71 = sin(qJ(4));
t73 = cos(qJ(4));
t94 = -t100 * t71 + t101 * t73 + mrSges(4,1);
t93 = t61 * mrSges(7,1) + t72 * mrSges(6,2) + t62 * mrSges(7,2) - mrSges(4,2) + mrSges(5,3) + (m(7) * pkin(5) + mrSges(6,1)) * t70 + t102 * pkin(7);
t92 = cos(qJ(1));
t91 = sin(qJ(1));
t90 = t92 * pkin(1) + t91 * qJ(2);
t89 = cos(pkin(10));
t88 = sin(pkin(10));
t87 = t92 * pkin(2) + t90;
t84 = t91 * pkin(1) - t92 * qJ(2);
t81 = t91 * pkin(2) + t84;
t56 = t92 * t88 - t91 * t89;
t55 = -t91 * t88 - t92 * t89;
t1 = (-mrSges(3,2) - mrSges(1,3) - mrSges(2,3) + mrSges(4,3) + (-m(2) - m(3)) * pkin(6) + (-m(4) - t102) * (-qJ(3) + pkin(6)) + t100 * t73 + t101 * t71) * g(3) + (-m(3) * t84 - m(4) * t81 - mrSges(1,2) - t98 * t92 + t99 * t91 - t102 * (-t56 * pkin(3) + t81) + t94 * t56 + t93 * t55) * g(2) + (-m(3) * t90 - m(4) * t87 - mrSges(1,1) + t99 * t92 + t98 * t91 - t102 * (-t55 * pkin(3) + t87) - t93 * t56 + t94 * t55) * g(1);
U  = t1;
