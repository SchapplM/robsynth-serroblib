% Calculate potential energy for
% S6RRRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 19:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPRR10_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR10_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR10_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR10_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR10_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR10_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:15:31
% EndTime: 2019-03-09 19:15:32
% DurationCPUTime: 0.66s
% Computational Cost: add. (206->81), mult. (371->79), div. (0->0), fcn. (388->10), ass. (0->38)
t118 = -m(5) - m(6);
t117 = mrSges(2,2) - mrSges(3,3);
t83 = sin(qJ(3));
t84 = sin(qJ(2));
t87 = cos(qJ(3));
t116 = (pkin(3) * t87 + qJ(4) * t83) * t84;
t115 = mrSges(5,2) + mrSges(4,3) - mrSges(6,3) - mrSges(7,3);
t88 = cos(qJ(2));
t90 = -pkin(10) - pkin(9);
t114 = -mrSges(3,1) * t88 - mrSges(2,1) + (-m(7) * t90 + mrSges(3,2)) * t84;
t113 = m(6) * pkin(9) - t115;
t81 = qJ(5) + qJ(6);
t75 = sin(t81);
t76 = cos(t81);
t86 = cos(qJ(5));
t112 = -t75 * mrSges(7,1) - t86 * mrSges(6,2) - t76 * mrSges(7,2) + mrSges(4,2) - mrSges(5,3);
t82 = sin(qJ(5));
t111 = -m(7) * (pkin(5) * t82 + qJ(4)) - t82 * mrSges(6,1) + t112;
t110 = -m(6) * pkin(4) - m(7) * (pkin(5) * t86 + pkin(4)) - t86 * mrSges(6,1) - t76 * mrSges(7,1) + t82 * mrSges(6,2) + t75 * mrSges(7,2) - mrSges(4,1) - mrSges(5,1);
t109 = t84 * pkin(2) + pkin(6);
t85 = sin(qJ(1));
t107 = t84 * t85;
t89 = cos(qJ(1));
t106 = t84 * t89;
t104 = t85 * t88;
t103 = t88 * t89;
t102 = t89 * pkin(1) + t85 * pkin(7);
t101 = t85 * pkin(1) - t89 * pkin(7);
t98 = pkin(2) * t103 + pkin(8) * t106 + t102;
t97 = -t88 * pkin(8) + t109;
t66 = t103 * t87 + t83 * t85;
t96 = t66 * pkin(3) + t98;
t94 = pkin(2) * t104 + pkin(8) * t107 + t101;
t64 = t104 * t87 - t83 * t89;
t93 = t64 * pkin(3) + t94;
t65 = t103 * t83 - t85 * t87;
t63 = t104 * t83 + t87 * t89;
t1 = (-mrSges(1,3) - mrSges(2,3) - m(4) * t97 - m(5) * (t97 + t116) + (-m(6) - m(7)) * (t109 + t116) + (-m(2) - m(3)) * pkin(6) + (-mrSges(3,2) - m(6) * (-pkin(8) + pkin(9)) - m(7) * (-pkin(8) - t90) + t115) * t88 + (t110 * t87 - mrSges(3,1) + ((-m(7) * pkin(5) - mrSges(6,1)) * t82 + t112) * t83) * t84) * g(3) + (-m(3) * t101 - m(4) * t94 - m(7) * t93 - mrSges(1,2) + t118 * (t63 * qJ(4) + t93) - t117 * t89 + t114 * t85 + t110 * t64 + t111 * t63 + t113 * t107) * g(2) + (-m(3) * t102 - m(4) * t98 - m(7) * t96 - mrSges(1,1) + t118 * (t65 * qJ(4) + t96) + t114 * t89 + t117 * t85 + t110 * t66 + t111 * t65 + t113 * t106) * g(1);
U  = t1;
