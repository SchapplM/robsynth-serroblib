% Calculate potential energy for
% S6RPRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
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
% Datum: 2019-03-09 03:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPRP4_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP4_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP4_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP4_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP4_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP4_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:11:02
% EndTime: 2019-03-09 03:11:03
% DurationCPUTime: 0.42s
% Computational Cost: add. (209->67), mult. (209->65), div. (0->0), fcn. (184->8), ass. (0->31)
t112 = mrSges(4,2) - mrSges(5,3);
t111 = -mrSges(4,1) + mrSges(5,2) - mrSges(6,3) - mrSges(7,2);
t110 = -m(6) - m(7);
t80 = qJ(1) + pkin(9);
t74 = sin(t80);
t86 = cos(qJ(3));
t104 = t74 * t86;
t83 = sin(qJ(3));
t97 = qJ(4) * t83;
t109 = pkin(3) * t104 + t74 * t97;
t108 = mrSges(3,2) - mrSges(4,3) - mrSges(5,1);
t107 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t106 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t105 = t111 * t86 + t112 * t83 - mrSges(3,1);
t84 = sin(qJ(1));
t78 = t84 * pkin(1);
t87 = cos(qJ(1));
t79 = t87 * pkin(1);
t75 = cos(t80);
t103 = t75 * t86;
t82 = sin(qJ(5));
t102 = t82 * t83;
t85 = cos(qJ(5));
t101 = t83 * t85;
t81 = qJ(2) + pkin(6);
t98 = t74 * pkin(2) + t78;
t96 = t83 * pkin(3) + t81;
t95 = t75 * pkin(2) + t74 * pkin(7) + t79;
t93 = -t75 * pkin(7) + t98;
t92 = pkin(3) * t103 + t75 * t97 + t95;
t1 = (-m(2) * pkin(6) - m(5) * t96 - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + t110 * (t83 * pkin(8) + t96) + (-m(3) - m(4)) * t81 + (-t106 * t85 + t107 * t82 + (m(5) - t110) * qJ(4) - t112) * t86 + t111 * t83) * g(3) + (-mrSges(1,2) - t84 * mrSges(2,1) - t87 * mrSges(2,2) - m(3) * t78 - m(4) * t93 - m(5) * (t93 + t109) + t110 * (pkin(8) * t104 + (-pkin(4) - pkin(7)) * t75 + t98 + t109) - t107 * (t74 * t102 - t75 * t85) + t106 * (t74 * t101 + t75 * t82) - t108 * t75 + t105 * t74) * g(2) + (-m(3) * t79 - m(4) * t95 - m(5) * t92 - t87 * mrSges(2,1) + t84 * mrSges(2,2) - mrSges(1,1) + t110 * (t74 * pkin(4) + pkin(8) * t103 + t92) - t107 * (t75 * t102 + t74 * t85) - t106 * (-t75 * t101 + t74 * t82) + t108 * t74 + t105 * t75) * g(1);
U  = t1;
