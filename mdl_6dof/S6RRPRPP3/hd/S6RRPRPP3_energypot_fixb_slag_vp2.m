% Calculate potential energy for
% S6RRPRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
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
% Datum: 2019-03-09 09:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRPP3_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP3_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP3_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP3_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP3_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP3_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:54:09
% EndTime: 2019-03-09 09:54:10
% DurationCPUTime: 0.54s
% Computational Cost: add. (224->81), mult. (289->74), div. (0->0), fcn. (278->8), ass. (0->40)
t87 = sin(pkin(9));
t88 = cos(pkin(9));
t124 = -m(4) * pkin(2) - t88 * mrSges(4,1) + t87 * mrSges(4,2) - mrSges(3,1);
t123 = -m(4) * qJ(3) - mrSges(7,1) + mrSges(3,2) - mrSges(4,3);
t122 = -m(6) - m(7);
t121 = -m(3) - m(4);
t120 = -mrSges(5,3) - mrSges(6,1);
t91 = sin(qJ(1));
t92 = cos(qJ(2));
t110 = t91 * t92;
t86 = pkin(9) + qJ(4);
t81 = sin(t86);
t82 = cos(t86);
t93 = cos(qJ(1));
t67 = t110 * t81 + t82 * t93;
t108 = t93 * t81;
t68 = t110 * t82 - t108;
t119 = t68 * pkin(4) + t67 * qJ(5);
t118 = -t87 * mrSges(4,1) - t88 * mrSges(4,2) + mrSges(2,2) - mrSges(3,3);
t117 = mrSges(5,2) - mrSges(7,2) - mrSges(6,3);
t89 = -pkin(8) - qJ(3);
t90 = sin(qJ(2));
t115 = -mrSges(2,1) + t124 * t92 + (-m(7) * (pkin(5) - t89) + t123) * t90;
t114 = -m(7) * qJ(6) - mrSges(5,1) + mrSges(6,2) - mrSges(7,3);
t113 = pkin(3) * t87;
t112 = t90 * t91;
t111 = t90 * t93;
t109 = t92 * t93;
t107 = t93 * pkin(1) + t91 * pkin(7);
t104 = t89 * t111;
t84 = t91 * pkin(1);
t103 = -pkin(7) * t93 + t84;
t79 = pkin(3) * t88 + pkin(2);
t102 = t79 * t109 + t91 * t113 + t107;
t69 = t108 * t92 - t91 * t82;
t70 = t109 * t82 + t91 * t81;
t96 = t70 * pkin(4) + t69 * qJ(5) + t102;
t71 = t79 * t110;
t95 = -t89 * t112 + t71 + t84 + (-pkin(7) - t113) * t93;
t1 = (-mrSges(1,3) - mrSges(2,3) + (-m(2) + t121) * pkin(6) + (m(7) * pkin(5) - t120 - t123) * t92 + (-m(5) + t122) * (t90 * t79 + t92 * t89 + pkin(6)) + (t122 * (pkin(4) * t82 + qJ(5) * t81) + t114 * t82 + t117 * t81 + t124) * t90) * g(3) + (-mrSges(1,2) - m(3) * t103 - m(4) * t84 - m(5) * t95 - m(6) * (t95 + t119) - m(7) * (t103 + t71 + t119) + t120 * t112 + (m(4) * pkin(7) + m(7) * t113 - t118) * t93 + t114 * t68 + t117 * t67 + t115 * t91) * g(2) + (-mrSges(1,1) - m(5) * (t102 - t104) - m(6) * (t96 - t104) - m(7) * t96 + t120 * t111 + t121 * t107 + t118 * t91 + t114 * t70 + t117 * t69 + t115 * t93) * g(1);
U  = t1;
