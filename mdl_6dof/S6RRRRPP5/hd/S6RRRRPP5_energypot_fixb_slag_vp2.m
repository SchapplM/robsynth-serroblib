% Calculate potential energy for
% S6RRRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
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
% Datum: 2019-03-09 21:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRPP5_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP5_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP5_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP5_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP5_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP5_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:05:01
% EndTime: 2019-03-09 21:05:01
% DurationCPUTime: 0.51s
% Computational Cost: add. (224->81), mult. (289->74), div. (0->0), fcn. (278->8), ass. (0->40)
t88 = sin(qJ(3));
t91 = cos(qJ(3));
t125 = -m(4) * pkin(2) - t91 * mrSges(4,1) + t88 * mrSges(4,2) - mrSges(3,1);
t124 = -m(4) * pkin(8) + mrSges(3,2) - mrSges(4,3) + mrSges(7,3);
t123 = -m(6) - m(7);
t122 = -m(3) - m(4);
t121 = -mrSges(5,3) - mrSges(6,2);
t87 = qJ(3) + qJ(4);
t83 = cos(t87);
t93 = cos(qJ(1));
t110 = t93 * t83;
t90 = sin(qJ(1));
t92 = cos(qJ(2));
t112 = t90 * t92;
t82 = sin(t87);
t69 = t112 * t82 + t110;
t111 = t93 * t82;
t70 = t112 * t83 - t111;
t120 = t70 * pkin(4) + t69 * qJ(5);
t119 = mrSges(5,2) - mrSges(7,2) - mrSges(6,3);
t118 = -t88 * mrSges(4,1) - t91 * mrSges(4,2) + mrSges(2,2) - mrSges(3,3);
t89 = sin(qJ(2));
t94 = -pkin(9) - pkin(8);
t116 = -mrSges(2,1) + t125 * t92 + (-m(7) * (-qJ(6) - t94) + t124) * t89;
t115 = -m(7) * pkin(5) - mrSges(5,1) - mrSges(6,1) - mrSges(7,1);
t114 = pkin(3) * t88;
t113 = t90 * t89;
t109 = t93 * t89;
t108 = t93 * pkin(1) + t90 * pkin(7);
t105 = t94 * t109;
t85 = t90 * pkin(1);
t104 = -t93 * pkin(7) + t85;
t80 = t91 * pkin(3) + pkin(2);
t103 = t93 * t92 * t80 + t90 * t114 + t108;
t71 = t111 * t92 - t90 * t83;
t72 = t110 * t92 + t90 * t82;
t97 = t72 * pkin(4) + t71 * qJ(5) + t103;
t73 = t80 * t112;
t96 = -t94 * t113 + t73 + t85 + (-pkin(7) - t114) * t93;
t1 = (-mrSges(1,3) - mrSges(2,3) + (-m(2) + t122) * pkin(6) + (-m(7) * qJ(6) - t121 - t124) * t92 + (-m(5) + t123) * (t89 * t80 + t92 * t94 + pkin(6)) + (t123 * (pkin(4) * t83 + qJ(5) * t82) + t115 * t83 + t119 * t82 + t125) * t89) * g(3) + (-mrSges(1,2) - m(3) * t104 - m(4) * t85 - m(5) * t96 - m(6) * (t96 + t120) - m(7) * (t104 + t73 + t120) + t121 * t113 + (m(4) * pkin(7) + m(7) * t114 - t118) * t93 + t115 * t70 + t119 * t69 + t116 * t90) * g(2) + (-mrSges(1,1) - m(5) * (t103 - t105) - m(6) * (t97 - t105) - m(7) * t97 + t121 * t109 + t122 * t108 + t118 * t90 + t115 * t72 + t119 * t71 + t116 * t93) * g(1);
U  = t1;
