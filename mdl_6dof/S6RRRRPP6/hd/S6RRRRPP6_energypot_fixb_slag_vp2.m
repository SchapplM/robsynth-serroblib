% Calculate potential energy for
% S6RRRRPP6
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
% Datum: 2019-03-09 21:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRPP6_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP6_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP6_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP6_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP6_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP6_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:11:04
% EndTime: 2019-03-09 21:11:05
% DurationCPUTime: 0.59s
% Computational Cost: add. (224->79), mult. (289->71), div. (0->0), fcn. (278->8), ass. (0->38)
t132 = -m(4) * pkin(8) - mrSges(6,1) - mrSges(7,1) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t91 = sin(qJ(3));
t94 = cos(qJ(3));
t131 = -m(4) * pkin(2) - t94 * mrSges(4,1) + t91 * mrSges(4,2) - mrSges(3,1);
t129 = -m(6) - m(7);
t92 = sin(qJ(2));
t95 = cos(qJ(2));
t97 = -pkin(9) - pkin(8);
t128 = t131 * t95 - mrSges(2,1) + (-m(7) * (pkin(5) - t97) + t132) * t92;
t127 = -m(3) - m(4);
t93 = sin(qJ(1));
t114 = t93 * t95;
t90 = qJ(3) + qJ(4);
t85 = sin(t90);
t86 = cos(t90);
t96 = cos(qJ(1));
t71 = t114 * t85 + t86 * t96;
t112 = t96 * t85;
t72 = t114 * t86 - t112;
t125 = t72 * pkin(4) + t71 * qJ(5);
t124 = mrSges(5,2) - mrSges(7,2) - mrSges(6,3);
t123 = -t91 * mrSges(4,1) - t94 * mrSges(4,2) + mrSges(2,2) - mrSges(3,3);
t119 = -m(7) * qJ(6) - mrSges(5,1) + mrSges(6,2) - mrSges(7,3);
t118 = pkin(3) * t91;
t115 = t92 * t97;
t113 = t95 * t96;
t111 = t96 * pkin(1) + t93 * pkin(7);
t108 = t96 * t115;
t88 = t93 * pkin(1);
t107 = -t96 * pkin(7) + t88;
t83 = pkin(3) * t94 + pkin(2);
t106 = t83 * t113 + t93 * t118 + t111;
t73 = t112 * t95 - t93 * t86;
t74 = t113 * t86 + t85 * t93;
t100 = t74 * pkin(4) + t73 * qJ(5) + t106;
t75 = t83 * t114;
t99 = -t93 * t115 + t75 + t88 + (-pkin(7) - t118) * t96;
t1 = (-mrSges(1,3) - mrSges(2,3) + (-m(2) + t127) * pkin(6) + (m(7) * pkin(5) - t132) * t95 + (-m(5) + t129) * (t92 * t83 + t95 * t97 + pkin(6)) + (t129 * (pkin(4) * t86 + qJ(5) * t85) + t119 * t86 + t124 * t85 + t131) * t92) * g(3) + (-mrSges(1,2) - m(3) * t107 - m(4) * t88 - m(5) * t99 - m(6) * (t99 + t125) - m(7) * (t107 + t75 + t125) + (m(4) * pkin(7) + m(7) * t118 - t123) * t96 + t119 * t72 + t124 * t71 + t128 * t93) * g(2) + (-mrSges(1,1) - m(5) * (t106 - t108) - m(6) * (t100 - t108) - m(7) * t100 + t127 * t111 + t123 * t93 + t119 * t74 + t124 * t73 + t128 * t96) * g(1);
U  = t1;
