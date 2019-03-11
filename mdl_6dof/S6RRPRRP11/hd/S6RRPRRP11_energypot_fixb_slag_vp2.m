% Calculate potential energy for
% S6RRPRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
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
% Datum: 2019-03-09 12:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRP11_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP11_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP11_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP11_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP11_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP11_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:45:18
% EndTime: 2019-03-09 12:45:19
% DurationCPUTime: 0.56s
% Computational Cost: add. (176->77), mult. (246->75), div. (0->0), fcn. (221->8), ass. (0->34)
t79 = sin(qJ(4));
t102 = pkin(4) * t79;
t78 = qJ(4) + qJ(5);
t70 = sin(t78);
t113 = -m(6) * t102 - m(7) * (pkin(5) * t70 + t102) + mrSges(3,2) - mrSges(4,3);
t85 = -pkin(9) - pkin(8);
t112 = m(6) * t85 + m(7) * (-qJ(6) + t85) - mrSges(3,1) + mrSges(4,2) - mrSges(6,3) - mrSges(7,3);
t80 = sin(qJ(2));
t111 = t80 * mrSges(5,2);
t110 = mrSges(6,1) + mrSges(7,1);
t109 = mrSges(6,2) + mrSges(7,2);
t81 = sin(qJ(1));
t93 = qJ(3) * t80;
t83 = cos(qJ(2));
t97 = t81 * t83;
t108 = pkin(2) * t97 + t81 * t93;
t107 = -m(4) - m(6) - m(7);
t106 = m(5) - t107;
t105 = -m(5) * pkin(8) - mrSges(5,3);
t104 = t112 * t83 + t113 * t80 - mrSges(2,1);
t82 = cos(qJ(4));
t69 = t82 * pkin(4) + pkin(3);
t71 = cos(t78);
t103 = m(6) * t69 + m(7) * (pkin(5) * t71 + t69) - t79 * mrSges(5,2) + mrSges(4,1) - mrSges(2,2) + mrSges(3,3);
t100 = t80 * t81;
t84 = cos(qJ(1));
t99 = t80 * t84;
t98 = t81 * t82;
t96 = t82 * t84;
t95 = t83 * t84;
t94 = t84 * pkin(1) + t81 * pkin(7);
t74 = t81 * pkin(1);
t92 = -t84 * pkin(7) + t74;
t1 = (-mrSges(1,3) - mrSges(2,3) + (-m(2) - m(3)) * pkin(6) - t106 * (t80 * pkin(2) + pkin(6)) + (t79 * mrSges(5,1) + t82 * mrSges(5,2) + t106 * qJ(3) + t109 * t71 + t110 * t70 - t113) * t83 + (t105 + t112) * t80) * g(3) + (-mrSges(1,2) - m(3) * t92 - m(5) * (pkin(8) * t97 + t108 + t74) - (t79 * t100 - t96) * mrSges(5,1) - t98 * t111 - mrSges(5,3) * t97 - t110 * (t70 * t100 - t71 * t84) - t109 * (t71 * t100 + t70 * t84) + t107 * (t92 + t108) + (-m(5) * (-pkin(3) - pkin(7)) + t103) * t84 + t104 * t81) * g(2) + (-mrSges(1,1) - m(3) * t94 - (t79 * t99 + t98) * mrSges(5,1) - t96 * t111 + t105 * t95 - t110 * (t70 * t99 + t71 * t81) - t109 * (-t70 * t81 + t71 * t99) - t106 * (pkin(2) * t95 + t84 * t93 + t94) + (-m(5) * pkin(3) - t103) * t81 + t104 * t84) * g(1);
U  = t1;
