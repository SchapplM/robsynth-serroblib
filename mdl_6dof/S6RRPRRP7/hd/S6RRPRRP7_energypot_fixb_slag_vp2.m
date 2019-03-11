% Calculate potential energy for
% S6RRPRRP7
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
% Datum: 2019-03-09 12:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRP7_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP7_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP7_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP7_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP7_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP7_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:16:40
% EndTime: 2019-03-09 12:16:40
% DurationCPUTime: 0.51s
% Computational Cost: add. (186->72), mult. (360->79), div. (0->0), fcn. (377->8), ass. (0->35)
t124 = -mrSges(3,1) - mrSges(4,1);
t123 = mrSges(3,2) - mrSges(4,3);
t122 = -m(6) - m(7);
t121 = mrSges(5,2) - mrSges(7,2) - mrSges(6,3);
t94 = sin(qJ(2));
t98 = cos(qJ(2));
t120 = t123 * t94 + t124 * t98 - mrSges(2,1);
t119 = mrSges(2,2) - mrSges(3,3) - mrSges(4,2) + mrSges(5,3);
t118 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t117 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t97 = cos(qJ(4));
t116 = t94 * t97;
t95 = sin(qJ(1));
t115 = t95 * t98;
t99 = cos(qJ(1));
t114 = t98 * t99;
t113 = t99 * pkin(1) + t95 * pkin(7);
t112 = qJ(3) * t94;
t111 = t95 * pkin(1) - pkin(7) * t99;
t110 = pkin(2) * t114 + t99 * t112 + t113;
t109 = t94 * pkin(2) - qJ(3) * t98 + pkin(6);
t93 = sin(qJ(4));
t75 = t93 * t94 + t97 * t98;
t106 = pkin(2) * t115 + t95 * t112 + t111;
t105 = t94 * pkin(3) + t109;
t104 = pkin(3) * t115 + t99 * pkin(8) + t106;
t103 = pkin(3) * t114 - pkin(8) * t95 + t110;
t96 = cos(qJ(5));
t92 = sin(qJ(5));
t76 = -t93 * t98 + t116;
t73 = t75 * t99;
t72 = t93 * t114 - t99 * t116;
t71 = t75 * t95;
t70 = t93 * t115 - t95 * t116;
t1 = (-m(4) * t109 - m(5) * t105 - mrSges(1,3) - mrSges(2,3) - t123 * t98 + t124 * t94 + t122 * (t76 * pkin(4) + pkin(9) * t75 + t105) + (-m(2) - m(3)) * pkin(6) + (t117 * t92 + t118 * t96 - mrSges(5,1)) * t76 + t121 * t75) * g(3) + (-m(3) * t111 - m(4) * t106 - m(5) * t104 - t71 * mrSges(5,1) - mrSges(1,2) + t122 * (t71 * pkin(4) + t70 * pkin(9) + t104) + t118 * (t71 * t96 + t92 * t99) + t117 * (t71 * t92 - t99 * t96) + t120 * t95 + t121 * t70 - t119 * t99) * g(2) + (-m(3) * t113 - m(4) * t110 - m(5) * t103 - t73 * mrSges(5,1) - mrSges(1,1) + t122 * (t73 * pkin(4) + pkin(9) * t72 + t103) + t118 * (t73 * t96 - t92 * t95) + t117 * (t73 * t92 + t95 * t96) + t120 * t99 + t121 * t72 + t119 * t95) * g(1);
U  = t1;
