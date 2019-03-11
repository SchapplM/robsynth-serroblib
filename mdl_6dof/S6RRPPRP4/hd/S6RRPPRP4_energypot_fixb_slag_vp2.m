% Calculate potential energy for
% S6RRPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
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
% Datum: 2019-03-09 08:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPPRP4_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP4_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP4_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP4_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP4_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP4_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:36:56
% EndTime: 2019-03-09 08:36:57
% DurationCPUTime: 0.56s
% Computational Cost: add. (203->78), mult. (401->87), div. (0->0), fcn. (430->8), ass. (0->35)
t131 = -m(6) - m(7);
t130 = -mrSges(4,1) - mrSges(5,1);
t129 = mrSges(2,2) - mrSges(3,3);
t128 = mrSges(4,2) - mrSges(5,3);
t101 = sin(qJ(2));
t104 = cos(qJ(2));
t127 = -t104 * mrSges(3,1) + t101 * mrSges(3,2) - mrSges(2,1);
t126 = -mrSges(4,3) - mrSges(5,2) + mrSges(6,3) + mrSges(7,2);
t125 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t124 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t102 = sin(qJ(1));
t105 = cos(qJ(1));
t123 = pkin(1) * t105 + pkin(7) * t102;
t98 = sin(pkin(9));
t122 = t101 * t98;
t99 = cos(pkin(9));
t121 = t101 * t99;
t120 = t101 * t102;
t119 = t101 * t105;
t118 = t102 * t104;
t117 = t104 * t105;
t116 = pkin(1) * t102 - pkin(7) * t105;
t115 = pkin(2) * t117 + qJ(3) * t119 + t123;
t114 = pkin(2) * t101 - qJ(3) * t104 + pkin(6);
t112 = pkin(2) * t118 + qJ(3) * t120 + t116;
t111 = pkin(3) * t121 + qJ(4) * t122 + t114;
t82 = -t102 * t99 + t117 * t98;
t83 = t102 * t98 + t117 * t99;
t110 = pkin(3) * t83 + qJ(4) * t82 + t115;
t80 = t105 * t99 + t118 * t98;
t81 = -t105 * t98 + t118 * t99;
t108 = pkin(3) * t81 + qJ(4) * t80 + t112;
t103 = cos(qJ(5));
t100 = sin(qJ(5));
t1 = (-m(4) * t114 - m(5) * t111 - mrSges(1,3) - mrSges(2,3) + t131 * (pkin(4) * t121 + pkin(8) * t104 + t111) + t124 * (t100 * t121 - t103 * t122) + (-m(2) - m(3)) * pkin(6) + (-mrSges(3,2) - t126) * t104 + (t125 * (t100 * t98 + t103 * t99) + t128 * t98 + t130 * t99 - mrSges(3,1)) * t101) * g(3) + (-m(3) * t116 - m(4) * t112 - m(5) * t108 - mrSges(1,2) + t130 * t81 + t128 * t80 + t131 * (pkin(4) * t81 - pkin(8) * t120 + t108) + t125 * (t100 * t80 + t103 * t81) + t124 * (t100 * t81 - t103 * t80) - t129 * t105 + t127 * t102 + t126 * t120) * g(2) + (-m(3) * t123 - m(4) * t115 - m(5) * t110 - mrSges(1,1) + t130 * t83 + t128 * t82 + t131 * (pkin(4) * t83 - pkin(8) * t119 + t110) + t125 * (t100 * t82 + t103 * t83) + t124 * (t100 * t83 - t103 * t82) + t127 * t105 + t129 * t102 + t126 * t119) * g(1);
U  = t1;
