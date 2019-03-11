% Calculate potential energy for
% S6PRRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-03-08 22:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRRPP3_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP3_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPP3_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP3_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP3_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPP3_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:54:21
% EndTime: 2019-03-08 22:54:22
% DurationCPUTime: 0.53s
% Computational Cost: add. (335->95), mult. (748->111), div. (0->0), fcn. (903->10), ass. (0->54)
t148 = -mrSges(4,3) + mrSges(3,2);
t147 = mrSges(5,2) - mrSges(6,3) - mrSges(7,2);
t146 = mrSges(4,2) - mrSges(5,3) - mrSges(6,1) - m(7) * (pkin(5) + pkin(9)) - mrSges(7,1);
t145 = -m(7) * qJ(6) - mrSges(5,1) + mrSges(6,2) - mrSges(7,3);
t113 = cos(pkin(10));
t115 = sin(qJ(3));
t112 = sin(pkin(6));
t142 = cos(qJ(3));
t132 = t112 * t142;
t111 = sin(pkin(10));
t118 = cos(qJ(2));
t116 = sin(qJ(2));
t140 = cos(pkin(6));
t131 = t116 * t140;
t97 = t111 * t118 + t113 * t131;
t85 = t113 * t132 + t97 * t115;
t144 = pkin(9) * t85;
t99 = -t111 * t131 + t113 * t118;
t87 = -t111 * t132 + t115 * t99;
t143 = pkin(9) * t87;
t136 = t112 * t116;
t100 = t115 * t136 - t140 * t142;
t141 = t100 * pkin(9);
t139 = t111 * t112;
t138 = t112 * t113;
t137 = t112 * t115;
t135 = t112 * t118;
t134 = t113 * pkin(1) + pkin(7) * t139;
t133 = t140 * pkin(7) + qJ(1);
t130 = t118 * t140;
t128 = t111 * pkin(1) - pkin(7) * t138;
t98 = t111 * t130 + t113 * t116;
t127 = t99 * pkin(2) + pkin(8) * t98 + t134;
t88 = t111 * t137 + t99 * t142;
t126 = t88 * pkin(3) + t127;
t125 = pkin(2) * t136 - pkin(8) * t135 + t133;
t101 = t140 * t115 + t116 * t132;
t124 = t101 * pkin(3) + t125;
t96 = t111 * t116 - t113 * t130;
t123 = t97 * pkin(2) + pkin(8) * t96 + t128;
t86 = -t113 * t137 + t97 * t142;
t122 = t86 * pkin(3) + t123;
t114 = sin(qJ(4));
t117 = cos(qJ(4));
t80 = t114 * t88 - t98 * t117;
t81 = t114 * t98 + t117 * t88;
t121 = t81 * pkin(4) + t80 * qJ(5) + t126;
t89 = t101 * t114 + t117 * t135;
t90 = t101 * t117 - t114 * t135;
t120 = t90 * pkin(4) + t89 * qJ(5) + t124;
t78 = t114 * t86 - t96 * t117;
t79 = t114 * t96 + t117 * t86;
t119 = t79 * pkin(4) + t78 * qJ(5) + t122;
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - m(3) * t133 - t140 * mrSges(3,3) - (t116 * mrSges(3,1) + t118 * mrSges(3,2)) * t112 - m(4) * t125 - t101 * mrSges(4,1) + mrSges(4,3) * t135 - m(5) * (t124 + t141) - m(6) * (t120 + t141) - m(7) * t120 + t145 * t90 + t147 * t89 + t146 * t100) * g(3) + (-mrSges(1,2) - t111 * mrSges(2,1) - t113 * mrSges(2,2) - m(3) * t128 - t97 * mrSges(3,1) + mrSges(3,3) * t138 - m(4) * t123 - t86 * mrSges(4,1) - m(5) * (t122 + t144) - m(6) * (t119 + t144) - m(7) * t119 + t148 * t96 + t145 * t79 + t147 * t78 + t146 * t85) * g(2) + (-mrSges(1,1) - t113 * mrSges(2,1) + t111 * mrSges(2,2) - m(3) * t134 - t99 * mrSges(3,1) - mrSges(3,3) * t139 - m(4) * t127 - t88 * mrSges(4,1) - m(5) * (t126 + t143) - m(6) * (t121 + t143) - m(7) * t121 + t148 * t98 + t145 * t81 + t147 * t80 + t146 * t87) * g(1);
U  = t1;
