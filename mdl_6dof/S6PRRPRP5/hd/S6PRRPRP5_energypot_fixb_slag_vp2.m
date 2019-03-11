% Calculate potential energy for
% S6PRRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
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
% Datum: 2019-03-08 21:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRPRP5_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP5_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP5_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP5_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP5_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP5_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:45:32
% EndTime: 2019-03-08 21:45:32
% DurationCPUTime: 0.54s
% Computational Cost: add. (301->90), mult. (659->107), div. (0->0), fcn. (780->10), ass. (0->45)
t145 = -m(6) - m(7);
t144 = mrSges(4,2) - mrSges(5,3);
t143 = mrSges(4,3) + mrSges(5,1);
t142 = mrSges(3,2) - t143;
t141 = -mrSges(4,1) + mrSges(5,2) - mrSges(6,3) - mrSges(7,2);
t140 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t139 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t138 = cos(qJ(3));
t137 = cos(pkin(6));
t109 = sin(pkin(10));
t110 = sin(pkin(6));
t136 = t109 * t110;
t114 = sin(qJ(2));
t135 = t110 * t114;
t116 = cos(qJ(2));
t134 = t110 * t116;
t111 = cos(pkin(10));
t133 = t111 * t110;
t132 = t111 * pkin(1) + pkin(7) * t136;
t131 = t137 * pkin(7) + qJ(1);
t130 = pkin(8) * t134;
t129 = pkin(2) * t135 + t131;
t128 = t110 * t138;
t127 = t114 * t137;
t126 = t116 * t137;
t125 = t109 * pkin(1) - pkin(7) * t133;
t96 = t109 * t126 + t111 * t114;
t97 = -t109 * t127 + t111 * t116;
t124 = t97 * pkin(2) + t96 * pkin(8) + t132;
t113 = sin(qJ(3));
t98 = t113 * t135 - t137 * t138;
t99 = t137 * t113 + t114 * t128;
t123 = t99 * pkin(3) + t98 * qJ(4) + t129;
t94 = t109 * t114 - t111 * t126;
t95 = t109 * t116 + t111 * t127;
t122 = t95 * pkin(2) + t94 * pkin(8) + t125;
t85 = -t109 * t128 + t97 * t113;
t86 = t113 * t136 + t97 * t138;
t121 = t86 * pkin(3) + t85 * qJ(4) + t124;
t83 = t111 * t128 + t95 * t113;
t84 = -t113 * t133 + t95 * t138;
t120 = t84 * pkin(3) + t83 * qJ(4) + t122;
t115 = cos(qJ(5));
t112 = sin(qJ(5));
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - m(3) * t131 - t137 * mrSges(3,3) - (t114 * mrSges(3,1) + t116 * mrSges(3,2)) * t110 - m(4) * (t129 - t130) - m(5) * (t123 - t130) + t144 * t98 + t145 * (t99 * pkin(9) + (-pkin(4) - pkin(8)) * t134 + t123) + t140 * (t98 * t112 - t115 * t134) + t139 * (t112 * t134 + t98 * t115) + t143 * t134 + t141 * t99) * g(3) + (-m(3) * t125 - m(4) * t122 - m(5) * t120 - t109 * mrSges(2,1) - t95 * mrSges(3,1) - t111 * mrSges(2,2) + mrSges(3,3) * t133 - mrSges(1,2) + t144 * t83 + t145 * (t94 * pkin(4) + t84 * pkin(9) + t120) + t140 * (t83 * t112 + t94 * t115) - t139 * (t94 * t112 - t83 * t115) + t142 * t94 + t141 * t84) * g(2) + (-m(3) * t132 - m(4) * t124 - m(5) * t121 - t111 * mrSges(2,1) - t97 * mrSges(3,1) + t109 * mrSges(2,2) - mrSges(3,3) * t136 - mrSges(1,1) + t144 * t85 + t145 * (t96 * pkin(4) + t86 * pkin(9) + t121) + t140 * (t85 * t112 + t96 * t115) - t139 * (t96 * t112 - t85 * t115) + t142 * t96 + t141 * t86) * g(1);
U  = t1;
