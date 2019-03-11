% Calculate potential energy for
% S6PRRPRP4
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
% Datum: 2019-03-08 21:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRPRP4_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP4_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP4_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP4_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:40:05
% EndTime: 2019-03-08 21:40:06
% DurationCPUTime: 0.58s
% Computational Cost: add. (289->96), mult. (621->109), div. (0->0), fcn. (727->10), ass. (0->50)
t143 = pkin(4) + pkin(8);
t149 = -mrSges(6,1) - mrSges(7,1);
t148 = -mrSges(6,2) - mrSges(7,2);
t112 = sin(qJ(5));
t147 = -m(7) * (pkin(5) * t112 + qJ(4)) + mrSges(4,2) - mrSges(5,3);
t115 = cos(qJ(5));
t146 = m(6) * t143 + m(7) * (t115 * pkin(5) + t143) + mrSges(5,1) + mrSges(4,3);
t145 = mrSges(3,2) - t146;
t144 = -m(6) * pkin(9) + m(7) * (-qJ(6) - pkin(9)) - mrSges(4,1) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t108 = sin(pkin(10));
t110 = cos(pkin(10));
t114 = sin(qJ(2));
t116 = cos(qJ(2));
t138 = cos(pkin(6));
t125 = t116 * t138;
t92 = t108 * t114 - t110 * t125;
t142 = t92 * pkin(8);
t94 = t108 * t125 + t110 * t114;
t141 = t94 * pkin(8);
t140 = cos(qJ(3));
t109 = sin(pkin(6));
t137 = t108 * t109;
t136 = t109 * t114;
t135 = t109 * t116;
t134 = t110 * t109;
t133 = t110 * pkin(1) + pkin(7) * t137;
t132 = t138 * pkin(7) + qJ(1);
t131 = pkin(8) * t135;
t126 = t114 * t138;
t95 = -t108 * t126 + t110 * t116;
t130 = t95 * pkin(2) + t133;
t129 = pkin(2) * t136 + t132;
t128 = t109 * t140;
t113 = sin(qJ(3));
t86 = t113 * t137 + t95 * t140;
t124 = t86 * pkin(3) + t130;
t97 = t138 * t113 + t114 * t128;
t123 = t97 * pkin(3) + t129;
t122 = t108 * pkin(1) - pkin(7) * t134;
t93 = t108 * t116 + t110 * t126;
t121 = t93 * pkin(2) + t122;
t84 = -t113 * t134 + t93 * t140;
t120 = t84 * pkin(3) + t121;
t85 = -t108 * t128 + t95 * t113;
t119 = t85 * qJ(4) + t124;
t96 = t113 * t136 - t138 * t140;
t118 = t96 * qJ(4) + t123;
t83 = t110 * t128 + t93 * t113;
t117 = t83 * qJ(4) + t120;
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - m(3) * t132 - t138 * mrSges(3,3) - (t114 * mrSges(3,1) + t116 * mrSges(3,2)) * t109 - m(4) * (t129 - t131) - m(5) * (t118 - t131) - m(6) * t118 - m(7) * t123 + t147 * t96 + t149 * (t96 * t112 - t115 * t135) + t148 * (t112 * t135 + t96 * t115) + t146 * t135 + t144 * t97) * g(3) + (-mrSges(1,2) - t108 * mrSges(2,1) - t110 * mrSges(2,2) - m(3) * t122 - t93 * mrSges(3,1) + mrSges(3,3) * t134 - m(4) * (t121 + t142) - m(5) * (t117 + t142) - m(6) * t117 - m(7) * t120 + t147 * t83 + t149 * (t83 * t112 + t92 * t115) + t148 * (-t92 * t112 + t83 * t115) + t145 * t92 + t144 * t84) * g(2) + (-mrSges(1,1) - t110 * mrSges(2,1) + t108 * mrSges(2,2) - m(3) * t133 - t95 * mrSges(3,1) - mrSges(3,3) * t137 - m(4) * (t130 + t141) - m(5) * (t119 + t141) - m(6) * t119 - m(7) * t124 + t147 * t85 + t149 * (t85 * t112 + t94 * t115) + t148 * (-t94 * t112 + t85 * t115) + t145 * t94 + t144 * t86) * g(1);
U  = t1;
