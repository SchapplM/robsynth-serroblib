% Calculate potential energy for
% S6RRPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
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
% Datum: 2019-03-09 09:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPPRR4_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR4_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_energypot_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR4_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR4_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:01:05
% EndTime: 2019-03-09 09:01:05
% DurationCPUTime: 0.67s
% Computational Cost: add. (336->90), mult. (733->111), div. (0->0), fcn. (881->12), ass. (0->48)
t148 = -m(6) - m(7);
t147 = mrSges(4,2) - mrSges(5,3);
t146 = -mrSges(3,3) - mrSges(4,3) - mrSges(5,1);
t145 = -m(3) * pkin(1) - mrSges(2,1);
t144 = m(7) * pkin(10) - mrSges(6,2) + mrSges(7,3);
t143 = m(3) * pkin(8) - t146;
t112 = sin(qJ(6));
t116 = cos(qJ(6));
t142 = -t112 * mrSges(7,1) - t116 * mrSges(7,2) - mrSges(4,1) + mrSges(5,2) - mrSges(6,3);
t141 = -m(7) * pkin(5) - t116 * mrSges(7,1) + t112 * mrSges(7,2) - mrSges(6,1);
t111 = cos(pkin(6));
t140 = t111 * pkin(8) + pkin(7);
t118 = cos(qJ(2));
t104 = pkin(2) * t118 + pkin(1);
t115 = sin(qJ(1));
t119 = cos(qJ(1));
t109 = sin(pkin(6));
t114 = sin(qJ(2));
t96 = pkin(2) * t111 * t114 + (-pkin(8) - qJ(3)) * t109;
t139 = t115 * t104 + t119 * t96;
t138 = t109 * t114;
t137 = t109 * t115;
t136 = t109 * t119;
t110 = cos(pkin(11));
t135 = t110 * t118;
t134 = t114 * t119;
t133 = t115 * t114;
t132 = t115 * t118;
t131 = t118 * t119;
t130 = t119 * t104 - t115 * t96;
t129 = pkin(2) * t138 + t111 * qJ(3) + t140;
t108 = sin(pkin(11));
t128 = t108 * t118 + t110 * t114;
t98 = -t108 * t114 + t135;
t125 = t98 * t111;
t83 = -t115 * t128 + t119 * t125;
t126 = t128 * t111;
t84 = t115 * t98 + t119 * t126;
t127 = t84 * pkin(3) - t83 * qJ(4) + t139;
t85 = -t115 * t125 - t119 * t128;
t86 = -t115 * t126 + t119 * t98;
t124 = t86 * pkin(3) - qJ(4) * t85 + t130;
t94 = t108 * t138 - t109 * t135;
t95 = t128 * t109;
t123 = t95 * pkin(3) + qJ(4) * t94 + t129;
t117 = cos(qJ(5));
t113 = sin(qJ(5));
t1 = (-mrSges(1,3) - m(2) * pkin(7) - mrSges(2,3) - m(3) * t140 - (t114 * mrSges(3,1) + t118 * mrSges(3,2)) * t109 - m(4) * t129 - m(5) * t123 + t147 * t94 + t148 * (t111 * pkin(4) + pkin(9) * t95 + t123) - t144 * (t111 * t113 - t94 * t117) + t141 * (t111 * t117 + t113 * t94) + t142 * t95 + t146 * t111) * g(3) + (-mrSges(1,2) - mrSges(2,2) * t119 - (t111 * t134 + t132) * mrSges(3,1) - (t111 * t131 - t133) * mrSges(3,2) - m(4) * t139 - m(5) * t127 - t147 * t83 + t148 * (-pkin(4) * t136 + t84 * pkin(9) + t127) + t144 * (t113 * t136 - t83 * t117) + t145 * t115 + t141 * (-t83 * t113 - t117 * t136) + t142 * t84 + t143 * t136) * g(2) + (-mrSges(1,1) + t115 * mrSges(2,2) - (-t111 * t133 + t131) * mrSges(3,1) - (-t111 * t132 - t134) * mrSges(3,2) - m(4) * t130 - m(5) * t124 - t147 * t85 + t148 * (pkin(4) * t137 + pkin(9) * t86 + t124) - t144 * (t113 * t137 + t85 * t117) + t145 * t119 + t141 * (-t113 * t85 + t117 * t137) + t142 * t86 - t143 * t137) * g(1);
U  = t1;
