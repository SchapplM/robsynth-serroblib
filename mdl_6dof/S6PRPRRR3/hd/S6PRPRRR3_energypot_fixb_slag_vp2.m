% Calculate potential energy for
% S6PRPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 20:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRRR3_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR3_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR3_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR3_energypot_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR3_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR3_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:31:35
% EndTime: 2019-03-08 20:31:36
% DurationCPUTime: 0.85s
% Computational Cost: add. (357->119), mult. (480->136), div. (0->0), fcn. (534->14), ass. (0->45)
t135 = -m(6) - m(7);
t134 = -m(7) * pkin(10) + mrSges(6,2) - mrSges(7,3);
t112 = -pkin(8) - qJ(3);
t133 = m(4) * qJ(3) - m(5) * t112 - mrSges(3,2) + mrSges(4,3) + mrSges(5,3);
t113 = sin(qJ(6));
t115 = cos(qJ(6));
t132 = -m(7) * pkin(5) - t115 * mrSges(7,1) + t113 * mrSges(7,2) - mrSges(6,1);
t131 = -t113 * mrSges(7,1) - t115 * mrSges(7,2) - mrSges(6,3) - t133;
t106 = sin(pkin(12));
t130 = pkin(3) * t106;
t109 = cos(pkin(12));
t96 = t109 * pkin(3) + pkin(2);
t110 = cos(pkin(11));
t107 = sin(pkin(11));
t108 = sin(pkin(6));
t128 = t107 * t108;
t129 = t110 * pkin(1) + pkin(7) * t128;
t127 = t108 * t110;
t114 = sin(qJ(2));
t126 = t108 * t114;
t116 = cos(qJ(2));
t125 = t108 * t116;
t111 = cos(pkin(6));
t124 = t111 * t114;
t123 = t111 * t116;
t122 = t111 * pkin(7) + qJ(1);
t105 = pkin(12) + qJ(4);
t121 = t106 * t128;
t100 = t107 * pkin(1);
t120 = -pkin(7) * t127 + t100;
t104 = -pkin(9) + t112;
t98 = cos(t105);
t88 = pkin(4) * t98 + t96;
t97 = sin(t105);
t89 = pkin(4) * t97 + t130;
t118 = t104 * t125 + t111 * t89 + t88 * t126 + t122;
t99 = qJ(5) + t105;
t95 = cos(t99);
t94 = sin(t99);
t85 = -t107 * t124 + t110 * t116;
t84 = t107 * t123 + t110 * t114;
t83 = t107 * t116 + t110 * t124;
t82 = t107 * t114 - t110 * t123;
t77 = t111 * t94 + t126 * t95;
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - m(6) * t118 - t77 * mrSges(6,1) + mrSges(6,3) * t125 - m(7) * (pkin(5) * t77 + t118) - (-t113 * t125 + t77 * t115) * mrSges(7,1) - (-t77 * t113 - t115 * t125) * mrSges(7,2) + t134 * (-t111 * t95 + t126 * t94) + (-m(3) - m(4) - m(5)) * t122 + (-m(5) * t130 - t106 * mrSges(4,1) - t97 * mrSges(5,1) - t109 * mrSges(4,2) - t98 * mrSges(5,2) - mrSges(3,3)) * t111 + (t133 * t116 + (-m(4) * pkin(2) - m(5) * t96 - t109 * mrSges(4,1) - t98 * mrSges(5,1) + t106 * mrSges(4,2) + t97 * mrSges(5,2) - mrSges(3,1)) * t114) * t108) * g(3) + (-mrSges(1,2) - t107 * mrSges(2,1) - t110 * mrSges(2,2) - m(3) * t120 - t83 * mrSges(3,1) + mrSges(3,3) * t127 - m(4) * (pkin(2) * t83 + t120) - (-t106 * t127 + t109 * t83) * mrSges(4,1) - (-t83 * t106 - t109 * t127) * mrSges(4,2) - m(5) * (t83 * t96 + t100 + (-pkin(7) - t130) * t127) - (-t127 * t97 + t83 * t98) * mrSges(5,1) - (-t127 * t98 - t83 * t97) * mrSges(5,2) + t135 * (t100 + t83 * t88 - t82 * t104 + (-pkin(7) - t89) * t127) + t134 * (t127 * t95 + t83 * t94) + t132 * (-t127 * t94 + t83 * t95) + t131 * t82) * g(2) + (-mrSges(1,1) - t110 * mrSges(2,1) + t107 * mrSges(2,2) - m(3) * t129 - t85 * mrSges(3,1) - mrSges(3,3) * t128 - m(4) * (pkin(2) * t85 + t129) - (t109 * t85 + t121) * mrSges(4,1) - (-t106 * t85 + t109 * t128) * mrSges(4,2) - m(5) * (pkin(3) * t121 + t85 * t96 + t129) - (t128 * t97 + t85 * t98) * mrSges(5,1) - (t128 * t98 - t85 * t97) * mrSges(5,2) + t135 * (-t84 * t104 + t89 * t128 + t85 * t88 + t129) + t134 * (-t128 * t95 + t85 * t94) + t132 * (t128 * t94 + t85 * t95) + t131 * t84) * g(1);
U  = t1;
