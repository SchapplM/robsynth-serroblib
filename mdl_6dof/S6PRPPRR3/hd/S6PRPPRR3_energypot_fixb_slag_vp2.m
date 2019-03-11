% Calculate potential energy for
% S6PRPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 19:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPPRR3_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR3_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_energypot_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR3_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPPRR3_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:21:21
% EndTime: 2019-03-08 19:21:22
% DurationCPUTime: 0.66s
% Computational Cost: add. (312->91), mult. (686->114), div. (0->0), fcn. (818->12), ass. (0->46)
t147 = -m(6) - m(7);
t146 = -mrSges(3,1) - mrSges(4,1);
t145 = -mrSges(4,3) + mrSges(3,2);
t144 = mrSges(3,3) + mrSges(4,2) - mrSges(5,3);
t143 = -m(7) * pkin(9) + mrSges(6,2) - mrSges(7,3);
t116 = sin(qJ(6));
t119 = cos(qJ(6));
t142 = -t116 * mrSges(7,1) - t119 * mrSges(7,2) + mrSges(5,2) - mrSges(6,3);
t141 = -m(7) * pkin(5) - t119 * mrSges(7,1) + t116 * mrSges(7,2) - mrSges(6,1);
t111 = sin(pkin(10));
t112 = sin(pkin(6));
t140 = t111 * t112;
t114 = cos(pkin(10));
t139 = t112 * t114;
t117 = sin(qJ(5));
t138 = t112 * t117;
t118 = sin(qJ(2));
t137 = t112 * t118;
t120 = cos(qJ(5));
t136 = t112 * t120;
t115 = cos(pkin(6));
t135 = t115 * t118;
t121 = cos(qJ(2));
t134 = t115 * t121;
t133 = t114 * pkin(1) + pkin(7) * t140;
t132 = t115 * pkin(7) + qJ(1);
t131 = pkin(2) * t137 + t132;
t130 = t111 * pkin(1) - pkin(7) * t139;
t100 = -t111 * t135 + t114 * t121;
t99 = t111 * t134 + t114 * t118;
t129 = t100 * pkin(2) + qJ(3) * t99 + t133;
t97 = t111 * t118 - t114 * t134;
t98 = t111 * t121 + t114 * t135;
t128 = t98 * pkin(2) + qJ(3) * t97 + t130;
t127 = t98 * pkin(3) + qJ(4) * t139 + t128;
t126 = t100 * pkin(3) - qJ(4) * t140 + t129;
t125 = -qJ(3) * t112 * t121 + pkin(3) * t137 - t115 * qJ(4) + t131;
t113 = cos(pkin(11));
t110 = sin(pkin(11));
t92 = (-t110 * t121 + t113 * t118) * t112;
t91 = (t110 * t118 + t113 * t121) * t112;
t85 = t100 * t113 + t110 * t99;
t84 = t100 * t110 - t99 * t113;
t83 = t110 * t97 + t113 * t98;
t82 = t110 * t98 - t97 * t113;
t1 = (-m(2) * qJ(1) - m(3) * t132 - m(4) * t131 - m(5) * t125 - t92 * mrSges(5,1) - mrSges(1,3) - mrSges(2,3) + t147 * (t92 * pkin(4) + t91 * pkin(8) + t125) + t141 * (-t115 * t117 + t120 * t92) + t142 * t91 + t143 * (t115 * t120 + t117 * t92) + ((m(4) * qJ(3) - t145) * t121 + t146 * t118) * t112 - t144 * t115) * g(3) + (-m(3) * t130 - m(4) * t128 - m(5) * t127 - mrSges(2,1) * t111 - t83 * mrSges(5,1) - mrSges(2,2) * t114 - mrSges(1,2) + t146 * t98 + t145 * t97 + t147 * (t83 * pkin(4) + pkin(8) * t82 + t127) + t141 * (t114 * t138 + t120 * t83) + t142 * t82 + t143 * (-t114 * t136 + t117 * t83) + t144 * t139) * g(2) + (-m(3) * t133 - m(4) * t129 - m(5) * t126 - mrSges(2,1) * t114 - t85 * mrSges(5,1) + mrSges(2,2) * t111 - mrSges(1,1) + t145 * t99 + t147 * (t85 * pkin(4) + pkin(8) * t84 + t126) + t141 * (-t111 * t138 + t120 * t85) + t142 * t84 + t143 * (t111 * t136 + t117 * t85) + t146 * t100 - t144 * t140) * g(1);
U  = t1;
