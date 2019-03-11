% Calculate potential energy for
% S6RPRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 06:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRP7_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP7_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP7_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP7_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP7_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP7_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:18:33
% EndTime: 2019-03-09 06:18:33
% DurationCPUTime: 0.50s
% Computational Cost: add. (239->68), mult. (236->60), div. (0->0), fcn. (215->10), ass. (0->31)
t133 = m(5) * pkin(8) - mrSges(4,2) + mrSges(7,2) + mrSges(5,3) + mrSges(6,3);
t94 = sin(qJ(4));
t96 = cos(qJ(4));
t132 = -m(5) * pkin(3) - t96 * mrSges(5,1) + t94 * mrSges(5,2) - mrSges(4,1);
t127 = -m(6) - m(7);
t130 = m(3) * qJ(2) + t96 * mrSges(5,2) - mrSges(2,2) + mrSges(3,3) + mrSges(4,3) + (-pkin(4) * t127 + mrSges(5,1)) * t94;
t83 = pkin(4) * t96 + pkin(3);
t89 = pkin(10) + qJ(3);
t84 = sin(t89);
t85 = cos(t89);
t91 = sin(pkin(10));
t92 = cos(pkin(10));
t98 = -pkin(9) - pkin(8);
t129 = -m(3) * pkin(1) - t92 * mrSges(3,1) + t91 * mrSges(3,2) - mrSges(2,1) + (t127 * t83 + t132) * t85 + (-t127 * t98 - t133) * t84;
t128 = -m(4) - m(5);
t121 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t120 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t118 = t91 * pkin(2) + pkin(6);
t90 = qJ(4) + qJ(5);
t86 = sin(t90);
t95 = sin(qJ(1));
t113 = t86 * t95;
t97 = cos(qJ(1));
t112 = t86 * t97;
t87 = cos(t90);
t111 = t87 * t97;
t110 = t95 * t87;
t93 = -pkin(7) - qJ(2);
t81 = pkin(2) * t92 + pkin(1);
t77 = t97 * t81;
t1 = (-t91 * mrSges(3,1) - t92 * mrSges(3,2) - mrSges(1,3) - mrSges(2,3) + t128 * t118 + t127 * (t84 * t83 + t85 * t98 + t118) + (-m(2) - m(3)) * pkin(6) + t133 * t85 + (t120 * t86 + t121 * t87 + t132) * t84) * g(3) + (-mrSges(1,2) + t121 * (t85 * t110 - t112) + t120 * (t85 * t113 + t111) + (t127 + t128) * (t95 * t81 + t97 * t93)) * g(2) + (-m(5) * t77 - mrSges(1,1) + t121 * (t85 * t111 + t113) + t120 * (t85 * t112 - t110) + (-m(4) + t127) * (-t95 * t93 + t77)) * g(1) + (t129 * g(1) + t130 * g(2)) * t97 + (t129 * g(2) + (m(5) * t93 - t130) * g(1)) * t95;
U  = t1;
