% Calculate potential energy for
% S6RPRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
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
% Datum: 2019-03-09 06:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRP8_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP8_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP8_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP8_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP8_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP8_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:22:46
% EndTime: 2019-03-09 06:22:46
% DurationCPUTime: 0.39s
% Computational Cost: add. (171->70), mult. (202->64), div. (0->0), fcn. (177->8), ass. (0->31)
t102 = m(3) + m(4);
t101 = -m(6) - m(7);
t100 = mrSges(6,3) + mrSges(7,2);
t70 = qJ(3) + qJ(4);
t64 = sin(t70);
t65 = cos(t70);
t72 = sin(qJ(3));
t75 = cos(qJ(3));
t99 = -t72 * mrSges(4,1) - t64 * mrSges(5,1) - t75 * mrSges(4,2) - t65 * mrSges(5,2) + mrSges(2,2) - mrSges(3,3);
t98 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t97 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t96 = -m(4) * pkin(7) - mrSges(2,1) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t95 = pkin(2) + pkin(6);
t94 = pkin(3) * t72;
t71 = sin(qJ(5));
t76 = cos(qJ(1));
t93 = t71 * t76;
t73 = sin(qJ(1));
t92 = t73 * t64;
t91 = t73 * t65;
t74 = cos(qJ(5));
t90 = t73 * t74;
t89 = t76 * t74;
t88 = t76 * pkin(1) + t73 * qJ(2);
t87 = t75 * pkin(3) + t95;
t86 = -qJ(2) - t94;
t67 = t73 * pkin(1);
t77 = -pkin(8) - pkin(7);
t85 = -t73 * t77 + t67;
t80 = t73 * t94 - t76 * t77 + t88;
t1 = (-m(4) * t95 - m(5) * t87 - t75 * mrSges(4,1) + t72 * mrSges(4,2) - mrSges(3,1) - mrSges(1,3) - mrSges(2,3) + t101 * (t65 * pkin(4) + t64 * pkin(9) + t87) + (-m(2) - m(3)) * pkin(6) + (-t97 * t71 + t98 * t74 - mrSges(5,1)) * t65 + (mrSges(5,2) - t100) * t64) * g(3) + (-m(5) * t85 - mrSges(1,2) + t101 * (t76 * t65 * pkin(9) + t85) - t102 * t67 + t98 * (-t64 * t89 + t71 * t73) + t97 * (t64 * t93 + t90) + t96 * t73 + (-m(5) * t86 + t101 * (-pkin(4) * t64 + t86) - t100 * t65 + t102 * qJ(2) - t99) * t76) * g(2) + (-m(5) * t80 - mrSges(1,1) + t100 * t91 - t102 * t88 + t101 * (pkin(4) * t92 - pkin(9) * t91 + t80) + t98 * (t64 * t90 + t93) - t97 * (t71 * t92 - t89) + t96 * t76 + t99 * t73) * g(1);
U  = t1;
