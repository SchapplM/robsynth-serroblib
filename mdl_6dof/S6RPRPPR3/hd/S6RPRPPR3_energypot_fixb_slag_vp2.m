% Calculate potential energy for
% S6RPRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2]';
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
% Datum: 2019-03-09 02:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPPR3_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR3_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR3_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR3_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR3_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR3_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:44:06
% EndTime: 2019-03-09 02:44:07
% DurationCPUTime: 0.36s
% Computational Cost: add. (192->69), mult. (186->57), div. (0->0), fcn. (153->8), ass. (0->31)
t102 = -m(7) * pkin(8) - mrSges(4,1) - mrSges(5,1) + mrSges(6,2) - mrSges(7,3);
t72 = sin(qJ(6));
t75 = cos(qJ(6));
t101 = t75 * mrSges(7,1) - t72 * mrSges(7,2) + mrSges(6,1) - mrSges(4,2) + mrSges(5,3);
t100 = m(6) + m(7);
t70 = qJ(1) + pkin(9);
t64 = sin(t70);
t73 = sin(qJ(3));
t91 = qJ(4) * t73;
t76 = cos(qJ(3));
t94 = t64 * t76;
t99 = pkin(3) * t94 + t64 * t91;
t65 = cos(t70);
t98 = pkin(4) * t94 + t65 * qJ(5);
t96 = -mrSges(3,1) + t102 * t76 + (-m(7) * pkin(5) - t101) * t73;
t95 = t72 * mrSges(7,1) + t75 * mrSges(7,2) + mrSges(3,2) - mrSges(5,2) - mrSges(4,3) + mrSges(6,3);
t74 = sin(qJ(1));
t68 = t74 * pkin(1);
t77 = cos(qJ(1));
t69 = t77 * pkin(1);
t93 = t65 * t76;
t71 = qJ(2) + pkin(6);
t92 = t64 * pkin(2) + t68;
t90 = t73 * pkin(3) + t71;
t89 = t65 * pkin(2) + t64 * pkin(7) + t69;
t88 = -t65 * pkin(7) + t92;
t87 = pkin(3) * t93 + t65 * t91 + t89;
t81 = -t76 * qJ(4) + t90;
t79 = t88 + t99;
t66 = t73 * pkin(4);
t1 = (-mrSges(1,3) - m(2) * pkin(6) - mrSges(2,3) - mrSges(3,3) - m(5) * t81 - m(6) * (t66 + t81) - m(7) * (t66 + t90) + (-m(3) - m(4)) * t71 + (-m(7) * (-pkin(5) - qJ(4)) + t101) * t76 + t102 * t73) * g(3) + (-mrSges(1,2) - t74 * mrSges(2,1) - t77 * mrSges(2,2) - m(3) * t68 - m(4) * t88 - m(5) * t79 - m(6) * (t79 + t98) - m(7) * (t92 + t98 + t99) + (m(7) * pkin(7) - t95) * t65 + t96 * t64) * g(2) + (-m(3) * t69 - m(4) * t89 - m(5) * t87 - t77 * mrSges(2,1) + t74 * mrSges(2,2) - mrSges(1,1) - t100 * (pkin(4) * t93 + t87) + t96 * t65 + (t100 * qJ(5) + t95) * t64) * g(1);
U  = t1;
