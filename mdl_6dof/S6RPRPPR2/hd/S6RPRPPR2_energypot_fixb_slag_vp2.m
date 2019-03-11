% Calculate potential energy for
% S6RPRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
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
% Datum: 2019-03-09 02:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPPR2_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR2_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR2_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR2_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR2_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR2_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:41:10
% EndTime: 2019-03-09 02:41:10
% DurationCPUTime: 0.41s
% Computational Cost: add. (214->77), mult. (168->71), div. (0->0), fcn. (135->10), ass. (0->37)
t102 = -mrSges(5,1) + mrSges(6,2);
t101 = mrSges(5,2) - mrSges(6,3);
t100 = -m(3) - m(4);
t99 = m(6) + m(7);
t71 = qJ(1) + pkin(9);
t66 = cos(t71);
t70 = qJ(3) + pkin(10);
t63 = sin(t70);
t88 = qJ(5) * t63;
t65 = cos(t70);
t92 = t65 * t66;
t98 = pkin(4) * t92 + t66 * t88;
t75 = sin(qJ(3));
t78 = cos(qJ(3));
t97 = -m(4) * pkin(2) - t78 * mrSges(4,1) + t75 * mrSges(4,2) + t101 * t63 + t102 * t65 - mrSges(3,1);
t96 = -m(4) * pkin(7) - mrSges(6,1) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t76 = sin(qJ(1));
t68 = t76 * pkin(1);
t79 = cos(qJ(1));
t69 = t79 * pkin(1);
t64 = sin(t71);
t95 = t64 * t65;
t74 = sin(qJ(6));
t94 = t64 * t74;
t77 = cos(qJ(6));
t93 = t64 * t77;
t91 = t66 * t74;
t90 = t66 * t77;
t73 = qJ(2) + pkin(6);
t62 = pkin(3) * t78 + pkin(2);
t89 = t66 * t62 + t69;
t72 = -qJ(4) - pkin(7);
t87 = t64 * t62 + t66 * t72 + t68;
t86 = t75 * pkin(3) + t73;
t84 = -t64 * t72 + t89;
t83 = pkin(4) * t95 + t64 * t88 + t87;
t1 = (-m(2) * pkin(6) - m(5) * t86 - t75 * mrSges(4,1) - t78 * mrSges(4,2) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) - t99 * (t63 * pkin(4) + t86) + t100 * t73 + (t74 * mrSges(7,1) + t77 * mrSges(7,2) + t99 * qJ(5) - t101) * t65 + (-m(7) * pkin(8) - mrSges(7,3) + t102) * t63) * g(3) + (-mrSges(1,2) - t76 * mrSges(2,1) - t79 * mrSges(2,2) - m(5) * t87 - m(6) * t83 - m(7) * (pkin(8) * t95 + t83) - (t63 * t94 - t90) * mrSges(7,1) - (t63 * t93 + t91) * mrSges(7,2) - mrSges(7,3) * t95 + t100 * t68 + (m(7) * pkin(5) - t96) * t66 + t97 * t64) * g(2) + (-mrSges(1,1) - t79 * mrSges(2,1) + t76 * mrSges(2,2) - m(5) * t84 - m(6) * (t84 + t98) - m(7) * (pkin(8) * t92 + t89 + t98) - (t63 * t91 + t93) * mrSges(7,1) - (t63 * t90 - t94) * mrSges(7,2) - mrSges(7,3) * t92 + t100 * t69 + t97 * t66 + (-m(7) * (pkin(5) - t72) + t96) * t64) * g(1);
U  = t1;
