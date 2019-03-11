% Calculate potential energy for
% S6RRPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-03-09 10:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRPR2_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR2_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR2_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR2_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR2_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR2_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:11:58
% EndTime: 2019-03-09 10:11:58
% DurationCPUTime: 0.37s
% Computational Cost: add. (217->77), mult. (182->71), div. (0->0), fcn. (149->10), ass. (0->36)
t106 = -mrSges(5,1) + mrSges(6,2);
t105 = mrSges(5,2) - mrSges(6,3);
t104 = m(6) + m(7);
t81 = cos(qJ(1));
t74 = qJ(2) + pkin(10);
t70 = qJ(4) + t74;
t65 = sin(t70);
t92 = qJ(5) * t65;
t66 = cos(t70);
t96 = t81 * t66;
t103 = pkin(4) * t96 + t81 * t92;
t80 = cos(qJ(2));
t67 = t80 * pkin(2) + pkin(1);
t68 = sin(t74);
t69 = cos(t74);
t77 = sin(qJ(2));
t102 = -m(3) * pkin(1) - m(4) * t67 - t80 * mrSges(3,1) - t69 * mrSges(4,1) + t77 * mrSges(3,2) + t68 * mrSges(4,2) + t105 * t65 + t106 * t66 - mrSges(2,1);
t75 = -qJ(3) - pkin(7);
t101 = -m(3) * pkin(7) + m(4) * t75 - mrSges(6,1) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t100 = t77 * pkin(2) + pkin(6);
t78 = sin(qJ(1));
t99 = t78 * t66;
t76 = sin(qJ(6));
t98 = t78 * t76;
t79 = cos(qJ(6));
t97 = t78 * t79;
t95 = t81 * t76;
t94 = t81 * t79;
t57 = pkin(3) * t69 + t67;
t73 = -pkin(8) + t75;
t93 = t78 * t57 + t81 * t73;
t91 = pkin(3) * t68 + t100;
t56 = t81 * t57;
t88 = -t78 * t73 + t56;
t86 = pkin(4) * t99 + t78 * t92 + t93;
t1 = (-m(4) * t100 - m(5) * t91 - t77 * mrSges(3,1) - t68 * mrSges(4,1) - t80 * mrSges(3,2) - t69 * mrSges(4,2) - mrSges(1,3) - mrSges(2,3) - t104 * (t65 * pkin(4) + t91) + (-m(2) - m(3)) * pkin(6) + (t76 * mrSges(7,1) + t79 * mrSges(7,2) + t104 * qJ(5) - t105) * t66 + (-m(7) * pkin(9) - mrSges(7,3) + t106) * t65) * g(3) + (-mrSges(1,2) - m(5) * t93 - m(6) * t86 - m(7) * (pkin(9) * t99 + t86) - (t65 * t98 - t94) * mrSges(7,1) - (t65 * t97 + t95) * mrSges(7,2) - mrSges(7,3) * t99 + (m(7) * pkin(5) - t101) * t81 + t102 * t78) * g(2) + (-mrSges(1,1) - m(5) * t88 - m(6) * (t88 + t103) - m(7) * (pkin(9) * t96 + t103 + t56) - (t65 * t95 + t97) * mrSges(7,1) - (t65 * t94 - t98) * mrSges(7,2) - mrSges(7,3) * t96 + t102 * t81 + (-m(7) * (pkin(5) - t73) + t101) * t78) * g(1);
U  = t1;
