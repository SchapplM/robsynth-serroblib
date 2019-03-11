% Calculate potential energy for
% S6RPPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2]';
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
% Datum: 2019-03-09 02:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPRRP3_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP3_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP3_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP3_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP3_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP3_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:02:38
% EndTime: 2019-03-09 02:02:39
% DurationCPUTime: 0.42s
% Computational Cost: add. (191->65), mult. (183->57), div. (0->0), fcn. (158->8), ass. (0->27)
t100 = mrSges(5,2) - mrSges(6,3) - mrSges(7,2);
t99 = -m(6) - m(7);
t98 = -mrSges(3,1) + mrSges(4,2) - mrSges(5,3);
t97 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t96 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t74 = sin(qJ(4));
t77 = cos(qJ(4));
t95 = -t74 * mrSges(5,1) - t100 * t77 + mrSges(3,2) - mrSges(4,3);
t94 = pkin(4) * t74;
t93 = pkin(8) * t77;
t75 = sin(qJ(1));
t68 = t75 * pkin(1);
t78 = cos(qJ(1));
t70 = t78 * pkin(1);
t73 = sin(qJ(5));
t92 = t73 * t74;
t76 = cos(qJ(5));
t91 = t74 * t76;
t72 = qJ(2) + pkin(6);
t71 = qJ(1) + pkin(9);
t65 = sin(t71);
t88 = t65 * pkin(2) + t68;
t87 = pkin(3) + t72;
t86 = t65 * pkin(7) + t88;
t66 = cos(t71);
t85 = t66 * pkin(2) + t65 * qJ(3) + t70;
t1 = (-m(2) * pkin(6) - m(5) * t87 - mrSges(4,1) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + t99 * (t77 * pkin(4) + t74 * pkin(8) + t87) + (-m(3) - m(4)) * t72 + (-t73 * t96 + t76 * t97 - mrSges(5,1)) * t77 + t100 * t74) * g(3) + (-m(3) * t68 - m(4) * t88 - m(5) * t86 - t75 * mrSges(2,1) - t78 * mrSges(2,2) - mrSges(1,2) + t99 * (t66 * t93 + t86) + t97 * (t65 * t73 - t66 * t91) + t96 * (t65 * t76 + t66 * t92) + t98 * t65 + (t99 * (-qJ(3) - t94) + (m(4) + m(5)) * qJ(3) - t95) * t66) * g(2) + (-m(3) * t70 - m(4) * t85 - t78 * mrSges(2,1) + t75 * mrSges(2,2) - mrSges(1,1) + (-m(5) + t99) * (t66 * pkin(7) + t85) + (t97 * t73 + t96 * t76 + t98) * t66 + (t99 * (-t93 + t94) + t97 * t91 - t96 * t92 + t95) * t65) * g(1);
U  = t1;
