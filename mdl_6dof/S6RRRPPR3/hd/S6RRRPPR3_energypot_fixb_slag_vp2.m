% Calculate potential energy for
% S6RRRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6]';
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
% Datum: 2019-03-09 15:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPPR3_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR3_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR3_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPPR3_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR3_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR3_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:28:18
% EndTime: 2019-03-09 15:28:19
% DurationCPUTime: 0.35s
% Computational Cost: add. (189->67), mult. (200->55), div. (0->0), fcn. (167->8), ass. (0->29)
t72 = sin(qJ(6));
t75 = cos(qJ(6));
t105 = t75 * mrSges(7,1) - t72 * mrSges(7,2) + mrSges(6,1) - mrSges(4,2) + mrSges(5,3);
t104 = -m(7) * pkin(9) - mrSges(4,1) - mrSges(5,1) + mrSges(6,2) - mrSges(7,3);
t103 = -m(6) - m(7);
t77 = cos(qJ(1));
t71 = qJ(2) + qJ(3);
t67 = sin(t71);
t94 = qJ(4) * t67;
t68 = cos(t71);
t96 = t77 * t68;
t102 = pkin(3) * t96 + t77 * t94;
t73 = sin(qJ(2));
t76 = cos(qJ(2));
t100 = -m(3) * pkin(1) - t76 * mrSges(3,1) + t73 * mrSges(3,2) - mrSges(2,1) + t104 * t68 + (-m(7) * pkin(5) - t105) * t67;
t99 = m(3) * pkin(7) - t72 * mrSges(7,1) - t75 * mrSges(7,2) - mrSges(2,2) + mrSges(5,2) + mrSges(3,3) + mrSges(4,3) - mrSges(6,3);
t98 = t73 * pkin(2) + pkin(6);
t74 = sin(qJ(1));
t97 = t74 * t68;
t65 = t76 * pkin(2) + pkin(1);
t78 = -pkin(8) - pkin(7);
t95 = t74 * t65 + t77 * t78;
t92 = t67 * pkin(3) + t98;
t58 = t77 * t65;
t90 = -t74 * t78 + t58;
t88 = pkin(3) * t97 + t74 * t94 + t95;
t82 = -t68 * qJ(4) + t92;
t63 = t67 * pkin(4);
t1 = (-mrSges(1,3) - mrSges(2,3) - t73 * mrSges(3,1) - t76 * mrSges(3,2) - m(4) * t98 - m(5) * t82 - m(6) * (t63 + t82) - m(7) * (t63 + t92) + (-m(2) - m(3)) * pkin(6) + (-m(7) * (-pkin(5) - qJ(4)) + t105) * t68 + t104 * t67) * g(3) + (-m(4) * t95 - m(5) * t88 - mrSges(1,2) + t103 * (pkin(4) * t97 + t77 * qJ(5) + t88) + t99 * t77 + t100 * t74) * g(2) + (-mrSges(1,1) - m(4) * t90 - m(5) * (t90 + t102) + t103 * (pkin(4) * t96 + t102 + t58) + t100 * t77 + (t103 * (-qJ(5) - t78) - t99) * t74) * g(1);
U  = t1;
