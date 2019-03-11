% Calculate potential energy for
% S6RPRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5]';
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
% Datum: 2019-03-09 03:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPRP10_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP10_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP10_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPRP10_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP10_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP10_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:30:52
% EndTime: 2019-03-09 03:30:53
% DurationCPUTime: 0.44s
% Computational Cost: add. (130->71), mult. (211->65), div. (0->0), fcn. (186->6), ass. (0->26)
t99 = -mrSges(4,1) + mrSges(5,2);
t98 = -m(6) - m(7);
t97 = -mrSges(6,3) - mrSges(7,2);
t96 = -mrSges(2,1) + mrSges(3,2) - mrSges(4,3) - mrSges(5,1);
t71 = sin(qJ(3));
t74 = cos(qJ(3));
t95 = -t74 * mrSges(4,2) + t99 * t71 + mrSges(2,2) - mrSges(3,3);
t94 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t93 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t92 = pkin(2) + pkin(6);
t72 = sin(qJ(1));
t89 = t72 * t71;
t88 = t72 * t74;
t75 = cos(qJ(1));
t87 = t75 * t74;
t64 = t72 * pkin(1);
t86 = t72 * pkin(7) + t64;
t85 = t75 * pkin(1) + t72 * qJ(2);
t84 = qJ(4) * t74;
t83 = t75 * t84 + t86;
t82 = t75 * pkin(7) + t85;
t81 = t74 * pkin(3) + t71 * qJ(4) + t92;
t79 = pkin(3) * t89 + t82;
t73 = cos(qJ(5));
t70 = sin(qJ(5));
t1 = (-m(4) * t92 - m(5) * t81 - mrSges(3,1) - mrSges(1,3) - mrSges(2,3) + t98 * (t74 * pkin(8) + t81) + (-m(2) - m(3)) * pkin(6) + (t97 + t99) * t74 + (t94 * t70 + t93 * t73 + mrSges(4,2) - mrSges(5,3)) * t71) * g(3) + (-m(3) * t64 - m(4) * t86 - m(5) * t83 - mrSges(1,2) + t98 * (t72 * pkin(4) + t83) + t94 * (t70 * t87 + t72 * t73) - t93 * (t72 * t70 - t73 * t87) + t96 * t72 + (-t74 * mrSges(5,3) + (m(5) * pkin(3) + t98 * (-pkin(3) - pkin(8)) - t97) * t71 + (m(3) + m(4) + m(5) - t98) * qJ(2) - t95) * t75) * g(2) + (-m(3) * t85 - m(4) * t82 - m(5) * t79 - mrSges(1,1) + t97 * t89 + t98 * (t75 * pkin(4) + pkin(8) * t89 - t72 * t84 + t79) + t94 * (-t70 * t88 + t75 * t73) - t93 * (t75 * t70 + t73 * t88) + t96 * t75 + (-(-m(5) * qJ(4) - mrSges(5,3)) * t74 + t95) * t72) * g(1);
U  = t1;
