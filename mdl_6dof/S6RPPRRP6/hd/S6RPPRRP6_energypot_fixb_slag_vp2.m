% Calculate potential energy for
% S6RPPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
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
% Datum: 2019-03-09 02:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPRRP6_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP6_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP6_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP6_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP6_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP6_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:10:07
% EndTime: 2019-03-09 02:10:08
% DurationCPUTime: 0.34s
% Computational Cost: add. (120->62), mult. (185->59), div. (0->0), fcn. (160->6), ass. (0->27)
t92 = -m(6) - m(7);
t91 = mrSges(6,3) + mrSges(7,2);
t64 = sin(qJ(4));
t67 = cos(qJ(4));
t90 = -t64 * mrSges(5,1) - t67 * mrSges(5,2) - mrSges(2,1) + mrSges(3,2) - mrSges(4,3);
t89 = mrSges(2,2) - mrSges(3,3) - mrSges(4,2) + mrSges(5,3);
t88 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t87 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t86 = pkin(2) + pkin(6);
t65 = sin(qJ(1));
t85 = t65 * t64;
t66 = cos(qJ(5));
t84 = t65 * t66;
t83 = t65 * t67;
t68 = cos(qJ(1));
t82 = t68 * t64;
t81 = t68 * t66;
t80 = t68 * t67;
t79 = t68 * pkin(1) + t65 * qJ(2);
t78 = pkin(3) + t86;
t77 = t68 * qJ(3) + t79;
t76 = t65 * pkin(1) - t68 * qJ(2);
t74 = t65 * qJ(3) + t76;
t72 = -t65 * pkin(7) + t77;
t71 = t68 * pkin(7) + t74;
t63 = sin(qJ(5));
t1 = (-m(4) * t86 - m(5) * t78 - mrSges(3,1) - mrSges(4,1) - mrSges(1,3) - mrSges(2,3) + t92 * (t67 * pkin(4) + t64 * pkin(8) + t78) + (-m(2) - m(3)) * pkin(6) + (t87 * t63 + t88 * t66 - mrSges(5,1)) * t67 + (mrSges(5,2) - t91) * t64) * g(3) + (-m(3) * t76 - m(4) * t74 - m(5) * t71 - mrSges(1,2) + t91 * t83 + t92 * (pkin(4) * t85 - pkin(8) * t83 + t71) + t88 * (t68 * t63 + t64 * t84) + t87 * (t63 * t85 - t81) - t89 * t68 + t90 * t65) * g(2) + (-m(3) * t79 - m(4) * t77 - m(5) * t72 - mrSges(1,1) + t91 * t80 + t92 * (pkin(4) * t82 - pkin(8) * t80 + t72) + t88 * (-t65 * t63 + t64 * t81) + t87 * (t63 * t82 + t84) + t90 * t68 + t89 * t65) * g(1);
U  = t1;
