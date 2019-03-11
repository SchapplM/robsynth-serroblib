% Calculate potential energy for
% S6RPPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
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
% Datum: 2019-03-09 02:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPRRP8_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP8_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP8_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP8_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP8_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP8_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:15:16
% EndTime: 2019-03-09 02:15:17
% DurationCPUTime: 0.38s
% Computational Cost: add. (171->70), mult. (202->63), div. (0->0), fcn. (177->8), ass. (0->32)
t103 = m(3) + m(4);
t102 = -m(6) - m(7);
t101 = mrSges(6,3) + mrSges(7,2);
t70 = pkin(9) + qJ(4);
t64 = sin(t70);
t65 = cos(t70);
t71 = sin(pkin(9));
t72 = cos(pkin(9));
t100 = -t71 * mrSges(4,1) - t64 * mrSges(5,1) - t72 * mrSges(4,2) - t65 * mrSges(5,2) + mrSges(2,2) - mrSges(3,3);
t99 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t98 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t97 = -m(4) * qJ(3) - mrSges(2,1) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t96 = pkin(2) + pkin(6);
t95 = pkin(3) * t71;
t94 = pkin(4) * t64;
t74 = sin(qJ(5));
t77 = cos(qJ(1));
t93 = t74 * t77;
t75 = sin(qJ(1));
t92 = t75 * t65;
t91 = t75 * t74;
t76 = cos(qJ(5));
t90 = t75 * t76;
t89 = t77 * t76;
t88 = t77 * pkin(1) + t75 * qJ(2);
t87 = t72 * pkin(3) + t96;
t86 = -qJ(2) - t95;
t68 = t75 * pkin(1);
t73 = -pkin(7) - qJ(3);
t85 = -t75 * t73 + t68;
t80 = -t77 * t73 + t75 * t95 + t88;
t1 = (-m(4) * t96 - m(5) * t87 - t72 * mrSges(4,1) + t71 * mrSges(4,2) - mrSges(3,1) - mrSges(1,3) - mrSges(2,3) + t102 * (t65 * pkin(4) + t64 * pkin(8) + t87) + (-m(2) - m(3)) * pkin(6) + (-t98 * t74 + t99 * t76 - mrSges(5,1)) * t65 + (mrSges(5,2) - t101) * t64) * g(3) + (-m(5) * t85 - mrSges(1,2) + t102 * (t77 * t65 * pkin(8) + t85) - t103 * t68 + t99 * (-t64 * t89 + t91) + t98 * (t64 * t93 + t90) + t97 * t75 + (-m(5) * t86 + t102 * (t86 - t94) - t101 * t65 + t103 * qJ(2) - t100) * t77) * g(2) + (-m(5) * t80 - mrSges(1,1) + t101 * t92 - t103 * t88 + t102 * (-pkin(8) * t92 + t75 * t94 + t80) + t99 * (t64 * t90 + t93) - t98 * (t64 * t91 - t89) + t97 * t77 + t100 * t75) * g(1);
U  = t1;
