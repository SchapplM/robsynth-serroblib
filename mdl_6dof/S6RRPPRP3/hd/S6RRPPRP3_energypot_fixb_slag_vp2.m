% Calculate potential energy for
% S6RRPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5]';
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
% Datum: 2019-03-09 08:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPPRP3_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP3_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP3_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPPRP3_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP3_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP3_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:33:29
% EndTime: 2019-03-09 08:33:29
% DurationCPUTime: 0.47s
% Computational Cost: add. (142->68), mult. (235->63), div. (0->0), fcn. (206->6), ass. (0->31)
t110 = mrSges(5,1) + mrSges(4,3) - mrSges(3,2);
t109 = -m(6) * pkin(8) + m(7) * (-qJ(6) - pkin(8)) - mrSges(3,1) - mrSges(4,1) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t108 = m(7) * pkin(5);
t107 = -m(6) - m(7);
t106 = mrSges(6,1) + mrSges(7,1);
t105 = -mrSges(6,2) - mrSges(7,2);
t104 = -m(5) + t107;
t103 = mrSges(2,2) - mrSges(3,3) - mrSges(4,2) + mrSges(5,3);
t78 = cos(qJ(5));
t67 = pkin(5) * t78 + pkin(4);
t76 = sin(qJ(2));
t79 = cos(qJ(2));
t102 = -mrSges(2,1) + (-m(6) * pkin(4) - m(7) * t67 - t110) * t76 + t109 * t79;
t101 = t76 * pkin(2) + pkin(6);
t75 = sin(qJ(5));
t80 = cos(qJ(1));
t100 = t75 * t80;
t77 = sin(qJ(1));
t99 = t77 * t75;
t98 = t77 * t78;
t97 = t77 * t79;
t96 = t78 * t80;
t95 = t79 * t80;
t94 = t80 * pkin(1) + t77 * pkin(7);
t93 = qJ(3) * t76;
t91 = t77 * pkin(1) - pkin(7) * t80;
t90 = pkin(2) * t95 + t80 * t93 + t94;
t89 = -qJ(3) * t79 + t101;
t85 = pkin(2) * t97 + t77 * t93 + t91;
t69 = t76 * pkin(3);
t1 = (-mrSges(1,3) - mrSges(2,3) - m(4) * t89 - m(5) * (t69 + t89) + t107 * (t69 + t101) + (-m(2) - m(3)) * pkin(6) + (-m(6) * (-pkin(4) - qJ(3)) - m(7) * (-qJ(3) - t67) + t106 * t78 + t105 * t75 + t110) * t79 + t109 * t76) * g(3) + (-t100 * t108 - m(3) * t91 - m(4) * t85 - mrSges(1,2) + t104 * (pkin(3) * t97 + t80 * qJ(4) + t85) - t106 * (t76 * t98 + t100) + t105 * (-t76 * t99 + t96) - t103 * t80 + t102 * t77) * g(2) + (t99 * t108 - m(3) * t94 - m(4) * t90 - mrSges(1,1) + t104 * (pkin(3) * t95 - t77 * qJ(4) + t90) - t106 * (t76 * t96 - t99) + t105 * (-t76 * t100 - t98) + t103 * t77 + t102 * t80) * g(1);
U  = t1;
