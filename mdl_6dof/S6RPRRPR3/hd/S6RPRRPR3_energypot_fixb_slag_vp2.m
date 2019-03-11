% Calculate potential energy for
% S6RPRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
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
% Datum: 2019-03-09 05:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRPR3_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR3_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR3_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR3_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR3_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR3_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:04:31
% EndTime: 2019-03-09 05:04:31
% DurationCPUTime: 0.47s
% Computational Cost: add. (250->74), mult. (269->76), div. (0->0), fcn. (264->10), ass. (0->38)
t111 = mrSges(4,2) - mrSges(6,2) - mrSges(5,3) + mrSges(7,3);
t110 = mrSges(3,2) - mrSges(4,3);
t82 = sin(qJ(4));
t83 = sin(qJ(3));
t86 = cos(qJ(4));
t109 = (pkin(4) * t86 + qJ(5) * t82) * t83;
t87 = cos(qJ(3));
t108 = -t87 * mrSges(4,1) + t111 * t83 - mrSges(3,1);
t81 = sin(qJ(6));
t85 = cos(qJ(6));
t107 = -t81 * mrSges(7,1) - t85 * mrSges(7,2) + mrSges(5,2) - mrSges(6,3);
t106 = -m(7) * pkin(5) - mrSges(7,1) * t85 + mrSges(7,2) * t81 - mrSges(5,1) - mrSges(6,1);
t105 = pkin(3) * t87;
t84 = sin(qJ(1));
t77 = t84 * pkin(1);
t88 = cos(qJ(1));
t78 = t88 * pkin(1);
t79 = qJ(1) + pkin(10);
t74 = sin(t79);
t104 = t74 * t83;
t75 = cos(t79);
t103 = t75 * t83;
t102 = t82 * t87;
t98 = t86 * t87;
t80 = qJ(2) + pkin(6);
t97 = pkin(3) * t83 + t80;
t96 = pkin(2) * t75 + pkin(7) * t74 + t78;
t95 = pkin(2) * t74 - pkin(7) * t75 + t77;
t94 = pkin(8) * t103 + t105 * t75 + t96;
t92 = -pkin(8) * t87 + t97;
t91 = pkin(8) * t104 + t105 * t74 + t95;
t62 = t102 * t75 - t74 * t86;
t63 = t74 * t82 + t75 * t98;
t90 = pkin(4) * t63 + qJ(5) * t62 + t94;
t60 = t102 * t74 + t75 * t86;
t61 = t74 * t98 - t75 * t82;
t89 = pkin(4) * t61 + qJ(5) * t60 + t91;
t1 = (-mrSges(1,3) - m(2) * pkin(6) - mrSges(2,3) - mrSges(3,3) - m(5) * t92 - m(6) * (t92 + t109) - m(7) * (t97 + t109) + (-m(3) - m(4)) * t80 + (-m(7) * (-pkin(8) + pkin(9)) - t111) * t87 + (t106 * t86 + t107 * t82 - mrSges(4,1)) * t83) * g(3) + (-mrSges(1,2) - t84 * mrSges(2,1) - t88 * mrSges(2,2) - m(3) * t77 - m(4) * t95 - m(5) * t91 - m(6) * t89 - m(7) * (-pkin(9) * t104 + t89) - t110 * t75 + t106 * t61 + t107 * t60 + t108 * t74) * g(2) + (-mrSges(1,1) - t88 * mrSges(2,1) + t84 * mrSges(2,2) - m(3) * t78 - m(4) * t96 - m(5) * t94 - m(6) * t90 - m(7) * (-pkin(9) * t103 + t90) + t110 * t74 + t106 * t63 + t107 * t62 + t108 * t75) * g(1);
U  = t1;
