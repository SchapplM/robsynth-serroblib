% Calculate potential energy for
% S6RPRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
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
% Datum: 2019-03-09 04:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRPP3_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP3_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP3_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP3_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP3_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP3_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:35:07
% EndTime: 2019-03-09 04:35:08
% DurationCPUTime: 0.43s
% Computational Cost: add. (232->71), mult. (243->72), div. (0->0), fcn. (228->8), ass. (0->36)
t109 = -mrSges(6,1) - mrSges(7,1) + mrSges(4,2) - mrSges(5,3);
t108 = mrSges(3,2) - mrSges(4,3);
t81 = sin(qJ(4));
t82 = sin(qJ(3));
t84 = cos(qJ(4));
t107 = (pkin(4) * t84 + qJ(5) * t81) * t82;
t106 = mrSges(5,2) - mrSges(7,2) - mrSges(6,3);
t85 = cos(qJ(3));
t105 = -t85 * mrSges(4,1) + t109 * t82 - mrSges(3,1);
t104 = -m(7) * qJ(6) - mrSges(5,1) + mrSges(6,2) - mrSges(7,3);
t103 = pkin(3) * t85;
t83 = sin(qJ(1));
t77 = t83 * pkin(1);
t86 = cos(qJ(1));
t78 = t86 * pkin(1);
t79 = qJ(1) + pkin(9);
t74 = sin(t79);
t102 = t74 * t82;
t75 = cos(t79);
t101 = t75 * t82;
t100 = t81 * t85;
t96 = t84 * t85;
t80 = qJ(2) + pkin(6);
t95 = pkin(3) * t82 + t80;
t94 = pkin(2) * t75 + pkin(7) * t74 + t78;
t93 = pkin(2) * t74 - pkin(7) * t75 + t77;
t92 = pkin(8) * t101 + t103 * t75 + t94;
t90 = -pkin(8) * t85 + t95;
t89 = pkin(8) * t102 + t103 * t74 + t93;
t61 = t100 * t75 - t74 * t84;
t62 = t74 * t81 + t75 * t96;
t88 = pkin(4) * t62 + qJ(5) * t61 + t92;
t59 = t100 * t74 + t75 * t84;
t60 = t74 * t96 - t75 * t81;
t87 = pkin(4) * t60 + qJ(5) * t59 + t89;
t1 = (-mrSges(1,3) - m(2) * pkin(6) - mrSges(2,3) - mrSges(3,3) - m(5) * t90 - m(6) * (t90 + t107) - m(7) * (t95 + t107) + (-m(3) - m(4)) * t80 + (-m(7) * (-pkin(5) - pkin(8)) - t109) * t85 + (t104 * t84 + t106 * t81 - mrSges(4,1)) * t82) * g(3) + (-mrSges(1,2) - t83 * mrSges(2,1) - t86 * mrSges(2,2) - m(3) * t77 - m(4) * t93 - m(5) * t89 - m(6) * t87 - m(7) * (pkin(5) * t102 + t87) - t108 * t75 + t104 * t60 + t106 * t59 + t105 * t74) * g(2) + (-mrSges(1,1) - t86 * mrSges(2,1) + t83 * mrSges(2,2) - m(3) * t78 - m(4) * t94 - m(5) * t92 - m(6) * t88 - m(7) * (pkin(5) * t101 + t88) + t108 * t74 + t104 * t62 + t106 * t61 + t105 * t75) * g(1);
U  = t1;
