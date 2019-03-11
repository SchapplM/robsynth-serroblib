% Calculate potential energy for
% S6RRPPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
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
% Datum: 2019-03-09 09:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPPRR10_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR10_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR10_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR10_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR10_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR10_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:34:10
% EndTime: 2019-03-09 09:34:10
% DurationCPUTime: 0.53s
% Computational Cost: add. (186->72), mult. (246->60), div. (0->0), fcn. (221->10), ass. (0->30)
t82 = sin(pkin(10));
t108 = pkin(4) * t82;
t81 = pkin(10) + qJ(5);
t74 = qJ(6) + t81;
t69 = sin(t74);
t70 = cos(t74);
t72 = sin(t81);
t73 = cos(t81);
t83 = cos(pkin(10));
t124 = -t82 * mrSges(5,1) - t83 * mrSges(5,2) - m(6) * t108 - t72 * mrSges(6,1) - t73 * mrSges(6,2) + mrSges(3,2) - mrSges(4,3) - m(7) * (pkin(5) * t72 + t108) - t69 * mrSges(7,1) - t70 * mrSges(7,2);
t84 = -pkin(8) - qJ(4);
t122 = -mrSges(3,1) + mrSges(4,2) + m(7) * (-pkin(9) + t84) - mrSges(7,3) + m(6) * t84 - mrSges(6,3) - m(5) * qJ(4);
t85 = sin(qJ(2));
t87 = cos(qJ(2));
t120 = t122 * t87 + t124 * t85 - mrSges(2,1);
t101 = qJ(3) * t85;
t86 = sin(qJ(1));
t104 = t86 * t87;
t119 = pkin(2) * t104 + t86 * t101;
t118 = -m(5) - m(6) - m(7);
t117 = m(4) - t118;
t109 = t83 * mrSges(5,1) + t73 * mrSges(6,1) + t70 * mrSges(7,1) - t82 * mrSges(5,2) - t72 * mrSges(6,2) - t69 * mrSges(7,2) + mrSges(4,1) - mrSges(2,2) + mrSges(3,3);
t71 = t83 * pkin(4) + pkin(3);
t88 = cos(qJ(1));
t103 = t87 * t88;
t102 = t88 * pkin(1) + t86 * pkin(7);
t78 = t86 * pkin(1);
t98 = -pkin(7) * t88 + t78;
t63 = pkin(5) * t73 + t71;
t1 = (-mrSges(1,3) - mrSges(2,3) + (-m(2) - m(3)) * pkin(6) - t117 * (t85 * pkin(2) + pkin(6)) + (t117 * qJ(3) - t124) * t87 + (-mrSges(5,3) + t122) * t85) * g(3) + (-mrSges(1,2) - m(3) * t98 - m(4) * (t98 + t119) - mrSges(5,3) * t104 + t118 * (t78 + t119) + (-m(5) * (-pkin(3) - pkin(7)) - m(6) * (-pkin(7) - t71) - m(7) * (-pkin(7) - t63) + t109) * t88 + t120 * t86) * g(2) + (-m(3) * t102 - mrSges(5,3) * t103 - mrSges(1,1) - t117 * (pkin(2) * t103 + t102) + (-m(5) * pkin(3) - m(6) * t71 - m(7) * t63 - t109) * t86 + (-t117 * t101 + t120) * t88) * g(1);
U  = t1;
