% Calculate potential energy for
% S6RRRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 01:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRRP5_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP5_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP5_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP5_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP5_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP5_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:18:02
% EndTime: 2019-03-10 01:18:03
% DurationCPUTime: 0.46s
% Computational Cost: add. (235->70), mult. (249->61), div. (0->0), fcn. (228->10), ass. (0->32)
t86 = cos(qJ(3));
t71 = t86 * pkin(3) + pkin(2);
t82 = qJ(3) + qJ(4);
t73 = cos(t82);
t65 = pkin(4) * t73 + t71;
t75 = qJ(5) + t82;
t70 = cos(t75);
t72 = sin(t82);
t83 = sin(qJ(3));
t114 = -m(6) * t65 - m(7) * (pkin(5) * t70 + t65) - mrSges(3,1) - m(5) * t71 - t73 * mrSges(5,1) + t72 * mrSges(5,2) - m(4) * pkin(2) - t86 * mrSges(4,1) + t83 * mrSges(4,2);
t89 = -pkin(9) - pkin(8);
t81 = -pkin(10) + t89;
t113 = m(6) * t81 + m(7) * (-qJ(6) + t81) + mrSges(3,2) - mrSges(6,3) - mrSges(7,3) + m(5) * t89 - mrSges(5,3) - m(4) * pkin(8) - mrSges(4,3);
t112 = -m(4) - m(5);
t111 = -mrSges(6,1) - mrSges(7,1);
t110 = mrSges(6,2) + mrSges(7,2);
t109 = -m(3) - m(6) - m(7);
t106 = t109 + t112;
t84 = sin(qJ(2));
t87 = cos(qJ(2));
t105 = t113 * t84 + t114 * t87 - mrSges(2,1);
t76 = t83 * pkin(3);
t66 = pkin(4) * t72 + t76;
t69 = sin(t75);
t104 = m(6) * t66 + m(7) * (pkin(5) * t69 + t66) - mrSges(2,2) + mrSges(3,3) + t72 * mrSges(5,1) + t73 * mrSges(5,2) + t83 * mrSges(4,1) + t86 * mrSges(4,2);
t85 = sin(qJ(1));
t103 = t85 * t87;
t88 = cos(qJ(1));
t102 = t88 * t69;
t101 = t88 * t70;
t78 = t85 * pkin(1);
t1 = (-mrSges(1,3) - mrSges(2,3) + (-m(2) + t106) * pkin(6) - t113 * t87 + (t110 * t69 + t111 * t70 + t114) * t84) * g(3) + (-mrSges(1,2) + t112 * t78 + t111 * (t70 * t103 - t102) - t110 * (-t69 * t103 - t101) + t109 * (-t88 * pkin(7) + t78) + (m(4) * pkin(7) - m(5) * (-pkin(7) - t76) + t104) * t88 + t105 * t85) * g(2) + (-mrSges(1,1) + t111 * (t87 * t101 + t85 * t69) - t110 * (-t87 * t102 + t85 * t70) + t106 * (t88 * pkin(1) + t85 * pkin(7)) + (-m(5) * t76 - t104) * t85 + t105 * t88) * g(1);
U  = t1;
