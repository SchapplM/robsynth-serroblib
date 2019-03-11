% Calculate potential energy for
% S6RPRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 06:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRP3_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP3_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP3_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP3_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP3_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP3_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:02:49
% EndTime: 2019-03-09 06:02:49
% DurationCPUTime: 0.48s
% Computational Cost: add. (250->69), mult. (222->61), div. (0->0), fcn. (201->10), ass. (0->31)
t84 = sin(qJ(4));
t87 = cos(qJ(4));
t116 = -m(5) * pkin(3) - t87 * mrSges(5,1) + t84 * mrSges(5,2) - mrSges(4,1);
t115 = -m(5) * pkin(8) + mrSges(4,2) - mrSges(7,2) - mrSges(5,3) - mrSges(6,3);
t73 = t87 * pkin(4) + pkin(3);
t85 = sin(qJ(3));
t88 = cos(qJ(3));
t90 = -pkin(9) - pkin(8);
t114 = t73 * t88 - t85 * t90;
t113 = m(4) + m(5);
t112 = -m(6) - m(7);
t111 = -t84 * mrSges(5,1) - t87 * mrSges(5,2) + mrSges(3,2) - mrSges(4,3);
t109 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t108 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t107 = t115 * t85 + t116 * t88 - mrSges(3,1);
t106 = pkin(4) * t84;
t86 = sin(qJ(1));
t79 = t86 * pkin(1);
t89 = cos(qJ(1));
t80 = t89 * pkin(1);
t82 = qJ(4) + qJ(5);
t77 = sin(t82);
t104 = t77 * t88;
t78 = cos(t82);
t103 = t78 * t88;
t83 = qJ(2) + pkin(6);
t81 = qJ(1) + pkin(10);
t75 = sin(t81);
t99 = t75 * pkin(2) + t79;
t76 = cos(t81);
t1 = (-m(2) * pkin(6) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + t112 * (t85 * t73 + t88 * t90 + t83) + (-m(3) - t113) * t83 - t115 * t88 + (t108 * t77 + t109 * t78 + t116) * t85) * g(3) + (-m(3) * t79 - t86 * mrSges(2,1) - t89 * mrSges(2,2) - mrSges(1,2) - t113 * t99 + t112 * ((-pkin(7) - t106) * t76 + t99 + t114 * t75) + t109 * (t75 * t103 - t76 * t77) + t108 * (t75 * t104 + t76 * t78) + (t113 * pkin(7) - t111) * t76 + t107 * t75) * g(2) + (-m(3) * t80 - t89 * mrSges(2,1) + t86 * mrSges(2,2) - mrSges(1,1) + (-t113 + t112) * (t76 * pkin(2) + t75 * pkin(7) + t80) + (t112 * t106 - t108 * t78 + t109 * t77 + t111) * t75 + (t109 * t103 + t108 * t104 + t112 * t114 + t107) * t76) * g(1);
U  = t1;
