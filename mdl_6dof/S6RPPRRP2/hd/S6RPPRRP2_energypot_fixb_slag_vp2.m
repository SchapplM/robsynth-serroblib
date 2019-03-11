% Calculate potential energy for
% S6RPPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-03-09 02:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPRRP2_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP2_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP2_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP2_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP2_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP2_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:00:11
% EndTime: 2019-03-09 02:00:11
% DurationCPUTime: 0.36s
% Computational Cost: add. (238->70), mult. (194->66), div. (0->0), fcn. (169->10), ass. (0->35)
t110 = -m(3) - m(4);
t109 = -m(6) - m(7);
t108 = -mrSges(6,3) - mrSges(7,2);
t79 = pkin(10) + qJ(4);
t72 = sin(t79);
t74 = cos(t79);
t81 = sin(pkin(10));
t82 = cos(pkin(10));
t107 = -m(4) * pkin(2) - t82 * mrSges(4,1) - t74 * mrSges(5,1) + t81 * mrSges(4,2) + t72 * mrSges(5,2) - mrSges(3,1);
t106 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t105 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t104 = m(4) * qJ(3) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3);
t103 = pkin(4) * t74;
t86 = sin(qJ(1));
t77 = t86 * pkin(1);
t88 = cos(qJ(1));
t78 = t88 * pkin(1);
t80 = qJ(1) + pkin(9);
t75 = cos(t80);
t102 = t72 * t75;
t73 = sin(t80);
t101 = t73 * t72;
t85 = sin(qJ(5));
t100 = t73 * t85;
t87 = cos(qJ(5));
t99 = t73 * t87;
t98 = t75 * t85;
t97 = t75 * t87;
t83 = qJ(2) + pkin(6);
t71 = pkin(3) * t82 + pkin(2);
t84 = -pkin(7) - qJ(3);
t96 = t73 * t71 + t75 * t84 + t77;
t95 = t81 * pkin(3) + t83;
t94 = t75 * t71 - t73 * t84 + t78;
t1 = (-m(2) * pkin(6) - m(5) * t95 - t81 * mrSges(4,1) - t82 * mrSges(4,2) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + t109 * (t72 * pkin(4) - pkin(8) * t74 + t95) + t110 * t83 + (-mrSges(5,2) - t108) * t74 + (t105 * t85 + t106 * t87 - mrSges(5,1)) * t72) * g(3) + (-m(5) * t96 - t86 * mrSges(2,1) - t88 * mrSges(2,2) - mrSges(1,2) + t109 * (pkin(8) * t101 + t73 * t103 + t96) + t110 * t77 + t106 * (t74 * t99 - t98) + t105 * (t74 * t100 + t97) + t108 * t101 + t104 * t75 + t107 * t73) * g(2) + (-m(5) * t94 - t88 * mrSges(2,1) + t86 * mrSges(2,2) - mrSges(1,1) + t109 * (pkin(8) * t102 + t75 * t103 + t94) + t110 * t78 + t106 * (t74 * t97 + t100) + t105 * (t74 * t98 - t99) + t108 * t102 + t107 * t75 - t104 * t73) * g(1);
U  = t1;
