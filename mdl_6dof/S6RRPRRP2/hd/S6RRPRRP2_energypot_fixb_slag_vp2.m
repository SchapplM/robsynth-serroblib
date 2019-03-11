% Calculate potential energy for
% S6RRPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 11:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRP2_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP2_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP2_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP2_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP2_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP2_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:43:00
% EndTime: 2019-03-09 11:43:00
% DurationCPUTime: 0.36s
% Computational Cost: add. (239->70), mult. (208->67), div. (0->0), fcn. (183->10), ass. (0->32)
t113 = -m(6) - m(7);
t112 = -mrSges(6,3) - mrSges(7,2);
t83 = qJ(2) + pkin(10);
t79 = qJ(4) + t83;
t73 = sin(t79);
t74 = cos(t79);
t89 = cos(qJ(2));
t76 = t89 * pkin(2) + pkin(1);
t77 = sin(t83);
t78 = cos(t83);
t86 = sin(qJ(2));
t111 = -m(3) * pkin(1) - m(4) * t76 - t89 * mrSges(3,1) - t78 * mrSges(4,1) - t74 * mrSges(5,1) + t86 * mrSges(3,2) + t77 * mrSges(4,2) + t73 * mrSges(5,2) - mrSges(2,1);
t84 = -qJ(3) - pkin(7);
t110 = -m(3) * pkin(7) + m(4) * t84 + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t109 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t108 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t107 = t86 * pkin(2) + pkin(6);
t87 = sin(qJ(1));
t106 = t87 * t73;
t85 = sin(qJ(5));
t105 = t87 * t85;
t88 = cos(qJ(5));
t104 = t87 * t88;
t90 = cos(qJ(1));
t103 = t90 * t73;
t102 = t90 * t74;
t65 = pkin(3) * t78 + t76;
t82 = -pkin(8) + t84;
t101 = t87 * t65 + t90 * t82;
t100 = pkin(3) * t77 + t107;
t97 = t90 * t65 - t82 * t87;
t1 = (-m(4) * t107 - m(5) * t100 - mrSges(3,1) * t86 - t77 * mrSges(4,1) - mrSges(3,2) * t89 - t78 * mrSges(4,2) - mrSges(1,3) - mrSges(2,3) + t113 * (t73 * pkin(4) - pkin(9) * t74 + t100) + (-m(2) - m(3)) * pkin(6) + (-mrSges(5,2) - t112) * t74 + (t108 * t85 + t109 * t88 - mrSges(5,1)) * t73) * g(3) + (-m(5) * t101 - mrSges(1,2) + t113 * (t87 * t74 * pkin(4) + pkin(9) * t106 + t101) + t109 * (t74 * t104 - t85 * t90) + t108 * (t74 * t105 + t88 * t90) + t112 * t106 - t110 * t90 + t111 * t87) * g(2) + (-m(5) * t97 - mrSges(1,1) + t113 * (pkin(4) * t102 + pkin(9) * t103 + t97) + t109 * (t88 * t102 + t105) + t108 * (t85 * t102 - t104) + t112 * t103 + t111 * t90 + t110 * t87) * g(1);
U  = t1;
