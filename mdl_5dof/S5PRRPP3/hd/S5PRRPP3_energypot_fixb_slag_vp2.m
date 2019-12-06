% Calculate potential energy for
% S5PRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:14
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRPP3_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP3_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP3_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP3_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP3_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPP3_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:11:55
% EndTime: 2019-12-05 16:11:55
% DurationCPUTime: 0.45s
% Computational Cost: add. (145->67), mult. (286->78), div. (0->0), fcn. (301->8), ass. (0->33)
t118 = mrSges(3,2) - mrSges(4,3);
t117 = -m(5) - m(6);
t116 = mrSges(2,2) - mrSges(3,3);
t115 = -mrSges(5,3) - mrSges(6,2);
t114 = mrSges(4,2) + t115;
t91 = sin(qJ(2));
t93 = cos(qJ(2));
t113 = -t93 * mrSges(3,1) + t118 * t91 - mrSges(2,1);
t112 = -m(6) * pkin(4) - mrSges(5,1) - mrSges(6,1);
t111 = -m(6) * qJ(5) + mrSges(5,2) - mrSges(6,3);
t110 = pkin(2) * t93;
t87 = sin(pkin(7));
t109 = t87 * t91;
t89 = cos(pkin(7));
t108 = t89 * t91;
t90 = sin(qJ(3));
t107 = t90 * t91;
t106 = t90 * t93;
t92 = cos(qJ(3));
t104 = t91 * t92;
t103 = t92 * t93;
t102 = t89 * pkin(1) + t87 * pkin(5);
t101 = t87 * pkin(1) - t89 * pkin(5);
t100 = pkin(6) * t108 + t89 * t110 + t102;
t99 = t91 * pkin(2) - t93 * pkin(6) + qJ(1);
t97 = pkin(6) * t109 + t87 * t110 + t101;
t88 = cos(pkin(8));
t86 = sin(pkin(8));
t72 = t89 * t103 + t87 * t90;
t71 = t89 * t106 - t87 * t92;
t68 = t87 * t103 - t89 * t90;
t67 = t87 * t106 + t89 * t92;
t1 = (-m(4) * t99 - mrSges(1,3) - mrSges(2,3) + t117 * (pkin(3) * t104 + qJ(4) * t107 + t99) - t118 * t93 + (-t92 * mrSges(4,1) + t90 * mrSges(4,2) - mrSges(3,1)) * t91 + t112 * (t88 * t104 - t93 * t86) + t111 * (t86 * t104 + t93 * t88) + t115 * t107 + (-m(2) - m(3)) * qJ(1)) * g(3) + (-m(3) * t101 - m(4) * t97 - t68 * mrSges(4,1) - mrSges(1,2) + t117 * (t68 * pkin(3) + t67 * qJ(4) + t97) - t116 * t89 + t112 * (t86 * t109 + t68 * t88) + t111 * (-t88 * t109 + t68 * t86) + t113 * t87 + t114 * t67) * g(2) + (-m(3) * t102 - m(4) * t100 - t72 * mrSges(4,1) - mrSges(1,1) + t117 * (t72 * pkin(3) + t71 * qJ(4) + t100) + t116 * t87 + t112 * (t86 * t108 + t72 * t88) + t111 * (-t88 * t108 + t72 * t86) + t113 * t89 + t114 * t71) * g(1);
U = t1;
