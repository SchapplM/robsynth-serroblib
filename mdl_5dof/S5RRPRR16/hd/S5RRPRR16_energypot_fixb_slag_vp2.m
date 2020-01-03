% Calculate potential energy for
% S5RRPRR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRR16_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR16_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR16_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR16_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR16_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR16_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:44:37
% EndTime: 2019-12-31 20:44:38
% DurationCPUTime: 0.53s
% Computational Cost: add. (175->79), mult. (365->90), div. (0->0), fcn. (405->10), ass. (0->39)
t117 = -m(5) - m(6);
t116 = -mrSges(3,1) + mrSges(4,2);
t115 = mrSges(3,2) - mrSges(4,3);
t114 = mrSges(3,3) + mrSges(4,1);
t113 = m(6) * pkin(9) - mrSges(5,2) + mrSges(6,3);
t86 = sin(qJ(5));
t90 = cos(qJ(5));
t112 = -mrSges(6,1) * t86 - mrSges(6,2) * t90 - mrSges(5,3) + t116;
t111 = -m(6) * pkin(4) - t90 * mrSges(6,1) + t86 * mrSges(6,2) - mrSges(5,1);
t85 = cos(pkin(5));
t110 = t85 * pkin(7) + pkin(6);
t84 = sin(pkin(5));
t88 = sin(qJ(2));
t109 = t84 * t88;
t89 = sin(qJ(1));
t108 = t84 * t89;
t92 = cos(qJ(2));
t107 = t84 * t92;
t106 = t88 * t89;
t93 = cos(qJ(1));
t105 = t88 * t93;
t104 = t89 * t92;
t103 = t92 * t93;
t102 = t93 * t84;
t101 = t93 * pkin(1) + pkin(7) * t108;
t100 = pkin(7) * t102;
t99 = pkin(2) * t109 + t110;
t70 = -t103 * t85 + t106;
t71 = t105 * t85 + t104;
t82 = t89 * pkin(1);
t98 = t71 * pkin(2) + t70 * qJ(3) + t82;
t72 = t104 * t85 + t105;
t73 = -t106 * t85 + t103;
t97 = t73 * pkin(2) + qJ(3) * t72 + t101;
t96 = t85 * pkin(3) + pkin(8) * t109 - qJ(3) * t107 + t99;
t91 = cos(qJ(4));
t87 = sin(qJ(4));
t69 = -t107 * t87 + t85 * t91;
t1 = (-mrSges(1,3) - m(2) * pkin(6) - mrSges(2,3) - m(3) * t110 - m(4) * t99 - m(5) * t96 - t69 * mrSges(5,1) - mrSges(5,3) * t109 - m(6) * (pkin(4) * t69 + t96) - (t109 * t86 + t69 * t90) * mrSges(6,1) - (t109 * t90 - t69 * t86) * mrSges(6,2) - t114 * t85 + ((m(4) * qJ(3) - t115) * t92 + t116 * t88) * t84 - t113 * (t107 * t91 + t85 * t87)) * g(3) + (-mrSges(1,2) - t89 * mrSges(2,1) - t93 * mrSges(2,2) - m(3) * (t82 - t100) - m(4) * (t98 - t100) + t117 * (t71 * pkin(8) + (-pkin(3) - pkin(7)) * t102 + t98) + t115 * t70 + t113 * (t102 * t87 + t70 * t91) + t114 * t102 + t111 * (-t102 * t91 + t70 * t87) + t112 * t71) * g(2) + (-m(3) * t101 - m(4) * t97 - t93 * mrSges(2,1) + t89 * mrSges(2,2) - mrSges(1,1) + t117 * (pkin(3) * t108 + pkin(8) * t73 + t97) + t115 * t72 - t113 * (t108 * t87 - t72 * t91) - t114 * t108 + t111 * (t108 * t91 + t72 * t87) + t112 * t73) * g(1);
U = t1;
