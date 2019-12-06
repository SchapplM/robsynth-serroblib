% Calculate potential energy for
% S5PRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1]';
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
% Datum: 2019-12-05 16:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRPRR8_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR8_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR8_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR8_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR8_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR8_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:02:20
% EndTime: 2019-12-05 16:02:21
% DurationCPUTime: 0.53s
% Computational Cost: add. (175->79), mult. (365->92), div. (0->0), fcn. (405->10), ass. (0->37)
t115 = -m(5) - m(6);
t114 = -mrSges(3,1) + mrSges(4,2);
t113 = mrSges(3,2) - mrSges(4,3);
t112 = mrSges(3,3) + mrSges(4,1);
t111 = m(6) * pkin(8) - mrSges(5,2) + mrSges(6,3);
t88 = sin(qJ(5));
t91 = cos(qJ(5));
t110 = -t88 * mrSges(6,1) - t91 * mrSges(6,2) - mrSges(5,3) + t114;
t109 = -m(6) * pkin(4) - t91 * mrSges(6,1) + t88 * mrSges(6,2) - mrSges(5,1);
t84 = sin(pkin(9));
t85 = sin(pkin(5));
t108 = t84 * t85;
t90 = sin(qJ(2));
t107 = t85 * t90;
t93 = cos(qJ(2));
t106 = t85 * t93;
t86 = cos(pkin(9));
t105 = t86 * t85;
t87 = cos(pkin(5));
t104 = t87 * t90;
t103 = t87 * t93;
t102 = t86 * pkin(1) + pkin(6) * t108;
t101 = t87 * pkin(6) + qJ(1);
t100 = pkin(6) * t105;
t99 = pkin(2) * t107 + t101;
t68 = -t86 * t103 + t84 * t90;
t69 = t86 * t104 + t84 * t93;
t80 = t84 * pkin(1);
t98 = t69 * pkin(2) + qJ(3) * t68 + t80;
t70 = t84 * t103 + t86 * t90;
t71 = -t84 * t104 + t86 * t93;
t97 = t71 * pkin(2) + qJ(3) * t70 + t102;
t96 = t87 * pkin(3) + pkin(7) * t107 - qJ(3) * t106 + t99;
t92 = cos(qJ(4));
t89 = sin(qJ(4));
t73 = -t89 * t106 + t87 * t92;
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - m(3) * t101 - m(4) * t99 - m(5) * t96 - t73 * mrSges(5,1) - mrSges(5,3) * t107 - m(6) * (t73 * pkin(4) + t96) - (t88 * t107 + t73 * t91) * mrSges(6,1) - (t91 * t107 - t73 * t88) * mrSges(6,2) - t112 * t87 + ((m(4) * qJ(3) - t113) * t93 + t114 * t90) * t85 - t111 * (t92 * t106 + t87 * t89)) * g(3) + (-mrSges(1,2) - t84 * mrSges(2,1) - t86 * mrSges(2,2) - m(3) * (t80 - t100) - m(4) * (t98 - t100) + t115 * (pkin(7) * t69 + (-pkin(3) - pkin(6)) * t105 + t98) + t113 * t68 + t111 * (t89 * t105 + t68 * t92) + t112 * t105 + t109 * (-t92 * t105 + t68 * t89) + t110 * t69) * g(2) + (-m(3) * t102 - m(4) * t97 - t86 * mrSges(2,1) + t84 * mrSges(2,2) - mrSges(1,1) + t115 * (pkin(3) * t108 + pkin(7) * t71 + t97) + t113 * t70 - t111 * (t89 * t108 - t70 * t92) - t112 * t108 + t109 * (t92 * t108 + t70 * t89) + t110 * t71) * g(1);
U = t1;
