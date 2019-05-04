% Calculate potential energy for
% S6RRRRRR10V2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,d6]';
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
% Datum: 2019-04-11 14:56
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRRR10V2_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(6,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10V2_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR10V2_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S6RRRRRR10V2_energypot_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR10V2_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR10V2_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-11 14:41:27
% EndTime: 2019-04-11 14:41:27
% DurationCPUTime: 0.31s
% Computational Cost: add. (223->67), mult. (286->71), div. (0->0), fcn. (295->12), ass. (0->35)
t115 = -m(5) - m(6);
t116 = -m(7) + t115;
t118 = t116 * pkin(5);
t106 = mrSges(4,2) - mrSges(5,3);
t96 = cos(qJ(2));
t83 = t96 * pkin(2) + pkin(1);
t87 = qJ(2) + qJ(3);
t84 = sin(t87);
t85 = cos(t87);
t91 = sin(qJ(2));
t117 = -m(3) * pkin(1) - m(4) * t83 - mrSges(3,1) * t96 - mrSges(4,1) * t85 + mrSges(3,2) * t91 + t106 * t84 - mrSges(2,1) + t116 * (pkin(3) * t85 + t83);
t92 = sin(qJ(1));
t113 = t84 * t92;
t95 = cos(qJ(4));
t112 = t84 * t95;
t97 = cos(qJ(1));
t111 = t84 * t97;
t90 = sin(qJ(4));
t110 = t92 * t90;
t109 = t92 * t95;
t108 = t97 * t90;
t107 = t97 * t95;
t86 = t91 * pkin(2);
t105 = t84 * pkin(3) + t86;
t104 = mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t101 = -m(7) * pkin(6) + mrSges(6,2) - mrSges(7,3);
t88 = sin(qJ(6));
t93 = cos(qJ(6));
t100 = -t93 * mrSges(7,1) + t88 * mrSges(7,2) - mrSges(6,1);
t99 = -t88 * mrSges(7,1) - t93 * mrSges(7,2) + mrSges(5,2) - mrSges(6,3);
t94 = cos(qJ(5));
t89 = sin(qJ(5));
t73 = t85 * t107 + t110;
t71 = t85 * t109 - t108;
t1 = (-mrSges(1,3) - mrSges(2,3) - t91 * mrSges(3,1) - t96 * mrSges(3,2) - m(4) * t86 - m(7) * t105 + (m(7) * pkin(5) - t106) * t85 + t115 * (-t85 * pkin(5) + pkin(4) + t105) + t101 * (t89 * t112 + t85 * t94) + (-m(2) - m(3) - m(4) - m(7)) * pkin(4) + t100 * (t94 * t112 - t85 * t89) + (-mrSges(5,1) * t95 + t99 * t90 - mrSges(4,1)) * t84) * g(3) + (-mrSges(1,2) - t71 * mrSges(5,1) + t100 * (t89 * t113 + t71 * t94) + t99 * (t85 * t110 + t107) + t101 * (-t94 * t113 + t71 * t89) - t104 * t97 + t113 * t118 + t117 * t92) * g(2) + (-mrSges(1,1) - t73 * mrSges(5,1) + t100 * (t89 * t111 + t73 * t94) + t99 * (t85 * t108 - t109) + t101 * (-t94 * t111 + t73 * t89) + t104 * t92 + t111 * t118 + t117 * t97) * g(1);
U  = t1;
