% Calculate potential energy for
% S6RPRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
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
% Datum: 2019-03-09 06:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRP10_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP10_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP10_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP10_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP10_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP10_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:30:12
% EndTime: 2019-03-09 06:30:13
% DurationCPUTime: 0.48s
% Computational Cost: add. (171->67), mult. (224->58), div. (0->0), fcn. (203->8), ass. (0->31)
t76 = sin(qJ(4));
t79 = cos(qJ(4));
t112 = -m(5) * pkin(3) - t79 * mrSges(5,1) + t76 * mrSges(5,2) - mrSges(4,1);
t111 = -m(5) * pkin(8) + mrSges(4,2) - mrSges(7,2) - mrSges(5,3) - mrSges(6,3);
t67 = t79 * pkin(4) + pkin(3);
t77 = sin(qJ(3));
t80 = cos(qJ(3));
t82 = -pkin(9) - pkin(8);
t110 = t67 * t77 + t80 * t82;
t109 = -m(4) - m(5);
t108 = -m(6) - m(7);
t106 = -t76 * mrSges(5,1) - t79 * mrSges(5,2) - mrSges(2,1) + mrSges(3,2) - mrSges(4,3);
t105 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t104 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t103 = -t111 * t80 + t112 * t77 + mrSges(2,2) - mrSges(3,3);
t102 = pkin(2) + pkin(6);
t101 = pkin(4) * t76;
t75 = qJ(4) + qJ(5);
t68 = sin(t75);
t78 = sin(qJ(1));
t99 = t78 * t68;
t69 = cos(t75);
t98 = t78 * t69;
t81 = cos(qJ(1));
t94 = t81 * t68;
t93 = t81 * t69;
t72 = t78 * pkin(1);
t92 = t78 * pkin(7) + t72;
t91 = t81 * pkin(1) + t78 * qJ(2);
t89 = t81 * pkin(7) + t91;
t1 = (-mrSges(3,1) - mrSges(1,3) - mrSges(2,3) + t108 * (t80 * t67 - t77 * t82 + t102) + t109 * t102 + (-m(2) - m(3)) * pkin(6) + (-t104 * t68 + t105 * t69 + t112) * t80 + t111 * t77) * g(3) + (-m(3) * t72 - mrSges(1,2) + t109 * t92 + t108 * (t78 * t101 + t92) + t105 * (-t77 * t93 + t99) + t104 * (t77 * t94 + t98) + t106 * t78 + (t108 * (-qJ(2) - t110) + (m(3) - t109) * qJ(2) - t103) * t81) * g(2) + (-m(3) * t91 - mrSges(1,1) + t109 * t89 + t108 * (t81 * t101 + t89) + t105 * (t77 * t98 + t94) - t104 * (t77 * t99 - t93) + t106 * t81 + (t108 * t110 + t103) * t78) * g(1);
U  = t1;
