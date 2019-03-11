% Calculate potential energy for
% S6RPRPRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
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
% Datum: 2019-03-09 03:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPRP9_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP9_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP9_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP9_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP9_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP9_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:27:29
% EndTime: 2019-03-09 03:27:29
% DurationCPUTime: 0.47s
% Computational Cost: add. (171->67), mult. (224->58), div. (0->0), fcn. (203->8), ass. (0->31)
t79 = sin(pkin(9));
t80 = cos(pkin(9));
t115 = -m(5) * pkin(3) - t80 * mrSges(5,1) + t79 * mrSges(5,2) - mrSges(4,1);
t114 = -m(5) * qJ(4) + mrSges(4,2) - mrSges(7,2) - mrSges(5,3) - mrSges(6,3);
t70 = t80 * pkin(4) + pkin(3);
t81 = -pkin(8) - qJ(4);
t82 = sin(qJ(3));
t84 = cos(qJ(3));
t113 = t70 * t82 + t81 * t84;
t112 = -m(4) - m(5);
t111 = -m(6) - m(7);
t109 = -t79 * mrSges(5,1) - t80 * mrSges(5,2) - mrSges(2,1) + mrSges(3,2) - mrSges(4,3);
t108 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t107 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t106 = -t114 * t84 + t115 * t82 + mrSges(2,2) - mrSges(3,3);
t105 = pkin(2) + pkin(6);
t104 = pkin(4) * t79;
t78 = pkin(9) + qJ(5);
t71 = sin(t78);
t83 = sin(qJ(1));
t101 = t83 * t71;
t72 = cos(t78);
t100 = t83 * t72;
t85 = cos(qJ(1));
t97 = t85 * t71;
t96 = t85 * t72;
t75 = t83 * pkin(1);
t95 = t83 * pkin(7) + t75;
t94 = t85 * pkin(1) + t83 * qJ(2);
t92 = t85 * pkin(7) + t94;
t1 = (-mrSges(3,1) - mrSges(1,3) - mrSges(2,3) + t111 * (t84 * t70 - t82 * t81 + t105) + t112 * t105 + (-m(2) - m(3)) * pkin(6) + (-t107 * t71 + t108 * t72 + t115) * t84 + t114 * t82) * g(3) + (-m(3) * t75 - mrSges(1,2) + t112 * t95 + t111 * (t83 * t104 + t95) + t108 * (-t82 * t96 + t101) + t107 * (t82 * t97 + t100) + t109 * t83 + (t111 * (-qJ(2) - t113) + (m(3) - t112) * qJ(2) - t106) * t85) * g(2) + (-m(3) * t94 - mrSges(1,1) + t112 * t92 + t111 * (t85 * t104 + t92) + t108 * (t82 * t100 + t97) - t107 * (t82 * t101 - t96) + t109 * t85 + (t111 * t113 + t106) * t83) * g(1);
U  = t1;
