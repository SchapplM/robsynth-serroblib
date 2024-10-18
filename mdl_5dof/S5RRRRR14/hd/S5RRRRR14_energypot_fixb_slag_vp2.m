% Calculate potential energy for
% S5RRRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,d3,d4,d5]';
% m [6x1]
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
% Datum: 2024-09-27 18:44
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRRR14_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR14_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR14_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR14_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR14_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR14_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 18:42:17
% EndTime: 2024-09-27 18:42:17
% DurationCPUTime: 0.12s
% Computational Cost: add. (242->73), mult. (185->71), div. (0->0), fcn. (154->22), ass. (0->37)
t101 = sin(qJ(3));
t103 = cos(qJ(3));
t124 = -m(5) * pkin(3) * t101 - m(6) * (pkin(4) * sin(qJ(4)) * t103 + t101 * (cos(qJ(4)) * pkin(4) + pkin(3)));
t105 = pkin(8) + pkin(9);
t123 = m(5) * t105 + m(6) * (pkin(10) + t105) + mrSges(4,3) + mrSges(5,3) + mrSges(6,3) + m(4) * pkin(8);
t120 = -m(3) - m(4) - m(5) - m(6);
t121 = pkin(1) * t120 - mrSges(2,1);
t84 = t103 * pkin(3) + pkin(2);
t96 = qJ(3) + qJ(4);
t89 = cos(t96);
t91 = qJ(5) + t96;
t118 = -m(4) * pkin(2) - m(5) * t84 - m(6) * (pkin(4) * t89 + t84) - t89 * mrSges(5,1) - cos(t91) * mrSges(6,1) + sin(t96) * mrSges(5,2) + sin(t91) * mrSges(6,2) - mrSges(3,1);
t85 = pkin(5) + t96;
t79 = qJ(5) + t85;
t71 = sin(t79) / 0.2e1;
t86 = pkin(5) - t96;
t80 = -qJ(5) + t86;
t72 = cos(t80) / 0.2e1;
t73 = sin(t85) / 0.2e1;
t74 = cos(t86) / 0.2e1;
t75 = sin(t80);
t76 = cos(t79);
t77 = sin(t86);
t78 = cos(t85);
t98 = sin(pkin(5));
t99 = cos(pkin(5));
t117 = -(t73 - t77 / 0.2e1) * mrSges(5,1) - (t71 - t75 / 0.2e1) * mrSges(6,1) - (t74 + t78 / 0.2e1) * mrSges(5,2) - (t72 + t76 / 0.2e1) * mrSges(6,2) - mrSges(3,2) + t124 * t99 + t123 * t98;
t97 = qJ(1) + qJ(2);
t88 = sin(t97);
t110 = t88 * t101;
t109 = t88 * t103;
t90 = cos(t97);
t108 = t90 * t101;
t107 = t90 * t103;
t104 = cos(qJ(1));
t102 = sin(qJ(1));
t1 = (-mrSges(1,3) - m(2) * pkin(6) - mrSges(2,3) - mrSges(3,3) - (t74 - t78 / 0.2e1) * mrSges(5,1) - (t73 + t77 / 0.2e1) * mrSges(5,2) - (t72 - t76 / 0.2e1) * mrSges(6,1) - (t71 + t75 / 0.2e1) * mrSges(6,2) + (-t101 * mrSges(4,1) - t103 * mrSges(4,2) + t124) * t98 + t120 * (pkin(7) + pkin(6)) - t123 * t99) * g(3) + (-mrSges(1,2) - t104 * mrSges(2,2) - (t108 * t99 + t109) * mrSges(4,1) - (t107 * t99 - t110) * mrSges(4,2) + t118 * t88 + t117 * t90 + t121 * t102) * g(2) + (-mrSges(1,1) + t102 * mrSges(2,2) - (-t110 * t99 + t107) * mrSges(4,1) - (-t109 * t99 - t108) * mrSges(4,2) + t118 * t90 - t117 * t88 + t121 * t104) * g(1);
U = t1;
