% Calculate potential energy for
% S5PRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d2,d3,d4,d5,theta1]';
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
% Datum: 2024-09-27 21:46
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRRR11_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR11_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR11_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR11_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR11_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR11_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 21:45:07
% EndTime: 2024-09-27 21:45:07
% DurationCPUTime: 0.09s
% Computational Cost: add. (242->67), mult. (185->60), div. (0->0), fcn. (154->22), ass. (0->33)
t103 = sin(qJ(3));
t104 = cos(qJ(3));
t124 = (-m(5) * pkin(3) - mrSges(4,1)) * t103 - m(6) * (pkin(4) * sin(qJ(4)) * t104 + t103 * (cos(qJ(4)) * pkin(4) + pkin(3))) - t104 * mrSges(4,2);
t105 = pkin(7) + pkin(8);
t121 = m(5) * t105 + m(6) * (pkin(9) + t105) + mrSges(4,3) + mrSges(5,3) + mrSges(6,3) + m(4) * pkin(7);
t117 = -m(3) - m(4) - m(5) - m(6);
t118 = pkin(1) * t117 - mrSges(2,1);
t84 = t104 * pkin(3) + pkin(2);
t97 = qJ(3) + qJ(4);
t90 = cos(t97);
t93 = qJ(5) + t97;
t115 = -m(4) * pkin(2) - m(5) * t84 - m(6) * (pkin(4) * t90 + t84) - t104 * mrSges(4,1) - t90 * mrSges(5,1) - cos(t93) * mrSges(6,1) + t103 * mrSges(4,2) + sin(t97) * mrSges(5,2) + sin(t93) * mrSges(6,2) - mrSges(3,1);
t101 = cos(pkin(5));
t87 = pkin(5) + t97;
t79 = qJ(5) + t87;
t71 = sin(t79) / 0.2e1;
t88 = pkin(5) - t97;
t80 = -qJ(5) + t88;
t72 = cos(t80) / 0.2e1;
t73 = sin(t87) / 0.2e1;
t74 = cos(t88) / 0.2e1;
t75 = sin(t80);
t76 = cos(t79);
t77 = sin(t88);
t78 = cos(t87);
t99 = sin(pkin(5));
t114 = -(t73 - t77 / 0.2e1) * mrSges(5,1) - (t71 - t75 / 0.2e1) * mrSges(6,1) - (t74 + t78 / 0.2e1) * mrSges(5,2) - (t72 + t76 / 0.2e1) * mrSges(6,2) - mrSges(3,2) + t121 * t99 + t124 * t101;
t100 = cos(pkin(10));
t98 = sin(pkin(10));
t95 = pkin(10) + qJ(2);
t86 = cos(t95);
t85 = sin(t95);
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - mrSges(3,3) - (t74 - t78 / 0.2e1) * mrSges(5,1) - (t73 + t77 / 0.2e1) * mrSges(5,2) - (t72 - t76 / 0.2e1) * mrSges(6,1) - (t71 + t75 / 0.2e1) * mrSges(6,2) + t124 * t99 + t117 * (pkin(6) + qJ(1)) - t121 * t101) * g(3) + (-t100 * mrSges(2,2) + t114 * t86 + t115 * t85 + t118 * t98 - mrSges(1,2)) * g(2) + (t98 * mrSges(2,2) + t118 * t100 - t114 * t85 + t115 * t86 - mrSges(1,1)) * g(1);
U = t1;
