% Calculate potential energy for
% S5PPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
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
% Datum: 2019-12-05 15:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PPRPR3_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR3_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR3_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR3_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR3_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR3_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:04:41
% EndTime: 2019-12-05 15:04:42
% DurationCPUTime: 0.48s
% Computational Cost: add. (162->60), mult. (223->60), div. (0->0), fcn. (216->10), ass. (0->28)
t79 = sin(qJ(3));
t81 = cos(qJ(3));
t107 = -m(4) * pkin(2) - t81 * mrSges(4,1) + t79 * mrSges(4,2) - mrSges(3,1);
t78 = sin(qJ(5));
t80 = cos(qJ(5));
t106 = -m(4) * pkin(5) - t78 * mrSges(6,1) - t80 * mrSges(6,2) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t103 = m(3) + m(4);
t102 = -m(5) - m(6);
t101 = -t79 * mrSges(4,1) - t81 * mrSges(4,2) + mrSges(2,2) - mrSges(3,3);
t99 = -m(6) * pkin(6) + mrSges(5,2) - mrSges(6,3);
t98 = -m(6) * pkin(4) - mrSges(6,1) * t80 + mrSges(6,2) * t78 - mrSges(5,1);
t73 = sin(pkin(8));
t75 = cos(pkin(8));
t97 = t106 * t73 + t107 * t75 - mrSges(2,1);
t96 = pkin(3) * t79;
t77 = -qJ(4) - pkin(5);
t94 = t73 * t77;
t74 = sin(pkin(7));
t91 = t74 * t75;
t76 = cos(pkin(7));
t90 = t75 * t76;
t89 = t76 * pkin(1) + t74 * qJ(2);
t72 = qJ(3) + pkin(9);
t70 = t74 * pkin(1);
t68 = cos(t72);
t67 = sin(t72);
t66 = pkin(3) * t81 + pkin(2);
t1 = (-mrSges(1,3) - mrSges(2,3) + t102 * (t73 * t66 + t75 * t77 + qJ(1)) + (-m(2) - t103) * qJ(1) - t106 * t75 + (t99 * t67 + t98 * t68 + t107) * t73) * g(3) + (-mrSges(1,2) + t102 * (-t74 * t94 + t66 * t91 + t70 + (-qJ(2) - t96) * t76) - t103 * t70 + t99 * (t67 * t91 + t68 * t76) + (t103 * qJ(2) - t101) * t76 + t98 * (-t67 * t76 + t68 * t91) + t97 * t74) * g(2) + (-mrSges(1,1) - t103 * t89 + t102 * (t66 * t90 + t74 * t96 - t76 * t94 + t89) + t99 * (t67 * t90 - t74 * t68) + t101 * t74 + t98 * (t67 * t74 + t68 * t90) + t97 * t76) * g(1);
U = t1;
