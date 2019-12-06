% Calculate potential energy for
% S5PPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
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
% Datum: 2019-12-05 14:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PPPRR1_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR1_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPPRR1_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR1_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPPRR1_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPPRR1_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:57:45
% EndTime: 2019-12-05 14:57:45
% DurationCPUTime: 0.47s
% Computational Cost: add. (162->60), mult. (223->59), div. (0->0), fcn. (216->10), ass. (0->29)
t73 = sin(pkin(9));
t76 = cos(pkin(9));
t108 = -m(4) * pkin(2) - t76 * mrSges(4,1) + t73 * mrSges(4,2) - mrSges(3,1);
t80 = sin(qJ(5));
t81 = cos(qJ(5));
t107 = -m(4) * qJ(3) - t80 * mrSges(6,1) - t81 * mrSges(6,2) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t104 = m(3) + m(4);
t103 = -m(5) - m(6);
t102 = -t73 * mrSges(4,1) - t76 * mrSges(4,2) + mrSges(2,2) - mrSges(3,3);
t100 = -m(6) * pkin(6) + mrSges(5,2) - mrSges(6,3);
t99 = -m(6) * pkin(4) - t81 * mrSges(6,1) + t80 * mrSges(6,2) - mrSges(5,1);
t74 = sin(pkin(8));
t77 = cos(pkin(8));
t98 = t107 * t74 + t108 * t77 - mrSges(2,1);
t97 = pkin(3) * t73;
t79 = -pkin(5) - qJ(3);
t95 = t74 * t79;
t75 = sin(pkin(7));
t92 = t75 * t77;
t72 = pkin(9) + qJ(4);
t67 = sin(t72);
t78 = cos(pkin(7));
t91 = t78 * t67;
t68 = cos(t72);
t90 = t78 * t68;
t89 = t78 * pkin(1) + t75 * qJ(2);
t70 = t75 * pkin(1);
t66 = t76 * pkin(3) + pkin(2);
t1 = (-mrSges(1,3) - mrSges(2,3) + t103 * (t74 * t66 + t77 * t79 + qJ(1)) + (-m(2) - t104) * qJ(1) - t107 * t77 + (t100 * t67 + t99 * t68 + t108) * t74) * g(3) + (-mrSges(1,2) + t103 * (-t75 * t95 + t66 * t92 + t70 + (-qJ(2) - t97) * t78) - t104 * t70 + t100 * (t67 * t92 + t90) + (t104 * qJ(2) - t102) * t78 + t99 * (t68 * t92 - t91) + t98 * t75) * g(2) + (-mrSges(1,1) - t104 * t89 + t103 * (t75 * t97 + t89) + t100 * (-t75 * t68 + t77 * t91) + t102 * t75 + t99 * (t75 * t67 + t77 * t90) + (t103 * (t66 * t77 - t95) + t98) * t78) * g(1);
U = t1;
