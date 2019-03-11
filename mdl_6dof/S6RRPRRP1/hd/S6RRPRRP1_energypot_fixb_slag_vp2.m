% Calculate potential energy for
% S6RRPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 11:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRP1_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP1_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP1_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP1_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP1_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP1_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:38:54
% EndTime: 2019-03-09 11:38:55
% DurationCPUTime: 0.36s
% Computational Cost: add. (225->65), mult. (195->58), div. (0->0), fcn. (166->10), ass. (0->30)
t78 = cos(qJ(5));
t103 = -m(6) * pkin(4) - m(7) * (t78 * pkin(5) + pkin(4)) - mrSges(5,1);
t102 = m(6) * pkin(9) - m(7) * (-qJ(6) - pkin(9)) - mrSges(5,2) + mrSges(6,3) + mrSges(7,3);
t101 = m(7) * pkin(5);
t100 = -mrSges(6,1) - mrSges(7,1);
t99 = mrSges(6,2) + mrSges(7,2);
t98 = -m(5) - m(6) - m(7);
t74 = -qJ(3) - pkin(7);
t97 = -m(3) * pkin(7) + m(4) * t74 + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t72 = qJ(2) + pkin(10);
t68 = qJ(4) + t72;
t62 = sin(t68);
t63 = cos(t68);
t79 = cos(qJ(2));
t65 = t79 * pkin(2) + pkin(1);
t66 = sin(t72);
t67 = cos(t72);
t76 = sin(qJ(2));
t96 = -m(3) * pkin(1) - m(4) * t65 - t79 * mrSges(3,1) - t67 * mrSges(4,1) + t76 * mrSges(3,2) + t66 * mrSges(4,2) - t102 * t62 + t103 * t63 - mrSges(2,1);
t95 = t76 * pkin(2) + pkin(6);
t75 = sin(qJ(5));
t77 = sin(qJ(1));
t94 = t77 * t75;
t93 = t77 * t78;
t80 = cos(qJ(1));
t92 = t80 * t75;
t91 = t80 * t78;
t71 = -pkin(8) + t74;
t59 = pkin(3) * t67 + t65;
t1 = (-m(4) * t95 - t76 * mrSges(3,1) - t66 * mrSges(4,1) - t79 * mrSges(3,2) - t67 * mrSges(4,2) - mrSges(1,3) - mrSges(2,3) + (-m(2) - m(3)) * pkin(6) + t98 * (pkin(3) * t66 + t95) + t102 * t63 + (t100 * t78 + t99 * t75 + t103) * t62) * g(3) + (t92 * t101 - mrSges(1,2) + t98 * (t77 * t59 + t80 * t71) + t100 * (t63 * t93 - t92) - t99 * (-t63 * t94 - t91) - t97 * t80 + t96 * t77) * g(2) + (-t94 * t101 - mrSges(1,1) + t98 * (t80 * t59 - t77 * t71) + t100 * (t63 * t91 + t94) - t99 * (-t63 * t92 + t93) + t97 * t77 + t96 * t80) * g(1);
U  = t1;
