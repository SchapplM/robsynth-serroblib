% Calculate potential energy for
% S6RPRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-03-09 03:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPRR3_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR3_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR3_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR3_energypot_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR3_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR3_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:40:37
% EndTime: 2019-03-09 03:40:38
% DurationCPUTime: 0.36s
% Computational Cost: add. (245->61), mult. (207->51), div. (0->0), fcn. (182->12), ass. (0->26)
t70 = pkin(11) + qJ(5);
t65 = qJ(6) + t70;
t58 = sin(t65);
t59 = cos(t65);
t73 = cos(pkin(11));
t60 = t73 * pkin(4) + pkin(3);
t61 = sin(t70);
t63 = cos(t70);
t72 = sin(pkin(11));
t100 = -mrSges(4,1) - m(5) * pkin(3) - mrSges(5,1) * t73 + mrSges(5,2) * t72 - m(6) * t60 - mrSges(6,1) * t63 + mrSges(6,2) * t61 - m(7) * (pkin(5) * t63 + t60) - mrSges(7,1) * t59 + mrSges(7,2) * t58;
t75 = -pkin(8) - qJ(4);
t99 = -m(5) * qJ(4) + mrSges(4,2) - mrSges(5,3) + m(7) * (-pkin(9) + t75) - mrSges(7,3) + m(6) * t75 - mrSges(6,3);
t98 = m(4) + m(5) + m(6) + m(7);
t76 = sin(qJ(3));
t78 = cos(qJ(3));
t94 = t100 * t78 + t99 * t76 - mrSges(3,1);
t92 = pkin(4) * t72;
t93 = m(6) * t92 + m(7) * (pkin(5) * t61 + t92) - mrSges(3,2) + mrSges(4,3) + t58 * mrSges(7,1) + t59 * mrSges(7,2) + t61 * mrSges(6,1) + t63 * mrSges(6,2) + t72 * mrSges(5,1) + t73 * mrSges(5,2);
t77 = sin(qJ(1));
t67 = t77 * pkin(1);
t79 = cos(qJ(1));
t68 = t79 * pkin(1);
t71 = qJ(1) + pkin(10);
t64 = cos(t71);
t62 = sin(t71);
t1 = (-m(2) * pkin(6) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + (-m(3) - t98) * (qJ(2) + pkin(6)) - t99 * t78 + t100 * t76) * g(3) + (-m(3) * t67 - t77 * mrSges(2,1) - t79 * mrSges(2,2) - mrSges(1,2) - t98 * (t62 * pkin(2) + t67) + (t98 * pkin(7) + t93) * t64 + t94 * t62) * g(2) + (-m(3) * t68 - t79 * mrSges(2,1) + t77 * mrSges(2,2) - mrSges(1,1) - t98 * (t64 * pkin(2) + t62 * pkin(7) + t68) + t94 * t64 - t93 * t62) * g(1);
U  = t1;
