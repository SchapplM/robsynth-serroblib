% Calculate potential energy for
% S6RRRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 22:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRPR2_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR2_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR2_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR2_energypot_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR2_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR2_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:57:18
% EndTime: 2019-03-09 21:57:19
% DurationCPUTime: 0.34s
% Computational Cost: add. (235->61), mult. (195->49), div. (0->0), fcn. (166->12), ass. (0->26)
t68 = pkin(11) + qJ(6);
t61 = sin(t68);
t62 = cos(t68);
t71 = sin(pkin(11));
t72 = cos(pkin(11));
t100 = -mrSges(5,1) - m(6) * pkin(4) - t72 * mrSges(6,1) + t71 * mrSges(6,2) - m(7) * (pkin(5) * t72 + pkin(4)) - t62 * mrSges(7,1) + t61 * mrSges(7,2);
t99 = -m(6) * qJ(5) + mrSges(5,2) - mrSges(6,3) + m(7) * (-pkin(10) - qJ(5)) - mrSges(7,3);
t98 = m(5) + m(6) + m(7);
t70 = qJ(2) + qJ(3);
t65 = qJ(4) + t70;
t58 = sin(t65);
t59 = cos(t65);
t76 = cos(qJ(2));
t60 = t76 * pkin(2) + pkin(1);
t63 = sin(t70);
t64 = cos(t70);
t74 = sin(qJ(2));
t95 = -m(3) * pkin(1) - m(4) * t60 - t76 * mrSges(3,1) - t64 * mrSges(4,1) + t74 * mrSges(3,2) + t63 * mrSges(4,2) + t100 * t59 + t99 * t58 - mrSges(2,1);
t78 = -pkin(8) - pkin(7);
t94 = m(3) * pkin(7) - m(4) * t78 + t61 * mrSges(7,1) + t72 * mrSges(6,2) + t62 * mrSges(7,2) - mrSges(2,2) + mrSges(3,3) + mrSges(4,3) + mrSges(5,3) + (m(7) * pkin(5) + mrSges(6,1)) * t71;
t92 = t74 * pkin(2) + pkin(6);
t77 = cos(qJ(1));
t75 = sin(qJ(1));
t69 = -pkin(9) + t78;
t54 = pkin(3) * t64 + t60;
t1 = (-m(4) * t92 - t74 * mrSges(3,1) - t63 * mrSges(4,1) - t76 * mrSges(3,2) - t64 * mrSges(4,2) - mrSges(1,3) - mrSges(2,3) + (-m(2) - m(3)) * pkin(6) - t98 * (pkin(3) * t63 + t92) - t99 * t59 + t100 * t58) * g(3) + (-mrSges(1,2) - t98 * (t75 * t54 + t77 * t69) + t94 * t77 + t95 * t75) * g(2) + (-mrSges(1,1) + (t98 * t69 - t94) * t75 + (-t98 * t54 + t95) * t77) * g(1);
U  = t1;
