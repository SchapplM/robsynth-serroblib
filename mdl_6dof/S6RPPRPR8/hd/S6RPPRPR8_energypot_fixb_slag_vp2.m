% Calculate potential energy for
% S6RPPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3]';
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
% Datum: 2019-03-09 01:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPRPR8_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR8_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR8_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR8_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR8_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR8_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:55:13
% EndTime: 2019-03-09 01:55:14
% DurationCPUTime: 0.46s
% Computational Cost: add. (155->73), mult. (176->56), div. (0->0), fcn. (143->8), ass. (0->29)
t85 = -mrSges(5,1) + mrSges(6,2);
t84 = m(3) + m(4);
t83 = m(6) + m(7);
t60 = sin(qJ(6));
t62 = cos(qJ(6));
t82 = -t60 * mrSges(7,1) - t62 * mrSges(7,2) - mrSges(6,3);
t81 = m(7) * pkin(8) + mrSges(7,3);
t56 = pkin(9) + qJ(4);
t50 = sin(t56);
t51 = cos(t56);
t57 = sin(pkin(9));
t58 = cos(pkin(9));
t80 = -mrSges(4,1) * t57 - mrSges(4,2) * t58 - t51 * mrSges(5,2) + t85 * t50 + mrSges(2,2) - mrSges(3,3);
t59 = -pkin(7) - qJ(3);
t79 = -m(4) * qJ(3) - mrSges(2,1) - mrSges(6,1) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3) - m(7) * (pkin(5) - t59) - t62 * mrSges(7,1) + t60 * mrSges(7,2);
t78 = pkin(2) + pkin(6);
t77 = pkin(3) * t57;
t76 = pkin(4) * t50;
t61 = sin(qJ(1));
t63 = cos(qJ(1));
t74 = t63 * pkin(1) + t61 * qJ(2);
t73 = t58 * pkin(3) + t78;
t72 = t61 * t77 + t74;
t54 = t61 * pkin(1);
t70 = -t61 * t59 + t54;
t65 = -t59 * t63 + t72;
t46 = t61 * t76;
t45 = t63 * t51 * qJ(5);
t1 = (-m(4) * t78 - m(5) * t73 - t58 * mrSges(4,1) + t57 * mrSges(4,2) - mrSges(3,1) - mrSges(1,3) - mrSges(2,3) - t83 * (t51 * pkin(4) + t50 * qJ(5) + t73) + (-m(2) - m(3)) * pkin(6) + (-t81 + t85) * t51 + (mrSges(5,2) + t82) * t50) * g(3) + (-mrSges(1,2) - m(5) * t70 - m(6) * (t45 + t70) - m(7) * t45 + (-m(7) - t84) * t54 + (m(6) * t76 - (m(7) * (-pkin(4) - pkin(8)) - mrSges(7,3)) * t50 + (-m(5) - t83) * (-qJ(2) - t77) + t82 * t51 + t84 * qJ(2) - t80) * t63 + t79 * t61) * g(2) + (-mrSges(1,1) - m(5) * t65 - m(6) * (t46 + t65) - m(7) * (t46 + t72) - t84 * t74 + t79 * t63 + (-t81 * t50 + (t83 * qJ(5) - t82) * t51 + t80) * t61) * g(1);
U  = t1;
