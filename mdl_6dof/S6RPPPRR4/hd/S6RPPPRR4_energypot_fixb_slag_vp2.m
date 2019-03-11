% Calculate potential energy for
% S6RPPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3]';
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
% Datum: 2019-03-09 01:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPPRR4_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR4_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR4_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR4_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR4_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR4_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:35:14
% EndTime: 2019-03-09 01:35:14
% DurationCPUTime: 0.38s
% Computational Cost: add. (149->55), mult. (232->44), div. (0->0), fcn. (243->8), ass. (0->24)
t87 = m(6) + m(7);
t88 = m(5) + t87;
t59 = sin(qJ(6));
t61 = cos(qJ(6));
t86 = m(7) * pkin(5) + t61 * mrSges(7,1) - t59 * mrSges(7,2) + mrSges(6,1);
t85 = m(7) * pkin(8) - mrSges(6,2) + mrSges(7,3);
t83 = -mrSges(2,1) - mrSges(3,1);
t82 = mrSges(2,2) - mrSges(3,3);
t80 = t59 * mrSges(7,1) + t61 * mrSges(7,2) + t87 * pkin(7) + mrSges(4,1) - mrSges(5,2) + mrSges(6,3);
t60 = sin(qJ(5));
t62 = cos(qJ(5));
t79 = t88 * qJ(4) + t86 * t60 - t85 * t62 - mrSges(4,2) + mrSges(5,3);
t78 = cos(qJ(1));
t77 = sin(qJ(1));
t58 = -qJ(3) + pkin(6);
t75 = t78 * pkin(1) + t77 * qJ(2);
t74 = cos(pkin(9));
t73 = sin(pkin(9));
t71 = t78 * pkin(2) + t75;
t68 = t77 * pkin(1) - qJ(2) * t78;
t67 = t77 * pkin(2) + t68;
t49 = t73 * t78 - t77 * t74;
t48 = -t77 * t73 - t78 * t74;
t1 = (mrSges(5,1) - mrSges(3,2) - mrSges(1,3) - mrSges(2,3) + mrSges(4,3) - t87 * (-pkin(4) + t58) + t86 * t62 + t85 * t60 + (-m(4) - m(5)) * t58 + (-m(2) - m(3)) * pkin(6)) * g(3) + (-m(3) * t68 - m(4) * t67 - mrSges(1,2) - t82 * t78 + t83 * t77 - t88 * (-t49 * pkin(3) + t67) + t80 * t49 + t79 * t48) * g(2) + (-m(3) * t75 - m(4) * t71 - mrSges(1,1) + t83 * t78 + t82 * t77 - t88 * (-t48 * pkin(3) + t71) - t79 * t49 + t80 * t48) * g(1);
U  = t1;
