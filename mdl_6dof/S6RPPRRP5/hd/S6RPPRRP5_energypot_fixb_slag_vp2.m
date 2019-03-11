% Calculate potential energy for
% S6RPPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
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
% Datum: 2019-03-09 02:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPRRP5_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP5_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP5_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP5_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP5_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP5_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:07:37
% EndTime: 2019-03-09 02:07:37
% DurationCPUTime: 0.33s
% Computational Cost: add. (118->57), mult. (172->50), div. (0->0), fcn. (143->6), ass. (0->24)
t59 = cos(qJ(5));
t84 = -m(6) * pkin(4) - m(7) * (t59 * pkin(5) + pkin(4)) - mrSges(5,1);
t83 = m(6) * pkin(8) - m(7) * (-qJ(6) - pkin(8)) - mrSges(5,2) + mrSges(6,3) + mrSges(7,3);
t82 = m(7) * pkin(5);
t81 = -mrSges(6,1) - mrSges(7,1);
t80 = mrSges(6,2) + mrSges(7,2);
t79 = -m(5) - m(6) - m(7);
t78 = mrSges(2,2) - mrSges(3,3) - mrSges(4,2) + mrSges(5,3);
t57 = sin(qJ(4));
t60 = cos(qJ(4));
t77 = t84 * t57 + t83 * t60 - mrSges(2,1) + mrSges(3,2) - mrSges(4,3);
t76 = pkin(2) + pkin(6);
t56 = sin(qJ(5));
t58 = sin(qJ(1));
t75 = t58 * t56;
t74 = t58 * t59;
t61 = cos(qJ(1));
t73 = t61 * t56;
t72 = t61 * t59;
t71 = t61 * pkin(1) + t58 * qJ(2);
t69 = t61 * qJ(3) + t71;
t68 = t58 * pkin(1) - t61 * qJ(2);
t67 = t58 * qJ(3) + t68;
t1 = (-m(4) * t76 - mrSges(3,1) - mrSges(4,1) - mrSges(1,3) - mrSges(2,3) + (-m(2) - m(3)) * pkin(6) + t79 * (pkin(3) + t76) + (t80 * t56 + t81 * t59 + t84) * t60 - t83 * t57) * g(3) + (-t73 * t82 - m(3) * t68 - m(4) * t67 - mrSges(1,2) + t79 * (t61 * pkin(7) + t67) + t81 * (t57 * t74 + t73) - t80 * (-t57 * t75 + t72) - t78 * t61 + t77 * t58) * g(2) + (t75 * t82 - m(3) * t71 - m(4) * t69 - mrSges(1,1) + t79 * (-t58 * pkin(7) + t69) + t81 * (t57 * t72 - t75) - t80 * (-t57 * t73 - t74) + t78 * t58 + t77 * t61) * g(1);
U  = t1;
