% Calculate potential energy for
% S6RPPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3,theta4]';
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
% Datum: 2019-03-09 01:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPPRR3_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR3_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR3_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR3_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR3_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR3_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:33:14
% EndTime: 2019-03-09 01:33:14
% DurationCPUTime: 0.34s
% Computational Cost: add. (175->61), mult. (238->52), div. (0->0), fcn. (247->10), ass. (0->29)
t74 = sin(qJ(6));
t75 = cos(qJ(6));
t100 = m(7) * pkin(5) + t75 * mrSges(7,1) - t74 * mrSges(7,2) + mrSges(6,1);
t99 = -m(7) * pkin(8) + mrSges(6,2) - mrSges(7,3);
t98 = -m(4) - m(5);
t97 = -m(6) - m(7);
t96 = -mrSges(2,1) - mrSges(3,1);
t95 = mrSges(2,2) - mrSges(3,3);
t69 = pkin(10) + qJ(5);
t62 = sin(t69);
t63 = cos(t69);
t70 = sin(pkin(10));
t71 = cos(pkin(10));
t93 = m(5) * pkin(3) + t71 * mrSges(5,1) - t70 * mrSges(5,2) + t100 * t63 - t99 * t62 + mrSges(4,1);
t92 = m(5) * qJ(4) + t74 * mrSges(7,1) + t75 * mrSges(7,2) - mrSges(4,2) + mrSges(5,3) + mrSges(6,3);
t91 = cos(qJ(1));
t90 = sin(qJ(1));
t72 = -qJ(3) + pkin(6);
t89 = t91 * pkin(1) + t90 * qJ(2);
t88 = cos(pkin(9));
t87 = sin(pkin(9));
t86 = t91 * pkin(2) + t89;
t83 = t90 * pkin(1) - t91 * qJ(2);
t80 = t90 * pkin(2) + t83;
t73 = -pkin(7) - qJ(4);
t61 = pkin(4) * t71 + pkin(3);
t57 = t91 * t87 - t90 * t88;
t56 = -t90 * t87 - t91 * t88;
t1 = (mrSges(5,1) * t70 + mrSges(5,2) * t71 - mrSges(3,2) - mrSges(1,3) - mrSges(2,3) + mrSges(4,3) + t97 * (-pkin(4) * t70 + t72) + t98 * t72 + t99 * t63 + t100 * t62 + (-m(2) - m(3)) * pkin(6)) * g(3) + (-m(3) * t83 - mrSges(1,2) - t95 * t91 + t96 * t90 + t98 * t80 + t97 * (t56 * t73 - t57 * t61 + t80) + t93 * t57 + t92 * t56) * g(2) + (-m(3) * t89 - mrSges(1,1) + t96 * t91 + t95 * t90 + t98 * t86 + t97 * (-t56 * t61 - t57 * t73 + t86) - t92 * t57 + t93 * t56) * g(1);
U  = t1;
