% Calculate potential energy for
% S4RRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,d2,d3,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RRRR6_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(8,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR6_energypot_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR6_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4RRRR6_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR6_energypot_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR6_energypot_fixb_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:29:19
% EndTime: 2019-12-31 17:29:19
% DurationCPUTime: 0.35s
% Computational Cost: add. (130->64), mult. (277->81), div. (0->0), fcn. (311->10), ass. (0->33)
t99 = -m(4) - m(5);
t98 = -m(5) * pkin(8) + mrSges(4,2) - mrSges(5,3);
t73 = sin(qJ(4));
t77 = cos(qJ(4));
t97 = -t73 * mrSges(5,1) - t77 * mrSges(5,2) + mrSges(3,2) - mrSges(4,3);
t96 = -m(5) * pkin(3) - t77 * mrSges(5,1) + t73 * mrSges(5,2) - mrSges(4,1);
t72 = cos(pkin(4));
t95 = t72 * pkin(6) + pkin(5);
t71 = sin(pkin(4));
t75 = sin(qJ(2));
t94 = t71 * t75;
t76 = sin(qJ(1));
t93 = t71 * t76;
t78 = cos(qJ(3));
t92 = t71 * t78;
t79 = cos(qJ(2));
t91 = t71 * t79;
t90 = t75 * t76;
t80 = cos(qJ(1));
t89 = t75 * t80;
t88 = t76 * t79;
t87 = t79 * t80;
t86 = t80 * t71;
t85 = t80 * pkin(1) + pkin(6) * t93;
t84 = t76 * pkin(1) - pkin(6) * t86;
t82 = pkin(2) * t94 - pkin(7) * t91 + t95;
t74 = sin(qJ(3));
t62 = -t72 * t90 + t87;
t61 = t72 * t88 + t89;
t60 = t72 * t89 + t88;
t59 = -t72 * t87 + t90;
t58 = t72 * t74 + t75 * t92;
t1 = (-mrSges(1,3) - m(2) * pkin(5) - mrSges(2,3) - m(3) * t95 - t72 * mrSges(3,3) - (t75 * mrSges(3,1) + t79 * mrSges(3,2)) * t71 - m(4) * t82 - t58 * mrSges(4,1) + mrSges(4,3) * t91 - m(5) * (pkin(3) * t58 + t82) - (t58 * t77 - t73 * t91) * mrSges(5,1) - (-t58 * t73 - t77 * t91) * mrSges(5,2) + t98 * (-t72 * t78 + t74 * t94)) * g(3) + (-m(3) * t84 - t76 * mrSges(2,1) - t60 * mrSges(3,1) - t80 * mrSges(2,2) + mrSges(3,3) * t86 - mrSges(1,2) + t99 * (t60 * pkin(2) + t59 * pkin(7) + t84) + t96 * (t60 * t78 - t74 * t86) + t97 * t59 + t98 * (t60 * t74 + t78 * t86)) * g(2) + (-m(3) * t85 - t80 * mrSges(2,1) - t62 * mrSges(3,1) + t76 * mrSges(2,2) - mrSges(3,3) * t93 - mrSges(1,1) + t99 * (t62 * pkin(2) + pkin(7) * t61 + t85) + t96 * (t62 * t78 + t74 * t93) + t97 * t61 + t98 * (t62 * t74 - t76 * t92)) * g(1);
U = t1;
