% Calculate potential energy for
% S6RPRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-03-09 03:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPRP1_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP1_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP1_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP1_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP1_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP1_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:01:20
% EndTime: 2019-03-09 03:01:20
% DurationCPUTime: 0.36s
% Computational Cost: add. (224->65), mult. (181->58), div. (0->0), fcn. (152->10), ass. (0->31)
t76 = cos(qJ(5));
t98 = -m(6) * pkin(4) - m(7) * (t76 * pkin(5) + pkin(4)) - mrSges(5,1);
t97 = m(6) * pkin(8) - m(7) * (-qJ(6) - pkin(8)) - mrSges(5,2) + mrSges(6,3) + mrSges(7,3);
t96 = m(7) * pkin(5);
t95 = -m(3) - m(4);
t94 = -mrSges(6,1) - mrSges(7,1);
t93 = mrSges(6,2) + mrSges(7,2);
t92 = -m(5) - m(6) - m(7);
t68 = qJ(3) + pkin(10);
t61 = sin(t68);
t63 = cos(t68);
t74 = sin(qJ(3));
t77 = cos(qJ(3));
t91 = -m(4) * pkin(2) - t77 * mrSges(4,1) + t74 * mrSges(4,2) - t97 * t61 + t98 * t63 - mrSges(3,1);
t90 = m(4) * pkin(7) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3);
t75 = sin(qJ(1));
t66 = t75 * pkin(1);
t78 = cos(qJ(1));
t67 = t78 * pkin(1);
t69 = qJ(1) + pkin(9);
t62 = sin(t69);
t73 = sin(qJ(5));
t89 = t62 * t73;
t88 = t62 * t76;
t64 = cos(t69);
t87 = t64 * t73;
t86 = t64 * t76;
t72 = qJ(2) + pkin(6);
t71 = -qJ(4) - pkin(7);
t60 = t77 * pkin(3) + pkin(2);
t1 = (-m(2) * pkin(6) - t74 * mrSges(4,1) - t77 * mrSges(4,2) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + t95 * t72 + t92 * (t74 * pkin(3) + t72) + t97 * t63 + (t93 * t73 + t94 * t76 + t98) * t61) * g(3) + (t87 * t96 - t75 * mrSges(2,1) - t78 * mrSges(2,2) - mrSges(1,2) + t92 * (t62 * t60 + t64 * t71 + t66) + t95 * t66 + t94 * (t63 * t88 - t87) - t93 * (-t63 * t89 - t86) + t90 * t64 + t91 * t62) * g(2) + (-t89 * t96 - t78 * mrSges(2,1) + t75 * mrSges(2,2) - mrSges(1,1) + t92 * (t64 * t60 - t62 * t71 + t67) + t95 * t67 + t94 * (t63 * t86 + t89) - t93 * (-t63 * t87 + t88) - t90 * t62 + t91 * t64) * g(1);
U  = t1;
