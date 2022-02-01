% Calculate potential energy for
% S5RPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% m [6x1]
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
% Datum: 2022-01-23 09:26
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRPR5_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR5_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR5_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR5_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:24:48
% EndTime: 2022-01-23 09:24:49
% DurationCPUTime: 0.52s
% Computational Cost: add. (152->71), mult. (180->67), div. (0->0), fcn. (161->10), ass. (0->32)
t79 = qJ(4) + pkin(6);
t113 = -mrSges(4,3) - mrSges(5,3) - m(4) * pkin(6) - m(5) * t79 + mrSges(3,2) + m(6) * (-pkin(7) - t79) - mrSges(6,3);
t112 = -m(4) - m(5);
t82 = cos(qJ(3));
t101 = pkin(3) * t82;
t76 = qJ(3) + pkin(9);
t69 = qJ(5) + t76;
t64 = sin(t69);
t65 = cos(t69);
t68 = cos(t76);
t111 = -m(5) * t101 - mrSges(3,1) - m(6) * (pkin(4) * t68 + pkin(2) + t101) - t65 * mrSges(6,1) + t64 * mrSges(6,2);
t77 = sin(pkin(8));
t78 = cos(pkin(8));
t109 = -mrSges(2,1) + t112 * (pkin(2) * t78 + pkin(1)) + t111 * t78 + t113 * t77;
t108 = -m(3) - m(6);
t80 = sin(qJ(3));
t102 = pkin(3) * t80;
t104 = m(5) * (qJ(2) + t102) - mrSges(2,2) + mrSges(3,3) + t64 * mrSges(6,1) + t65 * mrSges(6,2);
t83 = cos(qJ(1));
t98 = t78 * t83;
t97 = t80 * t83;
t67 = sin(t76);
t81 = sin(qJ(1));
t96 = t81 * t67;
t95 = t81 * t68;
t93 = t81 * t80;
t92 = t81 * t82;
t91 = t82 * t83;
t88 = qJ(2) * t83;
t73 = t81 * pkin(1);
t62 = pkin(4) * t67 + t102;
t1 = (-mrSges(1,3) - mrSges(2,3) + t112 * (t77 * pkin(2) + pkin(5)) + (-m(2) + t108) * pkin(5) - t113 * t78 + (-t82 * mrSges(4,1) - t68 * mrSges(5,1) + t80 * mrSges(4,2) + t67 * mrSges(5,2) + t111) * t77) * g(3) + (-mrSges(1,2) - m(3) * (t73 - t88) + m(4) * t88 + t97 * mrSges(4,1) + t91 * mrSges(4,2) - m(6) * t73 + (-t92 * mrSges(4,1) - t95 * mrSges(5,1) + t93 * mrSges(4,2) + t96 * mrSges(5,2)) * t78 + (t67 * mrSges(5,1) + t68 * mrSges(5,2) - m(6) * (-qJ(2) - t62) + t104) * t83 + t109 * t81) * g(2) + (-mrSges(1,1) - (t78 * t91 + t93) * mrSges(4,1) - (-t78 * t97 + t92) * mrSges(4,2) - (t68 * t98 + t96) * mrSges(5,1) - (-t67 * t98 + t95) * mrSges(5,2) + (t108 * pkin(1) + t109) * t83 + (-m(6) * t62 - t104 + (-m(4) + t108) * qJ(2)) * t81) * g(1);
U = t1;
