% Calculate potential energy for
% S6RPRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPPR5_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR5_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR5_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR5_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR5_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPPR5_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:49:57
% EndTime: 2019-03-09 02:49:58
% DurationCPUTime: 0.47s
% Computational Cost: add. (206->98), mult. (204->120), div. (0->0), fcn. (184->10), ass. (0->40)
t107 = rSges(7,3) + pkin(8) + qJ(5);
t106 = rSges(6,3) + qJ(5);
t79 = sin(qJ(1));
t80 = cos(qJ(1));
t105 = g(1) * t80 + g(2) * t79;
t74 = sin(pkin(9));
t102 = t74 * pkin(2) + pkin(6);
t72 = pkin(9) + qJ(3);
t67 = sin(t72);
t101 = t67 * t80;
t71 = pkin(10) + qJ(6);
t68 = cos(t71);
t100 = t68 * t80;
t69 = cos(t72);
t99 = t69 * t80;
t73 = sin(pkin(10));
t98 = t73 * t80;
t75 = cos(pkin(10));
t97 = t75 * t80;
t66 = sin(t71);
t96 = t79 * t66;
t95 = t79 * t68;
t94 = t79 * t73;
t93 = t79 * t75;
t76 = cos(pkin(9));
t64 = pkin(2) * t76 + pkin(1);
t78 = -pkin(7) - qJ(2);
t91 = t79 * t64 + t80 * t78;
t90 = qJ(4) * t67;
t89 = rSges(3,3) + qJ(2);
t87 = t67 * pkin(3) + t102;
t86 = t67 * t98;
t85 = t67 * t94;
t59 = t80 * t64;
t84 = pkin(3) * t99 + t80 * t90 + t59;
t83 = t91 + (pkin(3) * t69 + t90) * t79;
t82 = -t79 * t78 + t84;
t81 = rSges(3,1) * t76 - rSges(3,2) * t74 + pkin(1);
t63 = pkin(5) * t75 + pkin(4);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t80 - t79 * rSges(2,2)) + g(2) * (t79 * rSges(2,1) + rSges(2,2) * t80) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(3) * (rSges(3,1) * t74 + rSges(3,2) * t76 + pkin(6)) + (g(1) * t81 - g(2) * t89) * t80 + (g(1) * t89 + g(2) * t81) * t79) - m(4) * (g(1) * (rSges(4,1) * t99 - rSges(4,2) * t101 + t59) + g(2) * (-rSges(4,3) * t80 + t91) + g(3) * (rSges(4,1) * t67 + rSges(4,2) * t69 + t102) + (g(1) * (rSges(4,3) - t78) + g(2) * (rSges(4,1) * t69 - rSges(4,2) * t67)) * t79) - m(5) * (g(1) * (-rSges(5,2) * t99 + rSges(5,3) * t101 + t84) + g(2) * (-rSges(5,1) * t80 + t83) + g(3) * (-rSges(5,2) * t67 + (-rSges(5,3) - qJ(4)) * t69 + t87) + (g(1) * (rSges(5,1) - t78) + g(2) * (-rSges(5,2) * t69 + rSges(5,3) * t67)) * t79) - m(6) * (g(1) * (t79 * pkin(4) + (t86 + t93) * rSges(6,1) + (t67 * t97 - t94) * rSges(6,2) + t82) + g(2) * (-t80 * pkin(4) + (t85 - t97) * rSges(6,1) + (t67 * t93 + t98) * rSges(6,2) + t83) + g(3) * (t106 * t67 + t87) + (g(3) * (-rSges(6,1) * t73 - rSges(6,2) * t75 - qJ(4)) + t105 * t106) * t69) - m(7) * (g(1) * (t79 * t63 + pkin(5) * t86 + (t66 * t101 + t95) * rSges(7,1) + (t67 * t100 - t96) * rSges(7,2) + t82) + g(2) * (-t80 * t63 + pkin(5) * t85 + (t67 * t96 - t100) * rSges(7,1) + (t66 * t80 + t67 * t95) * rSges(7,2) + t83) + g(3) * (t107 * t67 + t87) + (g(3) * (-rSges(7,1) * t66 - rSges(7,2) * t68 - pkin(5) * t73 - qJ(4)) + t105 * t107) * t69);
U  = t1;
