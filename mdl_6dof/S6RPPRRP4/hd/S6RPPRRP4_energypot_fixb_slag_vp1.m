% Calculate potential energy for
% S6RPPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
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
% Datum: 2019-03-09 02:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPRRP4_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP4_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP4_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP4_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP4_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRP4_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:05:05
% EndTime: 2019-03-09 02:05:05
% DurationCPUTime: 0.37s
% Computational Cost: add. (181->82), mult. (298->101), div. (0->0), fcn. (340->8), ass. (0->33)
t77 = sin(qJ(4));
t79 = cos(qJ(4));
t102 = -pkin(4) * t79 - pkin(8) * t77;
t101 = rSges(7,1) + pkin(5);
t100 = rSges(7,3) + qJ(6);
t97 = rSges(5,3) + pkin(7);
t96 = sin(qJ(1));
t76 = sin(qJ(5));
t95 = t76 * t79;
t78 = cos(qJ(5));
t94 = t78 * t79;
t93 = pkin(6) - qJ(3);
t80 = cos(qJ(1));
t92 = t80 * pkin(1) + t96 * qJ(2);
t91 = cos(pkin(9));
t90 = sin(pkin(9));
t89 = t79 * pkin(8) + t93;
t88 = t80 * pkin(2) + t92;
t65 = -t80 * t91 - t90 * t96;
t87 = -t65 * pkin(3) + t88;
t72 = t96 * pkin(1);
t86 = t96 * pkin(2) - t80 * qJ(2) + t72;
t66 = t80 * t90 - t91 * t96;
t85 = -g(1) * t65 - g(2) * t66;
t84 = -rSges(5,1) * t79 + rSges(5,2) * t77;
t83 = -t66 * pkin(3) + t86;
t82 = pkin(7) * t66 + t102 * t65 + t87;
t81 = -t65 * pkin(7) + t102 * t66 + t83;
t56 = -t65 * t94 + t66 * t76;
t55 = -t65 * t95 - t66 * t78;
t54 = -t65 * t76 - t66 * t94;
t53 = t65 * t78 - t66 * t95;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t80 * rSges(2,1) - rSges(2,2) * t96) + g(2) * (rSges(2,1) * t96 + t80 * rSges(2,2)) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (t80 * rSges(3,1) + rSges(3,3) * t96 + t92) + g(2) * (t96 * rSges(3,1) + t72 + (-rSges(3,3) - qJ(2)) * t80) + g(3) * (pkin(6) + rSges(3,2))) - m(4) * (g(1) * (-rSges(4,1) * t65 - rSges(4,2) * t66 + t88) + g(2) * (-t66 * rSges(4,1) + t65 * rSges(4,2) + t86) + g(3) * (-rSges(4,3) + t93)) - m(5) * (g(1) * t87 + g(2) * t83 + g(3) * (-rSges(5,1) * t77 - rSges(5,2) * t79 + t93) + (g(1) * t97 + g(2) * t84) * t66 + (g(1) * t84 - g(2) * t97) * t65) - m(6) * (g(1) * (rSges(6,1) * t56 - rSges(6,2) * t55 + t82) + g(2) * (t54 * rSges(6,1) - t53 * rSges(6,2) + t81) + g(3) * (rSges(6,3) * t79 + t89) + (g(3) * (-rSges(6,1) * t78 + rSges(6,2) * t76 - pkin(4)) + t85 * rSges(6,3)) * t77) - m(7) * (g(1) * (t100 * t55 + t101 * t56 + t82) + g(2) * (t100 * t53 + t101 * t54 + t81) + g(3) * (rSges(7,2) * t79 + t89) + (g(3) * (-t100 * t76 - t101 * t78 - pkin(4)) + t85 * rSges(7,2)) * t77);
U  = t1;
