% Calculate potential energy for
% S6RPRPRP5
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
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPRP5_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP5_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP5_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP5_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP5_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP5_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:14:21
% EndTime: 2019-03-09 03:14:22
% DurationCPUTime: 0.43s
% Computational Cost: add. (239->95), mult. (227->115), div. (0->0), fcn. (215->10), ass. (0->41)
t110 = rSges(7,1) + pkin(5);
t109 = rSges(5,3) + qJ(4);
t108 = rSges(7,3) + qJ(6);
t85 = sin(qJ(1));
t86 = cos(qJ(1));
t107 = g(1) * t86 + g(2) * t85;
t80 = sin(pkin(9));
t104 = t80 * pkin(2) + pkin(6);
t78 = pkin(9) + qJ(3);
t73 = sin(t78);
t103 = rSges(4,2) * t73;
t75 = cos(t78);
t102 = t85 * t75;
t79 = sin(pkin(10));
t101 = t85 * t79;
t81 = cos(pkin(10));
t100 = t85 * t81;
t99 = t86 * t75;
t98 = t86 * t79;
t97 = t86 * t81;
t82 = cos(pkin(9));
t70 = t82 * pkin(2) + pkin(1);
t84 = -pkin(7) - qJ(2);
t94 = t85 * t70 + t86 * t84;
t93 = rSges(3,3) + qJ(2);
t65 = t86 * t70;
t91 = -t85 * t84 + t65;
t69 = t81 * pkin(4) + pkin(3);
t83 = -pkin(8) - qJ(4);
t90 = t73 * t69 + t75 * t83 + t104;
t89 = pkin(4) * t101 + t69 * t99 + t91;
t88 = rSges(3,1) * t82 - rSges(3,2) * t80 + pkin(1);
t87 = -pkin(4) * t98 + t69 * t102 + t94;
t77 = pkin(10) + qJ(5);
t74 = cos(t77);
t72 = sin(t77);
t60 = t85 * t72 + t74 * t99;
t59 = t72 * t99 - t85 * t74;
t58 = t74 * t102 - t86 * t72;
t57 = t72 * t102 + t86 * t74;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t86 * rSges(2,1) - t85 * rSges(2,2)) + g(2) * (t85 * rSges(2,1) + t86 * rSges(2,2)) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(3) * (t80 * rSges(3,1) + t82 * rSges(3,2) + pkin(6)) + (g(1) * t88 - g(2) * t93) * t86 + (g(1) * t93 + g(2) * t88) * t85) - m(4) * (g(1) * (rSges(4,1) * t99 - t86 * t103 + t65) + g(2) * (-t86 * rSges(4,3) + t94) + g(3) * (t73 * rSges(4,1) + t75 * rSges(4,2) + t104) + (g(1) * (rSges(4,3) - t84) + g(2) * (rSges(4,1) * t75 - t103)) * t85) - m(5) * (g(1) * (pkin(3) * t99 + (t75 * t97 + t101) * rSges(5,1) + (-t75 * t98 + t100) * rSges(5,2) + t91) + g(2) * (pkin(3) * t102 + (t75 * t100 - t98) * rSges(5,1) + (-t75 * t101 - t97) * rSges(5,2) + t94) + g(3) * (-t109 * t75 + t104) + (g(3) * (rSges(5,1) * t81 - rSges(5,2) * t79 + pkin(3)) + t107 * t109) * t73) - m(6) * (g(1) * (t60 * rSges(6,1) - t59 * rSges(6,2) + t89) + g(2) * (t58 * rSges(6,1) - t57 * rSges(6,2) + t87) + g(3) * (-t75 * rSges(6,3) + t90) + (g(3) * (rSges(6,1) * t74 - rSges(6,2) * t72) + t107 * (rSges(6,3) - t83)) * t73) - m(7) * (g(1) * (t108 * t59 + t110 * t60 + t89) + g(2) * (t108 * t57 + t110 * t58 + t87) + g(3) * (-t75 * rSges(7,2) + t90) + (g(3) * (t108 * t72 + t110 * t74) + t107 * (rSges(7,2) - t83)) * t73);
U  = t1;
