% Calculate potential energy for
% S5RRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRRR11_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR11_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR11_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR11_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR11_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR11_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:38:43
% EndTime: 2019-12-31 22:38:43
% DurationCPUTime: 0.34s
% Computational Cost: add. (219->95), mult. (435->125), div. (0->0), fcn. (511->12), ass. (0->44)
t110 = rSges(4,3) + pkin(8);
t82 = cos(pkin(5));
t109 = t82 * pkin(7) + pkin(6);
t108 = pkin(9) + rSges(5,3);
t81 = sin(pkin(5));
t85 = sin(qJ(2));
t107 = t81 * t85;
t86 = sin(qJ(1));
t106 = t81 * t86;
t88 = cos(qJ(3));
t105 = t81 * t88;
t89 = cos(qJ(2));
t104 = t81 * t89;
t90 = cos(qJ(1));
t103 = t81 * t90;
t102 = t86 * t85;
t101 = t86 * t89;
t100 = t90 * t85;
t99 = t90 * t89;
t98 = t90 * pkin(1) + pkin(7) * t106;
t97 = pkin(10) + pkin(9) + rSges(6,3);
t96 = pkin(2) * t107 + t109;
t68 = -t82 * t102 + t99;
t95 = t68 * pkin(2) + t98;
t66 = t82 * t100 + t101;
t78 = t86 * pkin(1);
t94 = t66 * pkin(2) - pkin(7) * t103 + t78;
t80 = qJ(4) + qJ(5);
t75 = sin(t80);
t76 = cos(t80);
t87 = cos(qJ(4));
t93 = t76 * rSges(6,1) - t75 * rSges(6,2) + t87 * pkin(4) + pkin(3);
t83 = sin(qJ(4));
t92 = t75 * rSges(6,1) + t76 * rSges(6,2) + t83 * pkin(4) + pkin(8);
t84 = sin(qJ(3));
t67 = t82 * t101 + t100;
t65 = -t82 * t99 + t102;
t64 = t85 * t105 + t82 * t84;
t63 = t84 * t107 - t82 * t88;
t60 = t84 * t106 + t68 * t88;
t59 = -t86 * t105 + t68 * t84;
t58 = -t84 * t103 + t66 * t88;
t57 = t88 * t103 + t66 * t84;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t90 * rSges(2,1) - t86 * rSges(2,2)) + g(2) * (t86 * rSges(2,1) + t90 * rSges(2,2)) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (t68 * rSges(3,1) - t67 * rSges(3,2) + t98) + g(2) * (t66 * rSges(3,1) - t65 * rSges(3,2) + t78) + g(3) * (t82 * rSges(3,3) + t109) + (g(1) * rSges(3,3) * t86 + g(3) * (rSges(3,1) * t85 + rSges(3,2) * t89) + g(2) * (-rSges(3,3) - pkin(7)) * t90) * t81) - m(4) * (g(1) * (t60 * rSges(4,1) - t59 * rSges(4,2) + t110 * t67 + t95) + g(2) * (t58 * rSges(4,1) - t57 * rSges(4,2) + t110 * t65 + t94) + g(3) * (t64 * rSges(4,1) - t63 * rSges(4,2) - t110 * t104 + t96)) - m(5) * (g(1) * (t60 * pkin(3) + t67 * pkin(8) + (t60 * t87 + t67 * t83) * rSges(5,1) + (-t60 * t83 + t67 * t87) * rSges(5,2) + t108 * t59 + t95) + g(2) * (t58 * pkin(3) + t65 * pkin(8) + (t58 * t87 + t65 * t83) * rSges(5,1) + (-t58 * t83 + t65 * t87) * rSges(5,2) + t108 * t57 + t94) + g(3) * (t64 * pkin(3) - pkin(8) * t104 + (-t83 * t104 + t64 * t87) * rSges(5,1) + (-t87 * t104 - t64 * t83) * rSges(5,2) + t108 * t63 + t96)) - m(6) * (g(1) * (t97 * t59 + t93 * t60 + t92 * t67 + t95) + g(2) * (t97 * t57 + t93 * t58 + t92 * t65 + t94) + g(3) * (-t92 * t104 + t97 * t63 + t93 * t64 + t96));
U = t1;
