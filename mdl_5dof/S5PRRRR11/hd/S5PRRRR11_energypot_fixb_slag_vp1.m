% Calculate potential energy for
% S5PRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d2,d3,d4,d5,theta1]';
% m [6x1]
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
% Datum: 2024-09-27 21:46
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRRR11_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR11_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR11_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR11_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR11_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR11_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 21:45:07
% EndTime: 2024-09-27 21:45:07
% DurationCPUTime: 0.06s
% Computational Cost: add. (242->90), mult. (178->105), div. (0->0), fcn. (154->22), ass. (0->40)
t99 = pkin(7) + pkin(8);
t111 = pkin(9) + t99 + rSges(6,3);
t110 = t99 + rSges(5,3);
t109 = rSges(4,3) + pkin(7);
t98 = cos(qJ(3));
t78 = t98 * pkin(3) + pkin(2);
t95 = cos(pkin(5));
t97 = sin(qJ(3));
t107 = t95 * t97;
t106 = t95 * t98;
t105 = pkin(6) + qJ(1);
t91 = qJ(3) + qJ(4);
t82 = pkin(5) - t91;
t81 = pkin(5) + t91;
t104 = t98 * sin(qJ(4)) * pkin(4) + t97 * (cos(qJ(4)) * pkin(4) + pkin(3));
t84 = cos(t91);
t103 = t84 * rSges(5,1) - sin(t91) * rSges(5,2) + t78;
t87 = qJ(5) + t91;
t102 = cos(t87) * rSges(6,1) - sin(t87) * rSges(6,2) + pkin(4) * t84 + t78;
t67 = sin(t81) / 0.2e1;
t68 = cos(t82) / 0.2e1;
t71 = sin(t82);
t72 = cos(t81);
t93 = sin(pkin(5));
t101 = (t67 - t71 / 0.2e1) * rSges(5,1) + (t68 + t72 / 0.2e1) * rSges(5,2) + pkin(3) * t107 - t110 * t93;
t73 = qJ(5) + t81;
t65 = sin(t73) / 0.2e1;
t74 = -qJ(5) + t82;
t66 = cos(t74) / 0.2e1;
t69 = sin(t74);
t70 = cos(t73);
t100 = (t65 - t69 / 0.2e1) * rSges(6,1) + (t66 + t70 / 0.2e1) * rSges(6,2) + t104 * t95 - t111 * t93;
t94 = cos(pkin(10));
t92 = sin(pkin(10));
t89 = pkin(10) + qJ(2);
t86 = t94 * pkin(1);
t85 = t92 * pkin(1);
t80 = cos(t89);
t79 = sin(t89);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t94 * rSges(2,1) - t92 * rSges(2,2)) + g(2) * (t92 * rSges(2,1) + t94 * rSges(2,2)) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (t80 * rSges(3,1) - t79 * rSges(3,2) + t86) + g(2) * (t79 * rSges(3,1) + t80 * rSges(3,2) + t85) + g(3) * (rSges(3,3) + t105)) - m(4) * (g(1) * (t80 * pkin(2) + t86 + (-t79 * t107 + t80 * t98) * rSges(4,1) + (-t79 * t106 - t80 * t97) * rSges(4,2)) + g(2) * (t79 * pkin(2) + t85 + (t80 * t107 + t79 * t98) * rSges(4,1) + (t80 * t106 - t79 * t97) * rSges(4,2)) + g(3) * (t109 * t95 + t105) + (g(3) * (rSges(4,1) * t97 + rSges(4,2) * t98) + (g(1) * t79 - g(2) * t80) * t109) * t93) - m(5) * (g(1) * (-t101 * t79 + t103 * t80 + t86) + g(2) * (t101 * t80 + t103 * t79 + t85) + g(3) * (t93 * t97 * pkin(3) + (t68 - t72 / 0.2e1) * rSges(5,1) + (t67 + t71 / 0.2e1) * rSges(5,2) + t110 * t95 + t105)) - m(6) * (g(1) * (-t100 * t79 + t102 * t80 + t86) + g(2) * (t100 * t80 + t102 * t79 + t85) + g(3) * ((t66 - t70 / 0.2e1) * rSges(6,1) + (t65 + t69 / 0.2e1) * rSges(6,2) + t111 * t95 + t104 * t93 + t105));
U = t1;
