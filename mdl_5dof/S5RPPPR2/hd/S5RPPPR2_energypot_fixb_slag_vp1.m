% Calculate potential energy for
% S5RPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% Datum: 2020-01-03 11:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPPR2_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR2_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR2_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR2_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:22:15
% EndTime: 2020-01-03 11:22:16
% DurationCPUTime: 0.46s
% Computational Cost: add. (158->92), mult. (313->122), div. (0->0), fcn. (348->10), ass. (0->43)
t78 = sin(pkin(7));
t81 = cos(pkin(7));
t108 = t81 * pkin(2) + t78 * qJ(3);
t107 = pkin(6) + rSges(6,3);
t105 = t78 * pkin(2) + pkin(5);
t104 = g(2) * qJ(2);
t77 = sin(pkin(8));
t103 = t77 * t78;
t80 = cos(pkin(8));
t102 = t78 * t80;
t83 = sin(qJ(1));
t101 = t78 * t83;
t100 = t83 * t77;
t99 = t83 * t80;
t85 = cos(qJ(1));
t98 = t85 * t77;
t97 = t85 * t78;
t96 = t85 * t80;
t94 = t83 * qJ(2);
t93 = -rSges(3,3) - qJ(2);
t75 = t83 * pkin(1);
t92 = t108 * t83 + t75;
t91 = rSges(3,1) * t81 - rSges(3,2) * t78;
t90 = -pkin(1) - t108;
t66 = t81 * t98 - t99;
t67 = -t81 * t96 - t100;
t89 = t67 * pkin(3) - t66 * qJ(4) - t94;
t88 = pkin(3) * t102 - t81 * qJ(3) + qJ(4) * t103 + t105;
t64 = t81 * t100 + t96;
t65 = t81 * t99 - t98;
t87 = t65 * pkin(3) + t64 * qJ(4) + t92;
t86 = (g(3) * t90 - t104) * t85;
t84 = cos(qJ(5));
t82 = sin(qJ(5));
t79 = cos(pkin(9));
t76 = sin(pkin(9));
t63 = t79 * t102 - t81 * t76;
t62 = t76 * t102 + t81 * t79;
t59 = t67 * t79 - t76 * t97;
t58 = t67 * t76 + t79 * t97;
t57 = t76 * t101 + t65 * t79;
t56 = -t79 * t101 + t65 * t76;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (pkin(5) + rSges(2,3)) + g(2) * (t83 * rSges(2,1) + t85 * rSges(2,2)) + g(3) * (-t85 * rSges(2,1) + t83 * rSges(2,2))) - m(3) * (g(1) * (t78 * rSges(3,1) + t81 * rSges(3,2) + pkin(5)) + g(2) * t75 + (g(2) * t91 + g(3) * t93) * t83 + (g(2) * t93 + g(3) * (-pkin(1) - t91)) * t85) - m(4) * (g(1) * ((-rSges(4,3) - qJ(3)) * t81 + (rSges(4,1) * t80 - rSges(4,2) * t77) * t78 + t105) + g(2) * (t65 * rSges(4,1) - t64 * rSges(4,2) + rSges(4,3) * t101 + t92) + g(3) * (t67 * rSges(4,1) + t66 * rSges(4,2) - t94) + (-t104 + g(3) * (-rSges(4,3) * t78 + t90)) * t85) - m(5) * (g(1) * (t63 * rSges(5,1) - t62 * rSges(5,2) + rSges(5,3) * t103 + t88) + g(2) * (t57 * rSges(5,1) - t56 * rSges(5,2) + t64 * rSges(5,3) + t87) + g(3) * (t59 * rSges(5,1) - t58 * rSges(5,2) - t66 * rSges(5,3) + t89) + t86) - m(6) * (g(1) * (t63 * pkin(4) + (t82 * t103 + t63 * t84) * rSges(6,1) + (t84 * t103 - t63 * t82) * rSges(6,2) + t107 * t62 + t88) + g(2) * (t57 * pkin(4) + (t57 * t84 + t64 * t82) * rSges(6,1) + (-t57 * t82 + t64 * t84) * rSges(6,2) + t87 + t107 * t56) + g(3) * (t59 * pkin(4) + (t59 * t84 - t66 * t82) * rSges(6,1) + (-t59 * t82 - t66 * t84) * rSges(6,2) + t89 + t107 * t58) + t86);
U = t1;
