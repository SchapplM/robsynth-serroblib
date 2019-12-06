% Calculate potential energy for
% S5PPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
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
% Datum: 2019-12-05 14:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PPPRR2_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR2_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPPRR2_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR2_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPPRR2_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPPRR2_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:59:22
% EndTime: 2019-12-05 14:59:22
% DurationCPUTime: 0.31s
% Computational Cost: add. (158->92), mult. (313->125), div. (0->0), fcn. (348->10), ass. (0->39)
t106 = rSges(5,3) + pkin(5);
t105 = pkin(6) + rSges(6,3);
t80 = sin(pkin(9));
t81 = sin(pkin(8));
t104 = t80 * t81;
t82 = sin(pkin(7));
t103 = t81 * t82;
t87 = sin(qJ(4));
t102 = t81 * t87;
t89 = cos(qJ(4));
t101 = t81 * t89;
t84 = cos(pkin(8));
t100 = t82 * t84;
t85 = cos(pkin(7));
t99 = t85 * t80;
t83 = cos(pkin(9));
t98 = t85 * t83;
t97 = t85 * pkin(1) + t82 * qJ(2);
t96 = qJ(3) * t81;
t95 = t81 * pkin(2) + qJ(1);
t94 = t97 + (pkin(2) * t84 + t96) * t85;
t64 = t82 * t80 + t84 * t98;
t93 = t64 * pkin(3) + t94;
t78 = t82 * pkin(1);
t92 = pkin(2) * t100 - t85 * qJ(2) + t82 * t96 + t78;
t62 = t83 * t100 - t99;
t91 = t62 * pkin(3) + t92;
t90 = t81 * t83 * pkin(3) + pkin(5) * t104 - t84 * qJ(3) + t95;
t88 = cos(qJ(5));
t86 = sin(qJ(5));
t66 = t83 * t101 - t84 * t87;
t65 = t83 * t102 + t84 * t89;
t63 = -t82 * t83 + t84 * t99;
t61 = t80 * t100 + t98;
t58 = t85 * t102 + t64 * t89;
t57 = -t85 * t101 + t64 * t87;
t56 = t82 * t102 + t62 * t89;
t55 = -t82 * t101 + t62 * t87;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t85 * rSges(2,1) - t82 * rSges(2,2)) + g(2) * (t82 * rSges(2,1) + t85 * rSges(2,2)) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (t82 * rSges(3,3) + t97) + g(2) * (rSges(3,1) * t100 - rSges(3,2) * t103 + t78) + g(3) * (t81 * rSges(3,1) + t84 * rSges(3,2) + qJ(1)) + (g(1) * (rSges(3,1) * t84 - rSges(3,2) * t81) + g(2) * (-rSges(3,3) - qJ(2))) * t85) - m(4) * (g(1) * (t85 * t81 * rSges(4,3) + t64 * rSges(4,1) - t63 * rSges(4,2) + t94) + g(2) * (t62 * rSges(4,1) - t61 * rSges(4,2) + rSges(4,3) * t103 + t92) + g(3) * ((-rSges(4,3) - qJ(3)) * t84 + (rSges(4,1) * t83 - rSges(4,2) * t80) * t81 + t95)) - m(5) * (g(1) * (t58 * rSges(5,1) - t57 * rSges(5,2) + t106 * t63 + t93) + g(2) * (t56 * rSges(5,1) - t55 * rSges(5,2) + t106 * t61 + t91) + g(3) * (t66 * rSges(5,1) - t65 * rSges(5,2) + rSges(5,3) * t104 + t90)) - m(6) * (g(1) * (t58 * pkin(4) + t63 * pkin(5) + (t58 * t88 + t63 * t86) * rSges(6,1) + (-t58 * t86 + t63 * t88) * rSges(6,2) + t105 * t57 + t93) + g(2) * (t56 * pkin(4) + t61 * pkin(5) + (t56 * t88 + t61 * t86) * rSges(6,1) + (-t56 * t86 + t61 * t88) * rSges(6,2) + t105 * t55 + t91) + g(3) * (t66 * pkin(4) + (t86 * t104 + t66 * t88) * rSges(6,1) + (t88 * t104 - t66 * t86) * rSges(6,2) + t105 * t65 + t90));
U = t1;
