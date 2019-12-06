% Calculate potential energy for
% S5PRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
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
% Datum: 2019-12-05 16:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRPRR7_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR7_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR7_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR7_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR7_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR7_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:59:39
% EndTime: 2019-12-05 15:59:39
% DurationCPUTime: 0.41s
% Computational Cost: add. (116->84), mult. (172->109), div. (0->0), fcn. (156->8), ass. (0->28)
t80 = rSges(5,3) + pkin(6);
t79 = rSges(6,3) + pkin(7) + pkin(6);
t54 = sin(pkin(8));
t55 = cos(pkin(8));
t78 = g(1) * t55 + g(2) * t54;
t57 = sin(qJ(2));
t74 = t54 * t57;
t59 = cos(qJ(2));
t73 = t54 * t59;
t72 = t55 * t57;
t56 = sin(qJ(4));
t71 = t56 * t57;
t58 = cos(qJ(4));
t70 = t57 * t58;
t68 = t55 * pkin(1) + t54 * pkin(5);
t67 = qJ(3) * t57;
t66 = t57 * pkin(2) + qJ(1);
t65 = t54 * t71;
t64 = t55 * t71;
t50 = t54 * pkin(1);
t63 = pkin(2) * t73 + t54 * t67 + t50;
t62 = t68 + (pkin(2) * t59 + t67) * t55;
t61 = -t55 * pkin(5) + t63;
t53 = qJ(4) + qJ(5);
t48 = cos(t53);
t47 = sin(t53);
t46 = pkin(4) * t58 + pkin(3);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t55 - rSges(2,2) * t54) + g(2) * (rSges(2,1) * t54 + rSges(2,2) * t55) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,3) * t54 + t68) + g(2) * (rSges(3,1) * t73 - rSges(3,2) * t74 + t50) + g(3) * (rSges(3,1) * t57 + rSges(3,2) * t59 + qJ(1)) + (g(1) * (rSges(3,1) * t59 - rSges(3,2) * t57) + g(2) * (-rSges(3,3) - pkin(5))) * t55) - m(4) * (g(1) * (rSges(4,1) * t54 + t62) + g(2) * (-rSges(4,2) * t73 + rSges(4,3) * t74 + t63) + g(3) * (-rSges(4,2) * t57 + (-rSges(4,3) - qJ(3)) * t59 + t66) + (g(1) * (-rSges(4,2) * t59 + rSges(4,3) * t57) + g(2) * (-rSges(4,1) - pkin(5))) * t55) - m(5) * (g(1) * (t54 * pkin(3) + (t54 * t58 + t64) * rSges(5,1) + (-t54 * t56 + t55 * t70) * rSges(5,2) + t62) + g(2) * (-t55 * pkin(3) + (-t55 * t58 + t65) * rSges(5,1) + (t54 * t70 + t55 * t56) * rSges(5,2) + t61) + g(3) * (t80 * t57 + t66) + (g(3) * (-rSges(5,1) * t56 - rSges(5,2) * t58 - qJ(3)) + t78 * t80) * t59) - m(6) * (g(1) * (t54 * t46 + pkin(4) * t64 + (t47 * t72 + t48 * t54) * rSges(6,1) + (-t47 * t54 + t48 * t72) * rSges(6,2) + t62) + g(2) * (-t55 * t46 + pkin(4) * t65 + (t47 * t74 - t48 * t55) * rSges(6,1) + (t47 * t55 + t48 * t74) * rSges(6,2) + t61) + g(3) * (t79 * t57 + t66) + (g(3) * (-rSges(6,1) * t47 - rSges(6,2) * t48 - pkin(4) * t56 - qJ(3)) + t78 * t79) * t59);
U = t1;
