% Calculate potential energy for
% S5RPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRR7_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR7_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR7_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR7_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR7_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR7_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:03:08
% EndTime: 2019-12-31 19:03:08
% DurationCPUTime: 0.36s
% Computational Cost: add. (157->77), mult. (141->100), div. (0->0), fcn. (125->10), ass. (0->32)
t81 = rSges(5,3) + pkin(7);
t80 = rSges(6,3) + pkin(8) + pkin(7);
t53 = qJ(1) + pkin(9);
t47 = sin(t53);
t48 = cos(t53);
t79 = g(1) * t48 + g(2) * t47;
t56 = sin(qJ(3));
t75 = rSges(4,2) * t56;
t55 = sin(qJ(4));
t74 = t47 * t55;
t59 = cos(qJ(3));
t73 = t47 * t59;
t72 = t48 * t55;
t54 = qJ(4) + qJ(5);
t49 = sin(t54);
t71 = t49 * t59;
t50 = cos(t54);
t70 = t50 * t59;
t69 = t55 * t59;
t58 = cos(qJ(4));
t68 = t58 * t59;
t46 = t58 * pkin(4) + pkin(3);
t67 = t59 * t46;
t65 = pkin(5) + qJ(2);
t57 = sin(qJ(1));
t51 = t57 * pkin(1);
t64 = t47 * pkin(2) + t51;
t60 = cos(qJ(1));
t52 = t60 * pkin(1);
t63 = t48 * pkin(2) + t47 * pkin(6) + t52;
t62 = -t48 * pkin(6) + t64;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t60 * rSges(2,1) - t57 * rSges(2,2)) + g(2) * (t57 * rSges(2,1) + t60 * rSges(2,2)) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(1) * (t48 * rSges(3,1) - t47 * rSges(3,2) + t52) + g(2) * (t47 * rSges(3,1) + t48 * rSges(3,2) + t51) + g(3) * (rSges(3,3) + t65)) - m(4) * (g(1) * (t47 * rSges(4,3) + t63) + g(2) * (rSges(4,1) * t73 - t47 * t75 + t64) + g(3) * (t56 * rSges(4,1) + t59 * rSges(4,2) + t65) + (g(1) * (rSges(4,1) * t59 - t75) + g(2) * (-rSges(4,3) - pkin(6))) * t48) - m(5) * (g(1) * (t48 * t59 * pkin(3) + (t48 * t68 + t74) * rSges(5,1) + (t47 * t58 - t48 * t69) * rSges(5,2) + t63) + g(2) * (pkin(3) * t73 + (t47 * t68 - t72) * rSges(5,1) + (-t47 * t69 - t48 * t58) * rSges(5,2) + t62) + g(3) * (-t81 * t59 + t65) + (g(3) * (rSges(5,1) * t58 - rSges(5,2) * t55 + pkin(3)) + t79 * t81) * t56) - m(6) * (g(1) * (t48 * t67 + pkin(4) * t74 + (t47 * t49 + t48 * t70) * rSges(6,1) + (t47 * t50 - t48 * t71) * rSges(6,2) + t63) + g(2) * (t47 * t67 - pkin(4) * t72 + (t47 * t70 - t48 * t49) * rSges(6,1) + (-t47 * t71 - t48 * t50) * rSges(6,2) + t62) + g(3) * (-t80 * t59 + t65) + (g(3) * (rSges(6,1) * t50 - rSges(6,2) * t49 + t46) + t79 * t80) * t56);
U = t1;
