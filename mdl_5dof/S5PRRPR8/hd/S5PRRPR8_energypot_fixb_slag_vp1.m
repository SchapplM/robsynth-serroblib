% Calculate potential energy for
% S5PRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-31 17:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRPR8_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR8_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR8_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR8_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR8_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPR8_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:42:15
% EndTime: 2019-12-31 17:42:15
% DurationCPUTime: 0.32s
% Computational Cost: add. (150->70), mult. (130->86), div. (0->0), fcn. (110->10), ass. (0->32)
t77 = rSges(6,3) + pkin(7);
t62 = -pkin(6) - pkin(5);
t76 = rSges(3,3) + pkin(5);
t61 = cos(qJ(2));
t48 = t61 * pkin(2) + pkin(1);
t55 = qJ(2) + qJ(3);
t49 = pkin(9) + t55;
t45 = sin(t49);
t74 = rSges(5,2) * t45;
t56 = sin(pkin(8));
t58 = sin(qJ(5));
t73 = t56 * t58;
t60 = cos(qJ(5));
t72 = t56 * t60;
t46 = cos(t49);
t57 = cos(pkin(8));
t71 = t57 * t46;
t70 = t57 * t58;
t69 = t57 * t60;
t68 = rSges(4,3) - t62;
t51 = cos(t55);
t43 = pkin(3) * t51 + t48;
t54 = -qJ(4) + t62;
t67 = t56 * t43 + t57 * t54;
t59 = sin(qJ(2));
t66 = t59 * pkin(2) + qJ(1);
t50 = sin(t55);
t65 = pkin(3) * t50 + t66;
t64 = rSges(3,1) * t61 - rSges(3,2) * t59 + pkin(1);
t63 = rSges(4,1) * t51 - rSges(4,2) * t50 + t48;
t42 = t57 * t43;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t57 * rSges(2,1) - t56 * rSges(2,2)) + g(2) * (t56 * rSges(2,1) + t57 * rSges(2,2)) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(3) * (t59 * rSges(3,1) + t61 * rSges(3,2) + qJ(1)) + (g(1) * t64 - g(2) * t76) * t57 + (g(1) * t76 + g(2) * t64) * t56) - m(4) * (g(3) * (t50 * rSges(4,1) + t51 * rSges(4,2) + t66) + (g(1) * t63 - g(2) * t68) * t57 + (g(1) * t68 + g(2) * t63) * t56) - m(5) * (g(1) * (rSges(5,1) * t71 - t57 * t74 + t42) + g(2) * (-t57 * rSges(5,3) + t67) + g(3) * (t45 * rSges(5,1) + t46 * rSges(5,2) + t65) + (g(1) * (rSges(5,3) - t54) + g(2) * (rSges(5,1) * t46 - t74)) * t56) - m(6) * (g(1) * (pkin(4) * t71 + t42 - t56 * t54 + (t46 * t69 + t73) * rSges(6,1) + (-t46 * t70 + t72) * rSges(6,2)) + g(2) * (t56 * t46 * pkin(4) + (t46 * t72 - t70) * rSges(6,1) + (-t46 * t73 - t69) * rSges(6,2) + t67) + g(3) * (-t77 * t46 + t65) + (g(3) * (rSges(6,1) * t60 - rSges(6,2) * t58 + pkin(4)) + (g(1) * t57 + g(2) * t56) * t77) * t45);
U = t1;
